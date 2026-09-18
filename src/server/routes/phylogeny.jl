# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/reference-trees and /api/v1/studies/{study}/placements
#
# Reference trees are a library shared by every study, each a directory
# reference_trees/{id}/ holding reference.json (name and settings that override
# the machine's phylogeny section), the reference sequences and the output of
# each step. A placement belongs to one study: projects/{study}/phylogeny/{id}/
# holds placement.json (name, reference tree, where the queries come from and
# settings), the queries and the output of each step. A finished placement copies
# its trees into projects/{study}/trees/ for the viewer.
using JSON3, Dates, UUIDs, YAML, SHA
const Phylo = MetaManifold.Phylogeny

const _FASTA_MAX_BYTES = 100 * 1024 * 1024
const _phylo_lock = ReentrantLock()
# The active job of each reference tree ("ref:<id>") and placement ("pl:<study>/<id>").
const _phylo_jobs = Dict{String,String}()

_valid_phylo_id(id::String) = occursin(r"^[0-9a-f]{8}$", id)
_new_phylo_id() = string(uuid4())[1:8]
_stamp() = string(now(UTC)) * "Z"

_refs_root() = joinpath(dirname(ServerState.data_dir()), "reference_trees")
_ref_dir(id::String) = joinpath(_refs_root(), id)
_ref_doc_path(id::String) = joinpath(_ref_dir(id), "reference.json")

_placements_dir(study::String) = joinpath(ServerState.projects_dir(), study, "phylogeny")
_placement_dir(study::String, id::String) = joinpath(_placements_dir(study), id)
_placement_doc_path(study::String, id::String) = joinpath(_placement_dir(study, id), "placement.json")

_read_doc(path) = try JSON3.read(read(path, String), Dict{String,Any}) catch; nothing end
function _write_doc(path, doc::AbstractDict)
    doc["modified"] = _stamp()
    _write_atomic(path, JSON3.write(doc))
end

# The machine's cascade, which reference trees read: factory defaults and config/pipeline.yml.
function _machine_config()
    config_dir = joinpath(dirname(ServerState.data_dir()), "config")
    load_merged_config([joinpath(config_dir, "defaults", "pipeline.yml"), joinpath(config_dir, "pipeline.yml")])
end

# A study's cascade: the machine's, then the study's pipeline.yml.
function _study_config(study::String)
    config_dir = joinpath(dirname(ServerState.data_dir()), "config")
    study_dir  = joinpath(ServerState.data_dir(), study)
    load_merged_config(config_dir, study_dir, study_dir)
end

_fasta_count(path) = isfile(path) ? count(l -> startswith(l, '>'), eachline(path)) : 0

function _active_job(key::String)
    lock(_phylo_lock) do
        jid = get(_phylo_jobs, key, nothing)
        isnothing(jid) && return nothing
        job = get_job(jid)
        (isnothing(job) || JobQueue.is_terminal(job.status)) ? nothing : job
    end
end

function _state(dir, key)
    job = _active_job(key)
    st  = Phylo.read_status(dir)
    !isnothing(job) ? "running" : isnothing(st) ? "new" : string(get(st, "state", "new"))
end

# Submit `f` as the one job of `key`, or return the job already running for it.
function _submit_phylo(f::Function, key::String; study=nothing)
    lock(_phylo_lock) do
        existing = let jid = get(_phylo_jobs, key, nothing)
            isnothing(jid) ? nothing : get_job(jid)
        end
        !isnothing(existing) && !JobQueue.is_terminal(existing.status) && return existing
        j = submit_job!(f, "phylogeny"; study)
        _phylo_jobs[key] = j.id
        j
    end
end

function _fail_status!(dir, e)
    Phylo._update_status!(dir, st -> begin
        st["state"] = "failed"
        st["error"] = sprint(showerror, e)
        st["finished"] = _stamp()
        st["steps"] = Dict{String,Any}()
    end)
end

# What the page shows of a workflow besides its own document.
function _workflow_detail(steps, dir, files, cfg, overrides, key)
    inherited = Phylo.placement_settings(cfg)
    settings  = try Phylo.placement_settings(cfg, overrides) catch; inherited end
    threads   = Phylo._threads(get(inherited, "threads", nothing), 4)
    remote = Dict(s.name => (t = try Phylo.remote_step_target(cfg, s.remote; threads) catch; nothing end;
                             isnothing(t) ? nothing : t.host) for s in steps)
    qc = [s.name for s in steps if isfile(joinpath(dir, "qc", "$(s.name).json"))]
    (; status = Phylo.read_status(dir), state = _state(dir, key), settings, inherited, remote, qc)
end

_step_named(steps, name) = let i = findfirst(s -> s.name == name, steps)
    isnothing(i) ? nothing : steps[i]
end

function _log_response(steps, dir, step)
    isnothing(_step_named(steps, step)) && return json_error(400, "invalid_step", "Unknown step '$step'")
    path = joinpath(dir, "logs", "$step.log")
    isfile(path) || return json_error(404, "not_found", "No log for $step yet")
    HTTP.Response(200, ["Content-Type" => "text/plain; charset=utf-8"];
                  body = join(last(readlines(path), 400), "\n"))
end

function _qc_response(steps, dir, step)
    isnothing(_step_named(steps, step)) && return json_error(400, "invalid_step", "Unknown step '$step'")
    path = joinpath(dir, "qc", "$step.json")
    isfile(path) || return json_error(404, "not_found", "No QC for $step yet")
    HTTP.Response(200, ["Content-Type" => "application/json"]; body = read(path))
end

# raw is the align step's output, trimmed the trim step's.
function _alignment_response(steps, files, which)
    which in ("raw", "trimmed") || return json_error(400, "invalid", "raw or trimmed")
    name = steps === Phylo.PLACEMENT_STEPS ?
           (which == "raw" ? "combined.aln.fasta" : "combined.trim.fasta") :
           (which == "raw" ? "reference.aln.fasta" : "reference.trim.fasta")
    isfile(files[name]) || return json_error(404, "not_found", "No $which alignment yet")
    HTTP.Response(200, ["Content-Type" => "text/plain; charset=utf-8"]; body = read(files[name]))
end

function _trim_preview_response(req, steps, dir, files)
    body = try JSON3.read(String(req.body), Dict{String,Any}) catch
        return json_error(400, "invalid_body", "Expected a JSON object")
    end
    try
        json(Phylo.trim_preview(steps, dir, files, body))
    catch e
        json_error(400, "trim_failed", sprint(showerror, e))
    end
end

function _check_overrides(cfg, ov, allowed::String)
    ov isa AbstractDict || return "overrides must be an object"
    all(k -> string(k) in (allowed, "threads"), keys(ov)) || return "overrides may set phylogeny.$allowed only"
    try
        Phylo.placement_settings(cfg, ov)
        nothing
    catch e
        sprint(showerror, e)
    end
end

function _put_fasta(req, path, what)
    length(req.body) <= _FASTA_MAX_BYTES || return json_error(413, "too_large", "FASTA uploads are limited to 100 MB")
    text = String(req.body)
    records = try Phylo.parse_fasta(text; what) catch e
        return json_error(400, "invalid_fasta", sprint(showerror, e))
    end
    renamed = count(l -> startswith(l, '>') && Phylo.safe_name(l[2:end]) != strip(l[2:end]),
                    eachline(IOBuffer(text)))
    Phylo.write_fasta(path, records)
    json((; count = length(records), renamed))
end

function _get_fasta(path, name)
    isfile(path) || return json_error(404, "not_found", "No $name yet")
    HTTP.Response(200, ["Content-Type" => "text/plain; charset=utf-8",
                        "Content-Disposition" => "attachment; filename=\"$name.fasta\""]; body=read(path))
end

## Reference trees
function _ref_request_error(id::Union{String,Nothing}=nothing)
    isnothing(id) && return nothing
    _valid_phylo_id(id) || return json_error(400, "invalid_id", "Reference tree ids are 8 hex characters")
    isfile(_ref_doc_path(id)) || return json_error(404, "not_found", "No reference tree '$id'")
    nothing
end

# The placements, in any study, built on reference tree `id`.
function _ref_users(id::String)
    users = NamedTuple[]
    isdir(ServerState.projects_dir()) || return users
    for study in readdir(ServerState.projects_dir())
        dir = _placements_dir(study)
        isdir(dir) || continue
        for pid in readdir(dir)
            doc = _read_doc(joinpath(dir, pid, "placement.json"))
            (doc isa AbstractDict && get(doc, "reference", nothing) == id) || continue
            push!(users, (; study, id = pid, name = string(get(doc, "name", pid))))
        end
    end
    users
end

function _ref_summary(id::String)
    doc = something(_read_doc(_ref_doc_path(id)), Dict{String,Any}())
    dir = _ref_dir(id)
    (; id,
       name       = string(get(doc, "name", id)),
       modified   = string(get(doc, "modified", "")),
       references = _fasta_count(joinpath(dir, "references.fasta")),
       built      = isfile(Phylo.reference_files(dir)["reference.treefile"]),
       state      = _state(dir, "ref:$id"))
end

@get "/api/v1/reference-trees" function(req)
    root = _refs_root()
    isdir(root) || return json([])
    ids = filter(id -> _valid_phylo_id(id) && isfile(_ref_doc_path(id)), readdir(root))
    json(sort([_ref_summary(id) for id in ids]; by = r -> lowercase(r.name)))
end

@post "/api/v1/reference-trees" function(req)
    body = try JSON3.read(String(req.body), Dict{String,Any}) catch; Dict{String,Any}() end
    name = strip(string(get(body, "name", "")))
    isempty(name) && return json_error(400, "invalid_name", "A reference tree needs a name")
    id = _new_phylo_id()
    doc = Dict{String,Any}("id" => id, "name" => name, "description" => "",
                           "overrides" => Dict{String,Any}(), "created" => _stamp())
    _write_doc(_ref_doc_path(id), doc)
    json(doc)
end

@get "/api/v1/reference-trees/{id}" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    doc = _read_doc(_ref_doc_path(id))
    dir = _ref_dir(id)
    d = _workflow_detail(Phylo.REFERENCE_STEPS, dir, Phylo.reference_files(dir), _machine_config(),
                         get(doc, "overrides", nothing), "ref:$id")
    json((; doc, summary = _ref_summary(id), d..., used_by = _ref_users(id)))
end

@put "/api/v1/reference-trees/{id}" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    body = try JSON3.read(String(req.body), Dict{String,Any}) catch
        return json_error(400, "invalid_body", "Expected a JSON object")
    end
    doc = _read_doc(_ref_doc_path(id))
    if haskey(body, "name")
        name = strip(string(body["name"]))
        isempty(name) && return json_error(400, "invalid_name", "A reference tree needs a name")
        doc["name"] = name
    end
    haskey(body, "description") && (doc["description"] = string(body["description"]))
    if haskey(body, "overrides")
        problem = _check_overrides(_machine_config(), body["overrides"], "reference")
        isnothing(problem) || return json_error(400, "invalid_overrides", problem)
        doc["overrides"] = body["overrides"]
    end
    _write_doc(_ref_doc_path(id), doc)
    json(doc)
end

@delete "/api/v1/reference-trees/{id}" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    isnothing(_active_job("ref:$id")) || return json_error(409, "running", "The reference tree is being built")
    users = _ref_users(id)
    isempty(users) || return json_error(409, "in_use",
        "Used by " * join(("$(u.study)/$(u.name)" for u in users), ", ") * "; delete those placements first")
    rm(_ref_dir(id); recursive=true, force=true)
    json((; deleted = id))
end

@put "/api/v1/reference-trees/{id}/fasta" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    isnothing(_active_job("ref:$id")) || return json_error(409, "running", "The reference tree is being built")
    _put_fasta(req, Phylo.reference_files(_ref_dir(id))["references.fasta"], "References")
end

@get "/api/v1/reference-trees/{id}/fasta" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    _get_fasta(Phylo.reference_files(_ref_dir(id))["references.fasta"], "references")
end

@post "/api/v1/reference-trees/{id}/run" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    dir   = _ref_dir(id)
    files = Phylo.reference_files(dir)
    isfile(files["references.fasta"]) || return json_error(400, "no_references", "Upload the reference sequences first")
    job = _submit_phylo("ref:$id") do
        doc = _read_doc(_ref_doc_path(id))
        try
            Phylo.check_reference(files)
        catch e
            _fail_status!(dir, e)
            rethrow()
        end
        Phylo.run_workflow(Phylo.REFERENCE_STEPS, dir, files, _machine_config();
            overrides = get(doc, "overrides", nothing),
            run = Dict{String,Any}("reference_tree" => id, "name" => string(get(doc, "name", id))))
    end
    json(_job_to_namedtuple(job))
end

@get "/api/v1/reference-trees/{id}/log/{step}" function(req, id::String, step::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    _log_response(Phylo.REFERENCE_STEPS, _ref_dir(id), step)
end

@get "/api/v1/reference-trees/{id}/qc/{step}" function(req, id::String, step::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    _qc_response(Phylo.REFERENCE_STEPS, _ref_dir(id), step)
end

@get "/api/v1/reference-trees/{id}/alignment/{which}" function(req, id::String, which::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    _alignment_response(Phylo.REFERENCE_STEPS, Phylo.reference_files(_ref_dir(id)), which)
end

@get "/api/v1/reference-trees/{id}/treefile" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    path = Phylo.reference_files(_ref_dir(id))["reference.treefile"]
    isfile(path) || return json_error(404, "not_found", "The tree has not been built yet")
    HTTP.Response(200, ["Content-Type" => "text/plain; charset=utf-8"]; body = read(path))
end

# The tree as the viewer reads a study's trees, with its saved view beside it.
_ref_view_path(id::String) = joinpath(_ref_dir(id), "tree", "view.json")

@get "/api/v1/reference-trees/{id}/tree" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    path = Phylo.reference_files(_ref_dir(id))["reference.treefile"]
    isfile(path) || return json_error(404, "not_found", "The tree has not been built yet")
    content = read(path, String)
    vpath = _ref_view_path(id)
    view = isfile(vpath) ? try JSON3.read(read(vpath, String)) catch; nothing end : nothing
    name = string(get(something(_read_doc(_ref_doc_path(id)), Dict()), "name", id))
    json((; file = "$(_phylo_slug(name)).treefile", format = "newick",
            sha256 = bytes2hex(sha256(content)),
            modified = string(unix2datetime(stat(path).mtime)), content, view))
end

@put "/api/v1/reference-trees/{id}/tree/view" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    text = String(req.body)
    sizeof(text) <= _VIEW_MAX_BYTES || return json_error(413, "too_large", "View state is limited to 5 MB")
    doc = try JSON3.read(text) catch
        return json_error(400, "invalid_body", "Expected a JSON object")
    end
    doc isa JSON3.Object || return json_error(400, "invalid_body", "Expected a JSON object")
    _write_atomic(_ref_view_path(id), text)
    json((; saved = true))
end

@post "/api/v1/reference-trees/{id}/trim-preview" function(req, id::String)
    err = _ref_request_error(id)
    isnothing(err) || return err
    dir = _ref_dir(id)
    _trim_preview_response(req, Phylo.REFERENCE_STEPS, dir, Phylo.reference_files(dir))
end

## Placements
function _placement_request_error(study::String, id::Union{String,Nothing}=nothing)
    _study_exists(study) || return json_error(404, "study_not_found", "Study '$study' not found")
    isnothing(id) && return nothing
    _valid_phylo_id(id) || return json_error(400, "invalid_id", "Placement ids are 8 hex characters")
    isfile(_placement_doc_path(study, id)) || return json_error(404, "not_found", "No placement '$id'")
    nothing
end

function _placement_files(study::String, id::String, doc)
    ref = string(something(get(doc, "reference", nothing), ""))
    Phylo.placement_files(_placement_dir(study, id), _valid_phylo_id(ref) ? _ref_dir(ref) : "")
end

function _placement_summary(study::String, id::String)
    doc = something(_read_doc(_placement_doc_path(study, id)), Dict{String,Any}())
    dir = _placement_dir(study, id)
    ref = get(doc, "reference", nothing)
    (; id,
       name      = string(get(doc, "name", id)),
       modified  = string(get(doc, "modified", "")),
       reference = ref,
       queries   = _fasta_count(joinpath(dir, "queries.fasta")),
       state     = _state(dir, "pl:$study/$id"))
end

## Taxon queries
# The ASVs of the chosen runs whose rank column holds one of the chosen values
# and which reach min_reads, counting only the samples of the chosen subgroups
# when a run names any. Tips are named by run so the same SeqName from two runs
# stays two tips.
function _taxon_queries(study::String, q::AbstractDict)
    runs   = get(q, "runs", Any[])
    table  = string(get(q, "table", ""))
    rank   = string(get(q, "rank", ""))
    values = string.(collect(get(q, "values", String[])))
    min_reads = get(q, "min_reads", 1)
    min_reads isa Real || (min_reads = 1)
    isempty(runs)   && error("Choose at least one run")
    isempty(table)  && error("Choose a results table")
    isempty(rank)   && error("Choose a rank")
    isempty(values) && error("Choose at least one taxon")

    records = Pair{String,String}[]
    per_run = NamedTuple[]
    for spec in runs
        run   = string(get(spec, "run", ""))
        group = let g = get(spec, "group", nothing); isnothing(g) ? nothing : string(g) end
        subgroups = string.(collect(something(get(spec, "subgroups", nothing), String[])))
        label = isnothing(group) ? run : "$group/$run"
        err = _validate_run_request(study, run)
        isnothing(err) || error("Run '$label' not found")
        dir = _require_duckdb(study, run; group)
        isnothing(dir) && error("Run '$label' has no results yet")
        rows = with_results_db(dir) do con
            cols = _duckdb_columns(con, table)
            isempty(cols) && error("Run '$label' has no table '$table'")
            rank in cols || error("Table '$table' of run '$label' has no column '$rank'")
            seq_col = findfirst(c -> lowercase(c) == "sequence", cols)
            isnothing(seq_col) && error("Table '$table' of run '$label' has no Sequence column")
            id_col = something(findfirst(c -> c in ("SeqName", "OTU", "ASV"), cols), 1)
            counts = _filter_by_prefix(_sample_count_columns(con, table), subgroups)
            isempty(counts) && !isempty(subgroups) &&
                error("Run '$label' has no samples in $(join(subgroups, ", "))")
            total  = isempty(counts) ? "0" : join(["COALESCE(\"$c\", 0)" for c in counts], " + ")
            marks  = join(fill("?", length(values)), ", ")
            sql = "SELECT CAST(\"$(cols[id_col])\" AS VARCHAR) AS id, CAST(\"$(cols[seq_col])\" AS VARCHAR) AS seq " *
                  "FROM \"$table\" WHERE CAST(\"$rank\" AS VARCHAR) IN ($marks) AND ($total) >= ? ORDER BY 1"
            DataFrame(DBInterface.execute(con, sql, Any[values..., min_reads]))
        end
        prefix = replace(label, r"[^A-Za-z0-9._-]+" => "_")
        n = 0
        for r in eachrow(rows)
            (ismissing(r.seq) || isempty(r.seq)) && continue
            push!(records, "$(prefix)_$(r.id)" => String(r.seq))
            n += 1
        end
        push!(per_run, (; run, group, subgroups, count = n))
    end
    (; records, per_run)
end

function _write_taxon_queries!(study::String, id::String, q::AbstractDict)
    found = _taxon_queries(study, q)
    isempty(found.records) && error("No ASVs in the chosen runs match $(get(q, "rank", "")) = " *
                                    join(string.(collect(get(q, "values", String[]))), ", "))
    text = join((">$(n)\n$(s)" for (n, s) in found.records), "\n")
    records = Phylo.parse_fasta(text; what="Queries")
    Phylo.write_fasta(joinpath(_placement_dir(study, id), "queries.fasta"), records)
    found
end

## Publishing
_phylo_slug(name) = let s = strip(first(replace(strip(name), r"[^A-Za-z0-9._-]+" => "_"), 60), ['_', '.'])
    isempty(s) ? "placement" : s
end

function _publish_placement(study::String, name::String, ref_name::String, files)
    slug = _phylo_slug(name)
    wanted = ["reference.treefile" => "$(_phylo_slug(ref_name)).treefile",
              "placement.jplace"   => "$slug.jplace",
              "accumulated.jplace" => "$slug-accumulated.jplace"]
    published = String[]
    for (k, file) in wanted
        isfile(files[k]) || continue
        dest = joinpath(_trees_dir(study), file)
        mkpath(dirname(dest))
        cp(files[k], dest * ".tmp"; force=true)
        mv(dest * ".tmp", dest; force=true)
        push!(published, file)
    end
    published
end

@get "/api/v1/studies/{study}/placements" function(req, study::String)
    err = _placement_request_error(study)
    isnothing(err) || return err
    dir = _placements_dir(study)
    isdir(dir) || return json([])
    ids = filter(id -> _valid_phylo_id(id) && isfile(_placement_doc_path(study, id)), readdir(dir))
    json(sort([_placement_summary(study, id) for id in ids]; by = p -> lowercase(p.name)))
end

@post "/api/v1/studies/{study}/placements" function(req, study::String)
    err = _placement_request_error(study)
    isnothing(err) || return err
    body = try JSON3.read(String(req.body), Dict{String,Any}) catch; Dict{String,Any}() end
    name = strip(string(get(body, "name", "")))
    isempty(name) && return json_error(400, "invalid_name", "A placement needs a name")
    ref = get(body, "reference", nothing)
    isnothing(ref) || (ref isa AbstractString && isnothing(_ref_request_error(String(ref)))) ||
        return json_error(400, "invalid_reference", "No reference tree '$ref'")
    id = _new_phylo_id()
    doc = Dict{String,Any}(
        "id" => id, "name" => name, "reference" => ref,
        "queries" => Dict{String,Any}("source" => "taxon", "runs" => Any[], "table" => "merged",
                                      "rank" => "", "values" => Any[], "min_reads" => 1),
        "overrides" => Dict{String,Any}(), "created" => _stamp())
    _write_doc(_placement_doc_path(study, id), doc)
    json(doc)
end

@get "/api/v1/studies/{study}/placements/{id}" function(req, study::String, id::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    doc = _read_doc(_placement_doc_path(study, id))
    isnothing(doc) && return json_error(500, "unreadable", "placement.json could not be read")
    dir = _placement_dir(study, id)
    d = _workflow_detail(Phylo.PLACEMENT_STEPS, dir, _placement_files(study, id, doc), _study_config(study),
                         get(doc, "overrides", nothing), "pl:$study/$id")
    json((; doc, summary = _placement_summary(study, id), d...))
end

@put "/api/v1/studies/{study}/placements/{id}" function(req, study::String, id::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    body = try JSON3.read(String(req.body), Dict{String,Any}) catch
        return json_error(400, "invalid_body", "Expected a JSON object")
    end
    doc = _read_doc(_placement_doc_path(study, id))
    if haskey(body, "name")
        name = strip(string(body["name"]))
        isempty(name) && return json_error(400, "invalid_name", "A placement needs a name")
        doc["name"] = name
    end
    if haskey(body, "reference")
        ref = body["reference"]
        isnothing(ref) || (ref isa AbstractString && isnothing(_ref_request_error(String(ref)))) ||
            return json_error(400, "invalid_reference", "No reference tree '$ref'")
        doc["reference"] = ref
    end
    if haskey(body, "queries")
        q = body["queries"]
        (q isa AbstractDict && get(q, "source", "") in ("taxon", "fasta")) ||
            return json_error(400, "invalid_queries", "queries.source must be taxon or fasta")
        doc["queries"] = q
    end
    if haskey(body, "overrides")
        problem = _check_overrides(_study_config(study), body["overrides"], "placement")
        isnothing(problem) || return json_error(400, "invalid_overrides", problem)
        doc["overrides"] = body["overrides"]
    end
    _write_doc(_placement_doc_path(study, id), doc)
    json(doc)
end

@delete "/api/v1/studies/{study}/placements/{id}" function(req, study::String, id::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    isnothing(_active_job("pl:$study/$id")) ||
        return json_error(409, "running", "The placement is running; cancel its job first")
    rm(_placement_dir(study, id); recursive=true, force=true)
    json((; deleted = id))
end

@put "/api/v1/studies/{study}/placements/{id}/fasta/queries" function(req, study::String, id::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    _put_fasta(req, joinpath(_placement_dir(study, id), "queries.fasta"), "Queries")
end

@get "/api/v1/studies/{study}/placements/{id}/fasta/queries" function(req, study::String, id::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    _get_fasta(joinpath(_placement_dir(study, id), "queries.fasta"), "queries")
end

# How many ASVs a taxon selection finds, per run, without writing anything.
@post "/api/v1/studies/{study}/placements/preview" function(req, study::String)
    err = _placement_request_error(study)
    isnothing(err) || return err
    q = try JSON3.read(String(req.body), Dict{String,Any}) catch
        return json_error(400, "invalid_body", "Expected a JSON object")
    end
    found = try _taxon_queries(study, q) catch e
        return json_error(400, "invalid_queries", sprint(showerror, e))
    end
    json((; count = length(found.records), per_run = found.per_run))
end

@get "/api/v1/studies/{study}/placements/{id}/log/{step}" function(req, study::String, id::String, step::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    _log_response(Phylo.PLACEMENT_STEPS, _placement_dir(study, id), step)
end

@get "/api/v1/studies/{study}/placements/{id}/qc/{step}" function(req, study::String, id::String, step::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    _qc_response(Phylo.PLACEMENT_STEPS, _placement_dir(study, id), step)
end

@get "/api/v1/studies/{study}/placements/{id}/alignment/{which}" function(req, study::String, id::String, which::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    doc = _read_doc(_placement_doc_path(study, id))
    _alignment_response(Phylo.PLACEMENT_STEPS, _placement_files(study, id, doc), which)
end

@post "/api/v1/studies/{study}/placements/{id}/trim-preview" function(req, study::String, id::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    doc = _read_doc(_placement_doc_path(study, id))
    _trim_preview_response(req, Phylo.PLACEMENT_STEPS, _placement_dir(study, id), _placement_files(study, id, doc))
end

@post "/api/v1/studies/{study}/placements/{id}/run" function(req, study::String, id::String)
    err = _placement_request_error(study, id)
    isnothing(err) || return err
    doc = _read_doc(_placement_doc_path(study, id))
    ref = get(doc, "reference", nothing)
    (ref isa AbstractString && isnothing(_ref_request_error(String(ref)))) ||
        return json_error(400, "no_reference", "Choose a reference tree first")
    isnothing(_active_job("ref:$ref")) ||
        return json_error(409, "reference_running", "The reference tree is being rebuilt; run the placement when it finishes")
    isfile(Phylo.reference_files(_ref_dir(String(ref)))["reference.treefile"]) ||
        return json_error(400, "reference_not_built", "Build the reference tree first")
    dir = _placement_dir(study, id)
    q = get(doc, "queries", Dict{String,Any}())
    get(q, "source", "taxon") == "fasta" && !isfile(joinpath(dir, "queries.fasta")) &&
        return json_error(400, "no_queries", "Upload the query sequences first")

    job = _submit_phylo("pl:$study/$id"; study) do
        files = _placement_files(study, id, doc)
        name = string(get(doc, "name", id))
        try
            get(q, "source", "taxon") == "taxon" && _write_taxon_queries!(study, id, q)
            Phylo.check_placement(files)
        catch e
            _fail_status!(dir, e)
            rethrow()
        end
        Phylo.run_workflow(Phylo.PLACEMENT_STEPS, dir, files, _study_config(study);
            overrides = get(doc, "overrides", nothing),
            run = Dict{String,Any}("study" => study, "placement" => id, "name" => name,
                                   "reference_tree" => ref))
        ref_doc = something(_read_doc(_ref_doc_path(String(ref))), Dict{String,Any}())
        d = _read_doc(_placement_doc_path(study, id))
        d["published"] = _publish_placement(study, name, string(get(ref_doc, "name", ref)), files)
        _write_doc(_placement_doc_path(study, id), d)
    end
    json(_job_to_namedtuple(job))
end
