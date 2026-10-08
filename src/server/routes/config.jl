# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/studies/{study}/config
#         /api/v1/studies/{study}/runs/{run}/config
using JSON3, YAML
using ..PrimersLibrary

const _config_file_lock = ReentrantLock()

## Config resolution
function _load_yml(path::String)
    isfile(path) ? something(YAML.load_file(path), Dict()) : Dict()
end

# Write a YAML document over `path` atomically. The temp file and rename mean a
# reader never observes a half-written file; the lock means two writes cannot
# interleave their bytes. It does NOT make a read-modify-write sequence atomic:
# two concurrent saves can still each validate against their own snapshot and
# the later one silently revert the earlier.
function _atomic_write_yaml(path::String, doc)
    mkpath(dirname(path))
    lock(_config_file_lock) do
        tmp = path * ".tmp"
        open(tmp, "w") do io
            YAML.write(io, Config.canonical_yaml_doc(doc))
        end
        mv(tmp, path; force=true)
    end
    nothing
end

## Recursively convert JSON3 values into plain Dict/Vector/scalars for YAML.write.
_to_plain(x::JSON3.Object) = Dict{String,Any}(string(k) => _to_plain(v) for (k, v) in x)
_to_plain(x::JSON3.Array)  = Any[_to_plain(v) for v in x]
_to_plain(x) = x

function _flatten(d::Dict, prefix::String="")
    out = Dict{String,Any}()
    for (k, v) in d
        key = prefix == "" ? string(k) : "$(prefix).$(k)"
        if v isa Dict
            merge!(out, _flatten(v, key))
        else
            out[key] = v
        end
    end
    out
end

function _default_cfg(root::String)
    factory = _flatten(_load_yml(joinpath(root, "config", "defaults", "pipeline.yml")))
    global_ = _flatten(_load_yml(joinpath(root, "config", "pipeline.yml")))
    merge(factory, global_)  # global overrides factory; factory provides all keys
end

# Returns the merged config cascade for a run as a flat dict of dotted keys,
# each annotated with its effective value and source level.
#
# Cascade order (lowest wins): default -> study -> group -> run
#
# "source" is the deepest level that explicitly sets the key.
# Keys not set at any user level have source "default".
function _resolve_config(study::String, run::Union{String,Nothing}=nothing,
                         group::Union{String,Nothing}=nothing)
    root        = dirname(ServerState.data_dir())
    default_cfg = _default_cfg(root)
    study_cfg   = _flatten(_load_yml(joinpath(ServerState.data_dir(), study, "pipeline.yml")))
    group_cfg   = isnothing(group) ? Dict{String,Any}() :
                  _flatten(_load_yml(joinpath(ServerState.data_dir(), study, group, "pipeline.yml")))
    run_cfg     = isnothing(run) ? Dict{String,Any}() :
                  _flatten(_load_yml(joinpath(ServerState.data_dir(), study,
                                              isnothing(group) ? run : joinpath(group, run),
                                              "pipeline.yml")))

    all_keys = union(keys(default_cfg), keys(study_cfg), keys(group_cfg), keys(run_cfg))

    Dict(k => begin
        if haskey(run_cfg, k)
            (; value=run_cfg[k],   source="run")
        elseif haskey(group_cfg, k)
            (; value=group_cfg[k], source="group")
        elseif haskey(study_cfg, k)
            (; value=study_cfg[k], source="study")
        elseif haskey(default_cfg, k)
            (; value=default_cfg[k], source="default")
        end
    end for k in all_keys)
end

## Allowed config keys - derived from defaults/pipeline.yml at load time.
# Only keys present in the factory defaults can be set via the API.
const _ALLOWED_CONFIG_KEYS = let
    factory_path = joinpath(@__DIR__, "..", "..", "..", "config", "defaults", "pipeline.yml")
    isfile(factory_path) ? Set(keys(_flatten(_load_yml(factory_path)))) : Set{String}()
end

const _FACTORY_VALUES = let
    factory_path = joinpath(@__DIR__, "..", "..", "..", "config", "defaults", "pipeline.yml")
    isfile(factory_path) ? _flatten(_load_yml(factory_path)) : Dict{String,Any}()
end

# Values reaching a command line, as the pipeline itself checks them.
const _ENUM_VALUES = Dict(
    "vsearch.strand"                 => (nothing, "plus", "both"),
    "tagging.source"                 => ("VSEARCH", "DADA2"),
    "analysis.alpha.normalisation"   => ("none", "rarefy", "srs"),
    "analysis.beta.normalisation"    => ("hellinger", "none", "rarefy", "srs"),
    "analysis.differential.method"   => Differential.METHODS,
    "analysis.differential.offset"   => Differential.OFFSET_METHODS,
    "phylogeny.reference.align.strategy" => Validation.PHYLOGENY_ALIGN_STRATEGIES,
    "phylogeny.placement.align.strategy" => Validation.PHYLOGENY_ALIGN_STRATEGIES,
    "phylogeny.reference.trim.method"    => Validation.PHYLOGENY_TRIM_METHODS,
    "phylogeny.placement.trim.method"    => Validation.PHYLOGENY_TRIM_METHODS,
    "phylogeny.reference.tree.bootstrap" => Validation.PHYLOGENY_BOOTSTRAPS,
)

# String keys that become part of a file name: output stems and read suffixes.
const _FILE_NAME_PARTS = ("dada2.output.fasta_prefix", "dada2.output.seq_table_prefix",
                          "dada2.output.taxa_prefix", "cutadapt.r1_suffix", "cutadapt.r2_suffix")

# Numeric keys whose factory default is a number but which also take ~.
const _NULLABLE_NUMBERS = ("phylogeny.placement.place.heuristic",)

"""
    _config_value_error(key, value) -> Union{String,Nothing}

Why a value cannot be saved under this key, or nothing. The type follows the
factory default; a key whose default is ~ accepts a string, number or boolean.
"""
function _config_value_error(key::String, value)
    v = value isa JSON3.Array ? collect(value) : value
    if haskey(_ENUM_VALUES, key)
        allowed = _ENUM_VALUES[key]
        v in allowed || return "$key must be one of $(join(filter(!isnothing, collect(allowed)), ", "))"
        return nothing
    end
    if endswith(key, "optional_args")
        (v isa AbstractString && (isempty(strip(v)) || occursin(Tools._SAFE_ARGS_RE, strip(v)))) ||
            return "$key may contain only letters, digits, spaces and - _ . = , / : +"
        return nothing
    end
    if key in _NULLABLE_NUMBERS
        (isnothing(v) || (v isa Real && !(v isa Bool))) || return "$key must be a number or empty"
        return nothing
    end
    if key in ("remote.host", "remote.rscript")
        isnothing(v) && return nothing
        (v isa AbstractString && Validation.is_shell_safe_arg(v) && !startswith(v, "-")) ||
            return "$key may not start with '-' or contain spaces or shell metacharacters"
        return nothing
    end
    # The rest of the remote block is refused here exactly as validation refuses it.
    if startswith(key, "remote.")
        bad = Validation.remote_value_error(key[length("remote.")+1:end], v)
        isnothing(bad) || return bad
    end
    # Output stems and read-file suffixes become parts of file names (and the
    # suffixes are spliced into a \Q...\E regex), so they are names, not paths.
    if key in _FILE_NAME_PARTS
        (v isa AbstractString && Validation.is_safe_name(v)) ||
            return "$key may contain only letters, digits and . _ -"
        return nothing
    end
    default = get(_FACTORY_VALUES, key, nothing)
    if default isa Bool
        v isa Bool || return "$key must be true or false"
    elseif default isa Real
        (v isa Real && !(v isa Bool)) || return "$key must be a number"
    elseif default isa AbstractString
        v isa AbstractString || return "$key must be text"
    elseif default isa AbstractVector
        v isa AbstractVector || return "$key must be a list"
    elseif default isa AbstractDict
        v isa AbstractDict || return "$key must be a section"
    end
    nothing
end

function _validate_config_key(dotted_key::String)
    dotted_key in _ALLOWED_CONFIG_KEYS || error(
        "Config key '$dotted_key' is not a recognised pipeline option")
end

function _write_override(path::String, dotted_key::String, value)
    _validate_config_key(dotted_key)
    lock(_config_file_lock) do
        cfg = _load_yml(path)
        parts = split(dotted_key, ".")
        d = cfg
        for p in parts[1:end-1]
            child = get(d, p, nothing)
            if child isa Dict
                d = child
            elseif isnothing(child)
                new_dict = Dict()
                d[p] = new_dict
                d = new_dict
            else
                error("Config key segment '$p' is not a section; cannot set sub-key")
            end
        end
        d[parts[end]] = value
        _atomic_write_yaml(path, cfg)
    end
end

function _delete_override(path::String, dotted_key::String)
    isfile(path) || return
    lock(_config_file_lock) do
        cfg   = _load_yml(path)
        parts = split(dotted_key, ".")
        d     = cfg
        for p in parts[1:end-1]
            haskey(d, p) || return
            d = d[p]
            d isa Dict || return
        end
        delete!(d, parts[end])
        _atomic_write_yaml(path, cfg)
    end
end

## Routes
# Default (global) config endpoints
@get "/api/v1/config" function(req)
    root = dirname(ServerState.data_dir())
    json(Dict(k => (; value=v, source="default") for (k, v) in _default_cfg(root)))
end

@patch "/api/v1/config" function(req)
    root      = dirname(ServerState.data_dir())
    user_path = joinpath(root, "config", "pipeline.yml")
    body = JSON3.read(req.body)
    for (k, v) in body
        string(k) in _ALLOWED_CONFIG_KEYS || return json_error(400, "invalid_config_key",
            "Config key '$(string(k))' is not a recognised pipeline option")
        bad = _config_value_error(string(k), v)
        isnothing(bad) || return json_error(400, "invalid_config_value", bad)
    end
    for (k, v) in body
        _write_override(user_path, string(k), v)
    end
    json(Dict(k => (; value=v, source="default") for (k, v) in _default_cfg(root)))
end

@delete "/api/v1/config/{key}" function(req, key::String)
    root      = dirname(ServerState.data_dir())
    user_path = joinpath(root, "config", "pipeline.yml")
    _delete_override(user_path, key)
    json(Dict(k => (; value=v, source="default") for (k, v) in _default_cfg(root)))
end

# List available primer pair names from config/primers.yml
@get "/api/v1/primers" function(req)
    primers_path = joinpath(dirname(ServerState.data_dir()), "config", "primers.yml")
    isfile(primers_path) || return json([])
    cfg = YAML.load_file(primers_path)
    cfg isa Dict || return json([])
    pairs = get(cfg, "Pairs", [])
    pairs isa Vector || return json([])
    # Each pair is a single-key Dict like Dict("EMP" => ["EMP515F", "EMP806R"])
    names = String[]
    for p in pairs
        p isa Dict || continue
        for k in keys(p)
            push!(names, string(k))
        end
    end
    json(names)
end

## Primers document (whole-file view and edit)
_primers_path() = joinpath(dirname(ServerState.data_dir()), "config", "primers.yml")

function _write_primers(native::AbstractDict)
    _atomic_write_yaml(_primers_path(), native)
end

# Validate then write; every primers write passes through here, so a dangling pair
# is caught before it reaches a run. Returns the saved canonical document, or an
# HTTP.Response error.
function _save_primers(doc::AbstractDict)
    native = PrimersLibrary.to_yaml_doc(doc)
    errors = PrimersLibrary.validate(doc; native)
    isempty(errors) || return json_error(400, "invalid_primers", join(errors, "; "))
    _write_primers(native)
    @info "Saved primers document ($(_primers_path()))"
    doc
end

# Read the document; an unparseable file becomes a 400 naming it. PrimersLibrary.load
# raises on a bad file because an empty document would pass validation and the next
# Save would overwrite the real primers. Returns the document, or an HTTP.Response error.
function _read_primers()
    try
        PrimersLibrary.load(_primers_path())
    catch e
        @error "primers.yml is not readable" path=_primers_path() exception=(e, catch_backtrace())
        json_error(400, "primers_unreadable",
            "config/primers.yml cannot be read as a primers document and must be repaired by hand: $(sprint(showerror, e))")
    end
end

## Cascade-resolved config references
# Read a pipeline.yml for the reference scan. A file anywhere in the tree may be
# malformed, or may not be a mapping at all. Neither is a save's business: the
# scan gathers advisory warnings only, so an unreadable file is skipped with a
# warning. Letting it throw would make one stray file block every save with a 500.
# Results are cached per call: an ancestor's pipeline.yml is on the cascade of
# every project beneath it and would otherwise be re-read once per descendant.
function _scan_yml(path::String, cache::Dict{String,Any})
    get!(cache, path) do
        cfg = try
            _load_yml(path)
        catch e
            @warn "Skipping an unreadable pipeline.yml while scanning for config references" path exception=e
            return Dict()
        end
        cfg isa AbstractDict ? cfg : Dict()
    end
end

# Sentinel for "this level does not set the key", distinct from a level that
# sets it to null. Only an absent key inherits: an explicit null overrides to
# null, which is what _resolve_config and Config.load_merged_config both do.
struct _Absent end
const _ABSENT = _Absent()

# The value at a dotted path, or _ABSENT when any segment is absent or is not a
# mapping. Distinguishing absent from present-and-null matters: only an absent
# key inherits, so _ABSENT (not nothing) is the sentinel for "this level does
# not set it"; a level that explicitly nulls the key returns nothing, which is
# a real value that overrides whatever an ancestor set.
function _getin(cfg, parts::Vector{String})
    cur = cfg
    for p in parts
        cur isa AbstractDict || return _ABSENT
        haskey(cur, p) || return _ABSENT
        cur = cur[p]
    end
    cur
end

# The names a resolved value contributes. A list key (cutadapt.primer_pairs)
# contributes each element; a scalar key (dada2.taxonomy.database) contributes
# itself. A level that is absent, or a value that is neither a list nor a
# scalar (a mapping, say), contributes no names. One scan therefore serves both.
_ref_names(v::AbstractVector) = String[string(x) for x in v]
_ref_names(::Nothing)         = String[]
_ref_names(::_Absent)         = String[]
_ref_names(::AbstractDict)    = String[]
_ref_names(v)                 = String[string(v)]

# Which projects effectively see each value of `dotted_key`, resolved through the
# real cascade: factory defaults, then the global config, then each ancestor
# pipeline.yml down to the project itself. Deepest wins. Returns
# value -> sorted project labels.
#
# Scanning raw files and reporting which ones NAME a value is wrong wherever the
# key has a factory default, because a project inherits the default without
# naming it anywhere. config/defaults/pipeline.yml sets both
# cutadapt.primer_pairs ("EMP") and dada2.taxonomy.database ("pr2"): on the real
# data tree three studies inherit the former and every project inherits the
# latter, so a name-scan reports nothing while the edit breaks all of them.
#
# The set of projects scanned is the union of every directory that owns a
# pipeline.yml and every run directory the pipeline route itself would submit
# (study/run and study/group/run, per _run_names/_group_names/_group_run_names).
# A run with no pipeline.yml of its own still inherits the factory default and
# would break identically, so it must be scanned too; a bare study or group
# directory that is not itself a run target and sets nothing of its own is not
# something that runs, so it is not forced in.
#
# Labels are project paths only. Neither the global config nor the factory
# defaults is a project, so neither is ever a label; they participate only as the
# base of each project's cascade.
function _cascade_refs(dotted_key::String)
    parts     = String[String(p) for p in split(dotted_key, ".")]
    root      = dirname(ServerState.data_dir())
    data_root = ServerState.data_dir()
    cache     = Dict{String,Any}()
    base = _getin(_scan_yml(joinpath(root, "config", "defaults", "pipeline.yml"), cache), parts)
    glob = _getin(_scan_yml(joinpath(root, "config", "pipeline.yml"), cache), parts)
    base_eff = glob isa _Absent ? base : glob
    refs = Dict{String,Vector{String}}()
    isdir(data_root) || return refs

    targets = Set{String}()
    # Skip unreadable directories; these warnings are advisory.
    for (dir, _, files) in walkdir(data_root; onerror = _ -> nothing)
        "pipeline.yml" in files || continue
        rel = relpath(dir, data_root)
        rel == "." || push!(targets, rel)
    end
    # _run_names reads pipeline.yml unguarded (via _is_pooled), so a malformed file
    # can throw; each study is guarded separately so one bad study is skipped with a warning.
    for study in _study_names()
        try
            for r in _run_names(study)
                push!(targets, joinpath(study, r))
            end
            for g in _group_names(study)
                for r in _group_run_names(study, g)
                    push!(targets, joinpath(study, g, r))
                end
            end
        catch e
            @warn "Skipping run/group enumeration for a study while scanning for config references" study exception=e
        end
    end

    for rel in targets
        segs = splitpath(rel)
        eff  = base_eff
        for i in 1:length(segs)
            v = _getin(_scan_yml(joinpath(data_root, segs[1:i]..., "pipeline.yml"), cache), parts)
            v isa _Absent || (eff = v)
        end
        for name in _ref_names(eff)
            push!(get!(refs, name, String[]), rel)
        end
    end
    for k in keys(refs)
        refs[k] = sort(unique(refs[k]))
    end
    refs
end

# Pairs present in `current` (the document as it stands on disk) but absent from
# the submitted `doc`, that some project's config cascade still resolves to. A
# rename reads as the old name disappearing, so this covers both deletion and
# rename. The write is not blocked; these are warnings only.
# `current` comes from the handler's guarded `_read_primers`, so a file that turns
# unparseable between reads cannot throw here.
function _pair_pipeline_refs_removed(doc::AbstractDict, current::AbstractDict)
    was = Set(PrimersLibrary.pair_names(current))
    now = Set(PrimersLibrary.pair_names(doc))
    removed = setdiff(was, now)
    isempty(removed) && return Any[]
    refs = _cascade_refs("cutadapt.primer_pairs")
    warnings = Any[]
    for name in sort(collect(removed))
        haskey(refs, name) || continue
        push!(warnings, Dict("pair" => name, "referenced_by" => refs[name]))
    end
    warnings
end

# Parse a request body into a primers document; a malformed body is a 400.
# Returns the document, or an HTTP.Response error.
function _primers_body(body::AbstractString)
    parsed = try
        JSON3.read(body)
    catch e
        return json_error(400, "invalid_json",
            "Request body is not valid JSON: $(sprint(showerror, e))")
    end
    parsed isa JSON3.Object || return json_error(400, "invalid_json",
        "Request body must be a JSON object with Forward, Reverse and Pairs")
    _to_plain(parsed)
end

# The whole primers document (Forward, Reverse, Pairs) in canonical shape.
@get "/api/v1/primers/document" function(req)
    doc = _read_primers()
    doc isa HTTP.Response && return doc
    json(doc)
end

# Replace the whole primers document. A pair naming an absent primer is a 400 with
# the file unchanged. The response warns about any dropped pair that a project's
# config cascade still resolves to. A corrupt file on disk refuses the write.
@put "/api/v1/primers" function(req)
    current = _read_primers()
    current isa HTTP.Response && return current
    doc = _primers_body(String(req.body))
    doc isa HTTP.Response && return doc
    warnings = _pair_pipeline_refs_removed(doc, current)
    result = _save_primers(doc)
    result isa HTTP.Response && return result
    json(Dict("document" => result, "warnings" => warnings))
end

@get "/api/v1/studies/{study}/config" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    json(_resolve_config(study))
end

@patch "/api/v1/studies/{study}/config" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    body = JSON3.read(req.body)
    for (k, v) in body
        string(k) in _ALLOWED_CONFIG_KEYS || return json_error(400, "invalid_config_key",
            "Config key '$(string(k))' is not a recognised pipeline option")
        bad = _config_value_error(string(k), v)
        isnothing(bad) || return json_error(400, "invalid_config_value", bad)
    end
    path = joinpath(ServerState.data_dir(), study, "pipeline.yml")
    for (k, v) in body
        _write_override(path, string(k), v)
    end
    json(_resolve_config(study))
end

@delete "/api/v1/studies/{study}/config/{key}" function(req, study::String, key::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    path = joinpath(ServerState.data_dir(), study, "pipeline.yml")
    _delete_override(path, key)
    json(_resolve_config(study))
end

## Override indicators - which downstream entities override a given level's keys
# Returns { "dotted.key": ["group1", "run1", ...] } for each study-level key
# that has a group or run override.
@get "/api/v1/studies/{study}/config/overrides" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    overrides = Dict{String,Vector{String}}()
    # Check group-level pipeline.yml files
    for group in _group_names(study)
        group_cfg = _flatten(_load_yml(joinpath(ServerState.data_dir(), study, group, "pipeline.yml")))
        for k in keys(group_cfg)
            push!(get!(overrides, k, String[]), group)
        end
        # Check runs within this group
        for run in _group_run_names(study, group)
            run_cfg = _flatten(_load_yml(joinpath(ServerState.data_dir(), study, group, run, "pipeline.yml")))
            for k in keys(run_cfg)
                push!(get!(overrides, k, String[]), "$group/$run")
            end
        end
    end
    # Check direct (ungrouped) runs
    for run in _run_names(study)
        run_cfg = _flatten(_load_yml(joinpath(ServerState.data_dir(), study, run, "pipeline.yml")))
        for k in keys(run_cfg)
            push!(get!(overrides, k, String[]), run)
        end
    end
    json(overrides)
end

# Returns { "dotted.key": ["run1", ...] } for each group-level key
# that has a run override.
@get "/api/v1/studies/{study}/groups/{group}/config/overrides" function(req, study::String, group::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    group in _group_names(study) || return json_error(404, "group_not_found",
                                                          "Group '$group' not found")
    overrides = Dict{String,Vector{String}}()
    for run in _group_run_names(study, group)
        run_cfg = _flatten(_load_yml(joinpath(ServerState.data_dir(), study, group, run, "pipeline.yml")))
        for k in keys(run_cfg)
            push!(get!(overrides, k, String[]), run)
        end
    end
    json(overrides)
end

## Group config endpoints
@get "/api/v1/studies/{study}/groups/{group}/config" function(req, study::String, group::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    group in _group_names(study) || return json_error(404, "group_not_found",
                                                          "Group '$group' not found")
    json(_resolve_config(study, nothing, group))
end

@patch "/api/v1/studies/{study}/groups/{group}/config" function(req, study::String, group::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    group in _group_names(study) || return json_error(404, "group_not_found",
                                                          "Group '$group' not found")
    body = JSON3.read(req.body)
    for (k, v) in body
        string(k) in _ALLOWED_CONFIG_KEYS || return json_error(400, "invalid_config_key",
            "Config key '$(string(k))' is not a recognised pipeline option")
        bad = _config_value_error(string(k), v)
        isnothing(bad) || return json_error(400, "invalid_config_value", bad)
    end
    path = joinpath(ServerState.data_dir(), study, group, "pipeline.yml")
    for (k, v) in body
        _write_override(path, string(k), v)
    end
    json(_resolve_config(study, nothing, group))
end

@delete "/api/v1/studies/{study}/groups/{group}/config/{key}" function(req, study::String,
                                                                            group::String,
                                                                            key::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    group in _group_names(study) || return json_error(404, "group_not_found",
                                                          "Group '$group' not found")
    path = joinpath(ServerState.data_dir(), study, group, "pipeline.yml")
    _delete_override(path, key)
    json(_resolve_config(study, nothing, group))
end

## Run config endpoints
# Helper: find all valid run names including those inside groups
function _all_run_names(study::String)
    runs = copy(_run_names(study))
    for group in _group_names(study)
        append!(runs, _group_run_names(study, group))
    end
    unique!(runs)
end

@get "/api/v1/studies/{study}/runs/{run}/config" function(req, study::String, run::String)
    err = _validate_run_request(study, run)
    isnothing(err) || return err
    group = let g = _req_group(req); isnothing(g) ? _run_group(study, run) : g end
    json(_resolve_config(study, run, group))
end

@patch "/api/v1/studies/{study}/runs/{run}/config" function(req, study::String, run::String)
    err = _validate_run_request(study, run)
    isnothing(err) || return err
    group = let g = _req_group(req); isnothing(g) ? _run_group(study, run) : g end
    run_path = isnothing(group) ? run : joinpath(group, run)
    body = JSON3.read(req.body)
    for (k, v) in body
        string(k) in _ALLOWED_CONFIG_KEYS || return json_error(400, "invalid_config_key",
            "Config key '$(string(k))' is not a recognised pipeline option")
        bad = _config_value_error(string(k), v)
        isnothing(bad) || return json_error(400, "invalid_config_value", bad)
    end
    path = joinpath(ServerState.data_dir(), study, run_path, "pipeline.yml")
    for (k, v) in body
        _write_override(path, string(k), v)
    end
    json(_resolve_config(study, run, group))
end

@delete "/api/v1/studies/{study}/runs/{run}/config/{key}" function(req, study::String,
                                                                        run::String,
                                                                        key::String)
    err = _validate_run_request(study, run)
    isnothing(err) || return err
    group = let g = _req_group(req); isnothing(g) ? _run_group(study, run) : g end
    run_path = isnothing(group) ? run : joinpath(group, run)
    path = joinpath(ServerState.data_dir(), study, run_path, "pipeline.yml")
    _delete_override(path, key)
    json(_resolve_config(study, run, group))
end
