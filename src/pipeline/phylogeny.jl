module Phylogeny

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Phylogenetic placement, as two workflows.
#
# A reference tree is built once and shared by any number of placements:
#   align       MAFFT aligns the reference sequences.
#   trim        trimAl removes columns (and optionally sequences) from it.
#   tree        IQ-TREE infers the tree with bootstrap support.
#
# A placement puts a set of queries on one reference tree:
#   align       MAFFT adds the queries to the trimmed reference alignment as
#               fragments (--addfragments).
#   trim        trimAl trims the combined alignment.
#   place       RAxML's EPA (-f v) places each query on the reference tree.
#   accumulate  gappa moves each query's placement mass onto the most basal
#               branch holding at least the threshold.
#
# Each workflow runs in a directory of its own. A step reruns when its inputs or
# its settings change; where it runs and how many threads it uses leave it
# current. Trimming and accumulation always run here, so a trimming setting can
# be tried without a round trip to the server.

export REFERENCE_STEPS, PLACEMENT_STEPS, Step, placement_settings, run_workflow,
       reference_files, placement_files, read_fasta, parse_fasta, safe_name,
       write_fasta, read_aligned, step_commands, trim_args, remote_step_target,
       read_status, trim_preview, step_qc, QC_STEPS, check_reference, check_placement

    using SHA, JSON3, Dates, Logging, OrderedCollections
    using ..PipelineLog
    using ..Validation
    using ..Tools: tool_bin, _sq, _run_logged, _safe_optional_args, _num
    using ..RemoteExec: run_remote, checked_target
    import ..Provenance

    """
    One step of a workflow. `inputs`, `outputs` and `extras` name files by the
    keys of the workflow's file map; `remote` is the `remote.stages` entry that
    sends it to the server, or nothing for a step that always runs here.
    """
    struct Step
        name    :: String
        section :: Tuple{String,String}
        remote  :: Union{String,Nothing}
        tools   :: Vector{String}
        inputs  :: Vector{String}
        outputs :: Vector{String}
        extras  :: Vector{String}
    end

    const REFERENCE_STEPS = [
        Step("align", ("reference", "align"), "phylogeny_align", ["mafft"],
             ["references.fasta"], ["reference.aln.fasta"], String[]),
        Step("trim", ("reference", "trim"), nothing, ["trimal"],
             ["reference.aln.fasta"], ["reference.trim.fasta", "reference.columns"], String[]),
        Step("tree", ("reference", "tree"), "phylogeny_tree", ["iqtree"],
             ["reference.trim.fasta"], ["reference.treefile"],
             ["reference.iqtree", "reference.log", "reference.contree"]),
    ]

    const PLACEMENT_STEPS = [
        Step("align", ("placement", "align"), "phylogeny_add", ["mafft"],
             ["queries.fasta", "reference.trim.fasta"], ["combined.aln.fasta"], String[]),
        Step("trim", ("placement", "trim"), nothing, ["trimal"],
             ["combined.aln.fasta"], ["combined.trim.fasta", "combined.columns"], String[]),
        Step("place", ("placement", "place"), "phylogeny_place", ["raxml"],
             ["combined.trim.fasta", "reference.treefile"], ["placement.jplace"],
             ["RAxML_info.epa", "RAxML_classification.epa",
              "RAxML_classificationLikelihoodWeights.epa", "RAxML_labelledTree.epa"]),
        Step("accumulate", ("placement", "accumulate"), nothing, ["gappa"],
             ["placement.jplace"], ["accumulated.jplace"], String[]),
    ]

    # Program names on the server when remote.tools does not name them.
    const REMOTE_TOOL_DEFAULTS = Dict("mafft" => "mafft", "iqtree" => "iqtree",
                                      "raxml" => "raxmlHPC-PTHREADS-SSE3")

    """
        reference_files(dir) -> Dict

    Where each file of a reference tree lives.
    """
    reference_files(dir::AbstractString) = Dict(
        "references.fasta"     => joinpath(dir, "references.fasta"),
        "reference.aln.fasta"  => joinpath(dir, "align", "reference.aln.fasta"),
        "reference.trim.fasta" => joinpath(dir, "trim", "reference.trim.fasta"),
        "reference.columns"    => joinpath(dir, "trim", "columns.txt"),
        "reference"            => joinpath(dir, "tree", "reference"),
        "reference.treefile"   => joinpath(dir, "tree", "reference.treefile"),
        "reference.iqtree"     => joinpath(dir, "tree", "reference.iqtree"),
        "reference.log"        => joinpath(dir, "tree", "reference.log"),
        "reference.contree"    => joinpath(dir, "tree", "reference.contree"),
    )

    """
        placement_files(dir, reference_dir) -> Dict

    Where each file of a placement lives; the reference alignment and tree are
    read from the reference tree's own directory.
    """
    function placement_files(dir::AbstractString, reference_dir::AbstractString)
        ref = reference_files(reference_dir)
        files = Dict(
            "queries.fasta"        => joinpath(dir, "queries.fasta"),
            "reference.trim.fasta" => ref["reference.trim.fasta"],
            "reference.treefile"   => ref["reference.treefile"],
            "combined.aln.fasta"   => joinpath(dir, "align", "combined.aln.fasta"),
            "combined.trim.fasta"  => joinpath(dir, "trim", "combined.trim.fasta"),
            "combined.columns"     => joinpath(dir, "trim", "columns.txt"),
            "placement.jplace"     => joinpath(dir, "place", "placement.jplace"),
            "accumulate"           => joinpath(dir, "accumulate"),
            "accumulated.jplace"   => joinpath(dir, "accumulate", "accumulated.jplace"),
        )
        for e in PLACEMENT_STEPS[3].extras
            files[e] = joinpath(dir, "place", e)
        end
        files
    end

    ## FASTA
    # RAxML and IQ-TREE refuse or silently rename taxa carrying these, and the
    # jplace tree would then name different tips from the alignment.
    const _BAD_NAME_CHARS = r"[\s:,();\[\]'\"]"
    const _SEQ_CHARS      = r"\A[ACGTURYKMSWBDHVN]*\z"

    safe_name(raw::AbstractString) = replace(strip(raw), _BAD_NAME_CHARS => "_")

    """
        parse_fasta(text; what) -> Vector{Pair{String,String}}

    Records from FASTA text, with names made safe for the tree tools and the
    sequences ungapped and upper-cased. Throws, naming the problem, on a record
    without sequence, a character outside the IUPAC nucleotide codes, or two
    records whose names coincide once made safe.
    """
    function parse_fasta(text::AbstractString; what::AbstractString="FASTA")
        records = Pair{String,String}[]
        name = nothing
        buf = IOBuffer()
        finish!() = if !isnothing(name)
            seq = uppercase(replace(String(take!(buf)), r"[\s\-.]" => ""))
            isempty(seq) && error("$what: '$name' has no sequence")
            occursin(_SEQ_CHARS, seq) ||
                error("$what: '$name' contains characters other than nucleotide codes")
            push!(records, name => seq)
        end
        for line in eachline(IOBuffer(text))
            if startswith(line, '>')
                finish!()
                name = safe_name(line[2:end])
                isempty(name) && error("$what: a record has no name")
            else
                isnothing(name) && !isempty(strip(line)) && error("$what: sequence before the first '>' header")
                print(buf, line)
            end
        end
        finish!()
        isempty(records) && error("$what: no sequences")
        seen = Dict{String,Int}()
        for (n, _) in records
            seen[n] = get(seen, n, 0) + 1
        end
        dups = sort([n for (n, c) in seen if c > 1])
        isempty(dups) || error("$what: duplicate names $(join(first(dups, 5), ", "))" *
                               (length(dups) > 5 ? " and $(length(dups) - 5) more" : ""))
        records
    end

    read_fasta(path::AbstractString; what::AbstractString=basename(path)) =
        parse_fasta(read(path, String); what)

    """
        read_aligned(path) -> Vector{Pair{String,String}}

    An alignment as written by MAFFT or trimAl, gaps kept, upper-cased.
    """
    function read_aligned(path::AbstractString)
        records = Pair{String,String}[]
        name = nothing
        buf = IOBuffer()
        for line in eachline(path)
            if startswith(line, '>')
                isnothing(name) || push!(records, name => uppercase(String(take!(buf))))
                name = String(strip(line[2:end]))
            else
                print(buf, strip(line))
            end
        end
        isnothing(name) || push!(records, name => uppercase(String(take!(buf))))
        records
    end

    function write_fasta(path::AbstractString, records)
        mkpath(dirname(path))
        tmp = path * ".tmp"
        open(tmp, "w") do io
            for (n, s) in records
                println(io, '>', n)
                println(io, s)
            end
        end
        mv(tmp, path; force=true)
        path
    end

    ## Settings
    _as_dict(x) = x isa AbstractDict ? Dict{String,Any}(string(k) => v for (k, v) in x) : Dict{String,Any}()

    function _merge(base::Dict{String,Any}, patch)
        out = copy(base)
        for (k, v) in _as_dict(patch)
            out[k] = (haskey(out, k) && out[k] isa AbstractDict && v isa AbstractDict) ?
                     _merge(_as_dict(out[k]), v) : v
        end
        out
    end

    """
        placement_settings(full_cfg, overrides=nothing) -> Dict

    The resolved `phylogeny` section with a reference tree's or placement's own
    overrides laid over it, validated as the config write gate would.
    """
    function placement_settings(full_cfg::AbstractDict, overrides=nothing)
        s = _merge(_as_dict(get(full_cfg, "phylogeny", nothing)), overrides)
        for w in ("reference", "placement")
            s[w] = _as_dict(get(s, w, nothing))
            for step in (w == "reference" ? REFERENCE_STEPS : PLACEMENT_STEPS)
                s[w][step.section[2]] = _as_dict(get(s[w], step.section[2], nothing))
            end
        end
        errors = Validation.ValidationError[]
        Validation._validate_phylogeny(errors, s, "phylogeny")
        isempty(errors) || error(join((e.message for e in errors), "; "))
        s
    end

    _section(s, step::Step) = s[step.section[1]][step.section[2]]

    _threads(v, fallback) = (v isa Integer && !(v isa Bool) && v >= 1) ? Int(v) : fallback

    ## Remote target
    """
        remote_step_target(full_cfg, stage; threads) -> NamedTuple or nothing

    Where a step with this `remote.stages` name runs: nothing for this machine,
    or the server when the stage is listed and a host is set. `remote.threads`
    applies there when it is a number; `threads` otherwise.
    """
    function remote_step_target(full_cfg::AbstractDict, stage::Union{String,Nothing}; threads::Integer)
        isnothing(stage) && return nothing
        rc = _as_dict(get(full_cfg, "remote", nothing))
        stage in string.(something(get(rc, "stages", nothing), String[])) || return nothing
        target = checked_target(rc, stage; threads = _threads(get(rc, "threads", nothing), threads))
        isnothing(target) && return nothing
        tools = merge(REMOTE_TOOL_DEFAULTS, Dict(string(k) => string(v)
                      for (k, v) in _as_dict(get(rc, "tools", nothing)) if !isnothing(v)))
        for (t, bin) in tools
            (Validation.is_shell_safe_arg(bin) && !startswith(bin, "-")) ||
                error("remote.tools.$t may not start with '-' or contain spaces or shell metacharacters (got: '$bin')")
        end
        (; target..., tools)
    end

    ## Commands
    _words(words...) = join(filter(!isempty, collect(words)), " ")

    function _mafft_mode(a::AbstractDict)
        strategy = string(get(a, "strategy", "auto"))
        strategy == "auto" && return "--auto"
        mi = _num(a, "maxiterate", 0)
        mi > 0 ? "--maxiterate $mi --$strategy" : "--$strategy"
    end

    """
        trim_args(t) -> String

    trimAl's selection flags for a trim section.
    """
    function trim_args(t::AbstractDict)
        method = string(get(t, "method", "manual"))
        parts = String[]
        if method == "manual"
            push!(parts, "-gt $(_num(t, "gap_threshold", 0.5))")
            c = _num(t, "conservation", nothing)
            isnothing(c) || push!(parts, "-cons $c")
            st = _num(t, "similarity_threshold", nothing)
            isnothing(st) || push!(parts, "-st $st")
        else
            push!(parts, "-$method")
        end
        ro = _num(t, "residue_overlap", nothing)
        so = _num(t, "sequence_overlap", nothing)
        (isnothing(ro) || isnothing(so)) || push!(parts, "-resoverlap $ro -seqoverlap $so")
        push!(parts, _safe_optional_args(t))
        _words(parts...)
    end

    """
        step_commands(step, s, bins, threads, seed; path, workdir) -> Vector{String}

    The shell commands of one step, in order. `path(name)` spells a file of the
    workflow's file map as the command line should see it; `workdir` is the
    absolute directory RAxML writes into, which it requires.
    """
    function step_commands(step::Step, s::AbstractDict, bins::AbstractDict,
                           threads::Integer, seed::Integer; path::Function, workdir::AbstractString="")
        c = _section(s, step)
        w = step.section[1]
        if w == "reference" && step.name == "align"
            return [_words(bins["mafft"], "--thread $threads", _mafft_mode(c), _safe_optional_args(c),
                           path("references.fasta"), ">", path("reference.aln.fasta"))]
        elseif step.name == "trim"
            src, dest, cols = w == "reference" ?
                ("reference.aln.fasta", "reference.trim.fasta", "reference.columns") :
                ("combined.aln.fasta", "combined.trim.fasta", "combined.columns")
            return [_words(bins["trimal"], "-in", path(src), "-out", path(dest), "-fasta",
                           trim_args(c), "-colnumbering", ">", path(cols))]
        elseif step.name == "tree"
            flag = string(get(c, "bootstrap", "standard")) == "ultrafast" ? "-bb" : "-b"
            return [_words(bins["iqtree"], "-s", path("reference.trim.fasta"),
                           "-m $(get(c, "model", "MFP"))", "$flag $(_num(c, "replicates", 100))",
                           "-nt $threads", "-seed $seed", "-pre", path("reference"), "-redo",
                           _safe_optional_args(c))]
        elseif step.name == "align"
            return [_words(bins["mafft"], "--thread $threads", _mafft_mode(c), _safe_optional_args(c),
                           "--addfragments", path("queries.fasta"), path("reference.trim.fasta"),
                           ">", path("combined.aln.fasta"))]
        elseif step.name == "place"
            h = _num(c, "heuristic", nothing)
            # The PTHREADS builds refuse fewer than two threads.
            return [_words(bins["raxml"], "-f v", "-m $(get(c, "model", "GTRCATI"))",
                           isnothing(h) ? "" : "-G $h", "-n epa",
                           "-s", path("combined.trim.fasta"), "-t", path("reference.treefile"),
                           "-T $(max(2, threads))", "-w", workdir, _safe_optional_args(c))]
        elseif step.name == "accumulate"
            return [_words(bins["gappa"], "edit accumulate",
                           "--jplace-path", path("placement.jplace"),
                           "--threshold $(_num(c, "threshold", 0.95))",
                           "--out-dir", path("accumulate"),
                           "--allow-file-overwriting", "--threads $threads")]
        end
        error("Unknown step: $(step.name)")
    end

    ## Status
    # status.json is what the page polls while a workflow runs.
    _status_path(dir) = joinpath(dir, "status.json")

    function read_status(dir::AbstractString)
        p = _status_path(dir)
        isfile(p) || return nothing
        try JSON3.read(read(p, String), Dict{String,Any}) catch; nothing end
    end

    const _status_lock = ReentrantLock()

    function _update_status!(dir, f::Function)
        lock(_status_lock) do
            mkpath(dir)
            st = something(read_status(dir), Dict{String,Any}())
            f(st)
            tmp = _status_path(dir) * ".tmp"
            write(tmp, JSON3.write(st))
            mv(tmp, _status_path(dir); force=true)
        end
    end

    _now() = string(now(UTC)) * "Z"

    function _set_step!(dir, step; kw...)
        _update_status!(dir, st -> begin
            steps = get!(st, "steps", Dict{String,Any}())
            entry = Dict{String,Any}(string(k) => v for (k, v) in get(steps, step, Dict{String,Any}()))
            for (k, v) in kw
                entry[string(k)] = v
            end
            steps[step] = entry
        end)
    end

    ## Step keys
    _canon(x::AbstractDict)   = "{" * join(["$k=$(_canon(x[k]))" for k in sort!(collect(keys(x)); by=string)], ",") * "}"
    _canon(x::AbstractVector) = "[" * join(_canon.(x), ",") * "]"
    _canon(x)                 = repr(x)

    function _step_key(files, step::Step, s, seed)
        io = IOBuffer()
        println(io, join(step.section, "."))
        println(io, _canon(_section(s, step)))
        step.name == "tree" && println(io, "seed=", seed)
        for name in step.inputs
            println(io, name, "=", bytes2hex(open(sha256, files[name])))
        end
        bytes2hex(sha256(take!(io)))
    end

    _step_dir(files, step::Step) = dirname(files[step.outputs[1]])
    _key_path(files, step::Step) = joinpath(_step_dir(files, step), "step.key")

    function _current(files, step::Step, key)
        kp = _key_path(files, step)
        isfile(kp) && strip(read(kp, String)) == key &&
            all(n -> isfile(files[n]) && filesize(files[n]) > 0, step.outputs)
    end

    ## Provenance
    function _local_env(step::Step, strict::Bool)
        tools = OrderedDict{String,Provenance.ToolRecord}()
        degraded = String[]
        for key in step.tools
            try
                tools[key] = Provenance.probe_tool(Provenance.TOOL_PROBES[key])
            catch err
                err isa Provenance.ProbeFailure || rethrow()
                strict && error("$key could not be run for step '$(step.name)': $(err.reason). " *
                                "Install it (install.sh) or set its path in config/tools.yml.")
                push!(degraded, key)
            end
        end
        Provenance.CapturedEnvironment(; tools, degraded_components=degraded)
    end

    # The version, path and checksum of a program on the server, read over the
    # stage's own connection.
    function _remote_record(ssh, target, key::AbstractString, bin::AbstractString)
        probe = Provenance.TOOL_PROBES[key]
        redirect = probe.stream === :stderr ? "2>&1 >/dev/null" : "2>/dev/null"
        script = "p=\$(command -v $bin) || exit 3; echo \"\$p\"; sha256sum \"\$p\" | cut -d' ' -f1; " *
                 "$bin $(join(probe.args, " ")) $redirect; true"
        out = try
            read(ssh(target.host, script), String)
        catch
            error("$key is not on the PATH of $(target.host) as '$bin'. Install it there or set remote.tools.$key.")
        end
        lines = split(out, '\n'; limit=3)
        length(lines) == 3 || error("Could not read the version of $bin on $(target.host)")
        Provenance.ToolRecord(probe.name, probe.parser(lines[3]), "$(target.host):$(strip(lines[1]))",
                              String(strip(lines[2])), Provenance.probed_by(probe))
    end

    ## Running
    function _clear_step!(files, step::Step)
        d = _step_dir(files, step)
        if isdir(d)
            if step.name == "tree"
                foreach(f -> startswith(f, "reference.") && rm(joinpath(d, f); force=true), readdir(d))
            elseif step.name == "place"
                foreach(f -> startswith(f, "RAxML_") && rm(joinpath(d, f); force=true), readdir(d))
            end
        end
        foreach(n -> rm(files[n]; force=true), step.outputs)
        rm(_key_path(files, step); force=true)
    end

    function _run_local_step(files, step::Step, s, seed, threads, log_path)
        bins = Dict(t => tool_bin(t) for t in step.tools)
        dir = _step_dir(files, step)
        mkpath(dir)
        cmds = step_commands(step, s, bins, threads, seed;
                             path = n -> _sq(files[n]), workdir = _sq(dir * "/"))
        for c in cmds
            _run_logged(c, log_path)
        end
        step.name == "place" &&
            mv(joinpath(dir, "RAxML_portableTree.epa.jplace"), files["placement.jplace"]; force=true)
        cmds
    end

    # Files go up and come back under their file-map names, flat in the staging
    # directory, so the commands can use the bare names.
    function _run_remote_step(files, step::Step, s, seed, target, log_path, emit, records)
        cmds = Ref(String[])
        inputs  = [files[n] => n for n in step.inputs]
        outputs = [(n == "placement.jplace" ? "RAxML_portableTree.epa.jplace" : n) => files[n]
                   for n in step.outputs]
        extras  = [n => files[n] for n in step.extras]
        before_run = function (ssh, _)
            for t in step.tools
                records[t] = _remote_record(ssh, target, t, target.tools[t])
            end
        end
        command_for = function (staging_dir)
            cmds[] = step_commands(step, s, target.tools, target.threads, seed;
                                   path = identity, workdir = "$staging_dir/")
            "cd $staging_dir && " * join(cmds[], " && ")
        end
        run_remote(emit, target, command_for; inputs, outputs, optional_outputs=extras,
                   before_run, log_path, stage=something(step.remote, step.name))
        cmds[]
    end

    # The tools read these files directly, so they must hold the safe names.
    function _normalise!(path, what)
        recs = read_fasta(path; what)
        text = join(">$n\n$q\n" for (n, q) in recs)
        read(path, String) == text || write_fasta(path, recs)
        recs
    end

    """
        check_reference(files)

    Throws unless the references are fit to build a tree from.
    """
    function check_reference(files)
        isfile(files["references.fasta"]) || error("Upload the reference sequences first")
        refs = _normalise!(files["references.fasta"], "References")
        length(refs) >= 4 || error("References: a tree needs at least 4 sequences (got $(length(refs)))")
        nothing
    end

    """
        check_placement(files)

    Throws unless the reference tree is built and the queries can go on it.
    """
    function check_placement(files)
        (isfile(files["reference.trim.fasta"]) && isfile(files["reference.treefile"])) ||
            error("The reference tree has not been built yet")
        isfile(files["queries.fasta"]) || error("There are no query sequences yet")
        qs = _normalise!(files["queries.fasta"], "Queries")
        refs = Set(first.(read_aligned(files["reference.trim.fasta"])))
        clash = sort([n for (n, _) in qs if n in refs])
        isempty(clash) || error("Queries share names with references: $(join(first(clash, 5), ", "))")
        nothing
    end

    """
        run_workflow(steps, dir, files, full_cfg; overrides, seed, run, emit)

    Run the steps that are not current, in order, writing the align and trim
    QC to `dir/qc/<step>.json` and each step's record to `dir/attestation.yml` whether it
    succeeds or not.
    """
    function run_workflow(steps::Vector{Step}, dir::AbstractString, files::AbstractDict,
                          full_cfg::AbstractDict;
                          overrides=nothing,
                          seed::Integer = Int(get(full_cfg, "seed", 123)),
                          run::AbstractDict = Dict{String,Any}(),
                          emit::Function = msg -> @info(msg))
        s = placement_settings(full_cfg, overrides)
        threads = _threads(get(s, "threads", nothing), 4)
        strict  = Provenance.strict_mode(full_cfg)
        section = Dict{String,Any}(steps[1].section[1] => s[steps[1].section[1]], "threads" => threads)
        att = Provenance.Attestation(;
            run = run,
            config = Dict{String,Any}("phylogeny" => section, "seed" => seed),
            config_sha256 = bytes2hex(sha256(_canon(Dict("phylogeny" => section, "seed" => seed)))))
        log_dir = joinpath(dir, "logs")
        mkpath(log_dir)
        _update_status!(dir, st -> begin
            st["state"] = "running"
            st["started"] = _now()
            delete!(st, "error")
            delete!(st, "finished")
            st["steps"] = Dict{String,Any}(step.name => Dict{String,Any}("state" => "pending") for step in steps)
        end)

        try
            for step in steps
                key = _step_key(files, step, s, seed)
                if _current(files, step, key)
                    emit("$(step.name): up to date")
                    isfile(joinpath(dir, "qc", "$(step.name).json")) || _write_qc(steps, dir, files, step)
                    _set_step!(dir, step.name; state="current")
                    continue
                end
                target = remote_step_target(full_cfg, step.remote; threads)
                where = isnothing(target) ? "local" : target.host
                emit("$(step.name): running ($where)")
                _set_step!(dir, step.name; state="running", where, started=_now())
                started  = Provenance.timestamp()
                log_path = joinpath(log_dir, "$(step.name).log")
                reset_tool_logs(log_path)
                _clear_step!(files, step)
                rm(joinpath(dir, "qc", "$(step.name).json"); force=true)
                env = Provenance.CapturedEnvironment()
                try
                    cmds = if isnothing(target)
                        env = _local_env(step, strict)
                        _run_local_step(files, step, s, seed, threads, log_path)
                    else
                        records = OrderedDict{String,Provenance.ToolRecord}()
                        c = _run_remote_step(files, step, s, seed, target, log_path, emit, records)
                        env = Provenance.CapturedEnvironment(; tools=records)
                        c
                    end
                    missing_out = filter(n -> !isfile(files[n]), step.outputs)
                    isempty(missing_out) || error("$(step.name) finished without writing $(join(missing_out, ", "))")
                    write(_key_path(files, step), key)
                    Provenance.record_stage!(att, step.name, env; started, commands=cmds,
                        outputs=[Provenance.output_record(files[n], dir) for n in step.outputs
                                 if startswith(files[n], dir)])
                    _write_qc(steps, dir, files, step)
                    _set_step!(dir, step.name; state="done", finished=_now())
                catch e
                    Provenance.record_stage!(att, step.name, env; started, status="failed")
                    tail = isfile(log_path) ? join(last(readlines(log_path), 15), "\n") : ""
                    msg = sprint(showerror, e)
                    _set_step!(dir, step.name; state="failed", finished=_now(),
                               error=isempty(tail) ? msg : "$msg\n$tail")
                    rethrow()
                end
            end
            _update_status!(dir, st -> (st["state"] = "done"; st["finished"] = _now()))
        catch e
            _update_status!(dir, st -> begin
                st["state"] = "failed"
                st["finished"] = _now()
                st["error"] = sprint(showerror, e)
            end)
            rethrow()
        finally
            Provenance.write_attestation(att, joinpath(dir, "attestation.yml"))
        end
        nothing
    end

    include("phylogeny/qc.jl")
end
