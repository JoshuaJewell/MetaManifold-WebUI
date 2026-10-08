module Validation

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export validate_environment, validate_project, ValidationError,
       DENOVO_METHODS, REMOTE_STAGES, DADA2_REMOTE_STAGES, PHYLOGENY_REMOTE_STAGES,
       PHYLOGENY_ALIGN_STRATEGIES, PHYLOGENY_BOOTSTRAPS, PHYLOGENY_TRIM_METHODS, SAFE_NAME_RE, is_safe_name, is_shell_safe,
       is_shell_safe_arg, remote_value_error, primer_document_errors,
       database_document_errors

    using YAML, Logging
    using ..PipelineTypes
    using ..Config

    ## Safe names
    # The guard for every user-supplied name that reaches a filesystem path or a
    # quoted SQL identifier: study, run, preset, filter, and category-set names.
    # `\z` because PCRE's `$` also matches before a trailing newline.
    const SAFE_NAME_RE = r"\A[A-Za-z0-9._-]+\z"

    is_safe_name(s::AbstractString) = occursin(SAFE_NAME_RE, s)

    # Allowed chimera detection methods accepted by DADA2's removeBimeraDenovo.
    # Kept here so configuration validation and the call-site guard in
    # pipeline/dada2/chimera.jl share a single source of truth.
    const DENOVO_METHODS = ("consensus", "pooled", "per-sample")

    # The stages that can be offloaded to a compute server. Defined here so the
    # config write gate can refuse any other name, which the pipeline would
    # silently run locally.
    const DADA2_REMOTE_STAGES     = ("learn_errors", "denoise", "chimera_removal", "assign_taxonomy")
    const PHYLOGENY_REMOTE_STAGES = ("phylogeny_align", "phylogeny_tree", "phylogeny_add", "phylogeny_place")
    const REMOTE_STAGES           = (DADA2_REMOTE_STAGES..., PHYLOGENY_REMOTE_STAGES...)

    # MAFFT strategies an align section may name; each but auto is the pairwise
    # flag of its method (--localpair with --maxiterate is L-INS-i).
    const PHYLOGENY_ALIGN_STRATEGIES = ("auto", "localpair", "genafpair", "globalpair", "6merpair")
    # trimAl's selection modes: manual reads the thresholds, the rest are its own flags.
    const PHYLOGENY_TRIM_METHODS = ("manual", "gappyout", "strict", "strictplus", "automated1",
                                    "nogaps", "noallgaps")
    const PHYLOGENY_BOOTSTRAPS       = ("standard", "ultrafast")

    # Characters that give a remote login shell something to do besides name a
    # file. These paths are interpolated into a command string that a remote
    # sshd runs through a login shell, and a path is never legitimately spelt
    # with any of them, so they are refused at the write gate rather than
    # allowed to reach that interpolation.
    const SHELL_METACHARACTERS = ['\'', '"', '`', '$', ';', '&', '|', '<', '>',
                                  '(', ')', '{', '}', '[', ']', '*', '?', '!',
                                  '\\', '\n', '\r', '\0']

    is_shell_safe(p::AbstractString) = !any(c -> c in SHELL_METACHARACTERS, p)

    """
        is_shell_safe_arg(p) -> Bool

    `is_shell_safe`, and additionally free of spaces.

    A space is not a metacharacter - it cannot make a shell run anything - so a
    path carrying one is safe to store, and databases.yml deliberately accepts
    it. It is still not usable as an unquoted word in a command line: the shell
    would split it into two arguments and the far side would be handed a
    truncated path. Anything that becomes such a word - the remote staging root,
    and every key=value parameter of a remote stage - is checked with this.
    """
    is_shell_safe_arg(p::AbstractString) = is_shell_safe(p) && !occursin(' ', p)

    # A thread count is DADA2's `multithread`: a Bool selecting all cores or one,
    # or a strict positive Integer. Floats and strings are refused so that YAML's
    # `true`, `4` and `"4"` cannot silently pick different code paths in R.
    _is_thread_count(v) = (v isa Bool) || (v isa Integer && v >= 1)

    struct ValidationError
        context::String
        message::String
    end

    function _err(errors::Vector{ValidationError}, ctx::String, msg::String)
        push!(errors, ValidationError(ctx, msg))
    end

    _is_number(val) = val isa Number && !isnan(Float64(val))

    # `get` returns its default only for an absent key, so a null section arrives
    # as `nothing`. Never throws: it backs a write gate that must return a 400.
    _seq(v) = v isa AbstractVector ? v : Any[]

    ## Tools Validation
    function _validate_tools(errors::Vector{ValidationError}, tools_config_path::String)
        ctx = "tools"
        isfile(tools_config_path) || begin
            _err(errors, ctx, "tools.yml not found at $tools_config_path - run install.jl first")
            return
        end
        cfg = YAML.load_file(tools_config_path)
        cfg isa Dict || begin
            _err(errors, ctx, "tools.yml is not a valid YAML mapping")
            return
        end

        required = ["cutadapt", "fastqc", "multiqc", "vsearch", "cd_hit_est"]
        for key in required
            entry = get(cfg, key, nothing)
            path  = entry isa Dict ? get(entry, "path", nothing) : entry
            if isnothing(path)
                _err(errors, ctx, "$key: not configured in tools.yml")
                continue
            end
            resolved = isfile(string(path)) ? string(path) : Sys.which(string(path))
            if isnothing(resolved)
                _err(errors, ctx, "$key: path '$path' not found or not executable")
            end
        end
    end

    ## Database document rules
    # The single source of truth for the STRUCTURE of a databases document, shared
    # by the environment validator (_validate_databases) and
    # DatabasesLibrary.validate. Operates on the native document shape (a
    # databases: mapping of dir plus one entry per database). Returns
    # human-readable errors; empty means valid. Never throws: every section and
    # entry is type-checked before use.
    #
    # This is deliberately thin, and it holds structure only. Two other kinds of
    # rule deliberately live elsewhere:
    #
    # Rules that would fail a databases.yml which validates today (an empty levels
    # list, a misspelt vsearch_format) are WRITE-time rules and live in
    # DatabasesLibrary.validate.
    #
    # Rules about the ENVIRONMENT the document is read in, of which "a local: path
    # names a file that exists right now" is the only one, live in
    # _validate_databases below. A document is not invalid for naming a file that
    # has yet to arrive: the user may legitimately save the path first, and the
    # pipeline reports the absence when it goes looking. Keeping such a rule here
    # put it in the write gate, where one stale local: anywhere in the file
    # rejected every save of the whole document, and where the 400 disclosed to
    # any client whether an arbitrary path exists on the server.
    function database_document_errors(cfg::AbstractDict)::Vector{String}
        errors = String[]
        dbs = get(cfg, "databases", nothing)
        dbs isa AbstractDict ||
            push!(errors, "databases.yml missing 'databases:' key")
        errors
    end

    ## Database Validation
    # isfile throws on a NUL byte and on an unreadable parent directory (EACCES),
    # both reachable from a hand-edited databases.yml; either is treated as not-a-file.
    function _is_named_file(p::AbstractString)
        occursin('\0', p) && return false
        try
            isfile(p)
        catch
            false
        end
    end

    # Whether each configured local: path names a file that exists. This is the
    # environment rule the validator exists to report: it says the config points
    # at something that is not there, which is true of the machine and not of the
    # document. It is deliberately not a write-time rule; see the note above.
    function _validate_database_files(errors::Vector{ValidationError}, ctx::String, cfg::AbstractDict)
        dbs = get(cfg, "databases", nothing)
        dbs isa AbstractDict || return
        for (db_name, db_cfg) in dbs
            # `dir` is the shared cache directory, so no database may be called dir.
            string(db_name) == "dir" && continue
            db_cfg isa AbstractDict || continue
            for method in ("dada2", "vsearch")
                mc = get(db_cfg, method, nothing)
                mc isa AbstractDict || continue
                local_path = get(mc, "local", nothing)
                isnothing(local_path) && continue
                _is_named_file(string(local_path)) ||
                    _err(errors, ctx, "$db_name.$method.local: file not found: $local_path")
            end
        end
    end

    function _validate_databases(errors::Vector{ValidationError}, databases_config_path::String)
        ctx = "databases"
        isfile(databases_config_path) || begin
            _err(errors, ctx, "databases.yml not found at $databases_config_path")
            return
        end
        cfg = YAML.load_file(databases_config_path)
        cfg isa Dict || begin
            _err(errors, ctx, "databases.yml is not a valid YAML mapping")
            return
        end
        for msg in database_document_errors(cfg)
            _err(errors, ctx, msg)
        end
        _validate_database_files(errors, ctx, cfg)
    end

    ## Primer document rules
    # The single source of truth for what a valid primers document is, shared by
    # the environment validator (_validate_primers) and PrimersLibrary.validate.
    # Operates on the native document shape (Forward/Reverse maps, Pairs a list of
    # single-key mappings name => [fwd, rev]). Returns human-readable errors; empty
    # means valid. Never throws: every section and entry is type-checked before use.
    #
    # Primer and pair names are deliberately NOT charset-checked. They are only
    # ever lookup keys: get_primer_args resolves a pair name to its sequences, and
    # only the sequence (checked against the IUPAC set below) reaches the cutadapt
    # command. No name reaches a filesystem path, a shell, or a SQL identifier, so
    # a charset rule here would reject a primers.yml that the pipeline runs
    # perfectly well, and would buy nothing.
    function primer_document_errors(cfg::AbstractDict)::Vector{String}
        errors = String[]
        fwd = get(cfg, "Forward", Dict())
        rev = get(cfg, "Reverse", Dict())

        fwd isa AbstractDict || push!(errors, "'Forward' must be a mapping of name -> sequence")
        rev isa AbstractDict || push!(errors, "'Reverse' must be a mapping of name -> sequence")

        valid_bases = Set("ACGTMRWSYKVHDBNacgtmrwsykvhdbn")
        for (name, seq) in merge(fwd isa AbstractDict ? fwd : Dict(),
                                 rev isa AbstractDict ? rev : Dict())
            if !(seq isa AbstractString)
                push!(errors, "primer '$name' sequence is not a string")
                continue
            end
            bad = filter(c -> c ∉ valid_bases, seq)
            isempty(bad) ||
                push!(errors, "primer '$name' contains invalid bases: $(join(unique(bad)))")
        end

        pairs = get(cfg, "Pairs", [])
        pairs isa AbstractVector || return errors
        for entry in pairs
            entry isa AbstractDict || continue
            for (pair_name, members) in entry
                if !(members isa AbstractVector && length(members) == 2)
                    push!(errors, "pair '$pair_name' must list exactly [ForwardName, ReverseName]")
                    continue
                end
                f_name, r_name = string(members[1]), string(members[2])
                fwd isa AbstractDict && haskey(fwd, f_name) ||
                    push!(errors, "pair '$pair_name' references unknown forward primer '$f_name'")
                rev isa AbstractDict && haskey(rev, r_name) ||
                    push!(errors, "pair '$pair_name' references unknown reverse primer '$r_name'")
            end
        end
        errors
    end

    ## Primers Validation
    function _validate_primers(errors::Vector{ValidationError}, primers_path::String)
        ctx = "primers ($primers_path)"
        isfile(primers_path) || begin
            _err(errors, ctx, "primers.yml not found")
            return
        end
        cfg = YAML.load_file(primers_path)
        cfg isa Dict || begin
            _err(errors, ctx, "not a valid YAML mapping")
            return
        end
        for m in primer_document_errors(cfg)
            _err(errors, ctx, m)
        end
    end

    _fraction(v) = _is_number(v) && 0 <= v <= 1
    _whole(v) = v isa Integer && !(v isa Bool)

    function _validate_align(errors, a::Dict, ctx, path)
        st = get(a, "strategy", "auto")
        st in PHYLOGENY_ALIGN_STRATEGIES ||
            _err(errors, ctx, "$path.strategy must be one of " *
                              "$(join(PHYLOGENY_ALIGN_STRATEGIES, ", ")) (got: $(repr(st)))")
        mi = get(a, "maxiterate", 0)
        (_whole(mi) && mi >= 0) ||
            _err(errors, ctx, "$path.maxiterate must be a whole number (got: $(repr(mi)))")
    end

    function _validate_trim(errors, t::Dict, ctx, path)
        m = get(t, "method", "manual")
        m in PHYLOGENY_TRIM_METHODS ||
            _err(errors, ctx, "$path.method must be one of $(join(PHYLOGENY_TRIM_METHODS, ", ")) (got: $(repr(m)))")
        gt = get(t, "gap_threshold", 0.5)
        _fraction(gt) || _err(errors, ctx, "$path.gap_threshold must be between 0 and 1 (got: $(repr(gt)))")
        c = get(t, "conservation", nothing)
        isnothing(c) || (_is_number(c) && 0 <= c <= 100) ||
            _err(errors, ctx, "$path.conservation must be between 0 and 100, or empty (got: $(repr(c)))")
        for k in ("similarity_threshold", "residue_overlap")
            v = get(t, k, nothing)
            isnothing(v) || _fraction(v) ||
                _err(errors, ctx, "$path.$k must be between 0 and 1, or empty (got: $(repr(v)))")
        end
        so = get(t, "sequence_overlap", nothing)
        isnothing(so) || (_is_number(so) && 0 <= so <= 100) ||
            _err(errors, ctx, "$path.sequence_overlap must be between 0 and 100, or empty (got: $(repr(so)))")
        isnothing(get(t, "residue_overlap", nothing)) == isnothing(so) ||
            _err(errors, ctx, "$path.residue_overlap and $path.sequence_overlap are set together")
    end

    _model_ok(m) = m isa AbstractString && occursin(r"\A[A-Za-z0-9+_.{},/-]+\z", m)

    function _validate_phylogeny(errors::Vector{ValidationError}, ph::Dict, ctx::String)
        th = get(ph, "threads", nothing)
        isnothing(th) || (_whole(th) && th >= 1) ||
            _err(errors, ctx, "phylogeny.threads must be a positive integer (got: $(repr(th)))")
        sub(d, k) = (v = get(d, k, Dict()); v isa Dict ? v : Dict())
        ref, pl = sub(ph, "reference"), sub(ph, "placement")

        _validate_align(errors, sub(ref, "align"), ctx, "phylogeny.reference.align")
        _validate_trim(errors, sub(ref, "trim"), ctx, "phylogeny.reference.trim")
        tr = sub(ref, "tree")
        bs = get(tr, "bootstrap", "standard")
        bs in PHYLOGENY_BOOTSTRAPS ||
            _err(errors, ctx, "phylogeny.reference.tree.bootstrap must be standard or ultrafast (got: $(repr(bs)))")
        reps = get(tr, "replicates", 100)
        if !(_whole(reps) && reps >= 1)
            _err(errors, ctx, "phylogeny.reference.tree.replicates must be a positive integer (got: $(repr(reps)))")
        elseif bs == "ultrafast" && reps < 1000
            _err(errors, ctx, "phylogeny.reference.tree.replicates must be at least 1000 for ultrafast bootstrap (got: $reps)")
        end
        _model_ok(get(tr, "model", "MFP")) ||
            _err(errors, ctx, "phylogeny.reference.tree.model must be a model name such as MFP (got: $(repr(get(tr, "model", nothing))))")

        _validate_align(errors, sub(pl, "align"), ctx, "phylogeny.placement.align")
        _validate_trim(errors, sub(pl, "trim"), ctx, "phylogeny.placement.trim")
        place = sub(pl, "place")
        _model_ok(get(place, "model", "GTRCATI")) ||
            _err(errors, ctx, "phylogeny.placement.place.model must be a model name such as GTRCATI (got: $(repr(get(place, "model", nothing))))")
        h = get(place, "heuristic", nothing)
        isnothing(h) || (_is_number(h) && 0 < h <= 1) ||
            _err(errors, ctx, "phylogeny.placement.place.heuristic must be above 0 and at most 1, or empty (got: $(repr(h)))")
        at = get(sub(pl, "accumulate"), "threshold", 0.95)
        (_is_number(at) && 0.5 <= at <= 1) ||
            _err(errors, ctx, "phylogeny.placement.accumulate.threshold must be between 0.5 and 1 (got: $(repr(at)))")
    end

    ## Pipeline Config Validation
    """
        _staging_dir_errors!(errors, sd, ctx)

    Record why `sd` cannot be `remote.staging_dir`: it must be an absolute path
    a remote login shell would pass through as one word.
    """
    function _staging_dir_errors!(errors::Vector{ValidationError}, sd, ctx::String)
        if isnothing(sd) || !(sd isa AbstractString)
            _err(errors, ctx, "remote.staging_dir must be set when remote.host is set")
        else
            startswith(sd, "/") ||
                _err(errors, ctx, "remote.staging_dir must be an absolute path on the " *
                                  "server (got: $(repr(sd)))")
            # staging_dir is interpolated into the command string the
            # remote sshd runs through a login shell, so it is gated here
            # exactly as databases.yml gates dada2.remote_path.
            is_shell_safe_arg(sd) ||
                _err(errors, ctx, "remote.staging_dir may not contain shell " *
                                  "metacharacters or spaces")
        end
    end

    """
        _validate_remote!(errors, rm_cfg, ctx)

    Record what is wrong with the `remote` block `rm_cfg`: its host, stages,
    threads, staging root (once a host is named) and tool names.
    """
    function _validate_remote!(errors::Vector{ValidationError}, rm_cfg, ctx::String)
        if rm_cfg isa Dict
            host = get(rm_cfg, "host", nothing)
            isnothing(host) || host isa AbstractString ||
                _err(errors, ctx, "remote.host must be a string (got: $(repr(host)))")

            stages = get(rm_cfg, "stages", nothing)
            if !isnothing(stages)
                if stages isa Vector
                    for st in stages
                        st isa AbstractString && st in REMOTE_STAGES ||
                            _err(errors, ctx, "remote.stages entries must be one of " *
                                              "$(join(REMOTE_STAGES, ", ")) (got: $(repr(st)))")
                    end
                else
                    _err(errors, ctx, "remote.stages must be a list (got: $(repr(stages)))")
                end
            end

            rth = get(rm_cfg, "threads", nothing)
            isnothing(rth) || _is_thread_count(rth) ||
                _err(errors, ctx, "remote.threads must be a Bool or positive integer (got: $(repr(rth)))")

            # Only checked once a host is named: the factory default carries a
            # placeholder staging_dir, and refusing that would make every config
            # invalid until the user configures a server they may never want.
            if !isnothing(host)
                _staging_dir_errors!(errors, get(rm_cfg, "staging_dir", nothing), ctx)
            end
        else
            _err(errors, ctx, "remote must be a mapping (got: $(repr(rm_cfg)))")
        end

        rt_tools = rm_cfg isa Dict ? get(rm_cfg, "tools", nothing) : nothing
        if rt_tools isa Dict
            for (tool, bin) in rt_tools
                # Each name becomes the first word of a command the remote login shell runs.
                (bin isa AbstractString && is_shell_safe_arg(bin) && !startswith(bin, "-")) ||
                    _err(errors, ctx, "remote.tools.$tool must be a program name or path " *
                                      "without spaces or shell metacharacters (got: $(repr(bin)))")
            end
        elseif !isnothing(rt_tools)
            _err(errors, ctx, "remote.tools must be a mapping (got: $(repr(rt_tools)))")
        end
    end

    """
        remote_value_error(key, value) -> Union{String,Nothing}

    Why `value` cannot be stored as `remote.<key>` (e.g. `key = "tools.mafft"`),
    or nothing. Runs the checks `validate_project` makes of the remote block on
    that one key, so the config write gate refuses what validation would;
    `staging_dir` is checked even when no host is named yet.
    """
    function remote_value_error(key::AbstractString, value)
        errors = ValidationError[]
        if key == "staging_dir"
            _staging_dir_errors!(errors, value, "")
        else
            nested = foldr((k, inner) -> Dict{String,Any}(String(k) => inner), split(key, '.'); init=value)
            _validate_remote!(errors, nested, "")
        end
        isempty(errors) ? nothing : join((e.message for e in errors), "; ")
    end

    """
        _validate_pipeline_cfg(errors, cfg, ctx)

    Record what is wrong with a merged pipeline config `cfg`, section by section.
    """
    function _validate_pipeline_cfg(errors::Vector{ValidationError}, cfg::Dict, ctx::String)
        rt = get(cfg, "r_threads", nothing)
        isnothing(rt) || _is_thread_count(rt) ||
            _err(errors, ctx, "r_threads must be a Bool or positive integer (got: $(repr(rt)))")

        _validate_remote!(errors, get(cfg, "remote", Dict()), ctx)

        ph = get(cfg, "phylogeny", nothing)
        ph isa Dict && _validate_phylogeny(errors, ph, ctx)

        ca = get(cfg, "cutadapt", Dict())
        if ca isa Dict
            pp = get(ca, "primer_pairs", nothing)
            pp isa Vector && !isempty(pp) ||
                _err(errors, ctx, "cutadapt.primer_pairs must be a non-empty list")
            ml = get(ca, "min_length", nothing)
            (_is_number(ml) && ml > 0) ||
                _err(errors, ctx, "cutadapt.min_length must be a positive number (got: $ml)")
        end

        da = get(cfg, "dada2", Dict())
        if da isa Dict
            ft = get(da, "filter_trim", Dict())
            if ft isa Dict
                tl = get(ft, "trunc_len", nothing)
                ml = get(ft, "min_len",   nothing)
                if tl isa Vector && length(tl) >= 1 && _is_number(tl[1]) && _is_number(ml)
                    tl[1] > ml ||
                        _err(errors, ctx, "dada2.filter_trim.trunc_len[1] ($(tl[1])) must be > min_len ($ml)")
                end
                ee = get(ft, "max_ee", nothing)
                if ee isa Vector
                    all(x -> _is_number(x) && x >= 0, ee) ||
                        _err(errors, ctx, "dada2.filter_trim.max_ee values must be non-negative numbers")
                end
            end

            tx = get(da, "taxonomy", Dict())
            if tx isa Dict
                db = get(tx, "database", nothing)
                isnothing(db) || db isa String ||
                    _err(errors, ctx, "dada2.taxonomy.database must be a string")
                # Deprecated in favour of the top-level r_threads, but still
                # honoured for assign_taxonomy, so it is still validated.
                mb = get(tx, "multithread", nothing)
                isnothing(mb) || _is_thread_count(mb) ||
                    _err(errors, ctx,
                         "dada2.taxonomy.multithread must be a Bool or positive integer (got: $(repr(mb)))")
            end
        end

        asv = get(cfg, "asv", Dict())
        if asv isa Dict
            dm = get(asv, "denovo_method", nothing)
            isnothing(dm) || (dm isa AbstractString && dm in DENOVO_METHODS) ||
                _err(errors, ctx,
                     "asv.denovo_method must be one of $(join(DENOVO_METHODS, ", ")) (got: $(repr(dm)))")
        end

        vs = get(cfg, "vsearch", Dict())
        if vs isa Dict
            id = get(vs, "identity", nothing)
            isnothing(id) || (_is_number(id) && 0 < id <= 1) ||
                _err(errors, ctx, "vsearch.identity must be between 0 and 1 (got: $id)")
            qc = get(vs, "query_cov", nothing)
            isnothing(qc) || (_is_number(qc) && 0 < qc <= 1) ||
                _err(errors, ctx, "vsearch.query_cov must be between 0 and 1 (got: $qc)")
        end

        cd = get(cfg, "cdhit", Dict())
        if cd isa Dict
            id = get(cd, "identity", nothing)
            isnothing(id) || (_is_number(id) && 0 < id <= 1) ||
                _err(errors, ctx, "cdhit.identity must be between 0 and 1 (got: $id)")
        end

        sw = get(cfg, "swarm", Dict())
        if sw isa Dict
            d = get(sw, "differences", nothing)
            isnothing(d) || (_is_number(d) && d >= 0) ||
                _err(errors, ctx, "swarm.differences must be a non-negative integer (got: $d)")
            id = get(sw, "identity", nothing)
            isnothing(id) || (_is_number(id) && 0 < id <= 1) ||
                _err(errors, ctx, "swarm.identity must be between 0 and 1 (got: $id)")
        end
    end

    ## Per-project Validation
    """
        validate_project(project, databases_config_path) -> Vector{ValidationError}

    Validate a single project: data files present, config coherent, primer pairs defined.
    """
    function validate_project(project::ProjectCtx,
                               databases_config_path::String)::Vector{ValidationError}
        errors = ValidationError[]
        ctx    = basename(project.dir)

        for d in project.data_dirs
            isdir(d) ||
                _err(errors, ctx, "data directory not found: $d")
        end
        fastqs = find_fastqs(project)
        isempty(fastqs) &&
            _err(errors, ctx, "no .fastq.gz files found in $(join(project.data_dirs, ", "))")

        primers_path = joinpath(project.config_dir, "primers.yml")
        _validate_primers(errors, primers_path)

        try
            config_path = write_run_config(project)
            cfg         = YAML.load_file(config_path)
            _validate_pipeline_cfg(errors, cfg, ctx)

            # Check cutadapt primer_pairs references exist in primers.yml
            ca = get(cfg, "cutadapt", Dict())
            pp = get(ca,  "primer_pairs", String[])
            if pp isa Vector && isfile(primers_path)
                pcfg  = YAML.load_file(primers_path)
                pairs = get(pcfg, "Pairs", [])
                defined_pairs = Set{String}()
                for entry in (pairs isa Vector ? pairs : [])
                    entry isa Dict && union!(defined_pairs, string.(keys(entry)))
                end
                for name in pp
                    string(name) in defined_pairs ||
                        _err(errors, ctx, "cutadapt.primer_pairs references '$name' which is not defined in primers.yml")
                end
            end

            da      = get(cfg, "dada2",    Dict())
            tx      = get(da,  "taxonomy", Dict())
            db_name = string(get(tx, "database", "pr2"))
            db_cfg  = YAML.load_file(databases_config_path)
            dbs     = get(db_cfg isa Dict ? db_cfg : Dict(), "databases", Dict())
            haskey(dbs isa Dict ? dbs : Dict(), db_name) ||
                _err(errors, ctx, "dada2.taxonomy.database '$db_name' not found in databases.yml")
        catch e
            _err(errors, ctx, "could not load merged config: $e")
        end

        return errors
    end

    ## Entry Point
    """
        validate_environment(projects, databases_config_path, tools_config_path)

    Validate tools, databases, and all projects. Logs all errors and returns
    the total count. Caller should abort if count > 0.
    """
    function validate_environment(projects::Vector{ProjectCtx},
                                   databases_config_path::String,
                                   tools_config_path::String)::Int
        all_errors = ValidationError[]

        _validate_tools(all_errors, tools_config_path)
        _validate_databases(all_errors, databases_config_path)

        for project in projects
            append!(all_errors, validate_project(project, databases_config_path))
        end

        if isempty(all_errors)
            @info "Validation: all checks passed"
            return 0
        end

        @error "Validation failed with $(length(all_errors)) error(s):"
        for e in all_errors
            @error "  [$(e.context)] $(e.message)"
        end
        return length(all_errors)
    end

end
