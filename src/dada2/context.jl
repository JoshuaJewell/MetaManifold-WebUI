# Shared pipeline context - helpers and setup called at the start of every stage.
# Included into module DADA2 by src/dada2.jl.

    # Emitter

    # Returns a function that routes progress messages to @info (CLI) or a
    # Channel{String} (Genie SSE). Pass the result as `emit` in stage functions.
    _emitter(::Nothing)           = msg -> @info msg
    _emitter(ch::Channel{String}) = msg -> put!(ch, msg)

    # Database helpers

    # Download uri to db_dir/basename(uri) if not already cached.
    # Respects the optional local: override in fmt_info.
    function _download_db_if_needed(key, fmt_info, db_dir, emit)
        local_p = get(fmt_info, "local", nothing)
        if !isnothing(local_p)
            local_p = string(local_p)
            if !isempty(local_p)
                isfile(local_p) && (emit("[$key] Using local file: $local_p"); return local_p)
                @warn "[$key] Configured local path not found: $local_p - falling back to uri"
            end
        end
        uri    = fmt_info["uri"]
        cached = joinpath(db_dir, basename(uri))
        if isfile(cached)
            emit("[$key] Using cached: $cached")
        else
            emit("[$key] Downloading: $uri")
            Downloads.download(uri, cached)
            emit("[$key] Saved to: $cached")
        end
        return cached
    end

    # Resolve the DADA2 taxonomy database from config when no external path is
    # provided. Reads taxonomy.database key and looks it up in the databases: section.
    function _resolve_taxonomy_db(cfg, emit)
        tax_cfg = cfg["taxonomy"]
        db_key  = string(tax_cfg["database"])
        db_cfg  = get(cfg, "databases", Dict())
        haskey(db_cfg, db_key) ||
            error("taxonomy.database = \"$db_key\" not found in databases: section")
        fmt_cfg = get(db_cfg[db_key], "dada2", nothing)
        isnothing(fmt_cfg) &&
            error("databases.$db_key.dada2 is not configured")
        db_dir = abspath(get(db_cfg, "dir", "./databases"))
        mkpath(db_dir)
        return _download_db_if_needed("$(db_key)_dada2", fmt_cfg, db_dir, emit)
    end

    # Config
    function validate_config(cfg)
        required = ["workspace", "file_patterns", "filter_trim", "dada",
                    "merge", "asv", "taxonomy", "output"]
        missing_secs = filter(k -> !haskey(cfg, k), required)
        isempty(missing_secs) ||
            error("Missing required config sections: $(join(missing_secs, ", "))")

        ws = cfg["workspace"]
        haskey(ws, "root")      || error("workspace.root is required")
        haskey(ws, "input_dir") || error("workspace.input_dir is required")
        isdir(ws["input_dir"])  ||
            error("workspace.input_dir does not exist: $(ws["input_dir"])")

        mode = get(cfg["file_patterns"], "mode", "paired")
        mode in ("paired", "forward", "reverse") ||
            error("file_patterns.mode must be one of: paired, forward, reverse")

        haskey(cfg["taxonomy"], "database") ||
            error("taxonomy.database is required")
    end

    # File discovery
    function find_fastq_files(input_dir, fwd_pattern, rev_pattern, mode)
        all_files = sort(readdir(input_dir, join=true))
        fwd = mode in ("paired", "forward") ?
            filter(f -> !isnothing(match(Regex(fwd_pattern), basename(f))), all_files) : String[]
        rev = mode in ("paired", "reverse") ?
            filter(f -> !isnothing(match(Regex(rev_pattern), basename(f))), all_files) : String[]
        fwd, rev
    end

    function validate_sample_files(fwd, rev, mode)
        if mode in ("paired", "forward")
            isempty(fwd) &&
                error("No forward FASTQ files found. " *
                    "Check workspace.input_dir and file_patterns.forward.")
            for f in fwd
                isfile(f) || error("Forward file not found: $f")
            end
        end
        if mode in ("paired", "reverse")
            isempty(rev) &&
                error("No reverse FASTQ files found. " *
                    "Check workspace.input_dir and file_patterns.reverse.")
            for f in rev
                isfile(f) || error("Reverse file not found: $f")
            end
        end
        if mode == "paired" && length(fwd) != length(rev)
            error("Forward file count ($(length(fwd))) does not match " *
                "reverse file count ($(length(rev))).")
        end
    end

    function extract_sample_names(files, split_char, split_index)
        [split(basename(f), split_char)[split_index] for f in files]
    end

    # Workspace
    function setup_workspace(root)
        dirs = Dict(
            "Tables"      => joinpath(root, "Tables"),
            "Checkpoints" => joinpath(root, "Checkpoints"),
            "Figures"     => joinpath(root, "Figures"),
            "Filtered"    => joinpath(root, "Filtered"),
        )
        for d in values(dirs)
            mkpath(d)
        end
        dirs
    end

    # Pipeline context
    # Shared setup called at the start of every stage: sources R functions, loads
    # and validates config, discovers files, handles the single-sample fallback, and
    # computes all path variables. Returns a NamedTuple so stage functions can
    # extract what they need without repeating boilerplate.
    #
    # Optional overrides (used when main.jl manages directory layout):
    #   input_dir      - overrides cfg["workspace"]["input_dir"]
    #   workspace_root - overrides cfg["workspace"]["root"]
    function _pipeline_context(config_path::String; input_dir=nothing, workspace_root=nothing)
        functions_r = joinpath(@__DIR__, "dada2_functions.r")
        R"source($functions_r)"

        cfg = YAML.load_file(config_path)

        if !isnothing(input_dir)
            cfg["workspace"]["input_dir"] = input_dir
        end
        if !isnothing(workspace_root)
            cfg["workspace"]["root"] = workspace_root
        end

        validate_config(cfg)
        verbose = get(cfg, "verbose", true)
        mode    = get(cfg["file_patterns"], "mode", "paired")
        root    = cfg["workspace"]["root"]
        dirs    = setup_workspace(root)

        fwd_files, rev_files = find_fastq_files(
            cfg["workspace"]["input_dir"],
            cfg["file_patterns"]["forward"],
            cfg["file_patterns"]["reverse"],
            mode)
        validate_sample_files(fwd_files, rev_files, mode)

        primary_files = isempty(fwd_files) ? rev_files : fwd_files
        sample_names  = extract_sample_names(
            primary_files,
            cfg["file_patterns"]["sample_name_split"],
            cfg["file_patterns"]["sample_name_index"])

        # Single-sample fallback: dada() returns a bare object (not a list) for a
        # single input file, breaking makeSequenceTable() and sapply() downstream.
        # Duplicate the paths so the pipeline sees 2 samples; the extra row is
        # dropped in chimera_removal(). ASV calls are unaffected.
        single_sample = length(sample_names) == 1
        if single_sample
            @warn "Only 1 sample found. Duplicating it to work around dada() " *
                "returning a bare object for single-file input. The duplicate " *
                "will be dropped from all outputs."
            fwd_files    = isempty(fwd_files) ? fwd_files : repeat(fwd_files, 2)
            rev_files    = isempty(rev_files) ? rev_files : repeat(rev_files, 2)
            sample_names = [sample_names[1], sample_names[1] * "_dup"]
        end

        filtered_dir = dirs["Filtered"]
        fwd_out = mode != "reverse" ?
            [joinpath(filtered_dir, s * "_R1_filt.fastq.gz") for s in sample_names] : String[]
        rev_out = mode != "forward" ?
            [joinpath(filtered_dir, s * "_R2_filt.fastq.gz") for s in sample_names] : String[]

        in_fwd      = mode != "reverse" ? fwd_files : rev_files
        out_fwd     = mode != "reverse" ? fwd_out   : rev_out
        in_rev_arg  = mode == "paired"  ? rev_files : nothing
        out_rev_arg = mode == "paired"  ? rev_out   : nothing

        ckpts = Dict(
            "filter"  => joinpath(dirs["Checkpoints"], "ckpt_filter.RData"),
            "errors"  => joinpath(dirs["Checkpoints"], "ckpt_errors.RData"),
            "denoise" => joinpath(dirs["Checkpoints"], "ckpt_denoise.RData"),
            "length"  => joinpath(dirs["Checkpoints"], "ckpt_length.RData"),
            "chimera" => joinpath(dirs["Checkpoints"], "ckpt_chimera.RData"),
        )

        (; cfg, verbose, mode, dirs, sample_names, single_sample,
        fwd_files, rev_files, fwd_out, rev_out,
        in_fwd, out_fwd, in_rev_arg, out_rev_arg, ckpts)
    end
