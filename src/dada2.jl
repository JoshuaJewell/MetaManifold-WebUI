module DADA2

# DADA2 amplicon sequencing pipeline — Julia orchestrator
#
# The pipeline is split into six independently callable stages so that
# intermediate outputs (quality profiles, error rate plots, length
# distributions, pipeline stats) can be reviewed and parameters adjusted
# before committing to the next step. Each stage saves its R objects to a
# per-stage checkpoint in Analysis/ and loads what it needs from the previous
# stage's checkpoint, making every stage re-runnable across Julia sessions.
#
# Stage order:
#   prefilter_qc     → filter_trim    → learn_errors →
#   denoise          → chimera_removal → assign_taxonomy
#
# Call dada2() to run all stages in sequence without stopping.
#
# Notice:
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).
#
# This work is based on the DADA2 tutorial by Benjamin J. Callahan, et al.,
# available at https://benjjneb.github.io/dada2/tutorial.html, with modification
# into a single module. The original material is licensed under the Creative
# Commons Attribution 4.0 International License (CC BY 4.0):
# https://creativecommons.org/licenses/by/4.0/.

export dada2, prefilter_qc, filter_trim, learn_errors, denoise,
       chimera_removal, assign_taxonomy

    import Downloads
    using Logging, RCall, YAML

    # Helpers

    # Returns a function that routes progress messages to @info (CLI) or a
    # Channel{String} (Genie SSE). Pass the result as `emit` in stage functions.
    _emitter(::Nothing)           = msg -> @info msg
    _emitter(ch::Channel{String}) = msg -> put!(ch, msg)

    # Database resolution helpers

    # Download `uri` to `db_dir/basename(uri)` if not already cached.
    # Respects the optional `local:` override in `fmt_info`.
    function _download_db_if_needed(key, fmt_info, db_dir, emit)
        local_p = get(fmt_info, "local", nothing)
        if !isnothing(local_p)
            local_p = string(local_p)
            if !isempty(local_p)
                isfile(local_p) && (emit("[$key] Using local file: $local_p"); return local_p)
                @warn "[$key] Configured local path not found: $local_p — falling back to uri"
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
    # provided. Supports the `databases: / taxonomy.database:` scheme (preferred)
    # and the legacy `taxonomy.uri:` key for backwards compatibility.
    function _resolve_taxonomy_db(cfg, emit)
        tax_cfg = cfg["taxonomy"]
        if haskey(tax_cfg, "database")
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
        # Legacy: taxonomy.uri
        tax_uri = tax_cfg["uri"]
        if isfile(tax_uri)
            emit("Using local taxonomy database: $tax_uri")
            return tax_uri
        end
        db_dir = abspath(get(cfg, "databases_dir", "./databases"))
        mkpath(db_dir)
        cached = joinpath(db_dir, basename(tax_uri))
        if !isfile(cached)
            emit("Downloading taxonomy database: $tax_uri")
            Downloads.download(tax_uri, cached)
            emit("Written: $cached")
        else
            emit("Using cached taxonomy database: $cached")
        end
        return cached
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

        boot_mode = get(cfg["output"], "bootstraps", "combined")
        boot_mode in ("none", "combined", "separate") ||
            error("output.bootstraps must be one of: none, combined, separate")

        if !get(cfg["taxonomy"], "skip", false)
            has_uri = haskey(cfg["taxonomy"], "uri")
            has_db  = haskey(cfg["taxonomy"], "database")
            has_uri || has_db ||
                error("taxonomy: configure `database` (referencing a databases: entry) " *
                      "or `uri` when skip is not true")
        end

        combined_mode = get(cfg["output"], "combined_mode", "regular")
        combined_mode in ("regular", "alternative") ||
            error("output.combined_mode must be one of: regular, alternative")
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

    # Pipeline context:
    # Shared setup called at the start of every stage: sources R, loads and
    # validates config, discovers files, handles the single-sample fallback, and
    # computes all path variables. Returns a NamedTuple so stage functions can
    # extract what they need without repeating boilerplate.
    #
    # Optional overrides (used when main.jl manages directory layout):
    #   input_dir      — overrides cfg["workspace"]["input_dir"]
    #   workspace_root — overrides cfg["workspace"]["root"]
    function _pipeline_context(config_path::String; input_dir=nothing, workspace_root=nothing)
        functions_r = joinpath(@__DIR__, "dada2_functions.r")
        R"source($functions_r)"

        cfg = YAML.load_file(config_path)

        # Apply caller-supplied overrides before validation so that validate_config
        # sees the final paths (including the isdir check on input_dir).
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

        # For single-end modes, route the relevant reads into the forward slots.
        in_fwd      = mode != "reverse" ? fwd_files : rev_files
        out_fwd     = mode != "reverse" ? fwd_out   : rev_out
        in_rev_arg  = mode == "paired"  ? rev_files : nothing
        out_rev_arg = mode == "paired"  ? rev_out   : nothing

        ckpts = Dict(
            "filter"  => joinpath(dirs["Checkpoints"], "ckpt_filter.RData"),
            "errors"  => joinpath(dirs["Checkpoints"], "ckpt_errors.RData"),
            "denoise" => joinpath(dirs["Checkpoints"], "ckpt_denoise.RData"),
            "chimera" => joinpath(dirs["Checkpoints"], "ckpt_chimera.RData"),
        )

        (; cfg, verbose, mode, dirs, sample_names, single_sample,
        fwd_files, rev_files, fwd_out, rev_out,
        in_fwd, out_fwd, in_rev_arg, out_rev_arg, ckpts)
    end

    # Pre-filter quality assessment
    """
        prefilter_qc(config_path; progress)

    **Stage 1** — Plot unfiltered quality profiles.

    Review `Figures/quality_unfiltered.pdf` to choose `truncLen` and `maxEE`
    values in config before running `filter_trim()`.
    """
    function prefilter_qc(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        emit = _emitter(progress)
        ctx  = _pipeline_context(config_path; input_dir, workspace_root)

        emit("Plotting unfiltered quality profiles")
        fwd_for_plot   = isempty(ctx.fwd_files) ? nothing : ctx.fwd_files
        rev_for_plot   = isempty(ctx.rev_files) ? nothing : ctx.rev_files
        unfiltered_pdf = joinpath(ctx.dirs["Figures"], "quality_unfiltered.pdf")
        R"plot_quality_profiles($fwd_for_plot, $rev_for_plot, $unfiltered_pdf)"
        emit("Written: $unfiltered_pdf")
        nothing
    end

    # Filter and trim
    """
        filter_trim(config_path; progress)

    **Stage 2** — Filter and trim reads; plot filtered quality profiles.

    Review `Figures/quality_filtered.pdf`. If filtering looks appropriate,
    proceed to `learn_errors()`; otherwise adjust `truncLen` / `maxEE` in
    config and re-run this stage.

    Saves: `Analysis/ckpt_filter.RData`
    """
    function filter_trim(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        ft      = ctx.cfg["filter_trim"]
        trunc_len = ft["trunc_len"]
        max_ee    = ft["max_ee"]
        verbose   = ctx.verbose
        in_fwd    = ctx.in_fwd
        out_fwd   = ctx.out_fwd
        in_rev    = ctx.in_rev_arg
        out_rev   = ctx.out_rev_arg
        fwd_out   = ctx.fwd_out
        rev_out   = ctx.rev_out

        emit("Filtering and trimming reads")
        if ctx.mode == "paired"
            R"""
            filter_stats <- filterAndTrim(
                $in_fwd,  $out_fwd,
                $in_rev,  $out_rev,
                truncQ   = $(ft["trunc_q"]),
                truncLen = $trunc_len,
                maxEE    = $max_ee,
                minLen   = $(ft["min_len"]),
                maxN     = $(ft["max_n"]),
                matchIDs = $(ft["match_ids"]),
                rm.phix  = $(ft["rm_phix"]),
                verbose  = $verbose
            )
            """
        else
            R"""
            filter_stats <- filterAndTrim(
                $in_fwd, $out_fwd,
                truncQ   = $(ft["trunc_q"]),
                truncLen = $(trunc_len[1]),
                maxEE    = $(max_ee[1]),
                minLen   = $(ft["min_len"]),
                maxN     = $(ft["max_n"]),
                rm.phix  = $(ft["rm_phix"]),
                verbose  = $verbose
            )
            """
        end

        emit("Plotting filtered quality profiles")
        fwd_filt = isempty(fwd_out) ? nothing : fwd_out
        rev_filt = isempty(rev_out) ? nothing : rev_out
        filtered_pdf = joinpath(ctx.dirs["Figures"], "quality_filtered.pdf")
        R"plot_quality_profiles($fwd_filt, $rev_filt, $filtered_pdf)"

        ckpt = ctx.ckpts["filter"]
        R"save(filter_stats, file=$ckpt)"
        emit("Written: $filtered_pdf")
        emit("Checkpoint: $ckpt")
        nothing
    end

    # Error model
    """
        learn_errors(config_path; progress)

    **Stage 3** — Learn substitution error rates and plot diagnostics.

    Review `Figures/error_rates.pdf`: the fitted line should closely follow the
    observed points. If not, increase `nbases` or `max_consist` in config and
    re-run this stage.

    Saves: `Analysis/ckpt_errors.RData`
    """
    function learn_errors(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        seed    = get(ctx.cfg["dada"], "seed", 123)
        nbases  = ctx.cfg["dada"]["nbases"]
        max_con = ctx.cfg["dada"]["max_consist"]
        verbose = ctx.verbose
        fwd_out = ctx.fwd_out
        rev_out = ctx.rev_out

        emit("Learning error rates")
        R"set.seed($seed)"

        ctx.mode != "reverse" ?
            R"fwd_errors <- learnErrors($fwd_out, nbases=$nbases, MAX_CONSIST=$max_con, verbose=$verbose)" :
            R"fwd_errors <- NULL"

        ctx.mode != "forward" ?
            R"rev_errors <- learnErrors($rev_out, nbases=$nbases, MAX_CONSIST=$max_con, verbose=$verbose)" :
            R"rev_errors <- NULL"

        error_pdf = joinpath(ctx.dirs["Figures"], "error_rates.pdf")
        R"plot_error_rates(fwd_errors, rev_errors, $error_pdf)"

        ckpt = ctx.ckpts["errors"]
        R"save(fwd_errors, rev_errors, file=$ckpt)"
        emit("Written: $error_pdf")
        emit("Checkpoint: $ckpt")
        nothing
    end

    # Denoise
    """
        denoise(config_path; progress)

    **Stage 4** — Denoise reads, merge pairs (paired mode), build the sequence
    table, and plot the raw ASV length distribution.

    Review `Figures/length_distribution.pdf` and set `band_size_min` /
    `band_size_max` in config to target the expected amplicon peak. Length
    filtering is applied in `chimera_removal()`, so this stage does not need
    to be re-run when adjusting the length cutoff.

    Requires: `Analysis/ckpt_errors.RData`
    Saves: `Analysis/ckpt_denoise.RData` (unfiltered seq_table)
    """
    function denoise(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        verbose = ctx.verbose
        fwd_out = ctx.fwd_out
        rev_out = ctx.rev_out

        isfile(ctx.ckpts["errors"]) ||
            error("Error model checkpoint not found. Run learn_errors() first.")
        errors_ckpt = ctx.ckpts["errors"]
        R"load($errors_ckpt)"

        emit("Denoising reads")
        pool_method = ctx.cfg["dada"]["pool_method"]

        ctx.mode != "reverse" ?
            R"dada_fwd <- dada($fwd_out, err=fwd_errors, pool=$pool_method, verbose=$verbose)" :
            R"dada_fwd <- NULL"

        ctx.mode != "forward" ?
            R"dada_rev <- dada($rev_out, err=rev_errors, pool=$pool_method, verbose=$verbose)" :
            R"dada_rev <- NULL"

        emit("Building sequence table")
        if ctx.mode == "paired"
            min_overlap   = ctx.cfg["merge"]["min_overlap"]
            max_mismatch  = ctx.cfg["merge"]["max_mismatch"]
            trim_overhang = ctx.cfg["merge"]["trim_overhang"]
            R"""
            merged <- mergePairs(
                dada_fwd, $fwd_out,
                dada_rev, $rev_out,
                minOverlap   = $min_overlap,
                maxMismatch  = $max_mismatch,
                trimOverhang = $trim_overhang,
                verbose      = $verbose
            )
            seq_table <- makeSequenceTable(merged)
            """
        else
            R"""
            merged    <- NULL
            seq_table <- makeSequenceTable(if (!is.null(dada_fwd)) dada_fwd else dada_rev)
            """
        end

        len_dist_pdf = joinpath(ctx.dirs["Figures"], "length_distribution.pdf")
        R"plot_length_distribution(seq_table, $len_dist_pdf)"

        ckpt = ctx.ckpts["denoise"]
        R"save(dada_fwd, dada_rev, merged, seq_table, file=$ckpt)"
        emit("Written: $len_dist_pdf")
        emit("Checkpoint: $ckpt")
        nothing
    end

    # Chimera removal
    """
        chimera_removal(config_path; progress)

    **Stage 5** — Optionally filter ASVs by length, remove chimeras, compute
    pipeline statistics, and write core output files.

    Length filtering (`band_size_min` / `band_size_max`) is applied here from
    the saved unfiltered seq_table, so this stage can be re-run with different
    length cutoffs without re-running `denoise()`.

    Review `Tables/pipeline_stats.csv` for unexpected read loss at any stage
    before committing to the (potentially long) `assign_taxonomy()` step.

    Requires: `Checkpoints/ckpt_filter.RData`, `Checkpoints/ckpt_denoise.RData`
    Saves: `Checkpoints/ckpt_chimera.RData`
    """
    function chimera_removal(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        verbose = ctx.verbose
        mode    = ctx.mode

        isfile(ctx.ckpts["filter"]) ||
            error("Filter checkpoint not found. Run filter_trim() first.")
        isfile(ctx.ckpts["denoise"]) ||
            error("Denoise checkpoint not found. Run denoise() first.")
        filter_ckpt  = ctx.ckpts["filter"]
        denoise_ckpt = ctx.ckpts["denoise"]
        R"load($filter_ckpt)"
        R"load($denoise_ckpt)"

        # Length filtering applied from the saved unfiltered seq_table.
        # Re-running this stage with new band sizes does not require re-denoising.
        band_min = get(ctx.cfg["asv"], "band_size_min", nothing)
        band_max = get(ctx.cfg["asv"], "band_size_max", nothing)
        if !isnothing(band_min) && !isnothing(band_max)
            emit("Filtering by length: $band_min-$band_max bp")
            R"seq_table <- filter_by_length(seq_table, $band_min, $band_max)"
            len_filt_pdf = joinpath(ctx.dirs["Figures"], "length_distribution_filtered.pdf")
            R"plot_length_distribution(seq_table, $len_filt_pdf)"
            emit("Written: $len_filt_pdf")
        end

        emit("Removing chimeras")
        denovo_method = ctx.cfg["asv"]["denovo_method"]
        R"""
        seq_table_nochim <- removeBimeraDenovo(seq_table, method=$denovo_method, verbose=$verbose)
        nochim_pct <- sum(seq_table_nochim) / sum(seq_table) * 100
        message("  Chimeric reads removed: ", round(100 - nochim_pct, 2),
                "% | Retained: ", round(nochim_pct, 2), "%")
        """

        # Drop duplicate sample row (single-sample fallback)
        final_names = ctx.sample_names
        if ctx.single_sample
            R"filter_stats     <- filter_stats[1, , drop=FALSE]"
            mode != "reverse" && R"dada_fwd <- dada_fwd[1]"
            mode != "forward" && R"dada_rev <- dada_rev[1]"
            mode == "paired"  && R"merged   <- merged[1]"
            R"seq_table_nochim <- seq_table_nochim[1, , drop=FALSE]"
            final_names = [ctx.sample_names[1]]
        end

        emit("Computing pipeline stats")
        stats_csv = joinpath(ctx.dirs["Tables"], "pipeline_stats.csv")
        R"""
        stats <- compute_pipeline_stats(filter_stats, dada_fwd, dada_rev, merged,
                                        seq_table_nochim, $final_names, $mode)
        write.csv(stats, $stats_csv, quote=FALSE)
        if ($verbose) print(stats)
        """

        emit("Writing core output tables")
        seq_prefix   = get(ctx.cfg["output"], "seq_table_prefix", "seqtab_nochim")
        fasta_prefix = get(ctx.cfg["output"], "fasta_prefix", "asvs")
        tables_dir   = ctx.dirs["Tables"]
        R"write_seq_table(seq_table_nochim, $tables_dir, $seq_prefix)"
        R"index <- write_fasta(seq_table_nochim, $tables_dir, $fasta_prefix)"

        ckpt = ctx.ckpts["chimera"]
        R"save(seq_table_nochim, index, file=$ckpt)"
        emit("Written: $stats_csv")
        emit("Written: $(joinpath(tables_dir, seq_prefix * ".csv"))")
        emit("Written: $(joinpath(tables_dir, fasta_prefix * ".fasta"))")
        emit("Checkpoint: $ckpt")
        nothing
    end

    # Taxonomy
    """
        assign_taxonomy(config_path; progress)

    **Stage 6** — Assign taxonomy to ASVs and write the combined output table.

    This is typically the longest step. Set `taxonomy.skip = true` in config to
    skip assignment and output sequence/count data only.

    Requires: `Checkpoints/ckpt_chimera.RData`
    Saves: `Checkpoints/checkpoint.RData` (full R environment snapshot)
    """
    function assign_taxonomy(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing, taxonomy_db=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        verbose = ctx.verbose

        isfile(ctx.ckpts["chimera"]) ||
            error("Chimera checkpoint not found. Run chimera_removal() first.")

        # Drop all data objects accumulated from prior stages before the
        # memory-intensive taxonomy subprocess runs. Named globals in R's
        # environment are reachable and gc() won't collect them; rm() them
        # explicitly so they don't inflate the subprocess's memory footprint.
        # lsf.str() returns function names; setdiff keeps those intact.
        R"""
        .data_objs <- setdiff(ls(), lsf.str())
        if (length(.data_objs) > 0L) rm(list = .data_objs)
        rm(.data_objs)
        gc()
        """

        chimera_ckpt = ctx.ckpts["chimera"]
        R"load($chimera_ckpt)"

        seq_prefix    = get(ctx.cfg["output"], "seq_table_prefix", "seqtab_nochim")
        fasta_prefix  = get(ctx.cfg["output"], "fasta_prefix", "asvs")
        taxa_prefix   = get(ctx.cfg["output"], "taxa_prefix", "taxonomy")
        combined_file = get(ctx.cfg["output"], "combined_filename", "tax_counts.xlsx")
        boot_mode     = get(ctx.cfg["output"], "bootstraps", "combined")
        comb_mode     = get(ctx.cfg["output"], "combined_mode", "regular")
        tables_dir    = ctx.dirs["Tables"]

        R"combined_input <- index"

        if !get(ctx.cfg["taxonomy"], "skip", false)
            emit("Assigning taxonomy")
            multithread = get(ctx.cfg["taxonomy"], "multithread", 4)
            min_boot    = get(ctx.cfg["taxonomy"], "min_boot", 0)
            tax_levels  = ctx.cfg["taxonomy"]["levels"]

            db_path = isnothing(taxonomy_db) ?
                _resolve_taxonomy_db(ctx.cfg, emit) : taxonomy_db

            R"""
            taxa_result <- run_assign_taxonomy(
                seq_table_nochim, $db_path,
                list(multithread=$multithread, min_boot=$min_boot, levels=$tax_levels),
                $verbose)
            taxa_df <- write_taxa_table(taxa_result$tax, taxa_result$boot, index,
                                        $tables_dir, $taxa_prefix, $boot_mode)
            """
            R"gc()"
            comb_mode == "regular" && R"combined_input <- taxa_df"
        else
            emit("Skipping taxonomy (taxonomy.skip = true)")
        end

        checkpoint = joinpath(ctx.dirs["Checkpoints"], "checkpoint.RData")
        R"save(seq_table_nochim, index, combined_input, file=$checkpoint)"
        emit("Checkpoint: $checkpoint")

        R"write_combined_table(combined_input, seq_table_nochim, $tables_dir, $combined_file)"

        emit("Pipeline complete. Outputs:")
        emit("  $(joinpath(tables_dir, seq_prefix * ".csv"))")
        emit("  $(joinpath(tables_dir, fasta_prefix * ".fasta"))")
        emit("  $(joinpath(tables_dir, fasta_prefix * ".csv"))")
        emit("  $(joinpath(tables_dir, taxa_prefix * ".csv"))")
        emit("  $(joinpath(tables_dir, combined_file))")
        emit("  $(joinpath(tables_dir, "pipeline_stats.csv"))")
        emit("  $checkpoint")
        nothing
    end

    # Run pipeline
    """
        dada2(config_path; progress)

    Run the complete DADA2 pipeline from raw reads to a taxonomy-annotated ASV
    table, calling all six stages in sequence:

        prefilter_qc → filter_trim → learn_errors →
        denoise → chimera_removal → assign_taxonomy

    For interactive use — reviewing intermediate outputs or adjusting parameters
    between steps — call the individual stage functions directly instead.

    ## Outputs
    Written to `workspace.root/Tables/`:
    - `seqtab_nochim.csv`       — chimera-free ASV count table
    - `asvs.fasta` / `asvs.csv` — ASV sequences with short identifiers
    - `taxonomy.csv`            — taxonomy assignments (with optional bootstraps)
    - `tax_counts.xlsx`         — combined taxonomy + per-sample counts
    - `pipeline_stats.csv`      — read counts at each pipeline stage

    Written to `workspace.root/Checkpoints/`:
    - `ckpt_filter.RData`   — filter_stats
    - `ckpt_errors.RData`   — fwd_errors, rev_errors
    - `ckpt_denoise.RData`  — dada objects, merged, unfiltered seq_table
    - `ckpt_chimera.RData`  — seq_table_nochim, index
    - `checkpoint.RData`    — full R environment snapshot
    """
    function dada2(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing, taxonomy_db=nothing)
        prefilter_qc(config_path;    progress, input_dir, workspace_root)
        filter_trim(config_path;     progress, input_dir, workspace_root);     R"gc()"
        learn_errors(config_path;    progress, input_dir, workspace_root);     R"gc()"
        denoise(config_path;         progress, input_dir, workspace_root);     R"gc()"
        chimera_removal(config_path; progress, input_dir, workspace_root);     R"gc()"
        assign_taxonomy(config_path; progress, input_dir, workspace_root, taxonomy_db)
    end

end