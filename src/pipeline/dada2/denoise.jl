# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Error Model
    """
        learn_errors(config_path; progress)

    **Stage 3** - Learn substitution error rates and plot diagnostics.

    Review `Figures/error_rates.pdf`: the fitted line should closely follow the
    observed points. If not, increase `nbases` or `max_consist` in config and
    re-run this stage.

    Saves: `Checkpoints/ckpt_errors.RData`
    """
    function learn_errors(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        lbl     = ctx.run_label
        emit    = _emitter(progress, lbl)
        @info "[$(lbl)] DADA2: Learn errors starting"

        errors_ckpt = ctx.ckpts["errors"]
        filter_ckpt = ctx.ckpts["filter"]
        hash_file   = joinpath(ctx.dirs["Checkpoints"], "learn_errors.hash")
        if isfile(errors_ckpt) && isfile(filter_ckpt) &&
           !_section_stale(config_path, stage_sections(:dada2_learn_errors), hash_file) &&
           mtime(errors_ckpt) > mtime(filter_ckpt)
            @info "[$(lbl)] DADA2: Skipping learn errors - checkpoint up to date"
            return nothing
        end

        R"rm(list=ls())"
        _source_r_functions(ctx)
        seed    = get(ctx.cfg["dada"], "seed", 123)
        nbases  = ctx.cfg["dada"]["nbases"]
        max_con = ctx.cfg["dada"]["max_consist"]
        verbose = ctx.verbose
        fwd_out = ctx.fwd_out
        rev_out = ctx.rev_out

        log_path = joinpath(ctx.dirs["Logs"], "learn_errors.log")
        open(log_path, "w") do io; println(io, "=== learn_errors ===\nconfig: $config_path") end
        R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
        try
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
        finally
            R"tryCatch({ sink(type='message'); sink(); close(con) }, error = function(e) NULL)"
        end
        _write_section_hash(config_path, stage_sections(:dada2_learn_errors), hash_file)
        emit("Written: $(joinpath(ctx.dirs["Figures"], "error_rates.pdf"))")
        emit("Checkpoint: $(ctx.ckpts["errors"])")
        emit("Log: $log_path")
        nothing
    end

    ## Denoise
    """
        denoise(config_path; progress)

    **Stage 4** - Denoise reads, merge pairs (paired mode), build the sequence
    table, and plot the raw ASV length distribution.

    Review `Figures/length_distribution.pdf` and set `band_size_min` /
    `band_size_max` in config to target the expected amplicon peak before
    running `filter_length()`.

    Requires: `Checkpoints/ckpt_errors.RData`
    Saves: `Checkpoints/ckpt_denoise.RData` (unfiltered seq_table)
    """
    function denoise(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        lbl     = ctx.run_label
        emit    = _emitter(progress, lbl)
        @info "[$(lbl)] DADA2: Denoise starting"
        verbose = ctx.verbose
        fwd_out = ctx.fwd_out
        rev_out = ctx.rev_out

        isfile(ctx.ckpts["errors"]) ||
            error("Error model checkpoint not found. Run learn_errors() first.")
        errors_ckpt = ctx.ckpts["errors"]

        denoise_ckpt = ctx.ckpts["denoise"]
        hash_file    = joinpath(ctx.dirs["Checkpoints"], "denoise.hash")
        if isfile(denoise_ckpt) &&
           !_section_stale(config_path, stage_sections(:dada2_denoise), hash_file) &&
           mtime(denoise_ckpt) > mtime(errors_ckpt)
            @info "[$(lbl)] DADA2: Skipping denoise - checkpoint up to date"
            return nothing
        end

        R"rm(list=ls())"
        _source_r_functions(ctx)
        log_path = joinpath(ctx.dirs["Logs"], "denoise.log")
        open(log_path, "w") do io; println(io, "=== denoise ===\nconfig: $config_path") end
        R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
        try
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
            R"""
            if (sum(seq_table) > 0) {
                plot_length_distribution(seq_table, $len_dist_pdf)
            } else {
                message("Skipping length distribution plot: seq_table is empty after merging")
            }
            """

            ckpt = ctx.ckpts["denoise"]
            R"save(dada_fwd, dada_rev, merged, seq_table, file=$ckpt)"
        finally
            R"tryCatch({ sink(type='message'); sink(); close(con) }, error = function(e) NULL)"
        end
        _write_section_hash(config_path, stage_sections(:dada2_denoise), hash_file)
        emit("Written: $(joinpath(ctx.dirs["Figures"], "length_distribution.pdf"))")
        emit("Checkpoint: $(ctx.ckpts["denoise"])")
        emit("Log: $log_path")
        nothing
    end

    ## Length Filter
    """
        filter_length(config_path; progress)

    **Stage 5** - Optionally filter the sequence table by amplicon length and
    plot the filtered length distribution.

    Set `asv.band_size_min` and `asv.band_size_max` in config to the expected
    amplicon length range. If both are null this stage is a passthrough. Re-run
    this stage alone to adjust length cutoffs without re-running `denoise()`.

    Requires: `Checkpoints/ckpt_denoise.RData`
    Saves: `Checkpoints/ckpt_length.RData`
    """
    function filter_length(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        ctx  = _pipeline_context(config_path; input_dir, workspace_root)
        lbl  = ctx.run_label
        emit = _emitter(progress, lbl)
        @info "[$(lbl)] DADA2: Filter length starting"

        isfile(ctx.ckpts["denoise"]) ||
            error("Denoise checkpoint not found. Run denoise() first.")
        denoise_ckpt = ctx.ckpts["denoise"]

        length_ckpt = ctx.ckpts["length"]
        hash_file   = joinpath(ctx.dirs["Checkpoints"], "filter_length.hash")
        if isfile(length_ckpt) &&
           !_section_stale(config_path, stage_sections(:dada2_filter_length), hash_file) &&
           mtime(length_ckpt) > mtime(denoise_ckpt)
            @info "[$(lbl)] DADA2: Skipping filter length - checkpoint up to date"
            return nothing
        end

        R"rm(list=ls())"
        _source_r_functions(ctx)
        log_path = joinpath(ctx.dirs["Logs"], "filter_length.log")
        open(log_path, "w") do io; println(io, "=== filter_length ===\nconfig: $config_path") end
        R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
        try
            R"load($denoise_ckpt)"

            band_min = get(ctx.cfg["asv"], "band_size_min", nothing)
            band_max = get(ctx.cfg["asv"], "band_size_max", nothing)
            if !isnothing(band_min) && !isnothing(band_max)
                emit("Filtering by length: $band_min-$band_max bp")
                R"seq_table <- filter_by_length(seq_table, $band_min, $band_max)"
                len_filt_pdf = joinpath(ctx.dirs["Figures"], "length_distribution_filtered.pdf")
                R"""
                if (sum(seq_table) > 0) {
                    plot_length_distribution(seq_table, $len_filt_pdf)
                } else {
                    message("Skipping filtered length distribution plot: seq_table is empty after length filter")
                }
                """
                emit("Written: $len_filt_pdf")
            else
                emit("No length filter configured (band_size_min/max not set) - passing through")
            end

            ckpt = ctx.ckpts["length"]
            R"save(dada_fwd, dada_rev, merged, seq_table, file=$ckpt)"
        finally
            R"tryCatch({ sink(type='message'); sink(); close(con) }, error = function(e) NULL)"
        end
        _write_section_hash(config_path, stage_sections(:dada2_filter_length), hash_file)
        emit("Checkpoint: $(ctx.ckpts["length"])")
        emit("Log: $log_path")
        nothing
    end
