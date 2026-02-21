# Web UI: QC page
# Stages: prefilter_qc, filter_trim

    # Pre-filter quality assessment
    """
        prefilter_qc(config_path; progress)

    **Stage 1** - Plot unfiltered quality profiles.

    Review `Figures/quality_unfiltered.pdf` to choose `truncLen` and `maxEE`
    values in config before running `filter_trim()`.
    """
    function prefilter_qc(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        emit = _emitter(progress)
        ctx  = _pipeline_context(config_path; input_dir, workspace_root)

        unfiltered_pdf = joinpath(ctx.dirs["Figures"], "quality_unfiltered.pdf")
        if isfile(unfiltered_pdf)
            all_inputs  = vcat(ctx.fwd_files, ctx.rev_files)
            input_mtime = isempty(all_inputs) ? mtime(config_path) :
                          max(mtime(config_path), maximum(mtime(f) for f in all_inputs))
            if mtime(unfiltered_pdf) > input_mtime
                @info "Skipping prefilter_qc: quality_unfiltered.pdf up to date"
                return nothing
            end
        end

        log_path = joinpath(ctx.dirs["Logs"], "prefilter_qc.log")
        open(log_path, "w") do io; println(io, "=== prefilter_qc ===\nconfig: $config_path") end
        R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
        try
            emit("Plotting unfiltered quality profiles")
            fwd_for_plot   = isempty(ctx.fwd_files) ? nothing : ctx.fwd_files
            rev_for_plot   = isempty(ctx.rev_files) ? nothing : ctx.rev_files
            unfiltered_pdf = joinpath(ctx.dirs["Figures"], "quality_unfiltered.pdf")
            R"plot_quality_profiles($fwd_for_plot, $rev_for_plot, $unfiltered_pdf)"

        finally
            R"tryCatch({ sink(type='message'); sink(); close(con) }, error = function(e) NULL)"
        end
        emit("Written: $(joinpath(ctx.dirs["Figures"], "quality_unfiltered.pdf"))")
        emit("Log: $log_path")
        nothing
    end

    # Filter and trim
    """
        filter_trim(config_path; progress)

    **Stage 2** - Filter and trim reads; plot filtered quality profiles.

    Review `Figures/quality_filtered.pdf`. If filtering looks appropriate,
    proceed to `learn_errors()`; otherwise adjust `truncLen` / `maxEE` in
    config and re-run this stage.

    Saves: `Checkpoints/ckpt_filter.RData`
    """
    function filter_trim(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)

        filter_ckpt = ctx.ckpts["filter"]
        if isfile(filter_ckpt)
            all_inputs  = vcat(ctx.fwd_files, ctx.rev_files)
            input_mtime = isempty(all_inputs) ? mtime(config_path) :
                          max(mtime(config_path), maximum(mtime(f) for f in all_inputs))
            if mtime(filter_ckpt) > input_mtime
                @info "Skipping filter_trim: checkpoint up to date"
                return nothing
            end
        end

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

        log_path = joinpath(ctx.dirs["Logs"], "filter_trim.log")
        open(log_path, "w") do io; println(io, "=== filter_trim ===\nconfig: $config_path") end
        R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
        try
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
        finally
            R"tryCatch({ sink(type='message'); sink(); close(con) }, error = function(e) NULL)"
        end
        emit("Written: $(joinpath(ctx.dirs["Figures"], "quality_filtered.pdf"))")
        emit("Checkpoint: $(ctx.ckpts["filter"])")
        emit("Log: $log_path")
        nothing
    end
