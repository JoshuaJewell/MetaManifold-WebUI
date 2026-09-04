# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Error Model
    # Bases in a filtered FASTQ: every fourth line from the second is a sequence.
    # derepFastq's abundance-weighted count - sum(uniques * nchar(names(uniques)))
    # - comes to the same total, since each read is counted once either way.
    function _fastq_bases(path::AbstractString)
        total = 0
        open(path) do raw
            stream = GzipDecompressorStream(raw)
            try
                for (i, line) in enumerate(eachline(stream))
                    (i - 2) % 4 == 0 && (total += length(line))
                end
            finally
                close(stream)
            end
        end
        return total
    end

    """
        _learn_errors_prefix(paths, nbases) -> Vector{String}

    The leading files `learnErrors` will actually consume.

    dada2's own loop dereplicates the files in the order given and breaks as soon
    as the cumulative base count exceeds `nbases`, then keeps `drps[1:i]`; with
    `randomize=FALSE`, which is what the local stage uses, that prefix is
    deterministic. Only those files contribute to the error model, so only those
    need to cross the wire - the model that comes back is identical either way.

    The counting stops exactly where the sending does, so no file is read that
    was not going to be sent, and a run whose reads never reach the budget sends
    all of them, as dada2 would use all of them.
    """
    function _learn_errors_prefix(paths::AbstractVector{<:AbstractString}, nbases::Real)
        cumulative = 0.0
        for (i, f) in enumerate(paths)
            cumulative += _fastq_bases(f)
            cumulative > nbases && return paths[1:i]
        end
        return paths
    end

    # The filtered reads and the manifests naming them, as upload pairs. Both
    # read-consuming remote stages send exactly this set, and each addresses the
    # files by basename under $staging/Filtered, so neither depends on this
    # machine's directory layout being reproduced on the server.
    #
    # `fwd`/`rev` default to the whole read set; learn_errors passes the prefix
    # its budget reaches instead.
    function _filtered_read_inputs(ctx, tmp::String; fwd=ctx.fwd_out, rev=ctx.rev_out)
        inputs = Pair{String,String}[]
        for (name, paths) in (("fwd", fwd), ("rev", rev))
            isempty(paths) && continue
            for f in paths
                push!(inputs, f => "Filtered/" * basename(f))
            end
            push!(inputs, _write_manifest(tmp, "$name.txt", basename.(paths)) => "$name.txt")
        end
        return inputs
    end

    # Learn the error model in the embedded R session. `threads` reaches
    # learnErrors, which forwards it to dada(): that path threads through
    # RcppParallel's in-process workers rather than forking, so raising it is
    # safe inside an R session hosted by a multithreaded Julia process. See the
    # BiocParallel note at the top of dada2_functions.r for the case that is not.
    function _learn_errors_local(ctx, emit, log_path::String, threads)
        R"rm(list=ls())"
        _source_r_functions(ctx)
        seed    = ctx.seed
        nbases  = ctx.cfg["dada"]["nbases"]
        max_con = ctx.cfg["dada"]["max_consist"]
        verbose = ctx.verbose
        fwd_out = ctx.fwd_out
        rev_out = ctx.rev_out

        R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
        try
            emit("Learning error rates")
            R"set.seed($seed, kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection')"

            if ctx.mode != "reverse"
                R"mm_fwd <- $fwd_out"
                _r_run_logged("fwd_errors <- learnErrors(mm_fwd, nbases=$(_r_lit(nbases)), " *
                              "MAX_CONSIST=$(_r_lit(max_con)), multithread=$(_r_lit(threads)), " *
                              "verbose=$(_r_lit(verbose)))")
            else
                R"fwd_errors <- NULL"
            end

            if ctx.mode != "forward"
                R"mm_rev <- $rev_out"
                _r_run_logged("rev_errors <- learnErrors(mm_rev, nbases=$(_r_lit(nbases)), " *
                              "MAX_CONSIST=$(_r_lit(max_con)), multithread=$(_r_lit(threads)), " *
                              "verbose=$(_r_lit(verbose)))")
            else
                R"rev_errors <- NULL"
            end

            error_pdf = joinpath(ctx.dirs["Figures"], "error_rates.pdf")
            R"plot_error_rates(fwd_errors, rev_errors, $error_pdf)"

            ckpt = ctx.ckpts["errors"]
            R"save(fwd_errors, rev_errors, file=$ckpt)"
        finally
            R"tryCatch({ sink(type='message'); sink(); close(con) }, error = function(e) NULL)"
        end
    end

    # Learn the error model on the bioserver. The filtered reads go up, the
    # error model and its diagnostic plot come back; nothing else crosses, and
    # nothing is assumed to be there already.
    function _learn_errors_remote(ctx, emit, target, log_path::String)
        nbases = ctx.cfg["dada"]["nbases"]
        mktempdir() do tmp
            # Ship only what the error model will read. On a study of any size
            # this is a fraction of the run: the budget is a fixed number of
            # bases, while the read set grows with every sample added.
            fwd = _learn_errors_prefix(ctx.fwd_out, nbases)
            rev = _learn_errors_prefix(ctx.rev_out, nbases)
            n_all = length(ctx.fwd_out) + length(ctx.rev_out)
            n_use = length(fwd) + length(rev)
            n_use < n_all &&
                emit("  Sending $n_use of $n_all read files - the rest are past " *
                     "the $(nbases)-base error-model budget")
            _run_remote_stage(emit, target, "learn_errors_remote.r",
                [
                    "nbases"      => string(nbases),
                    "max_consist" => string(ctx.cfg["dada"]["max_consist"]),
                    "seed"        => string(ctx.seed),
                    "multithread" => _mt_str(target.threads),
                    "verbose"     => string(ctx.verbose),
                ];
                inputs  = _filtered_read_inputs(ctx, tmp; fwd, rev),
                outputs = [
                    "error_rates.pdf"   => joinpath(ctx.dirs["Figures"], "error_rates.pdf"),
                    "ckpt_errors.RData" => ctx.ckpts["errors"],
                ],
                log_path, stage = "learn_errors")
        end
    end

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
        snap = _begin_section(config_path, stage_sections(:dada2_learn_errors), hash_file)

        log_path = joinpath(ctx.dirs["Logs"], "learn_errors.log")
        open(log_path, "w") do io; println(io, "=== learn_errors ===\nconfig: $config_path") end

        target = _remote_target(ctx.full_cfg, "learn_errors")
        if isnothing(target)
            _learn_errors_local(ctx, emit, log_path,
                                _r_threads(ctx.full_cfg; stage="learn_errors"))
        else
            emit("Learning error rates (remote: $(target.host))")
            _learn_errors_remote(ctx, emit, target, log_path)
        end
        _write_section_hash(config_path, stage_sections(:dada2_learn_errors), hash_file; snapshot=snap)
        emit("Written: $(joinpath(ctx.dirs["Figures"], "error_rates.pdf"))")
        emit("Checkpoint: $(ctx.ckpts["errors"])")
        emit("Log: $log_path")
        nothing
    end

    # Denoise in the embedded R session. dada() threads through RcppParallel,
    # not by forking, so `threads` is safe here for the same reason it is in
    # _learn_errors_local.
    function _denoise_local(ctx, emit, log_path::String, threads)
        verbose     = ctx.verbose
        fwd_out     = ctx.fwd_out
        rev_out     = ctx.rev_out
        errors_ckpt = ctx.ckpts["errors"]
        R"rm(list=ls())"
        _source_r_functions(ctx)
        R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
        try
            R"load($errors_ckpt)"

            emit("Denoising reads")
            R"set.seed($(ctx.seed), kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection')"
            pool_method = ctx.cfg["dada"]["pool_method"]

            # Bind the read-path vectors once; the dada and mergePairs commands
            # below reference these R names rather than inlining every path.
            ctx.mode != "reverse" && R"mm_fwd <- $fwd_out"
            ctx.mode != "forward" && R"mm_rev <- $rev_out"

            if ctx.mode != "reverse"
                _r_run_logged("dada_fwd <- dada(mm_fwd, err=fwd_errors, " *
                              "pool=$(_r_lit(pool_method)), multithread=$(_r_lit(threads)), " *
                              "verbose=$(_r_lit(verbose)))")
            else
                R"dada_fwd <- NULL"
            end

            if ctx.mode != "forward"
                _r_run_logged("dada_rev <- dada(mm_rev, err=rev_errors, " *
                              "pool=$(_r_lit(pool_method)), multithread=$(_r_lit(threads)), " *
                              "verbose=$(_r_lit(verbose)))")
            else
                R"dada_rev <- NULL"
            end

            # A single input comes back as a bare object; name it so the
            # sequence table carries the sample name.
            R"""
            name_one <- function(x, files) if (inherits(x, "dada") || is.data.frame(x)) setNames(list(x), basename(files)) else x
            if (!is.null(dada_fwd)) dada_fwd <- name_one(dada_fwd, mm_fwd)
            if (!is.null(dada_rev)) dada_rev <- name_one(dada_rev, mm_rev)
            """

            emit("Building sequence table")
            if ctx.mode == "paired"
                min_overlap   = ctx.cfg["merge"]["min_overlap"]
                max_mismatch  = ctx.cfg["merge"]["max_mismatch"]
                trim_overhang = ctx.cfg["merge"]["trim_overhang"]
                _r_run_logged("merged <- mergePairs(dada_fwd, mm_fwd, dada_rev, mm_rev" *
                              ", minOverlap=$(_r_lit(min_overlap))" *
                              ", maxMismatch=$(_r_lit(max_mismatch))" *
                              ", trimOverhang=$(_r_lit(trim_overhang))" *
                              ", verbose=$(_r_lit(verbose)))")
                R"merged <- name_one(merged, mm_fwd)"
                R"seq_table <- makeSequenceTable(merged)"
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
    end

    # Denoise on the bioserver. The error model is uploaded alongside the reads
    # rather than assumed to be there, so this stage runs on the server whether
    # or not learn_errors did.
    function _denoise_remote(ctx, emit, target, log_path::String)
        merge_cfg   = ctx.cfg["merge"]
        pool_method = ctx.cfg["dada"]["pool_method"]
        mktempdir() do tmp
            inputs = _filtered_read_inputs(ctx, tmp)
            push!(inputs, ctx.ckpts["errors"] => "ckpt_errors.RData")
            _run_remote_stage(emit, target, "denoise_remote.r",
                [
                    "mode"          => ctx.mode,
                    # dada(pool=) reads TRUE and "pseudo" differently, and the
                    # command line cannot carry that distinction on its own.
                    "pool_method"   => pool_method isa Bool ?
                                       (pool_method ? "true" : "false") : string(pool_method),
                    "pool_bool"     => string(pool_method isa Bool),
                    "min_overlap"   => string(merge_cfg["min_overlap"]),
                    "max_mismatch"  => string(merge_cfg["max_mismatch"]),
                    "trim_overhang" => string(merge_cfg["trim_overhang"]),
                    "seed"          => string(ctx.seed),
                    "multithread"   => _mt_str(target.threads),
                    "verbose"       => string(ctx.verbose),
                ];
                inputs,
                outputs = ["ckpt_denoise.RData" => ctx.ckpts["denoise"]],
                # Plotted only when the merge left something to plot, exactly as
                # the local stage does.
                optional_outputs = ["length_distribution.pdf" =>
                                    joinpath(ctx.dirs["Figures"], "length_distribution.pdf")],
                log_path, stage = "denoise")
        end
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
        snap = _begin_section(config_path, stage_sections(:dada2_denoise), hash_file)

        log_path = joinpath(ctx.dirs["Logs"], "denoise.log")
        open(log_path, "w") do io; println(io, "=== denoise ===\nconfig: $config_path") end

        target = _remote_target(ctx.full_cfg, "denoise")
        if isnothing(target)
            _denoise_local(ctx, emit, log_path, _r_threads(ctx.full_cfg; stage="denoise"))
        else
            emit("Denoising reads (remote: $(target.host))")
            _denoise_remote(ctx, emit, target, log_path)
        end
        _write_section_hash(config_path, stage_sections(:dada2_denoise), hash_file; snapshot=snap)
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
        snap = _begin_section(config_path, stage_sections(:dada2_filter_length), hash_file)

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
                _r_run_logged("seq_table <- filter_by_length(seq_table, " *
                              "$(_r_lit(band_min)), $(_r_lit(band_max)))")
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
        _write_section_hash(config_path, stage_sections(:dada2_filter_length), hash_file; snapshot=snap)
        emit("Checkpoint: $(ctx.ckpts["length"])")
        emit("Log: $log_path")
        nothing
    end
