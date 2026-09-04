# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Chimera Removal
    # Remove chimeras in the embedded R session. removeBimeraDenovo forwards
    # `multithread` to isBimeraDenovo/isBimeraDenovoTable, both of which thread
    # through RcppParallel rather than forking, so this is safe inside the
    # shared R session.
    function _chimera_removal_local(ctx, emit, log_path::String, threads, denovo_method)
        lbl          = ctx.run_label
        verbose      = ctx.verbose
        mode         = ctx.mode
        seq_prefix   = get(ctx.cfg["output"], "seq_table_prefix", "seqtab_nochim")
        fasta_prefix = get(ctx.cfg["output"], "fasta_prefix", "asvs")
        tables_dir   = ctx.dirs["Tables"]
        filter_ckpt  = ctx.ckpts["filter"]
        length_ckpt  = ctx.ckpts["length"]
        R"rm(list=ls())"
        _source_r_functions(ctx)

        R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
        try
            R"load($filter_ckpt)"
            R"load($length_ckpt)"

            has_data = rcopy(R"isTRUE(sum(seq_table, na.rm=TRUE) > 0)")
            emit("Removing chimeras")
            if has_data
                _r_run_logged("seq_table_nochim <- removeBimeraDenovo(seq_table, " *
                              "method=$(_r_lit(denovo_method)), " *
                              "multithread=$(_r_lit(threads)), verbose=$(_r_lit(verbose)))")
                R"""
                nochim_pct <- sum(seq_table_nochim) / sum(seq_table) * 100
                message("  Chimeric reads removed: ", round(100 - nochim_pct, 2),
                        "% | Retained: ", round(nochim_pct, 2), "%")
                """
            else
                R"seq_table_nochim <- seq_table"
                @info "[$(lbl)] DADA2: Chimera removal skipped - seq_table is empty"
            end

            final_names = ctx.sample_names
            kept        = ctx.kept

            emit("Computing pipeline stats")
            stats_csv = joinpath(ctx.dirs["Tables"], "pipeline_stats.csv")
            R"""
            stats <- compute_pipeline_stats(filter_stats, dada_fwd, dada_rev, merged,
                                            seq_table_nochim, $final_names, $mode, $kept)
            write.csv(stats, $stats_csv, quote=FALSE)
            if ($verbose) print(stats)
            """

            emit("Writing core output tables")
            if has_data
                R"write_seq_table(seq_table_nochim, $tables_dir, $seq_prefix)"
                R"index <- write_fasta(seq_table_nochim, $tables_dir, $fasta_prefix)"
            else
                R"index <- data.frame(SeqName=character(0), sequence=character(0))"
                touch(joinpath(tables_dir, seq_prefix   * ".csv"))
                touch(joinpath(tables_dir, fasta_prefix * ".fasta"))
                touch(joinpath(tables_dir, fasta_prefix * ".csv"))
            end

            ckpt = ctx.ckpts["chimera"]
            R"save(seq_table_nochim, index, file=$ckpt)"
        finally
            R"tryCatch({ sink(type='message'); sink(); close(con) }, error = function(e) NULL)"
        end
    end

    # Remove chimeras on the bioserver. Two checkpoints go up and the tables
    # come back: no read ever crosses, which makes this the cheapest of the four
    # stages to offload however many samples the run carries.
    function _chimera_removal_remote(ctx, emit, target, log_path::String, denovo_method)
        seq_prefix   = get(ctx.cfg["output"], "seq_table_prefix", "seqtab_nochim")
        fasta_prefix = get(ctx.cfg["output"], "fasta_prefix", "asvs")
        tables_dir   = ctx.dirs["Tables"]
        mktempdir() do tmp
            _run_remote_stage(emit, target, "chimera_remote.r",
                [
                    "mode"          => ctx.mode,
                    "denovo_method" => denovo_method,
                    "seq_prefix"    => seq_prefix,
                    "fasta_prefix"  => fasta_prefix,
                    "multithread"   => _mt_str(target.threads),
                    "verbose"       => string(ctx.verbose),
                ];
                inputs = [
                    ctx.ckpts["filter"] => "ckpt_filter.RData",
                    ctx.ckpts["length"] => "ckpt_length.RData",
                    _write_manifest(tmp, "samples.txt", ctx.sample_names) => "samples.txt",
                    _write_manifest(tmp, "kept.txt", string.(ctx.kept)) => "kept.txt",
                ],
                # The checkpoint comes back last, so a failed fetch never leaves a
                # new checkpoint beside old tables.
                outputs = [
                    "Tables/pipeline_stats.csv"     => joinpath(tables_dir, "pipeline_stats.csv"),
                    "Tables/$seq_prefix.csv"        => joinpath(tables_dir, "$seq_prefix.csv"),
                    "Tables/$fasta_prefix.fasta"    => joinpath(tables_dir, "$fasta_prefix.fasta"),
                    "Tables/$fasta_prefix.csv"      => joinpath(tables_dir, "$fasta_prefix.csv"),
                    "ckpt_chimera.RData"            => ctx.ckpts["chimera"],
                ],
                log_path, stage = "chimera_removal")
        end
    end

    """
        chimera_removal(config_path; progress)

    **Stage 6** - Remove chimeric sequences, drop the single-sample duplicate
    row if present, and write core output files (seq table, FASTA, pipeline
    stats).

    Review `Tables/pipeline_stats.csv` for unexpected read loss before
    committing to the (potentially long) `assign_taxonomy()` step.

    The `asv.denovo_method` config key selects DADA2's chimera detection mode.
    Allowed values are `"consensus"`, `"pooled"`, and `"per-sample"`; any
    other value raises before R is invoked.

    Requires: `Checkpoints/ckpt_filter.RData`, `Checkpoints/ckpt_length.RData`
    Saves: `Checkpoints/ckpt_chimera.RData`
    """
    function chimera_removal(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        lbl     = ctx.run_label
        emit    = _emitter(progress, lbl)
        @info "[$(lbl)] DADA2: Chimera removal starting"
        seq_prefix   = get(ctx.cfg["output"], "seq_table_prefix", "seqtab_nochim")
        fasta_prefix = get(ctx.cfg["output"], "fasta_prefix", "asvs")
        tables_dir   = ctx.dirs["Tables"]

        isfile(ctx.ckpts["filter"]) ||
            error("Filter checkpoint not found. Run filter_trim() first.")
        isfile(ctx.ckpts["length"]) ||
            error("Length filter checkpoint not found. Run filter_length() first.")
        filter_ckpt = ctx.ckpts["filter"]
        length_ckpt = ctx.ckpts["length"]

        chimera_ckpt = ctx.ckpts["chimera"]
        hash_file    = joinpath(ctx.dirs["Checkpoints"], "chimera_removal.hash")
        if isfile(chimera_ckpt) &&
           !_section_stale(config_path, stage_sections(:dada2_chimera_removal), hash_file) &&
           mtime(chimera_ckpt) > mtime(filter_ckpt) &&
           mtime(chimera_ckpt) > mtime(length_ckpt)
            @info "[$(lbl)] DADA2: Skipping chimera removal - checkpoint up to date"
            return nothing
        end
        snap = _begin_section(config_path, stage_sections(:dada2_chimera_removal), hash_file)

        log_path = joinpath(ctx.dirs["Logs"], "chimera_removal.log")
        open(log_path, "w") do io; println(io, "=== chimera_removal ===\nconfig: $config_path") end

        # Guard at the Julia/R boundary: an unrecognised method silently
        # degrades downstream chimera calls and produces an opaque R error
        # rather than a clear Julia validation failure. Checked before either
        # path so a remote run fails here too, not eight transfers later.
        denovo_method = ctx.cfg["asv"]["denovo_method"]
        denovo_method isa AbstractString && denovo_method in Validation.DENOVO_METHODS ||
            error("asv.denovo_method must be one of $(join(Validation.DENOVO_METHODS, ", ")) " *
                  "(got: $(repr(denovo_method)))")

        target = _remote_target(ctx.full_cfg, "chimera_removal")
        if isnothing(target)
            _chimera_removal_local(ctx, emit, log_path,
                                   _r_threads(ctx.full_cfg; stage="chimera_removal"),
                                   denovo_method)
        else
            emit("Removing chimeras (remote: $(target.host))")
            _chimera_removal_remote(ctx, emit, target, log_path, denovo_method)
        end
        _write_section_hash(config_path, stage_sections(:dada2_chimera_removal), hash_file; snapshot=snap)
        emit("Written: $(joinpath(tables_dir, "pipeline_stats.csv"))")
        emit("Written: $(joinpath(tables_dir, seq_prefix * ".csv"))")
        emit("Written: $(joinpath(tables_dir, fasta_prefix * ".fasta"))")
        emit("Checkpoint: $(ctx.ckpts["chimera"])")
        emit("Log: $log_path")
        nothing
    end
