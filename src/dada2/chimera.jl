# Web UI: Denoising page
# Stages: chimera_removal

    # Chimera removal
    """
        chimera_removal(config_path; progress)

    **Stage 6** - Remove chimeric sequences, drop the single-sample duplicate
    row if present, and write core output files (seq table, FASTA, pipeline
    stats).

    Review `Tables/pipeline_stats.csv` for unexpected read loss at any stage
    before committing to the (potentially long) `assign_taxonomy()` step.

    Requires: `Checkpoints/ckpt_filter.RData`, `Checkpoints/ckpt_length.RData`
    Saves: `Checkpoints/ckpt_chimera.RData`
    """
    function chimera_removal(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        verbose = ctx.verbose
        mode    = ctx.mode

        isfile(ctx.ckpts["filter"]) ||
            error("Filter checkpoint not found. Run filter_trim() first.")
        isfile(ctx.ckpts["length"]) ||
            error("Length filter checkpoint not found. Run filter_length() first.")
        filter_ckpt = ctx.ckpts["filter"]
        length_ckpt = ctx.ckpts["length"]
        R"load($filter_ckpt)"
        R"load($length_ckpt)"

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
