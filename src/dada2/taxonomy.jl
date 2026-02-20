# Web UI: Taxonomy page
# Stages: assign_taxonomy

    # Taxonomy assignment
    """
        assign_taxonomy(config_path; progress)

    **Stage 7** - Assign taxonomy to ASVs and write the combined output table.

    This is typically the longest step. Set `taxonomy.skip = true` in config to
    skip assignment and output sequence/count data only.

    Requires: `Checkpoints/ckpt_chimera.RData`
    Saves: `Checkpoints/checkpoint.RData`
    """
    function assign_taxonomy(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing, taxonomy_db=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        verbose = ctx.verbose

        isfile(ctx.ckpts["chimera"]) ||
            error("Chimera checkpoint not found. Run chimera_removal() first.")

        # Drop all data objects accumulated from prior stages before the
        # memory-intensive taxonomy assignment runs. Named globals in R's
        # environment are reachable and gc() won't collect them; rm() them
        # explicitly so they don't inflate the memory footprint.
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
        combined_file = get(ctx.cfg["output"], "combined_filename", "tax_counts.csv")
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
