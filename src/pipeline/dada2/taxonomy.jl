# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

    ## Taxonomy Assignment
    """
        assign_taxonomy(config_path; progress)

    **Stage 7** - Assign taxonomy to ASVs and write the combined output table.

    Requires: `Checkpoints/ckpt_chimera.RData`
    Saves: `Checkpoints/checkpoint.RData`
    """
    # Assign taxonomy on the bioserver. The chimera-free checkpoint goes up and
    # the taxonomy tables come back; the reference database goes with it only
    # when databases.yml does not already name a copy on the server.
    function _assign_taxonomy_remote(emit, target, chimera_ckpt, db_path, db_remote_path,
                                      tables_dir, checkpoint, taxa_prefix,
                                      min_boot, tax_levels, seed, verbose, log_path;
                                      db_sha256=nothing)
        !isnothing(db_remote_path) && !startswith(string(db_remote_path), "/") &&
            error("databases.yml dada2.remote_path must be an absolute path on the server " *
                  "(got: '$db_remote_path'). Do not include the hostname.")

        inputs = Pair{String,String}[chimera_ckpt => "ckpt_chimera.RData"]
        # A database already on the server costs nothing to reuse; one that is
        # not has to travel, and at ~100MB compressed it dominates this stage's
        # transfer. Set databases.yml dada2.remote_path to skip it.
        db_arg = if !isnothing(db_remote_path)
            emit("  Using remote database: $db_remote_path")
            string(db_remote_path)
        elseif !isnothing(db_path)
            emit("  Database will be transferred ($(basename(db_path)))")
            push!(inputs, db_path => basename(db_path))
            "REMOTE_STAGING/" * basename(db_path)
        else
            error("No taxonomy database: set databases.yml dada2.remote_path or provide a local database")
        end

        _run_remote_stage(emit, target, "taxonomy_remote.r",
            [
                "db"          => db_arg,
                "prefix"      => taxa_prefix,
                "multithread" => _mt_str(target.threads),
                "min_boot"    => string(min_boot),
                "levels"      => join(tax_levels, ","),
                "seed"        => string(seed),
                "verbose"     => string(verbose),
            ];
            inputs,
            # Whether the database was uploaded or was already on the server,
            # test its gzip stream there before assignTaxonomy reads it.
            verify_remote_gzip = endswith(lowercase(db_arg), ".gz") ? [db_arg] : String[],
            # A copy already on the server is checked against the configured release checksum.
            verify_remote_sha256 = (!isnothing(db_remote_path) && !isnothing(db_sha256)) ?
                Dict(db_arg => string(db_sha256)) : Dict{String,String}(),
            outputs = [
                "Tables/$taxa_prefix.csv"            => joinpath(tables_dir, "$taxa_prefix.csv"),
                "Tables/$(taxa_prefix)_bootstraps.csv" => joinpath(tables_dir, "$(taxa_prefix)_bootstraps.csv"),
                "Tables/$(taxa_prefix)_combined.csv"   => joinpath(tables_dir, "$(taxa_prefix)_combined.csv"),
                "checkpoint.RData"                   => checkpoint,
            ],
            log_path, stage = "assign_taxonomy")
    end

    function assign_taxonomy(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing, taxonomy_db=nothing, skip_classification=false)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        lbl     = ctx.run_label
        emit    = _emitter(progress, lbl)
        @info "[$(lbl)] DADA2: Assign taxonomy starting"
        verbose = ctx.verbose

        isfile(ctx.ckpts["chimera"]) ||
            error("Chimera checkpoint not found. Run chimera_removal() first.")

        chimera_ckpt = ctx.ckpts["chimera"]
        checkpoint   = joinpath(ctx.dirs["Checkpoints"], "checkpoint.RData")
        hash_file    = joinpath(ctx.dirs["Checkpoints"], "assign_taxonomy.hash")
        if isfile(checkpoint) &&
           !_section_stale(config_path, stage_sections(:dada2_assign_taxonomy), hash_file) &&
           mtime(checkpoint) > mtime(chimera_ckpt)
            @info "[$(lbl)] DADA2: Skipping assign taxonomy - checkpoint up to date"
            return nothing
        end
        snap = _begin_section(config_path, stage_sections(:dada2_assign_taxonomy), hash_file)

        # Free R data objects from prior stages to reduce memory before taxonomy
        R"""
        .data_objs <- setdiff(ls(), lsf.str())
        if (length(.data_objs) > 0L) rm(list = .data_objs)
        rm(.data_objs)
        gc()
        """
        _source_r_functions(ctx)

        seq_prefix    = get(ctx.cfg["output"], "seq_table_prefix", "seqtab_nochim")
        fasta_prefix  = get(ctx.cfg["output"], "fasta_prefix", "asvs")
        taxa_prefix   = get(ctx.cfg["output"], "taxa_prefix", "taxonomy")
        combined_file = "tax_counts.csv"
        asv_file      = "asv_counts.csv"
        tables_dir    = ctx.dirs["Tables"]
        target        = _remote_target(ctx.full_cfg, "assign_taxonomy")
        multithread   = isnothing(target) ?
                        _r_threads(ctx.full_cfg; stage="assign_taxonomy") : target.threads
        # Only a Bool or a positive Integer is valid; YAML can yield `"4"`.
        (multithread isa Bool) ||
            (multithread isa Integer && multithread >= 1) ||
            error("r_threads must be a Bool or positive integer (got: $(repr(multithread)))")
        min_boot      = get(ctx.cfg["taxonomy"], "min_boot", 0)
        db_key    = string(ctx.cfg["taxonomy"]["database"])
        dbs_path  = joinpath(@__DIR__, "..", "..", "..", "config", "databases.yml")
        dbs_cfg   = get(YAML.load_file(dbs_path), "databases", Dict())
        tax_levels = String[string(l) for l in get(get(dbs_cfg, db_key, Dict()), "levels", String[])]
        isempty(tax_levels) && error("No levels defined for database '$db_key' in $dbs_path")

        # Priority: remote_path > taxonomy_db arg > _resolve_taxonomy_db fallback
        db_remote_path = get(get(get(dbs_cfg, db_key, Dict()), "dada2", Dict()), "remote_path", nothing)
        if !isnothing(db_remote_path)
            db_remote_path = string(db_remote_path)
        end
        db_sha256 = get(get(get(dbs_cfg, db_key, Dict()), "dada2", Dict()), "sha256", nothing)
        db_path = isnothing(taxonomy_db) ? _resolve_taxonomy_db(ctx.cfg, emit) : taxonomy_db

        R"load($chimera_ckpt)"
        has_data = rcopy(R"isTRUE(sum(seq_table_nochim, na.rm=TRUE) > 0)")

        log_path = joinpath(ctx.dirs["Logs"], "assign_taxonomy.log")
        open(log_path, "w") do io; println(io, "=== assign_taxonomy ===\nconfig: $config_path") end

        if skip_classification
            @info "[$(lbl)] DADA2: Taxonomy classification skipped - writing count tables without assignments"
            R"""
            taxa_df <- data.frame(SeqName = index$SeqName, Sequence = index$Sequence,
                                  stringsAsFactors = FALSE)
            """
            R"write_combined_table(taxa_df, index, seq_table_nochim, $tables_dir, $combined_file, $asv_file)"
            R"save(seq_table_nochim, index, taxa_df, file=$checkpoint)"
            for suffix in (taxa_prefix * ".csv", taxa_prefix * "_bootstraps.csv",
                           taxa_prefix * "_combined.csv")
                touch(joinpath(tables_dir, suffix))
            end
        elseif !has_data
            @info "[$(lbl)] DADA2: Taxonomy assignment skipped - no ASVs in seq_table_nochim"
            R"taxa_df <- data.frame()"
            R"save(seq_table_nochim, index, taxa_df, file=$checkpoint)"
            for suffix in (taxa_prefix * ".csv", taxa_prefix * "_bootstraps.csv",
                           taxa_prefix * "_combined.csv", combined_file, asv_file)
                touch(joinpath(tables_dir, suffix))
            end
        elseif !isnothing(target)
            emit("Assigning taxonomy (remote: $(target.host))")
            _assign_taxonomy_remote(emit, target, chimera_ckpt, db_path, db_remote_path,
                                    tables_dir, checkpoint, taxa_prefix,
                                    min_boot, tax_levels, ctx.seed, verbose, log_path; db_sha256)
            R"load($checkpoint)"
            R"write_combined_table(taxa_df, index, seq_table_nochim, $tables_dir, $combined_file, $asv_file)"
            emit("Log: $log_path")
        else
            R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
            try
                emit("Assigning taxonomy")
                # assignTaxonomy bootstraps by randomly subsampling kmers, so
                # without a seed both the bootstrap values and the winning genus
                # vary between runs on identical input.
                R"set.seed($(ctx.seed), kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection')"
                _r_run_logged("taxa_result <- run_assign_taxonomy(seq_table_nochim, " *
                              "$(_r_lit(db_path)), list(multithread=$(_r_lit(multithread)), " *
                              "min_boot=$(_r_lit(min_boot)), levels=$(_r_lit(tax_levels))), " *
                              "$(_r_lit(verbose)))")
                R"""
                taxa_df <- write_taxa_table(taxa_result$tax, taxa_result$boot, index,
                                            $tables_dir, $taxa_prefix)
                """
                R"gc()"
                R"save(seq_table_nochim, index, taxa_df, file=$checkpoint)"
            finally
                R"tryCatch({ sink(type='message'); sink(); close(con) }, error = function(e) NULL)"
            end
            R"write_combined_table(taxa_df, index, seq_table_nochim, $tables_dir, $combined_file, $asv_file)"
            emit("Log: $log_path")
        end

        _write_section_hash(config_path, stage_sections(:dada2_assign_taxonomy), hash_file; snapshot=snap)
        emit("Checkpoint: $checkpoint")

        emit("Pipeline complete. Outputs:")
        emit("  $(joinpath(tables_dir, seq_prefix * ".csv"))")
        emit("  $(joinpath(tables_dir, fasta_prefix * ".fasta"))")
        emit("  $(joinpath(tables_dir, fasta_prefix * ".csv"))")
        emit("  $(joinpath(tables_dir, taxa_prefix * ".csv"))")
        emit("  $(joinpath(tables_dir, taxa_prefix * "_bootstraps.csv"))")
        emit("  $(joinpath(tables_dir, taxa_prefix * "_combined.csv"))")
        emit("  $(joinpath(tables_dir, combined_file))")
        emit("  $(joinpath(tables_dir, asv_file))")
        emit("  $(joinpath(tables_dir, "pipeline_stats.csv"))")
        emit("  $checkpoint")
        nothing
    end
