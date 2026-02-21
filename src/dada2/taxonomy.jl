# Web UI: Taxonomy page
# Stages: assign_taxonomy

    # Taxonomy assignment
    """
        assign_taxonomy(config_path; progress)

    **Stage 7** - Assign taxonomy to ASVs and write the combined output table.

    Requires: `Checkpoints/ckpt_chimera.RData`
    Saves: `Checkpoints/checkpoint.RData`
    """
    # DISCLAIMER: Remote taxonomy execution connects to a user-configured server via SSH,
    # transfers files via SCP, executes Rscript, and deletes the staging directory on
    # completion. The user is solely responsible for ensuring they have authorisation to
    # use the configured remote host, that the staging_dir is a safe path to write to
    # and delete from, and that the remote server has sufficient resources. The authors
    # of this software accept no liability for unintended data loss, unauthorised access,
    # or any other consequences arising from misconfiguration or misuse of this feature.
    function _assign_taxonomy_remote(emit, chimera_ckpt, db_path, tables_dir,
                                      checkpoint, taxa_prefix, multithread,
                                      min_boot, tax_levels, verbose, remote_cfg,
                                      log_path)
        host      = remote_cfg["host"]
        rscript   = get(remote_cfg, "rscript", "Rscript")
        base_dir  = get(remote_cfg, "staging_dir", nothing)
        remote_db = get(remote_cfg, "db_path", nothing)

        isnothing(base_dir) &&
            error("taxonomy.remote.staging_dir must be set explicitly in config")
        !isnothing(remote_db) && !startswith(string(remote_db), "/") &&
            error("taxonomy.remote.db_path must be an absolute path on the server " *
                  "(got: '$remote_db'). Do not include the hostname.")

        # Append a unique run ID so cleanup only ever touches this specific run's dir.
        run_id      = string(floor(Int, time()))
        staging_dir = "$base_dir/run_$run_id"

        scripts_dir   = @__DIR__
        functions_r   = joinpath(scripts_dir, "dada2_functions.r")
        remote_r      = joinpath(scripts_dir, "taxonomy_remote.r")
        remote_tables = "$staging_dir/Tables"
        remote_ckpt   = "$staging_dir/checkpoint.RData"

        # Use a ControlMaster socket so the password is entered once and all
        # subsequent ssh/scp calls reuse the existing connection silently.
        # ConnectTimeout: fail fast instead of hanging on unreachable hosts.
        # NumberOfPasswordPrompts=1: fail immediately on wrong password.
        # ServerAlive*: detect stale connections during long Rscript runs.
        ctl      = "/tmp/ssh_mux_$run_id"
        ssh_opts = `-o ControlMaster=auto -o ControlPath=$ctl -o ControlPersist=yes -o ConnectTimeout=15 -o NumberOfPasswordPrompts=1 -o ServerAliveInterval=30 -o ServerAliveCountMax=3`
        ssh = (args...) -> `ssh $ssh_opts $args`
        scp = (args...) -> `scp $ssh_opts $args`

        emit("  Connecting to $host (enter SSH password if prompted)...")
        emit("  Setting up staging directory on $host")
        run(ssh(host, "mkdir -p $remote_tables"))

        emit("  Transferring files to $host")
        run(scp(chimera_ckpt, "$host:$staging_dir/ckpt_chimera.RData"))
        run(scp(functions_r,  "$host:$staging_dir/dada2_functions.r"))
        run(scp(remote_r,     "$host:$staging_dir/taxonomy_remote.r"))

        remote_db_path = if !isnothing(remote_db)
            emit("  Using remote database: $remote_db")
            string(remote_db)
        else
            db_basename = basename(db_path)
            emit("  Transferring database ($db_basename) to $host")
            run(scp(db_path, "$host:$staging_dir/$db_basename"))
            "$staging_dir/$db_basename"
        end

        levels_str  = join(tax_levels, ",")
        verbose_str = verbose ? "true" : "false"
        # Julia Bool true/false -> R TRUE/FALSE; integers pass through as-is.
        mt_str      = multithread isa Bool ? (multithread ? "TRUE" : "FALSE") : string(multithread)
        remote_cmd  = "$rscript $staging_dir/taxonomy_remote.r " *
                      "functions=$staging_dir/dada2_functions.r " *
                      "ckpt=$staging_dir/ckpt_chimera.RData " *
                      "db=$remote_db_path " *
                      "tables=$remote_tables " *
                      "save=$remote_ckpt " *
                      "prefix=$taxa_prefix " *
                      "multithread=$mt_str " *
                      "min_boot=$min_boot " *
                      "levels=$levels_str " *
                      "verbose=$verbose_str"

        emit("  Running Rscript on $host:$staging_dir")
        open(log_path, "a") do io
            run(pipeline(ssh(host, remote_cmd); stdout=io, stderr=io))
        end

        emit("  Retrieving results from $host")
        run(scp("$host:$remote_tables/$taxa_prefix.csv",              "$tables_dir/"))
        run(scp("$host:$remote_tables/$(taxa_prefix)_bootstraps.csv", "$tables_dir/"))
        run(scp("$host:$remote_tables/$(taxa_prefix)_combined.csv",   "$tables_dir/"))
        run(scp("$host:$remote_ckpt",                                 checkpoint))

        # Safety check before cleanup: staging_dir must be at least 3 components
        # deep to guard against dangerous paths.
        parts = filter(!isempty, split(staging_dir, '/'))
        if length(parts) >= 3
            emit("  Cleaning up $host:$staging_dir")
            run(ssh(host, "rm -rf $staging_dir"))
        else
            @warn "Skipping remote cleanup: staging path '$staging_dir' looks too shallow to delete safely"
        end

        run(`ssh -o ControlPath=$ctl -O exit $host`)
    end

    function assign_taxonomy(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing, taxonomy_db=nothing)
        emit    = _emitter(progress)
        ctx     = _pipeline_context(config_path; input_dir, workspace_root)
        verbose = ctx.verbose

        isfile(ctx.ckpts["chimera"]) ||
            error("Chimera checkpoint not found. Run chimera_removal() first.")

        chimera_ckpt = ctx.ckpts["chimera"]
        checkpoint   = joinpath(ctx.dirs["Checkpoints"], "checkpoint.RData")
        hash_file    = joinpath(ctx.dirs["Checkpoints"], "config.hash")
        if isfile(checkpoint) &&
           !_section_stale(config_path, "dada2", hash_file) &&
           mtime(checkpoint) > mtime(chimera_ckpt)
            @info "Skipping assign_taxonomy: checkpoint up to date"
            return nothing
        end

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
        _source_r_functions()

        seq_prefix    = get(ctx.cfg["output"], "seq_table_prefix", "seqtab_nochim")
        fasta_prefix  = get(ctx.cfg["output"], "fasta_prefix", "asvs")
        taxa_prefix   = get(ctx.cfg["output"], "taxa_prefix", "taxonomy")
        combined_file = "tax_counts.csv"
        asv_file      = "asv_counts.csv"
        tables_dir    = ctx.dirs["Tables"]
        multithread   = get(ctx.cfg["taxonomy"], "multithread", 4)
        min_boot      = get(ctx.cfg["taxonomy"], "min_boot", 0)
        tax_levels    = ctx.cfg["taxonomy"]["levels"]

        remote_cfg = get(get(ctx.cfg, "taxonomy", Dict()), "remote", nothing)
        use_remote = !isnothing(remote_cfg) &&
                     !isnothing(get(remote_cfg, "host", nothing))

        # Skip local DB resolution when remote has its own db_path configured.
        skip_local_db = use_remote && !isnothing(get(remote_cfg, "db_path", nothing))
        db_path = skip_local_db ? nothing :
                  isnothing(taxonomy_db) ? _resolve_taxonomy_db(ctx.cfg, emit) : taxonomy_db

        R"load($chimera_ckpt)"
        has_data = rcopy(R"isTRUE(sum(seq_table_nochim, na.rm=TRUE) > 0)")

        log_path = joinpath(ctx.dirs["Logs"], "assign_taxonomy.log")
        open(log_path, "w") do io; println(io, "=== assign_taxonomy ===\nconfig: $config_path") end

        if !has_data
            @info "Taxonomy assignment skipped: no ASVs in seq_table_nochim"
            R"taxa_df <- data.frame()"
            R"save(seq_table_nochim, index, taxa_df, file=$checkpoint)"
            for suffix in (taxa_prefix * ".csv", taxa_prefix * "_bootstraps.csv",
                           taxa_prefix * "_combined.csv", combined_file, asv_file)
                touch(joinpath(tables_dir, suffix))
            end
        elseif use_remote
            emit("Assigning taxonomy (remote: $(remote_cfg["host"]))")
            _assign_taxonomy_remote(emit, chimera_ckpt, db_path, tables_dir,
                                    checkpoint, taxa_prefix, multithread,
                                    min_boot, tax_levels, verbose, remote_cfg,
                                    log_path)
            R"load($checkpoint)"
            R"write_combined_table(taxa_df, index, seq_table_nochim, $tables_dir, $combined_file, $asv_file)"
            emit("Log: $log_path")
        else
            R"con <- file($log_path, open='at'); sink(con); sink(con, type='message')"
            try
                emit("Assigning taxonomy")
                R"""
                taxa_result <- run_assign_taxonomy(
                    seq_table_nochim, $db_path,
                    list(multithread=$multithread, min_boot=$min_boot, levels=$tax_levels),
                    $verbose)
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

        _write_section_hash(config_path, "dada2", hash_file)
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
