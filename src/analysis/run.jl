    ## Level 1 - Per-run analysis
    """
        analyse_run(project, merged, asvs, db_meta; plot_lock=nothing)

    Produce per-run analysis outputs under `{project.dir}/analysis/`.
    When `plot_lock` is provided, CairoMakie calls are serialized behind
    it so multiple runs can prepare data in parallel.
    """
    function analyse_run(project::ProjectCtx, merged::MergedTables,
                         asvs::ASVResult, db_meta::DatabaseMeta;
                         plot_lock::Union{Nothing,ReentrantLock}=nothing)
        analysis_dir = joinpath(project.dir, "analysis")
        figures_dir  = joinpath(analysis_dir, "Figures")
        merged_csv   = merged.tables["merged"]

        filter_keys  = sort([k for k in keys(merged.tables) if k != "merged"])
        has_filters  = !isempty(filter_keys)
        all_source_keys = vcat(filter_keys, ["merged"])

        # Output paths.
        summary_path  = joinpath(analysis_dir, "pipeline_summary.csv")
        stages_pdf    = joinpath(figures_dir,  "pipeline_stages.pdf")
        alpha_pdf     = joinpath(figures_dir,  "alpha_diversity.pdf")
        filter_pdf    = joinpath(figures_dir,  "filter_composition.pdf")

        # Read analysis config from the cascade (needed for skip guard and charts).
        config_path  = write_run_config(project)
        analysis_cfg = get(YAML.load_file(config_path), "analysis", Dict())
        taxa_cfg     = get(analysis_cfg, "taxa_bar", Dict())
        top_n        = get(taxa_cfg, "top_n", 15)
        taxa_ranks   = _taxa_ranks(db_meta.levels, analysis_cfg)
        report_ranks = _report_ranks(db_meta.levels, analysis_cfg)

        required_outputs = String[summary_path, stages_pdf, alpha_pdf]
        has_filters && push!(required_outputs, filter_pdf)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(merged, method)
            for src in msrc_keys
                sd = _method_source_dirname(src)
                for (rankdir, _) in taxa_ranks
                    push!(required_outputs, joinpath(figures_dir, method, sd, rankdir, "taxa_bar.pdf"))
                    push!(required_outputs, joinpath(figures_dir, method, sd, rankdir, "taxa_bar_absolute.pdf"))
                end
            end
            for src in msrc_keys
                push!(required_outputs, joinpath(analysis_dir, "analysis_report_$(method)_$(_method_source_dirname(src)).txt"))
            end
        end

        # skip guard
        merged_mtime = isfile(merged_csv) ? mtime(merged_csv) : time()
        if all(isfile, required_outputs) &&
           all(f -> mtime(f) > merged_mtime, required_outputs)
            @info "Skipping analyse_run: outputs up to date in $analysis_dir"
            return
        end

        mkpath(figures_dir)

        # Helper: serialize CairoMakie calls if plot_lock is provided.
        _plot(f) = isnothing(plot_lock) ? f() : lock(f, plot_lock)

        # CSV cache
        cache = _CSVCache()

        src_key   = _source_key(merged)
        src_label = _source_label(src_key)
        subtitle  = "Source: $src_label"

        # pipeline summary
        stats_df = _pipeline_summary(project, merged, db_meta)

        # pipeline stages plot
        if !isnothing(stats_df)
            _plot() do
                pipeline_stats_plot(stats_df, stages_pdf; subtitle)
            end
            @info "Written: $stages_pdf"
            log_written(project, stages_pdf)
        end

        # taxa bar charts for all sources at multiple ranks (dual method)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(merged, method)
            for src in msrc_keys
                # For dada2 method with "merged" key, use merged.csv;
                # for _dada2 suffixed keys, use the corresponding CSV.
                csv_key = src
                src_csv_path = get(merged.tables, csv_key, "")
                if isempty(src_csv_path) || !isfile(src_csv_path) || filesize(src_csv_path) == 0
                    # For "merged" key under dada2 method, fall back to merged.csv
                    src == "merged" || continue
                    src_csv_path = merged.tables["merged"]
                    (!isfile(src_csv_path) || filesize(src_csv_path) == 0) && continue
                end
                raw_df, src_scols = _cached_read(cache, src_csv_path, db_meta)
                isempty(src_scols) && continue
                # Skip dada2 method if no dada2 taxonomy columns
                method == "dada2" && !_has_dada2(raw_df, db_meta.levels) && continue
                view_df = _method_df(raw_df, method, db_meta.levels)
                _total_seqs(view_df, src_scols) == 0 && continue
                sd = _method_source_dirname(src)
                src_sub = "Source: $sd ($method)"
                method_fig_dir = joinpath(figures_dir, method)
                _plot() do
                    _generate_taxa_charts(view_df, src_scols, method_fig_dir,
                                           top_n, src_sub, src == "merged" ? "merged" :
                                           (endswith(src, "_dada2") ? src[1:end-6] : src);
                                           ranks=taxa_ranks, rank_order=db_meta.levels)
                end
            end
        end

        # filter composition (priority-based)
        if has_filters
            filter_totals, comp_scols = _priority_filter_composition(merged, cache, db_meta)
            if !isnothing(filter_totals) && !isempty(filter_totals)
                _plot() do
                    filter_composition_plot(filter_totals, comp_scols, filter_pdf;
                                            subtitle, colour_overrides=merged.filter_colours)
                end
                @info "Written: $filter_pdf"
                log_written(project, filter_pdf)
            end
        end

        # alpha diversity (from primary source)
        primary_csv = merged.tables[src_key]
        if isfile(primary_csv) && filesize(primary_csv) > 0
            prim_df, prim_scols = _cached_read(cache, primary_csv, db_meta)
            if !isempty(prim_scols)
                alpha = _compute_alpha(prim_df, prim_scols)
                _plot() do
                    alpha_diversity_plot(alpha, alpha_pdf; subtitle)
                end
                @info "Written: $alpha_pdf"
                log_written(project, alpha_pdf)
            end
        end

        # analysis reports (dual method, one per source per method)
        run_name = basename(project.dir)
        merged_df_full, _ = _cached_read(cache, merged_csv, db_meta)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(merged, method)
            for src in msrc_keys
                csv_key = src
                src_csv_path = get(merged.tables, csv_key, "")
                if isempty(src_csv_path) || !isfile(src_csv_path) || filesize(src_csv_path) == 0
                    src == "merged" || continue
                    src_csv_path = merged.tables["merged"]
                    (!isfile(src_csv_path) || filesize(src_csv_path) == 0) && continue
                end
                raw_df, src_scols = _cached_read(cache, src_csv_path, db_meta)
                isempty(src_scols) && continue
                method == "dada2" && !_has_dada2(raw_df, db_meta.levels) && continue
                view_df = _method_df(raw_df, method, db_meta.levels)
                _total_seqs(view_df, src_scols) == 0 && continue
                sd = _method_source_dirname(src)
                report_path = joinpath(analysis_dir, "analysis_report_$(method)_$(sd).txt")
                clean_src = src == "merged" ? "merged" : (endswith(src, "_dada2") ? src[1:end-6] : src)
                _generate_report(view_df, src_scols, stats_df, merged_df_full,
                                 clean_src, report_path, run_name;
                                 stats_key=src, report_ranks=report_ranks)
                log_written(project, report_path)
            end
        end
    end
