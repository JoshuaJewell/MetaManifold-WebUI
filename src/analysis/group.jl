# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

    ## Level 2 - Per-group analysis
    function _analyse_group(group_dir::String,
                            members::Vector{Tuple{ProjectCtx, MergedTables}},
                            r_lock::ReentrantLock,
                            db_meta::DatabaseMeta;
                            plot_lock::Union{Nothing,ReentrantLock}=nothing)
        analysis_dir = joinpath(group_dir, "analysis")
        plots_dir    = joinpath(analysis_dir, "Plots")

        alpha_json  = joinpath(plots_dir, "alpha_comparison.json")
        nmds_json   = joinpath(plots_dir, "nmds.json")
        filter_json = joinpath(plots_dir, "filter_composition.json")

        # Determine source keys.
        src_key   = _source_key(members[1][2])
        src_label = _source_label(src_key)
        subtitle  = "Source: $src_label"

        # Derive taxa/report ranks from db_meta.
        analysis_cfg = Dict()
        try
            config_path = write_run_config(members[1][1])
            analysis_cfg = get(YAML.load_file(config_path), "analysis", Dict())
        catch; end
        taxa_ranks   = _taxa_ranks(db_meta.levels, analysis_cfg)
        report_ranks = _report_ranks(db_meta.levels, analysis_cfg)

        filter_keys     = sort([k for k in keys(members[1][2].tables) if k != "merged"])
        has_any_filters = any(pair -> length(pair[2].tables) > 1, members)

        # Skip guard.
        newest_merged = maximum(
            mtime(m.tables["merged"])
            for (_, m) in members if isfile(m.tables["merged"]);
            init=0.0
        )

        cache = _CSVCache()

        function _collect_source_data_sg(source::String)
            dfs    = DataFrame[]
            scols  = Vector{String}[]
            names_ = String[]
            labels = String[]
            all_s  = String[]
            for (proj, merged) in members
                haskey(merged.tables, source) || continue
                csv = merged.tables[source]
                isfile(csv) || continue
                df, sc = _cached_read(cache, csv, db_meta)
                isempty(sc) && continue
                rn = basename(proj.dir)
                push!(dfs, df); push!(scols, sc); push!(names_, rn)
                append!(all_s, sc); append!(labels, fill(rn, length(sc)))
            end
            return dfs, scols, names_, labels, all_s
        end

        required_outputs = String[alpha_json, nmds_json]
        has_any_filters && push!(required_outputs, filter_json)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(members[1][2], method)
            for src in msrc_keys
                s_dfs, _, _, _, s_all = _collect_source_data_sg(
                    endswith(src, "_dada2") ? src : src)
                if isempty(s_dfs) && method == "dada2" && src == "merged"
                    s_dfs, _, _, _, s_all = _collect_source_data_sg("merged")
                end
                isempty(s_dfs) && continue
                method == "dada2" && !_has_dada2(s_dfs[1], db_meta.levels) && continue
                view_dfs = [_method_df(df, method, db_meta.levels) for df in s_dfs]
                combined = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                _total_seqs(combined, s_all) == 0 && continue
                sd = _method_source_dirname(src)
                for (rankdir, _) in taxa_ranks
                    stem = "taxa_bar_$(method)_$(sd)_$(rankdir)"
                    push!(required_outputs, joinpath(plots_dir, "$(stem).json"))
                    push!(required_outputs, joinpath(plots_dir, "$(stem)_absolute.json"))
                    stem2 = "group_comparison_$(method)_$(sd)_$(rankdir)"
                    push!(required_outputs, joinpath(plots_dir, "$(stem2).json"))
                    push!(required_outputs, joinpath(plots_dir, "$(stem2)_absolute.json"))
                end
                push!(required_outputs, joinpath(analysis_dir, "analysis_report_$(method)_$(sd).txt"))
            end
        end

        if all(isfile, required_outputs) &&
           all(f -> mtime(f) > newest_merged, required_outputs)
            @info "Skipping group analysis: outputs up to date in $group_dir"
            return
        end

        reset_log(group_dir)
        mkpath(plots_dir)

        # Primary source data (for NMDS, alpha, filter composition).
        run_dfs, run_scols, run_names, run_labels, all_scols =
            _collect_source_data_sg(src_key)
        isempty(run_dfs) && return

        # Collect (method, src) work items for parallel dispatch.
        chart_work = Tuple{String,String,Vector{DataFrame},Vector{Vector{String}},Vector{String},Vector{String}}[]
        for method in _TAX_METHODS
            for src in _method_source_keys(members[1][2], method)
                s_dfs, s_scols, s_names, _, s_all = _collect_source_data_sg(
                    endswith(src, "_dada2") ? src : src)
                if isempty(s_dfs) && method == "dada2" && src == "merged"
                    s_dfs, s_scols, s_names, _, s_all = _collect_source_data_sg("merged")
                end
                isempty(s_dfs) && continue
                method == "dada2" && !_has_dada2(s_dfs[1], db_meta.levels) && continue
                view_dfs = [_method_df(df, method, db_meta.levels) for df in s_dfs]
                combined = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                _total_seqs(combined, s_all) == 0 && continue
                push!(chart_work, (method, src, view_dfs, s_scols, s_names, s_all))
            end
        end

        # Taxa bar + group comparison for all sources - parallel.
        chart_tasks = map(chart_work) do (method, src, view_dfs, s_scols, s_names, s_all)
            Threads.@spawn begin
                combined  = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                sd        = _method_source_dirname(src)
                src_sub   = "Source: $sd ($method)"
                clean_src = src == "merged" ? "merged" :
                            (endswith(src, "_dada2") ? src[1:end-6] : src)

                _generate_taxa_charts(combined, s_all, plots_dir, 15, src_sub,
                                      "$(method)_$(sd)";
                                      ranks=taxa_ranks, rank_order=db_meta.levels)

                for (rankdir, rank) in taxa_ranks
                    stem = "group_comparison_$(method)_$(sd)_$(rankdir)"
                    PipelinePlotsPlotly.group_comparison_chart(
                        view_dfs, s_scols, s_names,
                        joinpath(plots_dir, "$(stem).json");
                        top_n=15, rank, relative=true,
                        subtitle=src_sub, rank_order=db_meta.levels)
                    PipelinePlotsPlotly.group_comparison_chart(
                        view_dfs, s_scols, s_names,
                        joinpath(plots_dir, "$(stem)_absolute.json");
                        top_n=15, rank, relative=false,
                        subtitle=src_sub, rank_order=db_meta.levels)
                end
            end
        end

        # Group alpha diversity (boxplot by run).
        all_alpha    = DataFrame[]
        alpha_labels = String[]
        for (df, scols, rname) in zip(run_dfs, run_scols, run_names)
            alpha = _compute_alpha(df, scols)
            push!(all_alpha, alpha)
            push!(alpha_labels, rname)
        end
        if !isempty(all_alpha)
            PipelinePlotsPlotly.alpha_boxplot(all_alpha, alpha_labels, alpha_json; subtitle)
            @info "Written: $alpha_json"
            log_written(group_dir, alpha_json)
        end

        # Group NMDS.
        mat, _, _ = _build_combined_counts(run_dfs, run_scols, db_meta.levels)
        nmds_stress = NaN
        if size(mat, 1) >= 3
            coords, stress = _run_nmds(mat, r_lock)
            if !any(isnan, coords)
                nmds_cfg = Dict()
                try
                    config_path = write_run_config(members[1][1])
                    nmds_cfg = get(get(YAML.load_file(config_path), "analysis", Dict()),
                                   "nmds", Dict())
                catch; end
                max_stress = get(nmds_cfg, "max_stress", 0.2)
                if !isnan(stress) && stress > max_stress
                    @warn "NMDS stress $(round(stress; digits=3)) exceeds threshold $max_stress for group $(basename(group_dir))"
                    pipeline_log(group_dir, "WARN: NMDS stress $(round(stress; digits=3)) exceeds threshold $max_stress")
                end
                nmds_stress = stress
                PipelinePlotsPlotly.nmds_plot(coords, all_scols, nmds_json;
                                              colour_by=run_labels, stress=stress, subtitle)
                @info "Written: $nmds_json"
                log_written(group_dir, nmds_json)
            else
                @warn "NMDS failed for group $(basename(group_dir))"
                pipeline_log(group_dir, "WARN: NMDS failed for group $(basename(group_dir))")
            end
        else
            @warn "Too few samples ($(size(mat, 1))) for NMDS in group $(basename(group_dir))"
        end

        # Group filter composition.
        if has_any_filters
            group_filter_totals  = Dict{String, Vector{Float64}}()
            group_sample_names   = String[]
            group_filter_colours = Dict{String,String}()
            for (proj, merged) in members
                ft, scols_m = _priority_filter_composition(merged, cache, db_meta)
                isnothing(ft) && continue
                append!(group_sample_names, scols_m)
                merge!(group_filter_colours, merged.filter_colours)
                for (fk, vals) in ft
                    existing = get(group_filter_totals, fk, Float64[])
                    group_filter_totals[fk] = vcat(existing, vals)
                end
            end

            if haskey(group_filter_totals, "Unclassified") &&
               !any(>(0), group_filter_totals["Unclassified"])
                delete!(group_filter_totals, "Unclassified")
            end

            if !isempty(group_filter_totals) && !isempty(group_sample_names)
                PipelinePlotsPlotly.filter_composition_plot(
                    group_filter_totals, group_sample_names, filter_json;
                    subtitle, colour_overrides=group_filter_colours)
                @info "Written: $filter_json"
                log_written(group_dir, filter_json)
            end
        end

        # Analysis reports (dual method, per source) - parallel.
        shared_sections = Pair{String, String}[]

        buf = IOBuffer()
        print(buf, rpad("Run", 20))
        println(buf, rpad("Samples", 10), rpad("Mean richness", 16),
                rpad("Mean Shannon", 16), "Mean Simpson")
        for (adf, rname) in zip(all_alpha, alpha_labels)
            n  = nrow(adf)
            mr = round(mean(adf.richness); digits=0)
            ms = round(mean(adf.shannon); digits=3)
            mp = round(mean(adf.simpson); digits=3)
            if n > 1
                sr = round(std(adf.richness); digits=0)
                ss = round(std(adf.shannon); digits=2)
                sp = round(std(adf.simpson); digits=2)
                print(buf, rpad(rname, 20))
                println(buf, rpad(string(n), 10),
                        rpad("$(Int(mr)) +/- $(Int(sr))", 16),
                        rpad("$ms +/- $ss", 16),
                        "$mp +/- $sp")
            else
                print(buf, rpad(rname, 20))
                println(buf, rpad(string(n), 10),
                        rpad(string(Int(mr)), 16),
                        rpad(string(ms), 16),
                        string(mp))
            end
        end
        push!(shared_sections, "Per-Run Summary" => String(take!(buf)))

        if !isnan(nmds_stress)
            quality = nmds_stress < 0.2 ? "Good: below 0.2 threshold" : "Poor: above 0.2 threshold"
            push!(shared_sections, "NMDS" =>
                  "Stress: $(round(nmds_stress; digits=3)) ($quality)")
        end

        report_tasks = map(chart_work) do (method, src, view_dfs, _, _, s_all)
            Threads.@spawn begin
                combined = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                sd       = _method_source_dirname(src)
                sl       = sd == "unfiltered" ? "unfiltered" : sd
                sections = copy(shared_sections)
                for (rank_name, rank_col) in report_ranks
                    section = _top_taxa_section(combined, s_all, rank_name, rank_col; n=20)
                    !isempty(section) && push!(sections, "Top 20 $(rank_name) ($(sl))" => section)
                end
                rpath = joinpath(analysis_dir, "analysis_report_$(method)_$(sd).txt")
                _write_report(rpath,
                              "Comparison Report\n  Runs: $(join(run_names, ", "))\n  Source: $sl ($method)",
                              sections)
                log_written(group_dir, rpath)
            end
        end

        foreach(fetch, chart_tasks)
        foreach(fetch, report_tasks)
    end
