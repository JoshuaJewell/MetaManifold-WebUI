    ## Level 3 - Study-level analysis
    function _analyse_study_level(study_dir::String,
                                  valid::Vector{Tuple{ProjectCtx, MergedTables, DatabaseMeta}},
                                  r_lock::ReentrantLock,
                                  default_db_meta::DatabaseMeta;
                                  plot_lock::Union{Nothing,ReentrantLock}=nothing)
        analysis_dir = joinpath(study_dir, "analysis")
        figures_dir  = joinpath(analysis_dir, "Figures")
        nmds_pdf     = joinpath(figures_dir, "nmds.pdf")
        alpha_pdf    = joinpath(figures_dir, "alpha_comparison.pdf")
        perm_txt     = joinpath(analysis_dir, "permanova.txt")

        src_key   = _source_key(valid[1][2])
        filter_keys = sort([k for k in keys(valid[1][2].tables) if k != "merged"])

        # Derive taxa/report ranks from the first project's db_meta.
        analysis_cfg = Dict()
        try
            config_path = write_run_config(valid[1][1])
            analysis_cfg = get(YAML.load_file(config_path), "analysis", Dict())
        catch; end
        taxa_ranks   = _taxa_ranks(default_db_meta.levels, analysis_cfg)
        report_ranks = _report_ranks(default_db_meta.levels, analysis_cfg)

        # skip guard
        newest_merged = maximum(
            mtime(m.tables["merged"])
            for (_, m) in valid if isfile(m.tables["merged"]);
            init=0.0
        )
        all_mergeds = [m for (_, m, _) in valid]

        # Initialise cache here so it is shared between the skip guard and execution.
        cache = _CSVCache()

        required_outputs = [nmds_pdf, alpha_pdf]
        for method in _TAX_METHODS
            msrc_keys = _all_method_source_keys(all_mergeds, method)
            for src in msrc_keys
                # Mirror the execution loop: only include the report path if at
                # least one run has a non-empty CSV for this source key.
                has_data = any(valid) do (proj, merged, dm)
                    csv = get(merged.tables, src, "")
                    if isempty(csv) && method == "dada2" && src == "merged"
                        csv = get(merged.tables, "merged", "")
                    end
                    (isempty(csv) || !isfile(csv) || filesize(csv) == 0) && return false
                    df, sc = _cached_read(cache, csv, dm)
                    isempty(sc) && return false
                    method == "dada2" && !_has_dada2(df, default_db_meta.levels) && return false
                    vdf = _method_df(df, method, default_db_meta.levels)
                    _total_seqs(vdf, sc) > 0
                end
                has_data || continue
                push!(required_outputs, joinpath(analysis_dir, "analysis_report_$(method)_$(_method_source_dirname(src)).txt"))
            end
        end
        if all(isfile, required_outputs) &&
           all(f -> mtime(f) > newest_merged, required_outputs)
            @info "Skipping study-level analysis: outputs up to date"
            return
        end

        reset_log(study_dir)
        mkpath(figures_dir)

        # Helper: serialize CairoMakie calls if plot_lock is provided.
        _plot(f) = isnothing(plot_lock) ? f() : lock(f, plot_lock)

        src_label = _source_label(src_key)
        subtitle  = "Source: $src_label"

        # group projects by parent dir
        groups = Dict{String, Vector{Tuple{ProjectCtx, MergedTables, DatabaseMeta}}}()
        for (proj, merged, dm) in valid
            gdir = dirname(proj.dir)
            push!(get!(groups, gdir, []), (proj, merged, dm))
        end

        # collect data across all groups
        run_dfs      = DataFrame[]
        run_scols    = Vector{String}[]
        all_scols    = String[]
        group_labels = String[]       # group name per sample
        run_labels   = String[]       # run name per sample (for shape)
        group_names  = String[]

        for (gdir, members) in groups
            gname = basename(gdir)
            push!(group_names, gname)
            for (proj, merged, dm) in members
                src_csv = merged.tables[_source_key(merged)]
                isfile(src_csv) || continue
                df, scols = _cached_read(cache, src_csv, dm)
                isempty(scols) && continue
                rname = basename(proj.dir)
                push!(run_dfs, df)
                push!(run_scols, scols)
                append!(all_scols, scols)
                append!(group_labels, fill(gname, length(scols)))
                append!(run_labels, fill(rname, length(scols)))
            end
        end

        isempty(run_dfs) && return

        # study NMDS
        mat, _, _ = _build_combined_counts(run_dfs, run_scols, default_db_meta.levels)
        nmds_stress = NaN
        if size(mat, 1) >= 3
            coords, stress = _run_nmds(mat, r_lock)
            if !any(isnan, coords)
                nmds_stress = stress
                _plot() do
                    nmds_plot(coords, all_scols, nmds_pdf;
                              colour_by=group_labels, shape_by=run_labels,
                              colour_label="Group", shape_label="Run",
                              stress=stress, subtitle)
                end
                @info "Written: $nmds_pdf"
                log_written(study_dir, nmds_pdf)
            else
                @warn "Study-level NMDS failed"
                pipeline_log(study_dir, "WARN: Study-level NMDS failed")
            end
        else
            @warn "Too few samples ($(size(mat, 1))) for study-level NMDS"
        end

        # study alpha boxplot by group
        all_alpha    = DataFrame[]
        alpha_labels = String[]
        for (gdir, members) in groups
            gname = basename(gdir)
            combined_alpha = DataFrame(sample=String[], richness=Int[],
                                       shannon=Float64[], simpson=Float64[])
            for (_, merged, dm) in members
                src_csv = merged.tables[_source_key(merged)]
                isfile(src_csv) || continue
                df, scols = _cached_read(cache, src_csv, dm)
                isempty(scols) && continue
                append!(combined_alpha, _compute_alpha(df, scols))
            end
            if nrow(combined_alpha) > 0
                push!(all_alpha, combined_alpha)
                push!(alpha_labels, gname)
            end
        end

        if !isempty(all_alpha)
            _plot() do
                alpha_boxplot(all_alpha, alpha_labels, alpha_pdf; subtitle)
            end
            @info "Written: $alpha_pdf"
            log_written(study_dir, alpha_pdf)
        end

        # PERMANOVA - always runs using group and run as covariates (derived from the
        # project structure).  Optionally, place a metadata.csv in the study directory
        # with a 'sample' column and extra covariate columns; those columns are merged
        # in when every sample name in all_scols is unique (i.e. no cross-database
        # duplication of the same biological sample name).
        n_perm_samples = size(mat, 1)
        if n_perm_samples in 3:5
            @warn "PERMANOVA: only $n_perm_samples samples — results have very low statistical power and p-values should not be interpreted."
        end
        if n_perm_samples >= 3
            # Base metadata: one row per entry in all_scols (may contain duplicates).
            base_meta = DataFrame(sample=all_scols, group=group_labels, run=run_labels)

            # Optionally merge user-supplied extra covariates.
            user_meta = load_metadata(study_dir, study_dir)
            if !isnothing(user_meta)
                sample_col = "sample" in names(user_meta) ? "sample" :
                             "Sample" in names(user_meta) ? "Sample" : nothing
                if isnothing(sample_col)
                    @warn "metadata.csv has no 'sample' column - using group/run only for PERMANOVA"
                elseif length(unique(all_scols)) < length(all_scols)
                    @warn "Duplicate sample names across runs - extra metadata.csv columns ignored; using group/run only"
                else
                    extra_cols = [c for c in names(user_meta) if c != sample_col]
                    if !isempty(extra_cols)
                        lookup = Dict(String(r[sample_col]) => r for r in eachrow(user_meta))
                        matched = [get(lookup, s, nothing) for s in all_scols]
                        if all(!isnothing, matched)
                            for col in extra_cols
                                base_meta[!, col] = [m[col] for m in matched]
                            end
                        else
                            @warn "Not all samples found in metadata.csv - extra columns ignored"
                        end
                    end
                end
            end

            perm_result = _run_permanova(mat, base_meta, r_lock)
            if !isnothing(perm_result)
                mkpath(dirname(perm_txt))
                write(perm_txt, perm_result)
                @info "Written: $perm_txt"
                log_written(study_dir, perm_txt)
            end
        end

        # dual analysis reports
        report_sections = Pair{String, String}[]

        # Group overview.
        buf = IOBuffer()
        println(buf, "Groups: ", join(group_names, ", "))
        println(buf, "Total samples: ", length(all_scols))
        push!(report_sections, "Overview" => String(take!(buf)))

        # NMDS info.
        if !isnan(nmds_stress)
            quality = nmds_stress < 0.2 ? "Good: below 0.2 threshold" : "Poor: above 0.2 threshold"
            push!(report_sections, "NMDS" =>
                  "Stress: $(round(nmds_stress; digits=3)) ($quality)")
        end

        # Alpha diversity by group.
        if !isempty(all_alpha)
            buf = IOBuffer()
            print(buf, rpad("Group", 20))
            println(buf, rpad("Samples", 10), rpad("Mean richness", 16),
                    rpad("Mean Shannon", 16), "Mean Simpson")
            for (adf, gname) in zip(all_alpha, alpha_labels)
                n = nrow(adf)
                mr = round(mean(adf.richness); digits=0)
                ms = round(mean(adf.shannon); digits=3)
                mp = round(mean(adf.simpson); digits=3)
                print(buf, rpad(gname, 20))
                if n > 1
                    sr = round(std(adf.richness); digits=0)
                    ss = round(std(adf.shannon); digits=2)
                    sp = round(std(adf.simpson); digits=2)
                    println(buf, rpad(string(n), 10),
                            rpad("$(Int(mr)) +/- $(Int(sr))", 16),
                            rpad("$ms +/- $ss", 16),
                            "$mp +/- $sp")
                else
                    println(buf, rpad(string(n), 10),
                            rpad(string(Int(mr)), 16),
                            rpad(string(ms), 16),
                            string(mp))
                end
            end
            push!(report_sections, "Alpha Diversity by Group" => String(take!(buf)))
        end

        # Collect (method, src) work items for parallel report writing.
        report_work = Tuple{String,String,Vector{DataFrame},Vector{String}}[]
        for method in _TAX_METHODS
            for src in _all_method_source_keys(all_mergeds, method)
                src_run_dfs   = DataFrame[]
                src_all_scols = String[]
                for (_, members_inner) in groups
                    for (_, merged, dm) in members_inner
                        csv = get(merged.tables, src, "")
                        if isempty(csv) && method == "dada2" && src == "merged"
                            csv = get(merged.tables, "merged", "")
                        end
                        (isempty(csv) || !isfile(csv)) && continue
                        df, sc = _cached_read(cache, csv, dm)
                        isempty(sc) && continue
                        push!(src_run_dfs, df)
                        append!(src_all_scols, sc)
                    end
                end
                isempty(src_run_dfs) && continue
                method == "dada2" && !_has_dada2(src_run_dfs[1], default_db_meta.levels) && continue
                view_dfs = [_method_df(df, method, default_db_meta.levels) for df in src_run_dfs]
                combined = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                _total_seqs(combined, src_all_scols) == 0 && continue
                push!(report_work, (method, src, view_dfs, src_all_scols))
            end
        end

        report_tasks = map(report_work) do (method, src, view_dfs, src_all_scols)
            Threads.@spawn begin
                combined = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                sd       = _method_source_dirname(src)
                sl       = sd == "unfiltered" ? "unfiltered" : sd
                sections = copy(report_sections)
                for (rank_name, rank_col) in report_ranks
                    section = _top_taxa_section(combined, src_all_scols, rank_name, rank_col; n=20)
                    !isempty(section) && push!(sections, "Top 20 $(rank_name) ($(sl))" => section)
                end
                rpath = joinpath(analysis_dir, "analysis_report_$(method)_$(sd).txt")
                _write_report(rpath, "Study Report\n  Source: $sl ($method)", sections)
                log_written(study_dir, rpath)
            end
        end

        foreach(fetch, report_tasks)
    end

    #  Public: analyse_study (entry point from main.jl)
    """
        analyse_study(projects, merged_results)

    Run group-level and study-level analysis after the per-run `@threads`
    loop has completed.

    Groups projects by `dirname(project.dir)` (= group directory).
    Calls `_analyse_group` for each group with ≥2 runs, then
    `_analyse_study_level` when there are ≥2 groups.
    """
    function analyse_study(projects::Vector{ProjectCtx},
                           merged_results::Vector{<:Union{MergedTables, Nothing}},
                           db_metas::Vector{DatabaseMeta};
                           plot_lock::Union{Nothing,ReentrantLock}=nothing)
        isempty(projects) && return

        r_lock    = ReentrantLock()
        study_dir = projects[1].study_dir

        # Pair projects with results and db_metas, filtering out failures.
        valid = Tuple{ProjectCtx, MergedTables, DatabaseMeta}[
            (projects[i], merged_results[i], db_metas[i])
            for i in eachindex(projects)
            if i <= length(merged_results) && !isnothing(merged_results[i])
        ]
        isempty(valid) && return

        # Group by parent directory.
        groups = Dict{String, Vector{Tuple{ProjectCtx, MergedTables, DatabaseMeta}}}()
        for (proj, merged, dm) in valid
            gdir = dirname(proj.dir)
            push!(get!(groups, gdir, []), (proj, merged, dm))
        end

        # Per-group analysis (skipped for single-run groups).
        # Each group uses the db_meta of its first member (all runs in a group share a DB).
        for (gdir, members) in groups
            if length(members) > 1
                group_db_meta = members[1][3]
                group_members = [(p, m) for (p, m, _) in members]
                _analyse_group(gdir, group_members, r_lock, group_db_meta; plot_lock)
            end
        end

        # Study-level analysis (only meaningful with ≥2 genuine groups).
        if length(groups) >= 2
            study_db_meta = valid[1][3]
            _analyse_study_level(study_dir, valid, r_lock, study_db_meta; plot_lock)
        end
    end
