# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: on-demand analysis (alpha diversity, taxa bars, pipeline stats, NMDS, PERMANOVA)

using CSV, DataFrames

using MetaManifold.Analysis: alpha_chart, taxa_bar_chart, pipeline_stats_chart,
                             sample_columns, taxonomy_levels, filtered_counts,
                             aggregate_by_taxon, combined_counts_across_runs,
                             alpha_boxplot, nmds_chart, run_nmds, run_permanova, r_available

function _validate_run_request(study::String, run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                 "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    nothing
end

function _require_run_results(study::String, run::String; group::Union{String,Nothing}=nothing)
    err = _validate_run_request(study, run)
    !isnothing(err) && return nothing, err

    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return nothing, json_error(404, "no_results",
                                                 "No results for run '$run' - run the pipeline first")
    dir, nothing
end

## Per-run alpha diversity
@post "/api/v1/studies/{study}/runs/{run}/analysis/alpha" function(req,
                                                                    study::String,
                                                                    run::String)
    group = _req_group(req)
    dir, err = _require_run_results(study, run; group)
    !isnothing(err) && return err

    body = JSON3.read(String(req.body))
    table = get(body, :table, "merged")
    params = _body_filter_params(body)

    with_results_db(dir) do con
        scols = sample_columns(con, table)
        isempty(scols) && return json_error(400, "no_samples", "No sample columns found in table '$table'")

        columns = _duckdb_columns(con, table)
        where_clause, where_params = _build_where(params, columns)
        mat = filtered_counts(con, table, scols, where_clause, where_params)

        r = Int[]; sh = Float64[]; si = Float64[]
        for i in eachindex(scols)
            counts = round.(Int, mat[i, :])
            push!(r, richness(counts))
            push!(sh, shannon(counts))
            push!(si, simpson(counts))
        end

        fig = alpha_chart(scols, r, sh, si)
        HTTP.Response(200, ["Content-Type" => "application/json"],
                      body=JSON3.write(fig))
    end
end

## Per-run taxa bar chart
@post "/api/v1/studies/{study}/runs/{run}/analysis/taxa-bar" function(req,
                                                                       study::String,
                                                                       run::String)
    group = _req_group(req)
    dir, err = _require_run_results(study, run; group)
    !isnothing(err) && return err

    body = JSON3.read(String(req.body))
    table = get(body, :table, "merged")
    rank = get(body, :rank, nothing)
    top_n = get(body, :top_n, 15)
    relative = get(body, :relative, true)
    params = _body_filter_params(body)

    with_results_db(dir) do con
        scols = sample_columns(con, table)
        isempty(scols) && return json_error(400, "no_samples", "No sample columns found")

        columns = _duckdb_columns(con, table)
        where_clause, where_params = _build_where(params, columns)

        if isnothing(rank)
            levels = taxonomy_levels(con, table)
            isempty(levels) && return json_error(400, "no_taxonomy", "No taxonomy columns found")
            rank = last(levels)
        end

        agg = aggregate_by_taxon(con, table, scols, string(rank), where_clause, where_params)
        nrow(agg) == 0 && return json_error(400, "no_data", "No data after filtering")

        taxon_labels = String.(agg.taxon)
        counts = Matrix{Float64}(agg[:, scols])

        fig = taxa_bar_chart(taxon_labels, scols, counts; top_n=Int(top_n), relative=Bool(relative))
        HTTP.Response(200, ["Content-Type" => "application/json"],
                      body=JSON3.write(fig))
    end
end

## Pipeline stats
@get "/api/v1/studies/{study}/runs/{run}/analysis/pipeline-stats" function(req,
                                                                            study::String,
                                                                            run::String)
    err = _validate_run_request(study, run)
    !isnothing(err) && return err
    group = _req_group(req)
    run_dir = _run_project_dir(study, run; group)
    stats_csv = joinpath(run_dir, "dada2", "Tables", "pipeline_stats.csv")
    isfile(stats_csv) || return json_error(404, "no_stats",
                                                "No pipeline stats - run DADA2 first")

    stats_df = CSV.read(stats_csv, DataFrame)
    first_col = names(stats_df)[1]
    first_col != "sample" && rename!(stats_df, first_col => "sample")

    fig = pipeline_stats_chart(stats_df)
    isnothing(fig) && return json_error(400, "no_data", "Empty stats table")
    HTTP.Response(200, ["Content-Type" => "application/json"],
                  body=JSON3.write(fig))
end

## Taxonomy ranks discovery
@get "/api/v1/studies/{study}/runs/{run}/analysis/ranks" function(req,
                                                                    study::String,
                                                                    run::String)
    err = _validate_run_request(study, run)
    !isnothing(err) && return err
    group = _req_group(req)
    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json([])

    table = try
        b = JSON3.read(String(req.body))
        get(b, :table, "merged")
    catch
        "merged"
    end

    with_results_db(dir) do con
        levels = taxonomy_levels(con, table)
        json(levels)
    end
end

## Cross-run helpers
_opt_string(spec, key) = let v = get(spec, key, nothing)
    isnothing(v) || isempty(string(v)) ? nothing : string(v)
end

function _resolve_run_duckdb(study::String, run_spec)
    run_name = string(get(run_spec, :run, ""))
    group  = _opt_string(run_spec, :group)
    prefix = _opt_string(run_spec, :prefix)
    dir = _require_duckdb(study, run_name; group)
    isnothing(dir) && return nothing
    label = if !isnothing(prefix)
        prefix
    elseif !isnothing(group)
        "$group/$run_name"
    else
        run_name
    end
    (; run=run_name, group, dir, label, prefix)
end

"""Filter sample columns to those matching a prefix (e.g. "SubgroupName_")."""
function _filter_by_prefix(scols::Vector{String}, prefix::Union{String,Nothing})
    isnothing(prefix) && return scols
    pfx = prefix * "_"
    filter(c -> startswith(c, pfx), scols)
end

function _resolved_group_label(resolved)
    if !isnothing(resolved.prefix)
        resolved.prefix
    elseif !isnothing(resolved.group)
        resolved.group
    else
        resolved.run
    end
end

function _resolved_run_label(resolved)
    if !isnothing(resolved.group)
        "$(resolved.group)/$(resolved.run)"
    else
        resolved.run
    end
end

function _study_seed(study::String)
    resolved = _resolve_config(study)
    Int(get(get(resolved, "seed", (; value=DEFAULT_SEED, source="default")), :value, DEFAULT_SEED))
end

function _study_alpha_plot_config(study::String)
    resolved = _resolve_config(study)
    get_value(key, default) = get(get(resolved, key, (; value=default, source="default")), :value, default)
    (;
        show_points = Bool(get_value("analysis.alpha.show_points", true)),
        annotate_significance = Bool(get_value("analysis.alpha.annotate_significance", false)),
        pairwise_brackets = Bool(get_value("analysis.alpha.pairwise_brackets", false)),
        paired_samples = Bool(get_value("analysis.alpha.paired_samples", false)),
        significance_test = String(get_value("analysis.alpha.significance_test", "kruskal_wallis")),
    )
end

function _comparison_sample_id(sample_col::String, resolved)
    sample_id = sample_col
    if !isnothing(resolved.prefix)
        prefix = resolved.prefix * "_"
        startswith(sample_id, prefix) && (sample_id = sample_id[length(prefix)+1:end])
    end
    # Normalise common naming schemes so paired alpha tests can match the same
    # biological sample across tissues/runs, e.g. 12c_v vs 12l_v -> 12.
    sample_id = replace(sample_id, r"_[^_]+$" => "")
    sample_id = replace(sample_id, r"([0-9])[A-Za-z]$" => s"\1")
    sample_id
end

function _comparison_pair_key(sample_col::String, resolved; mode::Symbol)
    if mode == :run
        # Study/group-level run comparison: preserve subgroup+tissue identity,
        # drop only the run/method suffix so Caecum_12c_m pairs with Caecum_12c_v.
        return replace(sample_col, r"_[^_]+$" => "")
    else
        # Within-run subgroup comparison: pair by the underlying animal/sample id,
        # ignoring subgroup prefix, run/method suffix, and tissue suffix.
        return _comparison_sample_id(sample_col, resolved)
    end
end

function _expand_comparison_run_specs(study::String, runs_spec)
    expanded = NamedTuple[]
    for spec in runs_spec
        run_name = string(get(spec, :run, ""))
        isempty(run_name) && continue
        group  = _opt_string(spec, :group)
        prefix = _opt_string(spec, :prefix)

        if !isnothing(prefix)
            push!(expanded, (; run=run_name, group, prefix))
            continue
        end

        data_rel = isnothing(group) ? run_name : joinpath(group, run_name)
        data_dir = joinpath(ServerState.data_dir(), study, data_rel)
        if _is_pooled(data_dir)
            for subgroup in _subgroup_names(data_dir)
                push!(expanded, (; run=run_name, group, prefix=subgroup))
            end
        else
            push!(expanded, (; run=run_name, group, prefix=nothing))
        end
    end
    expanded
end

function _collect_cross_run_asvs(study::String, runs_spec, table, params)
    run_data = Tuple{String, Vector{String}, DataFrame}[]
    group_labels = String[]
    run_labels = String[]
    missing_seq = false

    for spec in _expand_comparison_run_specs(study, runs_spec)
        resolved = _resolve_run_duckdb(study, spec)
        isnothing(resolved) && continue
        with_results_db(resolved.dir) do con
            scols = _filter_by_prefix(sample_columns(con, table), resolved.prefix)
            isempty(scols) && return

            columns = _duckdb_columns(con, table)
            if !("sequence" in columns)
                missing_seq = true
                return
            end

            where_clause, where_params = _build_where(params, columns)
            seq_filter = isempty(where_clause) ?
                "WHERE \"sequence\" IS NOT NULL AND TRIM(\"sequence\") != ''" :
                "$where_clause AND \"sequence\" IS NOT NULL AND TRIM(\"sequence\") != ''"

            sum_exprs = join(["SUM(COALESCE(\"$c\", 0)) AS \"$c\"" for c in scols], ", ")
            sql = """
                SELECT "sequence", $sum_exprs
                FROM "$table" $seq_filter
                GROUP BY "sequence"
            """
            df = DataFrame(DBInterface.execute(con, sql, where_params))
            nrow(df) == 0 && return

            push!(run_data, (resolved.label, scols, df))
            group_label = _resolved_group_label(resolved)
            run_label = _resolved_run_label(resolved)
            append!(group_labels, fill(group_label, length(scols)))
            append!(run_labels, fill(run_label, length(scols)))
        end
    end

    run_data, group_labels, run_labels, missing_seq
end

function _drop_empty_samples(mat::Matrix{Float64}, aligned::AbstractVector...)
    keep = vec(sum(mat; dims=2) .> 0)
    filtered = Any[mat[keep, :]]
    for values in aligned
        push!(filtered, values[keep])
    end
    Tuple(filtered), count(keep), length(keep)
end

## Cross-run alpha comparison
@post "/api/v1/studies/{study}/analysis/alpha" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    body = JSON3.read(String(req.body))
    runs_spec = get(body, :runs, [])
    isempty(runs_spec) && return json_error(400, "no_runs", "Provide at least one run")
    table = get(body, :table, "merged")
    params = _body_filter_params(body)

    expanded_specs = _expand_comparison_run_specs(study, runs_spec)
    resolved_specs = filter(!isnothing, [_resolve_run_duckdb(study, spec) for spec in expanded_specs])
    isempty(resolved_specs) && return json_error(400, "no_data", "No data found for any run")

    compare_mode = length(unique(_resolved_run_label(r) for r in resolved_specs)) > 1 ? :run : :subgroup
    groups = Dict{String, NamedTuple{(:sample_ids, :richness, :shannon, :simpson), Tuple{Vector{String}, Vector{Int}, Vector{Float64}, Vector{Float64}}}}()

    for resolved in resolved_specs
        with_results_db(resolved.dir) do con
            scols = _filter_by_prefix(sample_columns(con, table), resolved.prefix)
            isempty(scols) && return
            columns = _duckdb_columns(con, table)
            where_clause, where_params = _build_where(params, columns)
            mat = filtered_counts(con, table, scols, where_clause, where_params)
            r = Int[]; sh = Float64[]; si = Float64[]
            for i in eachindex(scols)
                counts = round.(Int, mat[i, :])
                push!(r, richness(counts))
                push!(sh, shannon(counts))
                push!(si, simpson(counts))
            end
            sample_ids = [_comparison_pair_key(s, resolved; mode=compare_mode) for s in scols]
            label = compare_mode == :run ? _resolved_run_label(resolved) : _resolved_group_label(resolved)
            bucket = get!(groups, label, (; sample_ids=String[], richness=Int[], shannon=Float64[], simpson=Float64[]))
            append!(bucket.sample_ids, sample_ids)
            append!(bucket.richness, r)
            append!(bucket.shannon, sh)
            append!(bucket.simpson, si)
        end
    end

    isempty(groups) && return json_error(400, "no_data", "No data found for any run")
    ordered = [(label, vals.sample_ids, vals.richness, vals.shannon, vals.simpson)
               for (label, vals) in sort(collect(groups); by = first)]
    alpha_cfg = _study_alpha_plot_config(study)
    fig = alpha_boxplot(ordered;
                        show_points=alpha_cfg.show_points,
                        annotate_significance=alpha_cfg.annotate_significance,
                        pairwise_brackets=alpha_cfg.pairwise_brackets,
                        paired_samples=alpha_cfg.paired_samples,
                        significance_test=alpha_cfg.significance_test)
    HTTP.Response(200, ["Content-Type" => "application/json"],
                  body=JSON3.write(fig))
end

## Server capabilities
@get "/api/v1/capabilities" function(req)
    json((; r_available=r_available()))
end

## Cross-run NMDS
@post "/api/v1/studies/{study}/analysis/nmds" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    !r_available() && return json_error(503, "r_unavailable",
                                             "R/vegan is not available - install R and the vegan package")
    body = JSON3.read(String(req.body))
    runs_spec = get(body, :runs, [])
    length(runs_spec) < 2 && return json_error(400, "too_few_runs",
                                                    "NMDS requires at least 2 runs")
    table = get(body, :table, "merged")
    params = _body_filter_params(body)

    run_data, group_labels, run_name_labels, missing_seq = _collect_cross_run_asvs(study, runs_spec, table, params)

    missing_seq && return json_error(400, "no_sequence_data",
                                         "Table '$table' does not contain sequence data — select an ASV-level table (e.g. merged or asv_counts)")
    length(run_data) < 2 && return json_error(400, "too_few_runs",
                                                   "Need data from at least 2 runs for NMDS")
    mat, all_samples, _, _ = combined_asv_counts_across_runs(run_data)
    (filtered, kept, total) = _drop_empty_samples(mat, all_samples, group_labels, run_name_labels)
    mat, all_samples, group_labels, run_name_labels = filtered
    kept < total && @info "NMDS: dropped $(total - kept) empty samples before ordination"
    size(mat, 1) < 3 && return json_error(400, "too_few_samples",
                                               "Need at least 3 non-empty samples for NMDS")

    coords, stress = run_nmds(mat; seed=_study_seed(study))
    any(isnan, coords) && return json_error(500, "nmds_failed", "NMDS computation failed")

    # Use a single visual dimension when only one factor varies; otherwise encode
    # run by colour and subgroup/group by shape.
    unique_groups = unique(group_labels)
    unique_runs = unique(run_name_labels)
    colour_by = if length(unique_runs) <= 1
        group_labels
    else
        run_name_labels
    end
    shape_by = if length(unique_runs) <= 1 || length(unique_groups) <= 1
        nothing
    else
        group_labels
    end
    fig = nmds_chart(coords, all_samples;
                     colour_by=colour_by, stress, shape_by=shape_by)
    HTTP.Response(200, ["Content-Type" => "application/json"],
                  body=JSON3.write(fig))
end

## Cross-run PERMANOVA
@post "/api/v1/studies/{study}/analysis/permanova" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    !r_available() && return json_error(503, "r_unavailable",
                                             "R/vegan is not available")
    body = JSON3.read(String(req.body))
    runs_spec = get(body, :runs, [])
    length(runs_spec) < 2 && return json_error(400, "too_few_runs",
                                                    "PERMANOVA requires at least 2 runs")
    table = get(body, :table, "merged")
    params = _body_filter_params(body)

    run_data, group_labels, run_labels, missing_seq = _collect_cross_run_asvs(study, runs_spec, table, params)

    missing_seq && return json_error(400, "no_sequence_data",
                                         "Table '$table' does not contain sequence data — select an ASV-level table (e.g. merged or asv_counts)")
    length(run_data) < 2 && return json_error(400, "too_few_runs",
                                                   "Need data from at least 2 runs")
    mat, all_samples, _, _ = combined_asv_counts_across_runs(run_data)
    (filtered, kept, total) = _drop_empty_samples(mat, all_samples, group_labels, run_labels)
    mat, all_samples, group_labels, run_labels = filtered
    kept < total && @info "PERMANOVA: dropped $(total - kept) empty samples before analysis"
    size(mat, 1) < 3 && return json_error(400, "too_few_samples",
                                               "Need at least 3 non-empty samples for PERMANOVA")

    metadata = DataFrame(sample=all_samples, group=group_labels, run=run_labels)
    result = run_permanova(mat, metadata; seed=_study_seed(study))
    isnothing(result) && return json_error(500, "permanova_failed", "PERMANOVA computation failed")
    hasproperty(result, :message) && return json_error(500, "permanova_failed",
                                                           "PERMANOVA computation failed: $(result.message)")
    json(result)
end
