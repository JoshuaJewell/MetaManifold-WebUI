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
function _resolve_run_duckdb(study::String, run_spec)
    run_name = string(get(run_spec, :run, ""))
    group = let g = get(run_spec, :group, nothing)
        isnothing(g) || isempty(string(g)) ? nothing : string(g)
    end
    prefix = let p = get(run_spec, :prefix, nothing)
        isnothing(p) || isempty(string(p)) ? nothing : string(p)
    end
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

function _collect_cross_run_taxa(study::String, runs_spec, table, params)
    run_data = Tuple{String, Vector{String}, DataFrame}[]
    group_labels = String[]
    run_labels = String[]

    for spec in runs_spec
        resolved = _resolve_run_duckdb(study, spec)
        isnothing(resolved) && continue
        with_results_db(resolved.dir) do con
            scols = _filter_by_prefix(sample_columns(con, table), resolved.prefix)
            isempty(scols) && return

            columns = _duckdb_columns(con, table)
            where_clause, where_params = _build_where(params, columns)
            levels = taxonomy_levels(con, table)
            isempty(levels) && return

            rank = last(levels)
            agg = aggregate_by_taxon(con, table, scols, rank, where_clause, where_params)
            nrow(agg) == 0 && return

            push!(run_data, (resolved.label, scols, agg))
            group_label = _resolved_group_label(resolved)
            append!(group_labels, fill(group_label, length(scols)))
            append!(run_labels, fill(resolved.label, length(scols)))
        end
    end

    run_data, group_labels, run_labels
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

    groups = Tuple{String, Vector{Int}, Vector{Float64}, Vector{Float64}}[]

    for spec in runs_spec
        resolved = _resolve_run_duckdb(study, spec)
        isnothing(resolved) && continue
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
            push!(groups, (resolved.label, r, sh, si))
        end
    end

    isempty(groups) && return json_error(400, "no_data", "No data found for any run")
    fig = alpha_boxplot(groups)
    HTTP.Response(200, ["Content-Type" => "application/json"],
                  body=JSON3.write(fig))
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

    run_data, group_labels, run_name_labels = _collect_cross_run_taxa(study, runs_spec, table, params)

    length(run_data) < 2 && return json_error(400, "too_few_runs",
                                                   "Need data from at least 2 runs for NMDS")
    mat, all_samples, _, _ = combined_counts_across_runs(run_data)
    size(mat, 1) < 3 && return json_error(400, "too_few_samples",
                                               "Need at least 3 samples for NMDS")

    coords, stress = run_nmds(mat)
    any(isnan, coords) && return json_error(500, "nmds_failed", "NMDS computation failed")

    # Colour by run name (same name = same colour), shape by group when multiple groups present
    has_groups = length(unique(group_labels)) > 1
    fig = nmds_chart(coords, all_samples;
                     colour_by=run_name_labels, stress,
                     shape_by=has_groups ? group_labels : nothing)
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

    run_data, group_labels, run_labels = _collect_cross_run_taxa(study, runs_spec, table, params)

    length(run_data) < 2 && return json_error(400, "too_few_runs",
                                                   "Need data from at least 2 runs")
    mat, all_samples, _, _ = combined_counts_across_runs(run_data)
    size(mat, 1) < 3 && return json_error(400, "too_few_samples",
                                               "Need at least 3 samples for PERMANOVA")

    metadata = DataFrame(sample=all_samples, group=group_labels, run=run_labels)
    result = run_permanova(mat, metadata)
    isnothing(result) && return json_error(500, "permanova_failed", "PERMANOVA computation failed")
    json(result)
end
