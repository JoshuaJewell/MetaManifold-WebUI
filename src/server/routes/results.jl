# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: results (plot catalog + Plotly JSON + table data)

using JSON3, CSV, DataFrames, OrderedCollections, YAML, XLSX, DuckDB, DBInterface

## Tables
function _merge_dir(study::String, run::String;
                    group::Union{String,Nothing}=nothing)
    joinpath(_run_project_dir(study, run; group), "merged")
end

function _require_duckdb(study::String, run::String;
                         group::Union{String,Nothing}=nothing)
    dir = _merge_dir(study, run; group)
    db_path = joinpath(dir, "results.duckdb")
    isfile(db_path) || return nothing
    dir
end

## DuckDB SQL helpers
function _duckdb_columns(con, table::String)
    result = DBInterface.execute(con, "SELECT column_name FROM information_schema.columns WHERE table_name = ?", [table])
    df = DataFrame(result)
    String[string(row.column_name) for row in eachrow(df)]
end

# Returns integer columns that represent sample read counts (excludes *_boot bootstrap columns).
function _sample_count_columns(con, table::String)
    result = DBInterface.execute(con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = ? AND data_type LIKE '%INT%'",
        [table])
    filter(c -> !endswith(c, "_boot"),
           String[string(r.column_name) for r in eachrow(DataFrame(result))])
end

function _sum_reads(con, table::String, count_cols::Vector{String}, where::String="", sql_params::Vector=Any[])
    isempty(count_cols) && return 0
    sum_expr = join(["COALESCE(SUM(\"$c\"), 0)" for c in count_cols], " + ")
    row = only(DataFrame(DBInterface.execute(con, "SELECT ($sum_expr) AS n FROM \"$table\" $where", sql_params)))
    coalesce(row.n, 0)
end

function _build_where(params::Dict{String,String}, columns::Vector{String})
    clauses = String[]
    sql_params = Any[]
    col_set = Set(columns)

    for (k, v) in params
        startswith(k, "col.") || continue
        col = k[5:end]
        col in col_set || continue
        isempty(v) && continue
        push!(clauses, "LOWER(CAST(\"$col\" AS VARCHAR)) LIKE ?")
        push!(sql_params, "%" * lowercase(v) * "%")
    end

    for (k, v) in params
        startswith(k, "col_in.") || continue
        col = k[8:end]
        col in col_set || continue
        vals = split(v, "|")
        if isempty(vals) || (length(vals) == 1 && isempty(vals[1]))
            push!(clauses, "1 = 0")  # empty include = match nothing
        else
            placeholders = join(["?" for _ in vals], ", ")
            push!(clauses, "CAST(\"$col\" AS VARCHAR) IN ($placeholders)")
            append!(sql_params, vals)
        end
    end

    for (k, v) in params
        startswith(k, "col_ex.") || continue
        col = k[8:end]
        col in col_set || continue
        isempty(v) && continue
        vals = split(v, "|")
        placeholders = join(["?" for _ in vals], ", ")
        push!(clauses, "CAST(\"$col\" AS VARCHAR) NOT IN ($placeholders)")
        append!(sql_params, vals)
    end

    for (k, v) in params
        startswith(k, "col_min.") || continue
        col = k[9:end]
        col in col_set || continue
        threshold = tryparse(Float64, v)
        isnothing(threshold) && continue
        push!(clauses, "TRY_CAST(\"$col\" AS DOUBLE) IS NOT NULL AND TRY_CAST(\"$col\" AS DOUBLE) >= ?")
        push!(sql_params, threshold)
    end

    for (k, v) in params
        startswith(k, "col_max.") || continue
        col = k[9:end]
        col in col_set || continue
        threshold = tryparse(Float64, v)
        isnothing(threshold) && continue
        push!(clauses, "TRY_CAST(\"$col\" AS DOUBLE) IS NOT NULL AND TRY_CAST(\"$col\" AS DOUBLE) <= ?")
        push!(sql_params, threshold)
    end

    filter_q = get(params, "filter", nothing)
    if !isnothing(filter_q) && !isempty(filter_q)
        col_checks = join(["LOWER(CAST(\"$c\" AS VARCHAR)) LIKE ?" for c in columns], " OR ")
        push!(clauses, "($col_checks)")
        for _ in columns
            push!(sql_params, "%" * lowercase(filter_q) * "%")
        end
    end

    where = isempty(clauses) ? "" : "WHERE " * join(clauses, " AND ")
    (where, sql_params)
end

function _order_clause(sort_by::Union{String,Nothing}, sort_dir::String, columns::Vector{String})
    isnothing(sort_by) && return ""
    sort_by in columns || return ""
    dir_str = sort_dir == "desc" ? "DESC" : "ASC"
    if sort_by == "SeqName"
        return """ORDER BY
            regexp_extract("SeqName", '^[A-Za-z_]*') $dir_str,
            CASE WHEN regexp_extract("SeqName", '(\\d+)') = '' THEN 0
                 ELSE CAST(regexp_extract("SeqName", '(\\d+)') AS INTEGER) END $dir_str"""
    end
    "ORDER BY \"$sort_by\" $dir_str"
end

function _duckdb_rows(con, sql::String, params::Vector=[])
    result = DBInterface.execute(con, sql, params)
    df = DataFrame(result)
    colnames = names(df)
    [OrderedDict(c => (ismissing(row[c]) ? nothing : row[c]) for c in colnames)
     for row in eachrow(df)]
end

## Shared body parser

function _body_filter_params(body)
    params = Dict{String,String}()

    filter_q = get(body, :filter, nothing)
    !isnothing(filter_q) && (params["filter"] = string(filter_q))
    sort_by = get(body, :sortBy, nothing)
    !isnothing(sort_by) && (params["sort"] = string(sort_by))
    sort_dir = get(body, :sortDir, nothing)
    !isnothing(sort_dir) && (params["sort_dir"] = string(sort_dir))

    col_filters = get(body, :colFilters, nothing)
    isnothing(col_filters) && return params

    for (col, f) in pairs(col_filters)
        col_str = string(col)
        include_vals = get(f, :include, nothing)
        if !isnothing(include_vals) && include_vals isa AbstractVector
            params["col_in.$col_str"] = join(string.(include_vals), "|")
        end
        exclude_vals = get(f, :exclude, nothing)
        if !isnothing(exclude_vals) && exclude_vals isa AbstractVector && !isempty(exclude_vals)
            params["col_ex.$col_str"] = join(string.(exclude_vals), "|")
        end
        text_val = get(f, :text, nothing)
        !isnothing(text_val) && !isempty(string(text_val)) && (params["col.$col_str"] = string(text_val))
        min_val = get(f, :min, nothing)
        !isnothing(min_val) && (params["col_min.$col_str"] = string(min_val))
        max_val = get(f, :max, nothing)
        !isnothing(max_val) && (params["col_max.$col_str"] = string(max_val))
    end

    params
end

## Table catalog
@get "/api/v1/studies/{study}/runs/{run}/results/tables" function(req,
                                                                    study::String,
                                                                    run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)
    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json([])

    with_results_db(dir) do con
        tables_result = DBInterface.execute(con, "SHOW TABLES")
        catalog = map(eachrow(DataFrame(tables_result))) do row
            tname = string(row.name)
            cnt = only(DataFrame(DBInterface.execute(con, "SELECT COUNT(*) AS n FROM \"$(tname)\"")))
            label = replace(tname, "_" => " ") |> titlecase

            count_cols  = _sample_count_columns(con, tname)
            n_samples   = length(count_cols)
            total_reads = _sum_reads(con, tname, count_cols)

            (; id=tname, label, rows=cnt.n, total_reads, n_samples)
        end
        hidden = Set(["cluster_membership"])
        json(collect(filter(t -> !(t.id in hidden), catalog)))
    end
end

## Distinct values
@post "/api/v1/studies/{study}/runs/{run}/results/tables/{table}/distinct/{column}" function(req,
                                                                            study::String,
                                                                            run::String,
                                                                            table::String,
                                                                            column::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)

    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json_error(404, "no_results",
                                             "No results database for run '$run' - run the pipeline first")

    body = length(req.body) > 0 ? JSON3.read(String(req.body)) : (;)
    params = _body_filter_params(body)

    with_results_db(dir) do con
        columns = _duckdb_columns(con, table)
        column in columns || return json_error(404, "column_not_found",
                                                    "Column '$column' not found")
        (where, sql_params) = _build_where(params, columns)

        num_sql = "SELECT MIN(TRY_CAST(\"$column\" AS DOUBLE)) AS mn,
                          MAX(TRY_CAST(\"$column\" AS DOUBLE)) AS mx,
                          COUNT(TRY_CAST(\"$column\" AS DOUBLE)) AS cnt,
                          SUM(TRY_CAST(\"$column\" AS DOUBLE)) AS sm,
                          AVG(TRY_CAST(\"$column\" AS DOUBLE)) AS avg,
                          MEDIAN(TRY_CAST(\"$column\" AS DOUBLE)) AS med
                   FROM \"$table\" $where"
        num_row = only(DataFrame(DBInterface.execute(con, num_sql, sql_params)))

        if !ismissing(num_row.mn) && num_row.cnt > 0
            json((; column, type="numeric", min=num_row.mn, max=num_row.mx, count=num_row.cnt,
                    sum=num_row.sm, mean=num_row.avg, median=num_row.med))
        else
            not_null = "\"$column\" IS NOT NULL"
            full_where = isempty(where) ? "WHERE $not_null" : "$where AND $not_null"
            vals_sql = "SELECT DISTINCT CAST(\"$column\" AS VARCHAR) AS val
                        FROM \"$table\" $full_where
                        ORDER BY val"
            vals_df = DataFrame(DBInterface.execute(con, vals_sql, sql_params))
            vals = String[string(row.val) for row in eachrow(vals_df)]
            json((; column, type="text", values=vals, count=length(vals)))
        end
    end
end

## Table query (paginated, filtered, sorted)
@post "/api/v1/studies/{study}/runs/{run}/results/tables/{table}/query" function(req,
                                                                            study::String,
                                                                            run::String,
                                                                            table::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)

    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json_error(404, "no_results",
                                             "No results database for run '$run' - run the pipeline first")

    body     = JSON3.read(String(req.body))
    params   = _body_filter_params(body)
    page     = get(body, :page, 1)
    per_page = get(body, :perPage, 100)
    sort_by  = get(params, "sort",     nothing)
    sort_dir = get(params, "sort_dir", "asc")

    with_results_db(dir) do con
        columns = _duckdb_columns(con, table)
        isempty(columns) && return json_error(404, "table_not_found",
                                                   "Table '$table' not found")

        total_unfiltered = only(DataFrame(DBInterface.execute(con,
            "SELECT COUNT(*) AS n FROM \"$table\""))).n

        (where, sql_params) = _build_where(params, columns)

        total = only(DataFrame(DBInterface.execute(con,
            "SELECT COUNT(*) AS n FROM \"$table\" $where", sql_params))).n

        count_cols              = _sample_count_columns(con, table)
        total_reads_unfiltered  = _sum_reads(con, table, count_cols)
        total_reads             = _sum_reads(con, table, count_cols, where, sql_params)

        order = _order_clause(sort_by, sort_dir, columns)

        offset = (page - 1) * per_page
        data_sql = "SELECT * FROM \"$table\" $where $order LIMIT $per_page OFFSET $offset"
        rows = _duckdb_rows(con, data_sql, sql_params)

        json((; total, total_unfiltered, total_reads, total_reads_unfiltered, page, per_page, columns, rows))
    end
end

## Filter presets
function _filters_dir()
    joinpath(dirname(ServerState.projects_dir()), "config", "filters")
end

@get "/api/v1/filter-presets" function(req)
    dir = _filters_dir()
    isdir(dir) || return json([])
    presets = map(filter(f -> endswith(f, ".yml"), readdir(dir))) do f
        path = joinpath(dir, f)
        name = splitext(f)[1]
        label = replace(name, "_" => " ") |> s -> replace(s, "." => " - ") |> titlecase
        desc = ""
        try
            lines = readlines(path)
            for line in lines
                m = match(r"^#\s*(.+)", line)
                if !isnothing(m)
                    desc = m.captures[1]
                    break
                end
            end
        catch; end
        (; name, label, file=f, description=desc)
    end
    json(presets)
end

@post "/api/v1/studies/{study}/runs/{run}/results/tables/{table}/apply-preset" function(req,
                                                                                         study::String,
                                                                                         run::String,
                                                                                         table::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)

    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json_error(404, "no_results",
                                             "No results database for run '$run' - run the pipeline first")

    body = JSON3.read(String(req.body))
    preset_file = get(body, :preset, nothing)
    isnothing(preset_file) && return json_error(400, "missing_preset", "Body must include 'preset'")

    filter_path = joinpath(_filters_dir(), string(preset_file))
    isfile(filter_path) || return json_error(404, "preset_not_found",
                                                  "Preset '$preset_file' not found")

    config = YAML.load_file(filter_path)
    raw_filters = get(config, "filters", [])

    has_pipeline     = any(f -> f isa Dict && haskey(f, "pattern"), raw_filters)
    has_remove_empty = !isempty(get(config, "remove_empty", []))

    col_filters = Dict{String, Any}()
    params      = Dict{String, String}()

    if has_pipeline || has_remove_empty
        with_results_db(dir) do con
            columns = _duckdb_columns(con, table)
            col_set = Set(columns)

            ref_cols = Set{String}()
            for item in raw_filters
                (item isa Dict && haskey(item, "pattern")) || continue
                col = string(get(item, "column", ""))
                !isempty(col) && col in col_set && push!(ref_cols, col)
            end
            for col in get(config, "remove_empty", [])
                col_str = string(col)
                col_str in col_set && push!(ref_cols, col_str)
            end

            col_include = Dict{String, Set{String}}()
            for col in ref_cols
                vals_df = DataFrame(DBInterface.execute(con,
                    "SELECT DISTINCT CAST(\"$col\" AS VARCHAR) AS val FROM \"$table\" WHERE \"$col\" IS NOT NULL ORDER BY val"))
                col_include[col] = Set(string(row.val) for row in eachrow(vals_df))
            end

            for col in get(config, "remove_empty", [])
                col_str = string(col)
                haskey(col_include, col_str) || continue
                filter!(v -> !isempty(strip(v)) && lowercase(strip(v)) != "blank",
                        col_include[col_str])
            end

            for item in raw_filters
                (item isa Dict && haskey(item, "pattern")) || continue
                col       = string(get(item, "column", ""))
                pat       = string(get(item, "pattern", ""))
                action    = lowercase(string(get(item, "action", "exclude")))
                use_regex = get(item, "regex", false) == true
                (isempty(col) || isempty(pat) || !haskey(col_include, col)) && continue
                matcher = use_regex ? (v -> occursin(Regex(pat), v)) :
                                      (v -> occursin(pat, v))
                if action == "keep"
                    filter!(v -> matcher(v), col_include[col])
                else  # exclude (default)
                    filter!(v -> !matcher(v), col_include[col])
                end
            end

            for (col, vals) in col_include
                vals_vec = sort!(collect(vals))
                params["col_in.$col"] = join(vals_vec, "|")
                col_filters[col] = Dict("include" => vals_vec)
            end

            for item in raw_filters
                (item isa Dict && haskey(item, "type")) || continue
                col_str = string(get(item, "column", ""))
                isempty(col_str) && continue
                typ = string(get(item, "type", ""))
                if typ == "include"
                    vals = get(item, "values", [])
                    vals_str = [string(v) for v in vals]
                    params["col_in.$col_str"] = join(vals_str, "|")
                    col_filters[col_str] = Dict("include" => vals_str)
                elseif typ == "min"
                    val = get(item, "value", nothing)
                    isnothing(val) && continue
                    params["col_min.$col_str"] = string(val)
                    get!(col_filters, col_str, Dict{String,Any}())["min"] = val
                elseif typ == "max"
                    val = get(item, "value", nothing)
                    isnothing(val) && continue
                    params["col_max.$col_str"] = string(val)
                    get!(col_filters, col_str, Dict{String,Any}())["max"] = val
                end
            end

            rows_before = only(DataFrame(DBInterface.execute(con,
                "SELECT COUNT(*) AS n FROM \"$table\""))).n
            (where, sql_params) = _build_where(params, columns)
            rows_after = only(DataFrame(DBInterface.execute(con,
                "SELECT COUNT(*) AS n FROM \"$table\" $where", sql_params))).n
            json((; preset=preset_file, rows_before, rows_after, filters=col_filters))
        end
    else
        for item in raw_filters
            item isa Dict || continue
            col_str = string(get(item, "column", ""))
            isempty(col_str) && continue
            typ = string(get(item, "type", ""))
            if typ == "include"
                vals = get(item, "values", [])
                vals_str = [string(v) for v in vals]
                params["col_in.$col_str"] = join(vals_str, "|")
                col_filters[col_str] = Dict("include" => vals_str)
            elseif typ == "min"
                val = get(item, "value", nothing)
                isnothing(val) && continue
                params["col_min.$col_str"] = string(val)
                get!(col_filters, col_str, Dict{String,Any}())["min"] = val
            elseif typ == "max"
                val = get(item, "value", nothing)
                isnothing(val) && continue
                params["col_max.$col_str"] = string(val)
                get!(col_filters, col_str, Dict{String,Any}())["max"] = val
            end
        end

        with_results_db(dir) do con
            columns = _duckdb_columns(con, table)
            rows_before = only(DataFrame(DBInterface.execute(con,
                "SELECT COUNT(*) AS n FROM \"$table\""))).n
            (where, sql_params) = _build_where(params, columns)
            rows_after = only(DataFrame(DBInterface.execute(con,
                "SELECT COUNT(*) AS n FROM \"$table\" $where", sql_params))).n
            json((; preset=preset_file, rows_before, rows_after, filters=col_filters))
        end
    end
end

## Save filters to config/filters/
@post "/api/v1/filter-presets/{name}" function(req, name::String)
    occursin(r"^[a-zA-Z0-9._-]+$", name) || return json_error(400, "invalid_name",
        "Filter name must contain only letters, numbers, dots, hyphens, and underscores")

    dir = _filters_dir()
    mkpath(dir)
    filename = endswith(name, ".yml") ? name : name * ".yml"
    path = joinpath(dir, filename)

    body = JSON3.read(String(req.body))
    filters = get(body, :filters, nothing)
    isnothing(filters) && return json_error(400, "missing_filters", "Body must include 'filters'")
    description = get(body, :description, "Saved from WebUI")

    lines = String["# $description", ""]
    filter_entries = collect(pairs(filters))
    if isempty(filter_entries)
        lines = vcat(lines, ["filters: []"])
    else
        push!(lines, "filters:")
        for (col, f) in filter_entries
            col_str = string(col)
            include_vals = get(f, :include, nothing)
            min_val = get(f, :min, nothing)
            max_val = get(f, :max, nothing)

            if !isnothing(include_vals) && include_vals isa AbstractVector
                push!(lines, "  - column: \"$col_str\"")
                push!(lines, "    type: include")
                push!(lines, "    values:")
                for v in include_vals
                    push!(lines, "      - \"$(string(v))\"")
                end
            end
            if !isnothing(min_val)
                push!(lines, "  - column: \"$col_str\"")
                push!(lines, "    type: min")
                push!(lines, "    value: $(string(min_val))")
            end
            if !isnothing(max_val)
                push!(lines, "  - column: \"$col_str\"")
                push!(lines, "    type: max")
                push!(lines, "    value: $(string(max_val))")
            end
        end
    end

    open(path, "w") do io
        println(io, join(lines, "\n"))
    end

    stem = splitext(filename)[1]
    label = replace(stem, "_" => " ") |> s -> replace(s, "." => " - ") |> titlecase
    json((; name=stem, label, file=filename, description=string(description)))
end

## Save filtered table to merged directory
@post "/api/v1/studies/{study}/runs/{run}/results/tables/{table}/save" function(req,
                                                                                study::String,
                                                                                run::String,
                                                                                table::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)

    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json_error(404, "no_results",
                                             "No results database for run '$run' - run the pipeline first")

    body = JSON3.read(String(req.body))
    save_name = get(body, :name, nothing)
    isnothing(save_name) && return json_error(400, "missing_name", "Body must include 'name'")
    save_name = string(save_name)
    occursin(r"^[a-zA-Z0-9._-]+$", save_name) || return json_error(400, "invalid_name",
        "Name must contain only letters, numbers, dots, hyphens, and underscores")

    params = _body_filter_params(body)
    sort_by  = get(params, "sort",     nothing)
    sort_dir = get(params, "sort_dir", "asc")

    with_results_db(dir) do con
        columns = _duckdb_columns(con, table)
        (where, sql_params) = _build_where(params, columns)

        order = _order_clause(sort_by, sort_dir, columns)

        data_sql = "SELECT * FROM \"$table\" $where $order"
        df = DataFrame(DBInterface.execute(con, data_sql, sql_params))

        output_path = joinpath(dir, save_name * ".csv")
        CSV.write(output_path, df)

        DBInterface.execute(con,
            "CREATE OR REPLACE TABLE \"$(save_name)\" AS SELECT * FROM read_csv_auto('$(output_path)', auto_detect=true)")

        @info "Saved filtered table: $output_path ($(nrow(df)) rows)"
        json((; name=save_name, path=output_path, rows=nrow(df)))
    end
end

## Export filtered table as Excel (.xlsx)
@post "/api/v1/studies/{study}/runs/{run}/results/tables/{table}/export" function(req,
                                                                                  study::String,
                                                                                  run::String,
                                                                                  table::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)

    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json_error(404, "no_results",
                                             "No results database for run '$run' - run the pipeline first")

    body   = JSON3.read(String(req.body))
    params = _body_filter_params(body)
    sort_by  = get(params, "sort",     nothing)
    sort_dir = get(params, "sort_dir", "asc")

    df = with_results_db(dir) do con
        columns = _duckdb_columns(con, table)
        (where, sql_params) = _build_where(params, columns)
        order = _order_clause(sort_by, sort_dir, columns)
        DataFrame(DBInterface.execute(con,
            "SELECT * FROM \"$table\" $where $order", sql_params))
    end

    tmp = tempname() * ".xlsx"
    try
        XLSX.writetable(tmp, df)
        data = read(tmp)
        filename = get(body, :filename, "$(study)_$(run)_filtered.xlsx")
        HTTP.Response(200, [
            "Content-Type" => "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            "Content-Disposition" => "attachment; filename=\"$filename\"",
        ]; body=data)
    finally
        isfile(tmp) && rm(tmp)
    end
end

## Delete a filter preset

@delete "/api/v1/filter-presets/{name}" function(req, name::String)
    filename = endswith(name, ".yml") ? name : name * ".yml"
    path = joinpath(_filters_dir(), filename)
    isfile(path) || return json_error(404, "preset_not_found",
                                          "Preset '$name' not found")
    rm(path)
    @info "Deleted filter preset: $path"
    json((; deleted=filename))
end

## Delete a table from merged directory
@delete "/api/v1/studies/{study}/runs/{run}/results/tables/{table}" function(req,
                                                                              study::String,
                                                                              run::String,
                                                                              table::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)
    dir = _merge_dir(study, run; group)
    path = joinpath(dir, table * ".csv")
    isfile(path) || return json_error(404, "table_not_found",
                                          "Table '$table' not found")
    rm(path)

    db_dir = _require_duckdb(study, run; group)
    if !isnothing(db_dir)
        try
            with_results_db(db_dir) do con
                DBInterface.execute(con, "DROP TABLE IF EXISTS \"$table\"")
            end
        catch; end
    end

    @info "Deleted table: $path"
    json((; deleted=table))
end

## OTU member lookup (JOIN cluster_membership to merged)
@get "/api/v1/studies/{study}/runs/{run}/results/otu-members/{otu}" function(req,
                                                                              study::String,
                                                                              run::String,
                                                                              otu::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)

    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json_error(404, "no_results",
                                             "No results database for run '$run' - run the pipeline first")

    with_results_db(dir) do con
        tables = Set(string.(DataFrame(DBInterface.execute(con, "SHOW TABLES")).name))
        "cluster_membership" in tables || return json_error(404, "no_cluster_data",
                                                                  "No cluster membership data available")
        "merged" in tables || return json_error(404, "no_merged_data",
                                                     "No merged table available")

        rows = _duckdb_rows(con,
            """SELECT m.*
               FROM cluster_membership cm
               JOIN merged m ON m."SeqName" = cm."ASV"
               WHERE cm."OTU" = ?
               ORDER BY
                 regexp_extract(m."SeqName", '^[A-Za-z_]*'),
                 CASE WHEN regexp_extract(m."SeqName", '(\\d+)') = '' THEN 0
                      ELSE CAST(regexp_extract(m."SeqName", '(\\d+)') AS INTEGER) END
            """, [otu])

        columns = _duckdb_columns(con, "merged")
        json((; otu, columns, rows))
    end
end

## QC and DADA2 output discovery
@get "/api/v1/studies/{study}/runs/{run}/results/qc" function(req,
                                                               study::String,
                                                               run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)
    run_dir = _run_project_dir(study, run; group)
    report = joinpath(run_dir, "QC", "multiqc_report.html")
    has_report = isfile(report)
    run_rel = isnothing(group) ? run : joinpath(group, run)
    json((; has_report,
           report_url = has_report ? "/files/$study/runs/$run_rel/QC/multiqc_report.html" : nothing))
end

@get "/api/v1/studies/{study}/runs/{run}/results/dada2" function(req,
                                                                   study::String,
                                                                   run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)
    run_dir = _run_project_dir(study, run; group)
    run_rel = isnothing(group) ? run : joinpath(group, run)

    fig_order = ["quality_unfiltered", "quality_filtered", "error_rates",
                 "length_distribution", "length_distribution_filtered"]
    fig_dir = joinpath(run_dir, "dada2", "Figures")
    figures = if isdir(fig_dir)
        available = Set(splitext(f)[1] for f in readdir(fig_dir) if endswith(f, ".pdf"))
        ordered = filter(n -> n in available, fig_order)
        extras = sort(collect(setdiff(available, Set(fig_order))))
        map(vcat(ordered, extras)) do name
            label = replace(name, "_" => " ") |> titlecase
            url = "/files/$study/runs/$run_rel/dada2/Figures/$name.pdf"
            (; name, label, url)
        end
    else
        []
    end

    stats_csv = joinpath(run_dir, "dada2", "Tables", "pipeline_stats.csv")
    has_stats = isfile(stats_csv)

    log_dir = joinpath(run_dir, "dada2", "Logs")
    logs = if isdir(log_dir)
        map(sort(readdir(log_dir))) do f
            url = "/files/$study/runs/$run_rel/dada2/Logs/$f"
            (; name=f, url)
        end
    else
        []
    end

    config_path = joinpath(run_dir, "run_config.yml")
    dada2_config = if isfile(config_path)
        cfg = YAML.load_file(config_path)
        d2 = get(cfg, "dada2", Dict())
        d2 isa Dict || (d2 = Dict())
        config = Dict{String,Any}()
        for key in ("file_patterns", "filter_trim", "dada", "merge", "asv", "taxonomy")
            haskey(d2, key) && (config[key] = d2[key])
        end
        if haskey(config, "dada")
            delete!(config["dada"], "seed")
        end
        if haskey(config, "taxonomy")
            delete!(config["taxonomy"], "multithread")
            delete!(config["taxonomy"], "remote")
        end
        config
    else
        Dict{String,Any}()
    end

    json((; figures, has_stats, logs, config=dada2_config))
end

## DADA2 pipeline stats as table data
@get "/api/v1/studies/{study}/runs/{run}/results/dada2/stats" function(req,
                                                                        study::String,
                                                                        run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)
    stats_path = joinpath(_run_project_dir(study, run; group), "dada2", "Tables", "pipeline_stats.csv")
    isfile(stats_path) || return json_error(404, "not_found", "No pipeline stats found")

    df = CSV.read(stats_path, DataFrame)
    columns = names(df)
    rows = [OrderedDict(c => (ismissing(row[c]) ? nothing : row[c]) for c in columns)
            for row in eachrow(df)]
    json((; columns, rows))
end

@get "/api/v1/studies/{study}/runs/{run}/results/otu-counts" function(req,
                                                                       study::String,
                                                                       run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    group = _req_group(req)

    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json_error(404, "no_results",
                                             "No results database for run '$run' - run the pipeline first")

    with_results_db(dir) do con
        tables = Set(string.(DataFrame(DBInterface.execute(con, "SHOW TABLES")).name))
        "cluster_membership" in tables || return json((; counts=Dict{String,Int}()))

        df = DataFrame(DBInterface.execute(con,
            """SELECT "OTU", COUNT(*) AS n FROM cluster_membership GROUP BY "OTU" """))
        counts = Dict(string(row.OTU) => row.n for row in eachrow(df))
        json((; counts))
    end
end
