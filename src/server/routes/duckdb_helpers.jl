# Shared DuckDB query helpers used by results.jl and annotations.jl

using JSON3, DataFrames, OrderedCollections, DuckDB, DBInterface

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
