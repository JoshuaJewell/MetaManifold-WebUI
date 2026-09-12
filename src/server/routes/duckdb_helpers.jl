# Shared DuckDB query helpers used by the route files
using JSON3, DataFrames, OrderedCollections, DuckDB, DBInterface

## DuckDB SQL helpers
function _duckdb_columns(con, table::String)
    result = DBInterface.execute(con, "SELECT column_name FROM information_schema.columns WHERE table_name = ? ORDER BY ordinal_position", [table])
    df = DataFrame(result)
    String[string(row.column_name) for row in eachrow(df)]
end

# Returns the columns that represent sample read counts. The data-table path and
# the chart path must not hold divergent notions of what a sample column is, or
# the two endpoints report different read totals for the same run.
# `Analysis.sample_columns` is the single definition: it admits every numeric
# type rather than integers alone, and it alone excludes the identifier columns
# (OTU, ASV, Pident) as well as the bootstraps and the derived totals written by
# `Annotation._add_totals!`.
function _sample_count_columns(con, table::String)
    Analysis.sample_columns(con, table)
end

# Largest value of each count column over the filtered rows, the scale for
# heatmap shading. Taken over the whole filtered table so shading does not
# change from page to page.
function _count_maxima(con, table::String, count_cols::Vector{String}, where::String="", sql_params::Vector=Any[])
    isempty(count_cols) && return Dict{String,Float64}()
    exprs = join(["COALESCE(MAX(TRY_CAST(\"$c\" AS DOUBLE)), 0) AS \"m$i\""
                  for (i, c) in enumerate(count_cols)], ", ")
    row = only(eachrow(DataFrame(DBInterface.execute(con, "SELECT $exprs FROM \"$table\" $where", sql_params))))
    Dict{String,Float64}(c => Float64(coalesce(row[Symbol("m$i")], 0.0)) for (i, c) in enumerate(count_cols))
end

function _sum_reads(con, table::String, count_cols::Vector{String}, where::String="", sql_params::Vector=Any[])
    isempty(count_cols) && return 0
    sum_expr = join(["COALESCE(SUM(\"$c\"), 0)" for c in count_cols], " + ")
    row = only(DataFrame(DBInterface.execute(con, "SELECT ($sum_expr) AS n FROM \"$table\" $where", sql_params)))
    coalesce(row.n, 0)
end

## Sample read-count filtering
#
# Samples are columns in these wide tables, so a threshold on a sample's total
# reads selects columns. It cannot be expressed through `_build_where` and is
# applied alongside it: the row filters narrow the rows, and this narrows the
# sample columns those rows are summed over.
#
# `basis` decides which rows the per-sample total is measured over:
#
#   "filtered" (default)  the rows surviving the row filters, so the threshold
#                         reads as "reads left in this sample after filtering".
#                         This is what composes correctly with the other
#                         filters: a sample carrying 10,000 contaminant reads
#                         and 3 real ones fails a 50-read floor, as it should.
#   "raw"                 every row in the table, i.e. the sample's library
#                         size before any filtering.
const _SAMPLE_READ_BASES = ("filtered", "raw")

# Reserved colFilters key the frontend uses to carry the bounds. Mirrored in
# frontend/src/api/types.ts as SAMPLE_READS_FILTER_KEY.
const SAMPLE_READS_FILTER_KEY = "__sample_reads__"

# The (min, max, basis) bounds a request asks for, or all-nothing when it asks
# for none. Parsed from the same params dict `_build_where` reads.
function _sample_read_bounds(params::Dict{String,String})
    parse_bound(key) = let v = get(params, key, "")
        isempty(v) ? nothing : tryparse(Float64, v)
    end
    basis = get(params, "sample_reads_basis", "filtered")
    basis in _SAMPLE_READ_BASES || (basis = "filtered")
    (; min=parse_bound("sample_min_reads"), max=parse_bound("sample_max_reads"), basis)
end

_has_sample_read_bounds(b) = !isnothing(b.min) || !isnothing(b.max)

"""
    _sample_read_totals(con, table, sample_cols, where, sql_params) -> Dict{String,Float64}

Total reads per sample column over the rows selected by `where`. One pass over
the table: a per-column SUM is far cheaper than a query each, and the tables
here run to hundreds of sample columns.
"""
function _sample_read_totals(con, table::String, sample_cols::Vector{String},
                             where::String="", sql_params::Vector=Any[])
    isempty(sample_cols) && return Dict{String,Float64}()
    # Positional aliases keep a sample column named `n`, or one with characters
    # the result-frame lookup would mangle, addressable; the position ties the
    # value back to the name.
    sel = join(["COALESCE(SUM(\"$c\"), 0) AS \"s$i\""
                for (i, c) in enumerate(sample_cols)], ", ")
    row = only(DataFrame(DBInterface.execute(con,
        "SELECT $sel FROM \"$table\" $where", sql_params)))
    totals = Dict{String,Float64}()
    for (i, c) in enumerate(sample_cols)
        v = row[Symbol("s$i")]
        totals[c] = ismissing(v) ? 0.0 : Float64(v)
    end
    totals
end

"""
    _retain_sample_columns(con, table, sample_cols, bounds, where, sql_params) -> Vector{String}

The subset of `sample_cols` whose read total falls within `bounds`, in the input
order. Returns `sample_cols` unchanged when no bound is set, so the common path
costs no query.
"""
function _retain_sample_columns(con, table::String, sample_cols::Vector{String},
                                bounds, where::String="", sql_params::Vector=Any[])
    (_has_sample_read_bounds(bounds) && !isempty(sample_cols)) || return sample_cols
    # "raw" measures the library size, so the row filters must not narrow it.
    (w, p) = bounds.basis == "raw" ? ("", Any[]) : (where, sql_params)
    totals = _sample_read_totals(con, table, sample_cols, w, p)
    filter(sample_cols) do c
        t = get(totals, c, 0.0)
        (isnothing(bounds.min) || t >= bounds.min) &&
        (isnothing(bounds.max) || t <= bounds.max)
    end
end

# Convenience wrapper for callers that hold the raw params dict.
_retain_sample_columns(con, table::String, sample_cols::Vector{String},
                       params::Dict{String,String}, where::String, sql_params::Vector) =
    _retain_sample_columns(con, table, sample_cols,
                           _sample_read_bounds(params), where, sql_params)

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

## Deterministic row order
# Identity columns, in the order they are preferred as a tiebreak/fallback sort
# key. Every merged or annotation table carries at least one of them.
const _IDENTITY_COLUMNS = ("SeqName", "OTU", "ASV", "sequence")

"""
    _stable_order_clause(columns) -> String

An ORDER BY on the table's identity column, for queries that would otherwise
return rows in whatever order the engine produced them. DuckDB parallelises
scans and aggregations, so an unordered SELECT is free to hand back rows in a
different order between two runs over the same data; anything that pages
through, exports, or indexes those rows positionally then varies run to run.
Falls back to `rowid` (every DuckDB base table has one) when the table carries
no identity column, and to no ordering at all for a relation that has neither.
"""
function _stable_order_clause(columns::Vector{String})
    for c in _IDENTITY_COLUMNS
        c in columns && return "ORDER BY \"$c\""
    end
    isempty(columns) ? "" : "ORDER BY rowid"
end

# The requested ordering, made total. A user sort on a column with duplicate
# values leaves the tied rows in engine order, so the identity column is
# appended as a final tiebreak: without it, two pages of a table sorted by, say,
# Genus can overlap or skip rows between requests.
function _total_order_clause(sort_by::Union{String,Nothing}, sort_dir::String,
                             sort_columns::Vector{String}, all_columns::Vector{String})
    primary = _order_clause(sort_by, sort_dir, sort_columns)
    stable  = _stable_order_clause(all_columns)
    isempty(primary) && return stable
    isempty(stable)  && return primary
    # `stable` is itself a complete ORDER BY; splice its key list onto the
    # primary clause as a secondary key.
    primary * ", " * stable[length("ORDER BY ") + 1:end]
end

function _duckdb_rows(con, sql::String, params::Vector=[])
    result = DBInterface.execute(con, sql, params)
    df = DataFrame(result)
    colnames = names(df)
    [OrderedDict(c => (ismissing(row[c]) ? nothing : row[c]) for c in colnames)
     for row in eachrow(df)]
end

## Distinct-value query (numeric or text)
function _duckdb_distinct(con, table::String, column::String, body;
                          extra_guard::Union{Function,Nothing}=nothing)
    params = _body_filter_params(body)
    columns = _duckdb_columns(con, table)

    # Optional pre-check (e.g. for BLAST Assignment column that may not exist yet)
    if !isnothing(extra_guard)
        result = extra_guard(column, columns)
        !isnothing(result) && return result
    end

    column in columns || return json_error(404, "column_not_found",
                                           "Column '$column' not found")

    (where, sql_params) = _build_where(params, columns)
    num_sql = "SELECT MIN(TRY_CAST(\"$column\" AS DOUBLE)) AS mn,
                      MAX(TRY_CAST(\"$column\" AS DOUBLE)) AS mx,
                      COUNT(TRY_CAST(\"$column\" AS DOUBLE)) AS cnt,
                      SUM(TRY_CAST(\"$column\" AS DOUBLE)) AS sm,
                      AVG(TRY_CAST(\"$column\" AS DOUBLE)) AS avg,
                      MEDIAN(TRY_CAST(\"$column\" AS DOUBLE)) AS med,
                      QUANTILE_CONT(TRY_CAST(\"$column\" AS DOUBLE), 0.25) AS q1,
                      QUANTILE_CONT(TRY_CAST(\"$column\" AS DOUBLE), 0.75) AS q3
               FROM \"$table\" $where"
    num_row = only(DataFrame(DBInterface.execute(con, num_sql, sql_params)))

    if !ismissing(num_row.mn) && num_row.cnt > 0
        return json((; column, type="numeric", min=num_row.mn, max=num_row.mx,
                      count=num_row.cnt, sum=num_row.sm, mean=num_row.avg,
                      median=num_row.med, q1=num_row.q1, q3=num_row.q3))
    end

    not_null = "\"$column\" IS NOT NULL"
    full_where = isempty(where) ? "WHERE $not_null" : "$where AND $not_null"
    limit_n = 1000
    vals_sql = "SELECT DISTINCT CAST(\"$column\" AS VARCHAR) AS val
                FROM \"$table\" $full_where
                ORDER BY val
                LIMIT $(limit_n + 1)"
    vals_df = DataFrame(DBInterface.execute(con, vals_sql, sql_params))
    truncated = nrow(vals_df) > limit_n
    vals = String[string(row.val) for row in eachrow(vals_df[1:min(nrow(vals_df), limit_n), :])]
    json((; column, type="text", values=vals, count=length(vals), truncated))
end

## Paginated table query
function _duckdb_paginated_query(con, table::String, body;
                                 select_expr::String="*",
                                 response_columns::Union{Vector{String},Nothing}=nothing)
    params = _body_filter_params(body)
    page = max(1, Int(get(body, :page, 1)))
    per_page = clamp(Int(get(body, :perPage, 100)), 1, 10_000)
    sort_by = get(params, "sort", nothing)
    sort_dir = get(params, "sort_dir", "asc")

    columns = _duckdb_columns(con, table)
    isempty(columns) && return json_error(404, "table_not_found",
                                          "Table '$table' not found")

    total_unfiltered = only(DataFrame(DBInterface.execute(
        con, "SELECT COUNT(*) AS n FROM \"$table\""))).n

    (where, sql_params) = _build_where(params, columns)
    total = only(DataFrame(DBInterface.execute(
        con, "SELECT COUNT(*) AS n FROM \"$table\" $where", sql_params))).n

    all_count_cols = _sample_count_columns(con, table)
    total_reads_unfiltered = _sum_reads(con, table, all_count_cols)

    # Sample read-count bounds drop whole sample columns, so they are resolved
    # before the reads are totalled and before the rows are selected: the page
    # must not carry a column the filter excluded, nor count its reads.
    count_cols = _retain_sample_columns(con, table, all_count_cols, params,
                                        where, sql_params)
    dropped = setdiff(all_count_cols, count_cols)
    total_reads = _sum_reads(con, table, count_cols, where, sql_params)
    count_max = _count_maxima(con, table, count_cols, where, sql_params)

    dropped_set = Set(dropped)
    out_columns = isnothing(response_columns) ? columns : response_columns
    isempty(dropped) || (out_columns = filter(c -> !(c in dropped_set), out_columns))

    # Both select expressions in use start with `*`; narrowing it to an explicit
    # column list keeps the excluded samples out of the result set rather than
    # fetching and discarding them.
    select_sql = if isempty(dropped) || !startswith(select_expr, "*")
        select_expr
    else
        kept = join(["\"$c\"" for c in columns if !(c in dropped_set)], ", ")
        kept * select_expr[2:end]
    end

    order = _total_order_clause(sort_by, sort_dir, out_columns, columns)
    offset = (page - 1) * per_page
    rows = _duckdb_rows(con,
        "SELECT $select_sql FROM \"$table\" $where $order LIMIT $per_page OFFSET $offset",
        sql_params)

    json((; total, total_unfiltered, total_reads, total_reads_unfiltered, page,
            per_page, columns=out_columns, sample_count_columns=count_cols, count_max,
            excluded_samples=dropped, rows))
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

    # Sample read-count bounds select columns, so they sit beside colFilters:
    # `{min, max, basis}` where basis is "filtered" (default) or "raw".
    # See `_sample_read_bounds`.
    sample_reads = get(body, :sampleReads, nothing)
    if !isnothing(sample_reads)
        smin = get(sample_reads, :min, nothing)
        !isnothing(smin) && (params["sample_min_reads"] = string(smin))
        smax = get(sample_reads, :max, nothing)
        !isnothing(smax) && (params["sample_max_reads"] = string(smax))
        basis = get(sample_reads, :basis, nothing)
        !isnothing(basis) && (params["sample_reads_basis"] = string(basis))
    end

    col_filters = get(body, :colFilters, nothing)
    isnothing(col_filters) && return params

    for (col, f) in pairs(col_filters)
        col_str = string(col)
        # The UI carries the sample read-count bounds as a reserved entry in
        # colFilters, so they travel with every consumer of that record
        # (presets, chart requests, save, export) without separate plumbing.
        # An explicit top-level `sampleReads` takes precedence.
        if col_str == SAMPLE_READS_FILTER_KEY
            for (field, key) in ((:min, "sample_min_reads"), (:max, "sample_max_reads"),
                                 (:basis, "sample_reads_basis"))
                v = get(f, field, nothing)
                (isnothing(v) || haskey(params, key)) && continue
                params[key] = string(v)
            end
            continue
        end
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
