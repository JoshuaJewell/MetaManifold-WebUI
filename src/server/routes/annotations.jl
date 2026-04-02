# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: functional annotation (FuncDB) for listing, generation, querying,
# distinct-value lookup, and contamination tagging.

using JSON3, CSV, DataFrames, OrderedCollections, DuckDB, DBInterface, Dates

include(joinpath(@__DIR__, "duckdb_helpers.jl"))

const FUNCDB_PATH = joinpath(dirname(dirname(dirname(@__DIR__))), "databases", "FuncDB_species.csv")

const _ANNOTATION_REQUIRED_COLS = Dict(
    "VSEARCH" => ("Species", "Genus"),
    "DADA2"   => ("Species_dada2", "Genus_dada2"),
)

function _validate_source(source::String)
    source in ("VSEARCH", "DADA2") || return false
    true
end

function _annotation_dir(study::String, run::String, source::String;
                         group::Union{String,Nothing}=nothing)
    joinpath(_run_project_dir(study, run; group), "annotated", source)
end

function _annotation_db_path(study::String, run::String, source::String;
                             group::Union{String,Nothing}=nothing)
    joinpath(_annotation_dir(study, run, source; group), "annotations.duckdb")
end

function _require_study_run_source(study::String, run::String, source::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                 "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    _validate_source(source) || return json_error(400, "invalid_source",
                                                  "Source must be 'VSEARCH' or 'DADA2'")
    nothing
end

function _with_annotation_db(f, study::String, run::String, source::String;
                             group::Union{String,Nothing}=nothing, readonly::Bool=true)
    ann_db_path = _annotation_db_path(study, run, source; group)
    isfile(ann_db_path) || return json_error(404, "no_annotations",
                                             "No annotation database found for $source")

    db = DuckDB.DB(ann_db_path; readonly)
    con = DBInterface.connect(db)
    try
        f(con)
    finally
        DBInterface.close!(con)
        close(db)
    end
end

function _annotation_status(csv_path::String, merged_csv_path::String, funcdb_path::String)
    isfile(csv_path) || return "missing"
    ann_mtime = mtime(csv_path)
    merged_mtime = isfile(merged_csv_path) ? mtime(merged_csv_path) : 0.0
    funcdb_mtime = isfile(funcdb_path) ? mtime(funcdb_path) : 0.0
    (ann_mtime < merged_mtime || ann_mtime < funcdb_mtime) ? "stale" : "fresh"
end

function _contamination_stats(con, table::String)
    count_cols = _sample_count_columns(con, table)
    sum_expr = isempty(count_cols) ? "0" :
               join(["COALESCE(SUM(\"$c\"), 0)" for c in count_cols], " + ")

    sql = """SELECT "Contamination" AS status,
                    COUNT(*) AS nrows,
                    ($(sum_expr)) AS nreads
             FROM "$(table)"
             GROUP BY "Contamination" """

    df = DataFrame(DBInterface.execute(con, sql))
    tally = Dict{String,NamedTuple}()

    for row in eachrow(df)
        status = string(coalesce(row.status, "unassigned"))
        tally[status] = (; rows=Int(coalesce(row.nrows, 0)), reads=Int(coalesce(row.nreads, 0)))
    end

    for status in ("yes", "no", "unassigned")
        haskey(tally, status) || (tally[status] = (; rows=0, reads=0))
    end

    total_rows = sum(bucket.rows for bucket in values(tally))
    total_reads = sum(bucket.reads for bucket in values(tally))

    (; yes=tally["yes"], no=tally["no"], unassigned=tally["unassigned"],
       total=(; rows=total_rows, reads=total_reads))
end

function _annotation_catalog(study::String, run::String, source::String;
                             group::Union{String,Nothing}=nothing)
    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return Any[]

    ann_dir = _annotation_dir(study, run, source; group)
    ann_db_path = _annotation_db_path(study, run, source; group)

    ann_table_rows = Dict{String,Int}()
    if isfile(ann_db_path)
        _with_annotation_db(study, run, source; group) do ann_con
            ann_tables = Set(string.(DataFrame(DBInterface.execute(ann_con, "SHOW TABLES")).name))
            for tname in ann_tables
                ann_table_rows[tname] = only(DataFrame(DBInterface.execute(
                    ann_con, "SELECT COUNT(*) AS n FROM \"$tname\""))).n
            end
        end
    end

    with_results_db(dir) do con
        table_names = string.(DataFrame(DBInterface.execute(con, "SHOW TABLES")).name)
        catalog = Any[]

        for tname in table_names
            cols = _duckdb_columns(con, tname)
            col_set = Set(cols)
            all(c -> c in col_set, _ANNOTATION_REQUIRED_COLS[source]) || continue

            csv_path = joinpath(ann_dir, tname * ".csv")
            merged_csv_path = joinpath(dir, tname * ".csv")
            status = _annotation_status(csv_path, merged_csv_path, FUNCDB_PATH)
            rows = get(ann_table_rows, tname, nothing)
            generated_at = status == "missing" ? nothing : unix2datetime(mtime(csv_path))

            push!(catalog, (; table=tname, source, status, rows, generated_at))
        end

        catalog
    end
end

function _load_source_table(study::String, run::String, table::String;
                            group::Union{String,Nothing}=nothing)
    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return nothing

    with_results_db(dir) do con
        columns = _duckdb_columns(con, table)
        isempty(columns) && return nothing
        DataFrame(DBInterface.execute(con, "SELECT * FROM \"$table\""))
    end
end

function _write_annotation_output!(study::String, run::String, source::String,
                                   table::String, df::DataFrame;
                                   group::Union{String,Nothing}=nothing)
    ann_dir = _annotation_dir(study, run, source; group)
    mkpath(ann_dir)

    csv_path = joinpath(ann_dir, table * ".csv")
    CSV.write(csv_path, df)

    db = DuckDB.DB(_annotation_db_path(study, run, source; group))
    con = DBInterface.connect(db)
    try
        DBInterface.execute(con,
            "CREATE OR REPLACE TABLE \"$(table)\" AS SELECT * FROM read_csv_auto('$(csv_path)', auto_detect=true)")
    finally
        DBInterface.close!(con)
        close(db)
    end

    csv_path
end

function _annotation_response(table::String, source::String, output_path::String, rows::Int;
                              contamination_stats=nothing)
    generated_at = unix2datetime(mtime(output_path))
    if isnothing(contamination_stats)
        (; table, source, status="fresh", rows, output_path, generated_at)
    else
        (; table, source, status="fresh", rows, output_path, generated_at,
           contamination_stats)
    end
end

function _body_list(body_obj, field::String, base_cfg, cfg_key::String)
    raw = get(body_obj, Symbol(field), nothing)
    if raw isa AbstractVector
        return String[string(v) for v in raw]
    end

    entry = get(base_cfg, cfg_key, nothing)
    if isnothing(entry) || !(entry.value isa AbstractVector)
        return String[]
    end

    String[string(v) for v in entry.value]
end

function _build_synthetic_contamination_cfg(body, base_cfg)
    bl_fields = get(body, :blacklist, (;))
    wl_fields = get(body, :whitelist, (;))

    Dict{String,Any}(
        "annotation.contamination.blacklist.function" =>
            (; value=_body_list(bl_fields, "function", base_cfg,
                                "annotation.contamination.blacklist.function"), source="run"),
        "annotation.contamination.blacklist.detailed_function" =>
            (; value=_body_list(bl_fields, "detailed_function", base_cfg,
                                "annotation.contamination.blacklist.detailed_function"), source="run"),
        "annotation.contamination.blacklist.associated_organism" =>
            (; value=_body_list(bl_fields, "associated_organism", base_cfg,
                                "annotation.contamination.blacklist.associated_organism"), source="run"),
        "annotation.contamination.blacklist.associated_material" =>
            (; value=_body_list(bl_fields, "associated_material", base_cfg,
                                "annotation.contamination.blacklist.associated_material"), source="run"),
        "annotation.contamination.blacklist.environment" =>
            (; value=_body_list(bl_fields, "environment", base_cfg,
                                "annotation.contamination.blacklist.environment"), source="run"),
        "annotation.contamination.whitelist.function" =>
            (; value=_body_list(wl_fields, "function", base_cfg,
                                "annotation.contamination.whitelist.function"), source="run"),
        "annotation.contamination.whitelist.detailed_function" =>
            (; value=_body_list(wl_fields, "detailed_function", base_cfg,
                                "annotation.contamination.whitelist.detailed_function"), source="run"),
        "annotation.contamination.whitelist.associated_organism" =>
            (; value=_body_list(wl_fields, "associated_organism", base_cfg,
                                "annotation.contamination.whitelist.associated_organism"), source="run"),
        "annotation.contamination.whitelist.associated_material" =>
            (; value=_body_list(wl_fields, "associated_material", base_cfg,
                                "annotation.contamination.whitelist.associated_material"), source="run"),
        "annotation.contamination.whitelist.environment" =>
            (; value=_body_list(wl_fields, "environment", base_cfg,
                                "annotation.contamination.whitelist.environment"), source="run"),
    )
end

function _persist_contamination_config!(study::String, run::String, group, synthetic_cfg)
    run_path = isnothing(group) ? run : joinpath(group, run)
    config_path = joinpath(ServerState.data_dir(), study, run_path, "pipeline.yml")

    for (cfg_key, entry) in synthetic_cfg
        _write_override(config_path, cfg_key, entry.value)
    end
end

function _apply_contamination_filter!(study::String, run::String, source::String, table::String,
                                      body; group::Union{String,Nothing}=nothing)
    base_cfg = _resolve_config(study, run, group)
    synthetic_cfg = _build_synthetic_contamination_cfg(body, base_cfg)
    rows = 0
    output_path = ""
    contam_stats = nothing

    result = _with_annotation_db(study, run, source; group, readonly=false) do con
        columns = _duckdb_columns(con, table)
        isempty(columns) && return json_error(404, "table_not_found",
                                              "Table '$table' not found in annotation database")
        "Contamination" in columns || return json_error(400, "no_contamination_column",
            "Table '$table' has no Contamination column - regenerate full annotation")

        df = DataFrame(DBInterface.execute(con, "SELECT * FROM \"$table\""))
        df[!, "Contamination"] .= "unassigned"
        FuncDBAnnotation.apply_contamination_filter!(df, synthetic_cfg)

        output_path = _write_annotation_output!(study, run, source, table, df; group)
        rows = nrow(df)
        contam_stats = _contamination_stats(con, table)
        nothing
    end
    result isa HTTP.Response && return result

    try
        _persist_contamination_config!(study, run, group, synthetic_cfg)
    catch e
        return json_error(500, "config_persist_failed",
            "Contamination filter applied, but saving pipeline.yml failed: $(sprint(showerror, e))")
    end

    @info "Contamination filter applied" table source rows
    json(_annotation_response(table, source, output_path, rows; contamination_stats=contam_stats))
end

function _generate_annotation!(study::String, run::String, source::String, table::String;
                               group::Union{String,Nothing}=nothing)
    source_df = _load_source_table(study, run, table; group)
    isnothing(source_df) && return json_error(404, "table_not_found",
                                              "Table '$table' not found in results database")

    col_set = Set(names(source_df))
    for col in _ANNOTATION_REQUIRED_COLS[source]
        col in col_set || return json_error(400, "missing_taxonomy_columns",
            "Table '$table' is missing required column '$col' for $source annotation")
    end

    annotated_df = try
        FuncDBAnnotation.annotate_table(source_df, source, FUNCDB_PATH)
    catch e
        @error "Annotation failed" exception=(e, catch_backtrace())
        return json_error(500, "annotation_failed", "Annotation failed: $(sprint(showerror, e))")
    end

    try
        cfg = _resolve_config(study, run, group)
        FuncDBAnnotation.apply_contamination_filter!(annotated_df, cfg)
    catch e
        @warn "Contamination filter failed, skipping" exception=(e, catch_backtrace())
    end

    output_path = _write_annotation_output!(study, run, source, table, annotated_df; group)
    rows = nrow(annotated_df)
    @info "Generated annotation" table source rows output_path
    json(_annotation_response(table, source, output_path, rows))
end

function _annotation_query(con, table::String, body)
    params = _body_filter_params(body)
    page = max(1, Int(get(body, :page, 1)))
    per_page = clamp(Int(get(body, :perPage, 100)), 1, 10_000)
    sort_by = get(params, "sort", nothing)
    sort_dir = get(params, "sort_dir", "asc")

    columns = _duckdb_columns(con, table)
    isempty(columns) && return json_error(404, "table_not_found",
                                          "Table '$table' not found in annotation database")

    total_unfiltered = only(DataFrame(DBInterface.execute(
        con, "SELECT COUNT(*) AS n FROM \"$table\""))).n

    (where, sql_params) = _build_where(params, columns)
    total = only(DataFrame(DBInterface.execute(
        con, "SELECT COUNT(*) AS n FROM \"$table\" $where", sql_params))).n

    count_cols = _sample_count_columns(con, table)
    total_reads_unfiltered = _sum_reads(con, table, count_cols)
    total_reads = _sum_reads(con, table, count_cols, where, sql_params)
    order = _order_clause(sort_by, sort_dir, columns)
    offset = (page - 1) * per_page
    rows = _duckdb_rows(con,
        "SELECT * FROM \"$table\" $where $order LIMIT $per_page OFFSET $offset",
        sql_params)

    json((; total, total_unfiltered, total_reads, total_reads_unfiltered, page,
            per_page, columns, sample_count_columns=count_cols, rows))
end

function _annotation_distinct(con, table::String, column::String, body)
    params = _body_filter_params(body)
    columns = _duckdb_columns(con, table)
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

function _matched_taxon_column(source::String, rank::String)
    rank_entry = findfirst(r -> r.rank == rank, FuncDBAnnotation.RANK_HIERARCHY)
    isnothing(rank_entry) && return nothing
    entry = FuncDBAnnotation.RANK_HIERARCHY[rank_entry]
    source == "VSEARCH" ? entry.vsearch : entry.dada2
end

function _sync_annotation_csv!(csv_path::String, taxon_col::String, rank::String,
                               taxon::String, status::String)
    isfile(csv_path) || return

    csv_df = CSV.read(csv_path, DataFrame; stringtype=String)
    if !("Contamination" in names(csv_df)) || !(taxon_col in names(csv_df))
        return
    end

    mask = csv_df.funcdb_match_rank .== rank .&&
           lowercase.(strip.(csv_df[!, taxon_col])) .== lowercase(strip(taxon))
    csv_df.Contamination[mask] .= status
    CSV.write(csv_path, csv_df)
end

## Annotation status
@get "/api/v1/studies/{study}/runs/{run}/annotations/{source}" function(req,
                                                                        study::String,
                                                                        run::String,
                                                                        source::String)
    err = _require_study_run_source(study, run, source)
    !isnothing(err) && return err
    json(_annotation_catalog(study, run, source; group=_req_group(req)))
end

## Annotation generation
@post "/api/v1/studies/{study}/runs/{run}/annotations/{source}/generate" function(req,
                                                                                   study::String,
                                                                                   run::String,
                                                                                   source::String)
    err = _require_study_run_source(study, run, source)
    !isnothing(err) && return err
    isfile(FUNCDB_PATH) || return json_error(400, "funcdb_missing",
        "FuncDB database not found at expected path")

    group = _req_group(req)
    dir = _require_duckdb(study, run; group)
    isnothing(dir) && return json_error(404, "no_results",
        "No results database for run '$run' - run the pipeline first")

    body = JSON3.read(String(req.body))
    table = get(body, :table, nothing)
    isnothing(table) && return json_error(400, "missing_table", "Body must include 'table'")
    table = string(table)

    if get(body, :contamination_only, false) === true
        return _apply_contamination_filter!(study, run, source, table, body; group)
    end

    _generate_annotation!(study, run, source, table; group)
end

## Annotation table query
@post "/api/v1/studies/{study}/runs/{run}/annotations/{source}/{table}/query" function(req,
                                                                                        study::String,
                                                                                        run::String,
                                                                                        source::String,
                                                                                        table::String)
    err = _require_study_run_source(study, run, source)
    !isnothing(err) && return err

    body = JSON3.read(String(req.body))
    _with_annotation_db(study, run, source; group=_req_group(req)) do con
        _annotation_query(con, table, body)
    end
end

## Annotation distinct values
@post "/api/v1/studies/{study}/runs/{run}/annotations/{source}/{table}/distinct/{column}" function(req,
                                                                                                    study::String,
                                                                                                    run::String,
                                                                                                    source::String,
                                                                                                    table::String,
                                                                                                    column::String)
    err = _require_study_run_source(study, run, source)
    !isnothing(err) && return err

    body = length(req.body) > 0 ? JSON3.read(String(req.body)) : (;)
    _with_annotation_db(study, run, source; group=_req_group(req)) do con
        _annotation_distinct(con, table, column, body)
    end
end

## FuncDB entry creation
@post "/api/v1/funcdb/entries" function(req)
    isfile(FUNCDB_PATH) || return json_error(400, "funcdb_missing",
        "FuncDB database not found at expected path")

    body = JSON3.read(String(req.body))
    entry = Dict{String, String}()
    for (k, v) in pairs(body)
        k == :modified_by && continue
        entry[String(k)] = string(v)
    end

    modified_by = string(get(body, :modified_by, ""))

    row = try
        FuncDBAnnotation.append_funcdb_entry(FUNCDB_PATH, entry; modified_by)
    catch e
        return json_error(400, "invalid_entry", sprint(showerror, e))
    end

    json(row)
end

## Contamination update
# Sets Contamination for every row that shares the same matched rank and taxon.
@patch "/api/v1/studies/{study}/runs/{run}/annotations/{source}/{table}/contamination" function(req,
                                                                                                 study::String,
                                                                                                 run::String,
                                                                                                 source::String,
                                                                                                 table::String)
    err = _require_study_run_source(study, run, source)
    !isnothing(err) && return err

    body = JSON3.read(String(req.body))
    rank = string(get(body, :rank, ""))
    taxon = string(get(body, :taxon, ""))
    status = string(get(body, :status, ""))

    isempty(rank) && return json_error(400, "missing_rank", "Body must include 'rank'")
    isempty(taxon) && return json_error(400, "missing_taxon", "Body must include 'taxon'")
    status in ("unassigned", "yes", "no") ||
        return json_error(400, "invalid_status",
                          "Status must be 'unassigned', 'yes', or 'no'")

    taxon_col = _matched_taxon_column(source, rank)
    isnothing(taxon_col) && return json_error(400, "invalid_rank",
                                              "Unknown rank '$rank'")

    group = _req_group(req)
    result = _with_annotation_db(study, run, source; group, readonly=false) do con
        columns = _duckdb_columns(con, table)
        "Contamination" in columns || return json_error(400, "no_contamination_column",
            "Table '$table' has no Contamination column - regenerate the annotation")
        taxon_col in columns || return json_error(400, "missing_taxon_column",
            "Column '$taxon_col' not found in table")

        DBInterface.execute(con,
            """UPDATE "$(table)"
               SET "Contamination" = \$1
               WHERE "funcdb_match_rank" = \$2
                 AND LOWER(TRIM("$(taxon_col)")) = LOWER(TRIM(\$3))""",
            [status, rank, taxon])

        cnt = only(DataFrame(DBInterface.execute(con,
            """SELECT COUNT(*) AS n FROM "$(table)"
               WHERE "funcdb_match_rank" = \$1
                 AND LOWER(TRIM("$(taxon_col)")) = LOWER(TRIM(\$2))""",
            [rank, taxon]))).n

        (; rows_affected=cnt)
    end
    result isa HTTP.Response && return result

    _sync_annotation_csv!(
        joinpath(_annotation_dir(study, run, source; group), table * ".csv"),
        taxon_col, rank, taxon, status)

    @info "Contamination updated" table rank taxon status rows_affected=result.rows_affected
    json((; table, rank, taxon, status, rows_affected=result.rows_affected))
end

## Contamination stats
@get "/api/v1/studies/{study}/runs/{run}/annotations/{source}/{table}/contamination/stats" function(req,
                                                                                                     study::String,
                                                                                                     run::String,
                                                                                                     source::String,
                                                                                                     table::String)
    err = _require_study_run_source(study, run, source)
    !isnothing(err) && return err

    _with_annotation_db(study, run, source; group=_req_group(req)) do con
        columns = _duckdb_columns(con, table)
        isempty(columns) && return json_error(404, "table_not_found",
                                              "Table '$table' not found")
        "Contamination" in columns || return json_error(400, "no_contamination_column",
            "Table '$table' has no Contamination column")
        json(_contamination_stats(con, table))
    end
end
