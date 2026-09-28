# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: a publication table built to the request (rows of taxa or categories,
# columns of runs, sub-groups or samples, and the chosen values), served as JSON
# for preview, as a formatted .xlsx workbook, or as CSV.

using CSV, DataFrames, XLSX, DuckDB, DBInterface
using ..Categories

## Table model
# A table is a Dict with:
#   title, notes
#   header_rows: rows of spanning group headers above the column labels,
#                each a vector of Dict("label", "span")
#   columns:     Dict("label", "kind") with kind in label | int | pct | pp
#   rows, footer: vectors of cell vectors; `nothing` marks an empty cell.
# Footer rows are set off from the body by a rule.

const _PUB_UNASSIGNED   = "Unassigned"
const _PUB_UNCLASSIFIED = "Unclassified"

_pub_col(label, kind) = Dict{String,Any}("label" => label, "kind" => kind)
_pub_span(label, span) = Dict{String,Any}("label" => label, "span" => span)
_pub_display(s::AbstractString) = replace(s, "_" => " ")
_pub_pct(part, total) = total > 0 ? 100 * part / total : nothing

_pub_total_reads(tally) = sum((v[1] for v in values(tally)); init=0)
_pub_total_asvs(tally)  = sum((v[2] for v in values(tally)); init=0)

const _PUB_VALUES = (asvs = ("ASVs", "int"), reads = ("Reads", "int"), pct = ("%", "pct"))

"""
    _pub_table(; title, first_col, groups, nested, order, values, diff, notes)

One row per label in `order`. `groups` are `(; label, columns)`, each column
`(; label, n, tally)` with `tally` mapping label => (reads, asvs). Each column
shows `values` (a subset of :asvs, :reads, :pct), percentages being of that
column's reads. With `nested`, group labels head their columns in a row of
their own; otherwise each group is one column headed by its label. With
`diff`, a group of exactly two columns gains their difference in percentage
points (second minus first).
"""
function _pub_table(; title, first_col, groups, nested, order, values, diff, notes)
    k = length(values)
    heading(c) = c.n == 1 ? _pub_display(c.label) : "$(_pub_display(c.label)) (n = $(c.n))"
    has_diff(g) = nested && diff && :pct in values && length(g.columns) == 2
    top, sub = [_pub_span("", 1)], [_pub_span("", 1)]
    cols = [_pub_col(first_col, "label")]
    for g in groups
        nested && push!(top, _pub_span(_pub_display(g.label), k * length(g.columns) + (has_diff(g) ? 1 : 0)))
        for c in g.columns
            push!(sub, _pub_span(heading(c), k))
            append!(cols, [_pub_col(_PUB_VALUES[v]...) for v in values])
        end
        if has_diff(g)
            push!(sub, _pub_span("", 1))
            push!(cols, _pub_col("Δ (pp)", "pp"))
        end
    end
    totals = Dict(c => _pub_total_reads(c.tally) for g in groups for c in g.columns)
    rows = map(order) do lab
        cells = Any[lab]
        for g in groups
            pcts = Any[]
            for c in g.columns
                reads, asvs = get(c.tally, lab, (0, 0))
                p = _pub_pct(reads, totals[c])
                push!(pcts, p)
                for v in values
                    push!(cells, v == :asvs ? asvs : v == :reads ? reads : p)
                end
            end
            has_diff(g) && push!(cells, any(isnothing, pcts) ? nothing : pcts[2] - pcts[1])
        end
        cells
    end
    footer = Any["Total"]
    for g in groups
        for c in g.columns, v in values
            t = totals[c]
            push!(footer, v == :asvs ? _pub_total_asvs(c.tally) : v == :reads ? t : (t > 0 ? 100.0 : nothing))
        end
        has_diff(g) && push!(footer, nothing)
    end
    Dict{String,Any}("title" => title,
                     "header_rows" => nested ? [top, sub] : [sub], "columns" => cols,
                     "rows" => rows, "footer" => [footer], "notes" => notes)
end

# Category set order first, then any other labels by total reads, Unassigned last.
function _pub_category_order(set_categories::Vector{String}, tallies)
    present = Dict{String,Int}()
    for t in tallies, (lab, (reads, _)) in t
        present[lab] = get(present, lab, 0) + reads
    end
    extra = sort([l for l in keys(present) if !(l in set_categories) && l != _PUB_UNASSIGNED];
                 by = l -> (-present[l], l))
    order = vcat(set_categories, extra)
    get(present, _PUB_UNASSIGNED, 0) > 0 && push!(order, _PUB_UNASSIGNED)
    order
end

# Taxa by mean relative abundance across units, Unclassified last. Labels with
# no reads in any unit are dropped.
function _pub_taxon_order(tallies)
    totals = [_pub_total_reads(t) for t in tallies]
    labels = Set{String}()
    for t in tallies, (lab, (reads, _)) in t
        reads > 0 && push!(labels, lab)
    end
    mean_pct(l) = sum((something(_pub_pct(get(t, l, (0, 0))[1], tot), 0.0)
                       for (t, tot) in zip(tallies, totals)); init=0.0) / max(length(tallies), 1)
    order = sort([l for l in labels if l != _PUB_UNCLASSIFIED]; by = l -> (-mean_pct(l), l))
    _PUB_UNCLASSIFIED in labels && push!(order, _PUB_UNCLASSIFIED)
    order
end

## Presentation: heatmap fills and hidden zeros
"""
    _pub_present(tbl; heatmap="none", hide_zeros=false) -> Dict

Copy of `tbl` with a `fills` matrix (hex colour or `nothing`, one per body
cell) and, with `hide_zeros`, exact zeros in the body blanked to "". Numeric
columns are graded from zero to their largest value: per column, or with
`heatmap = "table"` across all columns of the same kind, since reads and
percentages cannot share a scale. Differences get a scale centred on zero.
Totals are never graded or hidden.
"""
function _pub_present(tbl; heatmap::String="none", hide_zeros::Bool=false)
    out = copy(tbl)
    cols = tbl["columns"]
    rows = tbl["rows"]
    graded(j) = cols[j]["kind"] in ("int", "pct", "pp")
    key(j) = heatmap == "table" ? cols[j]["kind"] : string(j)
    scale = Dict{String,Float64}()
    for row in rows, (j, v) in enumerate(row)
        (graded(j) && v isa Real) || continue
        scale[key(j)] = max(get(scale, key(j), 0.0), abs(v))
    end
    fills = map(rows) do row
        map(enumerate(row)) do (j, v)
            (heatmap == "none" || !graded(j)) && return nothing
            _heat_colour(v, get(scale, key(j), 0.0); diverging=cols[j]["kind"] == "pp")
        end
    end
    out["fills"] = fills
    out["rows"] = hide_zeros ?
        [Any[(graded(j) && v isa Real && v == 0) ? "" : v for (j, v) in enumerate(row)] for row in rows] :
        rows
    out
end

## Flat and formatted renderings
# Column headings with their spanning group labels folded in, e.g. "Multiplex (n = 66) Reads".
function _pub_flat_headers(tbl)
    ncol = length(tbl["columns"])
    parts = [String[] for _ in 1:ncol]
    for hrow in tbl["header_rows"]
        c = 1
        for h in hrow
            for j in c:(c + h["span"] - 1)
                isempty(h["label"]) || push!(parts[j], h["label"])
            end
            c += h["span"]
        end
    end
    [join(vcat(parts[j], [tbl["columns"][j]["label"]]), " ") for j in 1:ncol]
end

# CSV keeps four decimals so that small but non-zero shares stay visible.
_pub_csv_cell(v, kind) = isnothing(v) ? "" :
    kind in ("pct", "pp") ? round(v; digits=4) : v

# The value a formatted cell shows. Shares that round to zero but are not zero
# show as "<0.01".
function _pub_display_value(v, kind)
    isnothing(v) && return "–"
    v isa AbstractString && return v
    kind == "pct" && 0 < v < 0.005 && return "<0.01"
    kind == "pp" && return (r = round(v; digits=2); r == 0 ? 0.0 : r)
    v
end

function _pub_csv(tbl)::String
    headers = _pub_flat_headers(tbl)
    kinds = [c["kind"] for c in tbl["columns"]]
    rows = vcat(tbl["rows"], tbl["footer"])
    df = DataFrame([h => [_pub_csv_cell(r[j], kinds[j]) for r in rows]
                    for (j, h) in enumerate(headers)]; makeunique=true)
    io = IOBuffer()
    CSV.write(io, df)
    String(take!(io))
end

const _PUB_FORMATS = Dict("int" => "#,##0", "pct" => "0.00", "pp" => "+0.00;-0.00;0.00")
const _PUB_FONT = "Times New Roman"

# Three-rule layout: a heavy rule above the headings, a light rule under them
# and above the totals, and a heavy rule closing the table.
function _pub_write_sheet!(sh, tbl)
    ncol = length(tbl["columns"])
    heavy = ["style" => "medium", "color" => "FF000000"]
    light = ["style" => "thin", "color" => "FF000000"]
    blank_row!(r) = for c in 1:ncol; sh[r, c] = ""; end

    blank_row!(1)
    sh[1, 1] = tbl["title"]
    XLSX.mergeCells(sh, 1:1, 1:ncol)
    XLSX.setFont(sh, 1, 1; bold=true, name=_PUB_FONT, size=11)
    XLSX.setAlignment(sh, 1, 1; horizontal="left", wrapText=true)
    XLSX.setRowHeight(sh, 1, 1; height=30)

    r = 2
    for hrow in tbl["header_rows"]
        blank_row!(r)
        c = 1
        for h in hrow
            sh[r, c] = h["label"]
            if h["span"] > 1
                XLSX.mergeCells(sh, r:r, c:(c + h["span"] - 1))
            end
            if !isempty(h["label"])
                # A short rule under each group label shows which columns it spans.
                XLSX.setBorder(sh, r:r, c:(c + h["span"] - 1); bottom=light)
            end
            c += h["span"]
        end
        XLSX.setAlignment(sh, r:r, 1:ncol; horizontal="center")
        r += 1
    end
    for (j, col) in enumerate(tbl["columns"])
        sh[r, j] = col["label"]
    end
    XLSX.setAlignment(sh, r:r, 2:ncol; horizontal="right")
    XLSX.setBorder(sh, 2:2, 1:ncol; top=heavy)
    XLSX.setBorder(sh, r:r, 1:ncol; bottom=light)
    head_end = r
    r += 1

    body_start = r
    for (i, row) in enumerate(vcat(tbl["rows"], tbl["footer"]))
        fills = i <= length(tbl["rows"]) ? get(tbl, "fills", nothing) : nothing
        for (j, v) in enumerate(row)
            sh[r, j] = _pub_display_value(v, tbl["columns"][j]["kind"])
            fill = isnothing(fills) ? nothing : fills[i][j]
            isnothing(fill) || _xlsx_fill!(sh, r, j, fill)
        end
        i == length(tbl["rows"]) + 1 && XLSX.setBorder(sh, r:r, 1:ncol; top=light)
        r += 1
    end
    last_row = r - 1
    if last_row >= body_start
        for (j, col) in enumerate(tbl["columns"])
            fmt = get(_PUB_FORMATS, col["kind"], nothing)
            isnothing(fmt) || XLSX.setFormat(sh, body_start:last_row, j:j; format=fmt)
        end
        XLSX.setAlignment(sh, body_start:last_row, 2:ncol; horizontal="right")
    end
    XLSX.setBorder(sh, last_row:last_row, 1:ncol; bottom=heavy)
    XLSX.setFont(sh, 2:last_row, 1:ncol; name=_PUB_FONT, size=10)
    XLSX.setFont(sh, 2:head_end, 1:ncol; bold=true, name=_PUB_FONT, size=10)

    for note in tbl["notes"]
        r += 1
        blank_row!(r)
        sh[r, 1] = note
        XLSX.mergeCells(sh, r:r, 1:ncol)
        XLSX.setFont(sh, r, 1; italic=true, name=_PUB_FONT, size=9)
        XLSX.setAlignment(sh, r, 1; horizontal="left", wrapText=true)
    end

    XLSX.setColumnWidth(sh, 1:1, 1:1; width=24)
    ncol > 1 && XLSX.setColumnWidth(sh, 1:1, 2:ncol; width=11)
    sh
end

function _pub_xlsx(tables)::Vector{UInt8}
    tmp = tempname() * ".xlsx"
    try
        XLSX.openxlsx(tmp; mode="w") do xf
            for (i, tbl) in enumerate(tables)
                name = "Table $i"
                sh = i == 1 ? xf[1] : XLSX.addsheet!(xf, name)
                i == 1 && XLSX.rename!(sh, name)
                _pub_write_sheet!(sh, tbl)
            end
        end
        read(tmp)
    finally
        isfile(tmp) && rm(tmp)
    end
end

## Data collection
# Per-label (reads, ASVs) over `scols`. An ASV counts once it has any read in
# the scope. Blank labels are folded into `missing_label`.
function _pub_tally(con, table::String, label_col::String, scols::Vector{String},
                    missing_label::String)
    isempty(scols) && return Dict{String,Tuple{Int,Int}}()
    s = join(["COALESCE(\"$c\", 0)" for c in scols], " + ")
    sql = """
        SELECT COALESCE(NULLIF(TRIM(CAST("$label_col" AS VARCHAR)), ''), ?) AS label,
               SUM($s) AS reads,
               COUNT(*) FILTER (WHERE ($s) > 0) AS asvs
        FROM "$table" GROUP BY 1
    """
    df = DataFrame(DBInterface.execute(con, sql, [missing_label]))
    Dict{String,Tuple{Int,Int}}(String(r.label) => (round(Int, coalesce(r.reads, 0)), Int(r.asvs))
                                for r in eachrow(df))
end

# Collapse the request's run specs into one unit per (group, run), keeping any
# sub-group prefixes in the order given.
function _pub_units(runs_spec)
    units = NamedTuple[]
    index = Dict{Tuple{Any,String},Int}()
    for spec in runs_spec
        run = string(get(spec, :run, ""))
        isempty(run) && continue
        group = _opt_string(spec, :group)
        prefix = _opt_string(spec, :prefix)
        i = get!(index, (group, run)) do
            push!(units, (; run, group, prefixes=String[]))
            length(units)
        end
        !isnothing(prefix) && !(prefix in units[i].prefixes) && push!(units[i].prefixes, prefix)
    end
    units
end

_pub_error(status, code, msg) = (; error=json_error(status, code, msg))

function _pub_category_column!(study, resolved, table, set_name)
    col = _category_column(set_name)
    has(cols) = col in cols
    present = _with_resolved_results_table(resolved, table) do _, columns
        has(columns)
    end
    isnothing(present) && return nothing
    present && return col
    # Older runs may predate this set's tags; backfill them as the charts do.
    _with_resolved_results_table(resolved, table; readonly=false) do con, _
        Categories.ensure_columns!(con, table,
                                   _tagging_source(study, resolved.run; group=resolved.group),
                                   [set_name]; library=_chart_library(), suffixed=_suffixed())
    end
    ok = _with_resolved_results_table(resolved, table) do _, columns
        has(columns)
    end
    ok === true ? col : ""
end

const _PUB_COLUMNS = ("run", "subgroup", "sample")

"""
    _publication_table(study, body) -> Dict or (; error)

The table `body` describes: `rows` ("rank" with `rank`, or "category" with
`category_set`), `columns` (one of $(join(_PUB_COLUMNS, ", "))), `values` (any
of asvs, reads, pct), `difference`, `table` and an optional `title`, over the
runs in `body`.
"""
function _publication_table(study::String, body)
    units = _pub_units(get(body, :runs, []))
    isempty(units) && return _pub_error(400, "no_runs", "Provide at least one run")
    table   = string(get(body, :table, "merged"))
    by_cat  = string(get(body, :rows, "rank")) == "category"
    rank    = string(get(body, :rank, "Genus"))
    columns = string(get(body, :columns, "run"))
    columns in _PUB_COLUMNS || return _pub_error(400, "bad_columns",
        "columns must be one of $(join(_PUB_COLUMNS, ", "))")
    asked  = Set(string.(get(body, :values, ["reads", "pct"])))
    values = [v for v in (:asvs, :reads, :pct) if string(v) in asked]
    isempty(values) && return _pub_error(400, "no_values", "Choose at least one of asvs, reads, pct")
    diff = Bool(get(body, :difference, false))

    set_name, set_label, set_categories = "", "", String[]
    if by_cat
        set_name = string(get(body, :category_set, "default"))
        set_cfg = get(_chart_library()["sets"], set_name, nothing)
        isnothing(set_cfg) && return _pub_error(404, "category_set_not_found",
                                                "Category set '$set_name' not found")
        set_label = string(get(set_cfg, "label", set_name))
        set_categories = [string(c["name"]) for c in get(set_cfg, "categories", []) if haskey(c, "name")]
    end
    missing_label = by_cat ? _PUB_UNASSIGNED : _PUB_UNCLASSIFIED

    groups = Any[]
    for u in units
        label = isnothing(u.group) ? u.run : "$(u.group)/$(u.run)"
        resolved = _resolve_run_duckdb(study, Dict(:run => u.run, :group => u.group))
        isnothing(resolved) && return _pub_error(404, "no_results", "No results for run '$label'")
        cat_col = nothing
        if by_cat
            cat_col = _pub_category_column!(study, resolved, table, set_name)
            isnothing(cat_col) && return _pub_error(404, "table_not_found",
                                                    "Table '$table' not found for run '$label'")
            isempty(cat_col) && return _pub_error(404, "category_set_not_found",
                                                  "Category set '$set_name' could not be applied to run '$label'")
        end
        before = length(groups)
        result = _with_resolved_results_table(resolved, table) do con, cols
            label_col = by_cat ? cat_col : _rank_column(cols, rank)
            isnothing(label_col) && return _pub_error(400, "bad_rank",
                                                      "Rank '$rank' not found in table '$table' for run '$label'")
            scols = _filter_by_prefix(sample_columns(con, table), u.prefixes)
            column(l, sc) = (; label=l, n=length(sc), tally=_pub_tally(con, table, label_col, sc, missing_label))
            cs = columns == "run" ? [column(label, scols)] :
                 columns == "sample" ? [column(c, [c]) for c in scols] :
                 isempty(u.prefixes) ? [column("All samples", scols)] :
                 [column(p, _filter_by_prefix(scols, p)) for p in u.prefixes]
            push!(groups, (; label, columns=cs))
            nothing
        end
        result isa NamedTuple && return result
        length(groups) == before &&
            return _pub_error(404, "table_not_found", "Table '$table' not found for run '$label'")
    end

    tallies = [c.tally for g in groups for c in g.columns]
    order = by_cat ? _pub_category_order(set_categories, tallies) : _pub_taxon_order(tallies)
    first_col = by_cat ? "Category" : rank
    per = columns == "run" ? "run" : columns == "subgroup" ? "sub-group" : "sample"
    default_title = (by_cat ? "Read composition by $(set_label) category" : "$(rank)-level read composition") *
                    " per $per, from table $(table)."
    title = strip(string(something(get(body, :title, nothing), "")))
    notes = [by_cat ?
        "Reads are assigned to categories of the $(set_label) set; reads matching no category are $(_PUB_UNASSIGNED)." :
        "Reads with no taxonomy at this rank are $(_PUB_UNCLASSIFIED); names ending _X are unresolved below the named parent."]
    :asvs in values && push!(notes, "ASVs: ASVs with at least one read in the column's samples.")
    columns != "sample" && push!(notes, "n: samples.")
    diff && :pct in values && columns != "run" && any(g -> length(g.columns) == 2, groups) &&
        push!(notes, "Δ (pp): difference in percentage points, second minus first.")
    _pub_table(; title=isempty(title) ? default_title : title, first_col, groups,
               nested=columns != "run", order, values, diff, notes)
end

_pub_filename(study, ext) = replace("$(study)_table.$ext", r"[^\w._-]" => "_")

## Publication table: JSON preview, .xlsx workbook, or CSV
@post "/api/v1/studies/{study}/analysis/publication-tables" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    body = JSON3.read(String(req.body))
    tbl = _publication_table(study, body)
    tbl isa NamedTuple && return tbl.error

    format = string(get(body, :format, "json"))
    heatmap = string(get(body, :heatmap, "none"))
    heatmap in HEATMAP_MODES || return json_error(400, "bad_heatmap",
        "heatmap must be one of $(join(HEATMAP_MODES, ", "))")
    # Fills and hidden zeros are presentation only; CSV stays plain data.
    presented = _pub_present(tbl; heatmap, hide_zeros=Bool(get(body, :hide_zeros, false)))
    if format == "json"
        return json((; table=presented))
    elseif format == "xlsx"
        return HTTP.Response(200, [
            "Content-Type" => "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            "Content-Disposition" => "attachment; filename=\"$(_pub_filename(study, "xlsx"))\"",
        ]; body=_pub_xlsx([presented]))
    elseif format == "csv"
        return HTTP.Response(200, [
            "Content-Type" => "text/csv",
            "Content-Disposition" => "attachment; filename=\"$(_pub_filename(study, "csv"))\"",
        ]; body=_pub_csv(tbl))
    end
    json_error(400, "bad_format", "format must be json, xlsx or csv")
end
