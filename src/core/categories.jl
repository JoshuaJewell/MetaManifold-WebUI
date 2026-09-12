# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Category-set loading and SQL-fragment generation for composition analysis.
# Self-contained: does not depend on the quarantined FuncDB annotation module
# or on any server-state globals. All directory paths are passed explicitly.
module Categories

using DuckDB, DBInterface, DataFrames
using ..Validation

## Rank-name translation between the VSEARCH-named merged columns and their
## DADA2 counterparts. Self-contained so categorisation does not depend on the
## quarantined FuncDB annotation module.
const _RANK_PAIRS = [
    ("Domain",     "Domain_dada2"),
    ("Supergroup", "Supergroup_dada2"),
    ("Division",   "Division_dada2"),
    ("Subdivision","Subdivision_dada2"),
    ("Class",      "Class_dada2"),
    ("Order",      "Order_dada2"),
    ("Family",     "Family_dada2"),
    ("Genus",      "Genus_dada2"),
    ("Species",    "Species_dada2"),
]

const _DADA2_SUFFIX = "_dada2"

# `Pident` is written by `merge_taxonomy_counts` for every vsearch hit and by
# nothing else, so its presence is the marker that vsearch taxonomy is in the
# table. Three layouts reach this code:
#
#   both enabled        plain ranks = VSEARCH, `<rank>_dada2` = DADA2, Pident, *_boot
#   vsearch disabled    plain ranks = DADA2,   no _dada2, no Pident, *_boot
#   dada2 classify off  plain ranks = VSEARCH, no _dada2, Pident,     no *_boot
#
# The plain rank columns therefore mean different things depending on whether
# vsearch ran, which is what these helpers exist to resolve.
const _VSEARCH_MARKER = "Pident"

"""
    vsearch_available(present) -> Bool

Whether the table carries vsearch taxonomy, i.e. whether the plain rank columns
are VSEARCH's rather than DADA2's.
"""
vsearch_available(present) = _VSEARCH_MARKER in present

"""
    rank_col(rank, source, present) -> Union{String,Nothing}

Resolve the column holding `source`'s taxonomy for `rank`, given the column
names actually `present`. Returns `nothing` when the table carries no taxonomy
from that source - callers must treat that as an error rather than reading the
other source's column, which would silently mislabel one method's calls as the
other's.
"""
function rank_col(rank::AbstractString, source::AbstractString, present)
    r = String(rank)
    if source == "VSEARCH"
        # vsearch always writes the plain names, and only writes them when it ran.
        return (vsearch_available(present) && r in present) ? r : nothing
    elseif source == "DADA2"
        d = r * _DADA2_SUFFIX
        d in present && return d
        # No suffixed column. The plain names are DADA2's only when vsearch did
        # not run; otherwise they are vsearch's and DADA2 simply is not here.
        return (!vsearch_available(present) && r in present) ? r : nothing
    end
    nothing
end

"""
    available_sources(present) -> Vector{String}

Which taxonomy sources the table actually carries, in ("VSEARCH", "DADA2")
order. Empty when the table has no taxonomy at all - which happens when vsearch
is disabled and DADA2 classification is skipped.
"""
function available_sources(present)
    ranks = first.(_RANK_PAIRS)
    filter(["VSEARCH", "DADA2"]) do src
        any(r -> !isnothing(rank_col(r, src, present)), ranks)
    end
end

# Return the column name used in the composed table for a category set.
# A named source gets its own column, so both classifiers' tags can sit in one table.
column_name(set_name::AbstractString, source::Union{Nothing,AbstractString}=nothing) =
    "Category__" * String(set_name) * (isnothing(source) ? "" : "__" * String(source))

## Filter to SQL translation

# Escape a value for embedding in a single-quoted SQL string literal: any
# embedded single-quote is doubled so the literal closes correctly.
_sql_str(s) = "'" * replace(string(s), "'" => "''") * "'"

# Build a mapping from VSEARCH column names (used in filter YAMLs) to the
# actual column names for the selected source. For VSEARCH this is identity;
# for DADA2 taxonomy columns become their _dada2 equivalents.
# Subdivision is not in RANK_HIERARCHY but exists in merged tables and is
# handled explicitly here.
#
# Pass `present` (the columns the table actually has) wherever it is known: it
# resolves each rank against the table's actual layout. Without it the mapping
# is the naive suffix rule, which names columns that need not exist.
function col_translate_map(source::String, present=nothing)
    isnothing(present) && return Dict{String,String}(
        vs => (source == "VSEARCH" ? vs : da) for (vs, da) in _RANK_PAIRS)

    m = Dict{String,String}()
    for (vs, da) in _RANK_PAIRS
        resolved = rank_col(vs, source, present)
        # Not in the table for this source: keep the naive name so the caller's
        # own "column missing" handling reports the rank it asked for, rather
        # than quietly substituting the other source's column.
        m[vs] = isnothing(resolved) ? (source == "VSEARCH" ? vs : da) : resolved
    end
    m
end

"""
    filter_to_sql_conditions(filter_config, col_set, table_alias; col_map) -> Vector{String}

Translate a filter YAML config (from config/filters/) into SQL condition
fragments. All conditions within a single filter are AND'd together.
`col_set` is a Set of available columns (for validation).
`table_alias` is the SQL alias prefix (e.g. "m" for the merged table).
Pass `""` to emit bare column references (for UPDATE SET contexts).
`col_map` translates filter column names to actual table column names
(for source-awareness).
"""
function filter_to_sql_conditions(filter_config::Dict, col_set::Set{String},
                                  table_alias::String;
                                  col_map::Dict{String,String}=Dict{String,String}())
    _xlate(c) = get(col_map, c, c)
    # Prefix a column name with the alias, or emit a bare quoted name when no alias.
    _ref(c) = isempty(table_alias) ? "\"$c\"" : "$table_alias.\"$c\""

    conditions = String[]
    raw_filters = get(filter_config, "filters", [])

    # remove_empty conditions
    for col in get(filter_config, "remove_empty", [])
        col_str = _xlate(string(col))
        col_str in col_set || continue
        r = _ref(col_str)
        push!(conditions,
            """($r IS NOT NULL AND TRIM(CAST($r AS VARCHAR)) != '' AND LOWER(TRIM(CAST($r AS VARCHAR))) != 'blank')""")
    end

    for item in raw_filters
        item isa Dict || continue

        # Pattern-based rules
        if haskey(item, "pattern")
            col = _xlate(string(get(item, "column", "")))
            isempty(col) && continue
            col in col_set || continue
            pat = string(get(item, "pattern", ""))
            isempty(pat) && continue
            action = lowercase(string(get(item, "action", "exclude")))
            use_regex = get(item, "regex", false) == true
            r = _ref(col)

            if use_regex
                match_expr = "REGEXP_MATCHES(CAST($r AS VARCHAR), $(_sql_str(pat)))"
            else
                # Escape the pattern before embedding in a LIKE literal so a
                # single-quote in `pat` cannot break out of the string.
                pat_esc = replace(pat, "'" => "''")
                match_expr = "CAST($r AS VARCHAR) LIKE '%$(pat_esc)%'"
            end

            # A NULL value counts as no match.
            if action == "keep"
                push!(conditions, "COALESCE($match_expr, FALSE)")
            else  # exclude
                push!(conditions, "NOT COALESCE($match_expr, FALSE)")
            end

        # Type-based rules (min/max/include)
        elseif haskey(item, "type")
            col = _xlate(string(get(item, "column", "")))
            isempty(col) && continue
            col in col_set || continue
            typ = string(get(item, "type", ""))
            r = _ref(col)
            if typ == "min"
                val = get(item, "value", nothing)
                isnothing(val) && continue
                # Validate numerically: an unvalidated string could inject SQL.
                parsed = tryparse(Float64, string(val))
                isnothing(parsed) && continue
                push!(conditions, "TRY_CAST($r AS DOUBLE) >= $parsed")
            elseif typ == "max"
                val = get(item, "value", nothing)
                isnothing(val) && continue
                parsed = tryparse(Float64, string(val))
                isnothing(parsed) && continue
                push!(conditions, "TRY_CAST($r AS DOUBLE) <= $parsed")
            elseif typ == "include"
                vals = get(item, "values", [])
                isempty(vals) && continue
                # Escape each value so single-quotes cannot break the IN literal.
                val_list = join([_sql_str(v) for v in vals], ", ")
                push!(conditions, "CAST($r AS VARCHAR) IN ($val_list)")
            end
        end
    end

    conditions
end

# Column names a filter config references, translated for the selected source
# via `col_map`. Covers `remove_empty` entries and any `filters` item carrying a
# `column`. Used to check a filter is fully realisable against a table's schema.
function filter_column_refs(filter_config::Dict;
                            col_map::Dict{String,String}=Dict{String,String}())
    _xlate(c) = get(col_map, c, c)
    refs = String[]
    for col in get(filter_config, "remove_empty", [])
        push!(refs, _xlate(string(col)))
    end
    for item in get(filter_config, "filters", [])
        item isa Dict || continue
        (haskey(item, "pattern") || haskey(item, "type")) || continue
        col = string(get(item, "column", ""))
        isempty(col) && continue
        push!(refs, _xlate(col))
    end
    refs
end

"""
    category_case_when(categories, merged_col_set, source; filters, table_alias, strict) -> Union{String,Nothing}

Build a SQL CASE WHEN expression that classifies each row into a category.
Filter column names are translated to the correct source columns (VSEARCH or
DADA2). The `funcdb_require` key on a category entry is silently ignored
(FuncDB is quarantined).
`table_alias` defaults to "m" (SELECT context with FROM ... m); pass "" for
UPDATE SET contexts where bare column references are required.

A category carrying no `filter` names the catch-all bucket: its name becomes the
CASE `ELSE` label (the last such category wins), so every row unmatched by the
filtered categories is labelled with it rather than the default `Unassigned`.
This lets a set model "everything else", contaminants included, as an explicit
named, colourable category.

By default (`strict=false`, for display) a filter that references an absent
column has that condition dropped and still contributes a branch, and the
result is always a string. With `strict=true` (for row-dropping figure
exclusion) the classification must be exactly realisable: if any category
names a filter missing from `filters`, or that filter references a column
absent from `merged_col_set`, or no category yields a branch, `nothing` is
returned so the caller drops nothing rather than the wrong rows.
"""
function category_case_when(categories::Vector, merged_col_set::Set{String},
                            source::String; filters::Dict,
                            table_alias::String="m", strict::Bool=false)
    col_map = col_translate_map(source, merged_col_set)
    branches = String[]
    catchall = "Unassigned"

    # A source the table does not carry must not silently borrow the other
    # source's columns: on a vsearch-free table the plain rank names hold DADA2
    # taxonomy, and a VSEARCH-sourced filter naming "Genus" would match against
    # it as though vsearch had run. Refuse instead, loudly.
    avail = available_sources(merged_col_set)
    if !(source in avail)
        @warn "category_case_when: the table carries no $source taxonomy, so nothing " *
              "can be categorised from it. Set the source to one this run produced." source available=avail
        return strict ? nothing : _sql_str(catchall)
    end

    # Resolve filter names against stringified keys. A caller that builds this
    # dict itself, rather than taking it from CompositionLibrary.load, may hold a
    # non-String key (YAML and JSON both parse a numeric name as an integer), and
    # a raw lookup would then dangle a filter that is in fact present. A dict
    # already keyed by String, which is what the library hands us, needs no
    # rebuild; this runs once per set, so the copy is not free.
    by_name = keytype(filters) === String ? filters :
              Dict{String,Any}(string(k) => v for (k, v) in filters)

    for cat in categories
        cat_name = get(cat, "name", "")
        isempty(cat_name) && continue

        filter_name = get(cat, "filter", nothing)
        # A filterless category names the catch-all bucket (the ELSE label).
        isnothing(filter_name) && (catchall = string(cat_name); continue)

        # Resolve the filter by name from the library. A dangling name is a
        # config error: dropped for display, fatal under strict. Warn in both
        # modes so the error is never silent (a bad migration that dangles
        # every category must not render as a plausible all-Unassigned figure).
        filter_config = get(by_name, string(filter_name), nothing)
        if isnothing(filter_config)
            @warn "category_case_when: category names a filter absent from the filter library" category=cat_name filter=filter_name
            strict && return nothing
            continue
        end

        # A strict caller drops rows on this classification, so a partially
        # degraded filter (any referenced column absent) is untrustworthy: bail
        # rather than silently over- or under-match. Display callers tolerate it.
        if strict && any(ref -> !(ref in merged_col_set),
                         filter_column_refs(filter_config; col_map))
            return nothing
        end

        conditions = filter_to_sql_conditions(filter_config, merged_col_set,
                                              table_alias; col_map)

        # funcdb_require is ignored: FuncDB is quarantined and its annotation
        # columns are not available in the pipeline-side composition path.
        isempty(conditions) && continue
        # Escape single-quotes in the category name so a crafted name cannot
        # break out of the THEN string literal.
        push!(branches, "WHEN (" * join(conditions, " AND ") * ") THEN $(_sql_str(cat_name))")
    end

    isempty(branches) && return strict ? nothing : _sql_str(catchall)
    "CASE\n    " * join(branches, "\n    ") * "\n    ELSE $(_sql_str(catchall))\nEND"
end

## DuckDB category-column writer and lazy backfill

"""
    write_category_columns!(con, table, source, set_names; library)

For each category set in `set_names`, add a `Category__<set>` VARCHAR column to
`table` (idempotent via ADD COLUMN IF NOT EXISTS) and UPDATE it with the
CASE WHEN classification expression for `source`. Columns already populated are
still refreshed by the UPDATE; use `ensure_columns!` for a pure backfill.
`library` is the whole composition library Dict (as returned by
`CompositionLibrary.load`): a set's categories are read from
`library["sets"][set_name]["categories"]`, and `library["filters"]` is passed
through as the `filters` a category's classification resolves against.
"""
function write_category_columns!(con, table::String, source::String,
                                 set_names::Vector{String};
                                 library::Dict, suffixed::Bool=false)
    cols = Set(string.(DataFrame(DBInterface.execute(con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = ?",
        [table])).column_name))
    sets = get(library, "sets", Dict())
    filters = get(library, "filters", Dict())
    for set_name in set_names
        # Guard: the set name is used as a quoted SQL identifier. Reject any
        # name that falls outside safe identifier characters so a double-quote
        # cannot escape the quoting and inject SQL.
        if !Validation.is_safe_name(set_name)
            @warn "write_category_columns!: skipping set with unsafe name" set_name
            continue
        end
        cfg = get(sets, set_name, nothing)
        if isnothing(cfg)
            @warn "write_category_columns!: skipping set absent from the library" set_name
            continue
        end
        cats = get(cfg, "categories", [])
        # Use empty alias: UPDATE runs against the bare table, no FROM alias.
        case = category_case_when(cats, cols, source;
                                  filters, table_alias="")
        colname = column_name(set_name, suffixed ? source : nothing)
        DBInterface.execute(con,
            "ALTER TABLE \"$table\" ADD COLUMN IF NOT EXISTS \"$colname\" VARCHAR")
        DBInterface.execute(con,
            "UPDATE \"$table\" SET \"$colname\" = $case")
    end
    nothing
end

"""
    apply_max_x!(con, table, rank_cols, max_x)

Delete rows from `table` whose count of unresolved `_X` placeholders across
`rank_cols` exceeds `max_x`. Does nothing when `max_x < 0`.

The count expression mirrors `_build_quality_where` in
`src/server/routes/composition.jl`: NULL/empty/NA columns each contribute 1;
non-empty columns contribute the number of `_X` substrings measured by the
LENGTH/REPLACE trick `(LENGTH(col) - LENGTH(REPLACE(col, '_X', ''))) / 2`.
"""
function apply_max_x!(con, table::String, rank_cols::Vector{String}, max_x::Int)
    max_x < 0 && return nothing
    isempty(rank_cols) && return nothing
    x_checks = String[]
    for c in rank_cols
        push!(x_checks, """CASE
            WHEN \"$c\" IS NULL OR TRIM(CAST(\"$c\" AS VARCHAR)) = '' OR LOWER(TRIM(CAST(\"$c\" AS VARCHAR))) = 'na'
                THEN 1
            ELSE (LENGTH(CAST(\"$c\" AS VARCHAR)) - LENGTH(REPLACE(CAST(\"$c\" AS VARCHAR), '_X', ''))) / 2
        END""")
    end
    count_expr = join(x_checks, " + ")
    DBInterface.execute(con,
        "DELETE FROM \"$table\" WHERE ($count_expr) > $max_x")
    nothing
end

"""
    ensure_columns!(con, table, source, set_names; library)

Lazy backfill: write only the `Category__<set>` columns that are absent from
`table`, leaving any already-present ones untouched. Idempotent. `library` is
the whole composition library Dict; see `write_category_columns!`.
"""
function ensure_columns!(con, table::String, source::String,
                         set_names::Vector{String};
                         library::Dict, suffixed::Bool=false)
    present = Set(string.(DataFrame(DBInterface.execute(con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = ?",
        [table])).column_name))
    missing_sets = filter(s -> !(column_name(s, suffixed ? source : nothing) in present), set_names)
    isempty(missing_sets) && return nothing
    write_category_columns!(con, table, source, missing_sets; library, suffixed)
    nothing
end

end # module Categories
