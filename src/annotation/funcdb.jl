module FuncDBAnnotation

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export load_funcdb, annotate_table, append_funcdb_entry,
       apply_contamination_filter!,
       CONTAM_BLACKLIST_MAP, CONTAM_WHITELIST_MAP

using DataFrames, CSV, Logging, Dates

## Constants

# FuncDB value columns to output column names.
const FUNCDB_VALUE_COLS = (
    :Function,
    :Detailed_function,
    :Associated_organism,
    :Associated_material,
    :Environment,
    :Potential_human_pathogen,
    :Comment,
    :Reference,
)

const FUNCDB_OUTPUT_COLS = [
    "funcdb_function",
    "funcdb_detailed_function",
    "funcdb_associated_organism",
    "funcdb_associated_material",
    "funcdb_environment",
    "funcdb_potential_human_pathogen",
    "funcdb_comment",
    "funcdb_reference",
]

# Rank hierarchy from finest to coarsest.
# Each entry: (rank, FuncDB column, VSEARCH merged column, DADA2 merged column).
# Subdivision is in merged tables but absent from FuncDB and is skipped.
const RANK_HIERARCHY = [
    (rank="species",    funcdb=:Species,    vsearch="Species",    dada2="Species_dada2"),
    (rank="genus",      funcdb=:Genus,      vsearch="Genus",      dada2="Genus_dada2"),
    (rank="family",     funcdb=:family,     vsearch="Family",     dada2="Family_dada2"),
    (rank="order",      funcdb=:order,      vsearch="Order",      dada2="Order_dada2"),
    (rank="class",      funcdb=:class,      vsearch="Class",      dada2="Class_dada2"),
    (rank="division",   funcdb=:division,   vsearch="Division",   dada2="Division_dada2"),
    (rank="supergroup", funcdb=:supergroup, vsearch="Supergroup", dada2="Supergroup_dada2"),
]

# VSEARCH-only columns dropped from DADA2-annotated output.
# Maps rank name to its index in RANK_HIERARCHY (1 = finest = species).
# Used in load_funcdb to enforce Assignment_level as a ceiling when building maps.
const _RANK_INDEX = Dict(r.rank => i for (i, r) in enumerate(RANK_HIERARCHY))

const VSEARCH_ONLY_COLS = Set([
    "Pident", "Accession", "rRNA", "Organellum", "specimen", "sequence",
    "Domain", "Supergroup", "Division", "Subdivision", "Class", "Order",
    "Family", "Genus", "Species",
])

const CONTAM_BLACKLIST_MAP = [
    ("funcdb_function",            "annotation.contamination.blacklist.function"),
    ("funcdb_detailed_function",   "annotation.contamination.blacklist.detailed_function"),
    ("funcdb_associated_organism", "annotation.contamination.blacklist.associated_organism"),
    ("funcdb_associated_material", "annotation.contamination.blacklist.associated_material"),
    ("funcdb_environment",         "annotation.contamination.blacklist.environment"),
]

const CONTAM_WHITELIST_MAP = [
    ("funcdb_function",            "annotation.contamination.whitelist.function"),
    ("funcdb_detailed_function",   "annotation.contamination.whitelist.detailed_function"),
    ("funcdb_associated_organism", "annotation.contamination.whitelist.associated_organism"),
    ("funcdb_associated_material", "annotation.contamination.whitelist.associated_material"),
    ("funcdb_environment",         "annotation.contamination.whitelist.environment"),
]

## Internal helpers

"""
    _normalize_key(val) -> String

Lowercase, trimmed string. Returns empty string for missing, blank, or "blank".
"""
function _normalize_key(val)::String
    (ismissing(val) || isnothing(val)) && return ""
    s = strip(lowercase(string(val)))
    s == "blank" && return ""
    return s
end

"""
    _extract_values(row) -> NamedTuple

Pull the FuncDB value columns from a DataFrame row into a NamedTuple of Strings.
"""
function _extract_values(row)
    vals = map(FUNCDB_VALUE_COLS) do col
        v = hasproperty(row, col) ? getproperty(row, col) : missing
        (ismissing(v) || isnothing(v)) ? "" : strip(string(v))
    end
    return NamedTuple{FUNCDB_VALUE_COLS}(vals)
end

## load_funcdb

"""
    load_funcdb(path) -> Dict{String, Dict{String, NamedTuple}}

Load `FuncDB_species.csv` and build a lookup map per taxonomic rank.

Returns a Dict mapping rank name (e.g. "species", "genus", "family") to a
Dict of normalized taxon key to NamedTuple of FuncDB values.

Duplicate species keys log a warning; duplicate keys at higher ranks are
expected (many species per genus, etc.) and silently keep the first entry.
"""
function load_funcdb(path::String)
    isfile(path) || error("FuncDB file not found: $path")

    df = CSV.read(path, DataFrame; stringtype=String)

    maps = Dict{String, Dict{String, NamedTuple}}(
        r.rank => Dict{String, NamedTuple}() for r in RANK_HIERARCHY
    )

    for row in eachrow(df)
        vals = _extract_values(row)

        # Assignment_level is the finest rank at which the annotation was made.
        # Only populate maps for ranks at or finer than this level; indexing into
        # coarser maps would let a genus-level entry match via a shared family name.
        assignment_level = _normalize_key(
            hasproperty(row, :Assignment_level) ? getproperty(row, :Assignment_level) : missing
        )
        max_idx = get(_RANK_INDEX, assignment_level, length(RANK_HIERARCHY))

        for (i, r) in enumerate(RANK_HIERARCHY)
            i > max_idx && break
            key = _normalize_key(hasproperty(row, r.funcdb) ? getproperty(row, r.funcdb) : missing)
            isempty(key) && continue
            if haskey(maps[r.rank], key)
                r.rank == "species" && @warn "Duplicate species key in FuncDB, keeping first" key=key
            else
                maps[r.rank][key] = vals
            end
        end
    end

    @info "FuncDB loaded" species=length(maps["species"]) genera=length(maps["genus"])
    return maps
end

## annotate_table

"""
    annotate_table(source_df, taxonomy_source, funcdb_path) -> DataFrame

Annotate a taxonomy DataFrame with FuncDB functional annotations.

- `taxonomy_source`: "VSEARCH" or "DADA2"
- Matching walks the rank hierarchy from species down to supergroup, using
  the first match found.
- Columns from the other taxonomy source are dropped from the output.
- Appends funcdb_* value columns plus `funcdb_match_rank`.
"""
function annotate_table(source_df::DataFrame, taxonomy_source::String, funcdb_path::String)::DataFrame
    taxonomy_source in ("VSEARCH", "DADA2") ||
        error("Invalid taxonomy_source: '$taxonomy_source'. Must be VSEARCH or DADA2.")

    merged_col_field = taxonomy_source == "VSEARCH" ? :vsearch : :dada2

    maps = load_funcdb(funcdb_path)

    nrows = nrow(source_df)

    out_vals    = [Vector{String}(undef, nrows) for _ in 1:length(FUNCDB_OUTPUT_COLS)]
    match_rank  = Vector{String}(undef, nrows)

    for i in 1:nrows
        row = source_df[i, :]

        matched_vals = nothing
        matched_rank = "unmatched"

        for r in RANK_HIERARCHY
            col = getfield(r, merged_col_field)
            # Skip if merged table doesn't have this rank column
            hasproperty(source_df, Symbol(col)) || continue

            key = _normalize_key(row[Symbol(col)])
            isempty(key) && continue

            if haskey(maps[r.rank], key)
                matched_vals = maps[r.rank][key]
                matched_rank = r.rank
                break
            end
        end

        for (j, col_sym) in enumerate(FUNCDB_VALUE_COLS)
            out_vals[j][i] = isnothing(matched_vals) ? "" : string(matched_vals[col_sym])
        end

        match_rank[i] = matched_rank
    end

    result = copy(source_df)

    # Drop columns from the other taxonomy source
    if taxonomy_source == "VSEARCH"
        drop = filter(c -> endswith(c, "_dada2") || endswith(c, "_boot"), names(result))
    else
        drop = filter(c -> c in VSEARCH_ONLY_COLS, names(result))
    end
    isempty(drop) || select!(result, Not(drop))

    # Append FuncDB value columns
    for (j, col_name) in enumerate(FUNCDB_OUTPUT_COLS)
        result[!, col_name] = out_vals[j]
    end

    result[!, "funcdb_match_rank"] = match_rank
    result[!, "Contamination"] = fill("unassigned", nrows)

    return result
end

"""
    apply_contamination_filter!(df, resolved_cfg)

Apply contamination and non-contamination auto-filters to `df` in-place.

`resolved_cfg` is the flat dotted-key dict returned by `_resolve_config`, where
each value is a NamedTuple with a `.value` field.
- Contamination filter: row's funcdb field value is in the set -> `"yes"`
- Non-contamination filter: row's funcdb field value is in the set -> `"no"`
- If both filters match the same row, leave it `"unassigned"`
Only rows currently `"unassigned"` are modified.
"""
function apply_contamination_filter!(df::DataFrame, resolved_cfg::Dict)
    _val(key, default) = begin
        entry = get(resolved_cfg, key, nothing)
        isnothing(entry) ? default : entry.value
    end

    function _build_set(field_map)
        Dict{String,Set{String}}(
            col => begin
                raw = _val(cfg_key, [])
                items = raw isa AbstractVector ? raw : []
                Set{String}(lowercase(strip(string(v))) for v in items
                            if !isempty(strip(string(v))))
            end
            for (col, cfg_key) in field_map
        )
    end

    blacklist_sets = _build_set(CONTAM_BLACKLIST_MAP)
    whitelist_sets = _build_set(CONTAM_WHITELIST_MAP)

    all(isempty(s) for s in values(blacklist_sets)) &&
    all(isempty(s) for s in values(whitelist_sets)) && return

    col_names = Set(names(df))
    "Contamination" in col_names || return

    for i in 1:nrow(df)
        df[i, "Contamination"] == "unassigned" || continue

        contamination_hit = false
        for (col, _) in CONTAM_BLACKLIST_MAP
            col in col_names || continue
            isempty(blacklist_sets[col]) && continue
            val = lowercase(strip(string(df[i, col])))
            isempty(val) && continue
            if val in blacklist_sets[col]
                contamination_hit = true
                break
            end
        end

        non_contamination_hit = false
        for (col, _) in CONTAM_WHITELIST_MAP
            col in col_names || continue
            isempty(whitelist_sets[col]) && continue
            val = lowercase(strip(string(df[i, col])))
            isempty(val) && continue
            if val in whitelist_sets[col]
                non_contamination_hit = true
                break
            end
        end

        if contamination_hit && non_contamination_hit
            continue
        elseif contamination_hit
            df[i, "Contamination"] = "yes"
        elseif non_contamination_hit
            df[i, "Contamination"] = "no"
        end
    end
end

## append_funcdb_entry

"""
    append_funcdb_entry(path, entry::Dict; modified_by="") -> NamedTuple

Append a single entry to the FuncDB CSV file.

`entry` must contain at least one taxonomy key (Species, Genus, etc.) and at
least `Main_function`. Missing columns default to empty strings.
Returns a NamedTuple of the written row.
"""
function append_funcdb_entry(path::String, entry::Dict; modified_by::String="")
    isfile(path) || error("FuncDB file not found: $path")

    # Read the existing header to match column order exactly.
    existing_cols = Symbol.(names(CSV.read(path, DataFrame; limit=0)))

    # Build the row, filling missing keys with ""
    row = Dict{Symbol, String}()
    for col in existing_cols
        row[col] = get(entry, String(col), get(entry, col, ""))
    end
    row[:Modified_by]   = modified_by
    row[:Modified_date] = string(today())

    # Validate: at least one taxonomy key and a main function
    has_taxon = any(!isempty(get(row, Symbol(r.rank == "species" ? "Species" :
                    r.rank == "genus" ? "Genus" : r.rank), "")) for r in RANK_HIERARCHY)
    has_taxon || error("Entry must have at least one taxonomy key")
    isempty(get(row, :Function, "")) && error("Entry must have Function")

    df = DataFrame(; (col => [get(row, col, "")] for col in existing_cols)...)
    CSV.write(path, df; append=true)

    @info "FuncDB entry appended" species=get(row, :Species, "") genus=get(row, :Genus, "") modified_by
    return NamedTuple{Tuple(existing_cols)}(Tuple(get(row, c, "") for c in existing_cols))
end

end
