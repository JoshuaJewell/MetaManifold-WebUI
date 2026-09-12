# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: organism composition - category-based classification of ASVs/OTUs
# for relative abundance analysis across broad organism groups.
using JSON3, CSV, DataFrames, OrderedCollections, DuckDB, DBInterface, YAML
using ..Categories, ..CompositionLibrary
using ..Analysis: pinned_segment_order

# Catch-all bucket for taxa matching no category, and its fixed legend colour.
const _UNASSIGNED_CATEGORY = "Unassigned"
const _UNASSIGNED_COLOUR = "#95a5a6"

## Composition library
function _library_path()
    joinpath(dirname(ServerState.projects_dir()), "config", "composition.yml")
end

_library() = CompositionLibrary.load(_library_path())

# Write the whole library back to config/composition.yml.
_write_library(lib::Dict) = _atomic_write_yaml(_library_path(), lib)

# List the category sets defined in the composition library.
function _list_category_sets()
    sets = Dict{String,Any}[]
    for (name, config) in sort(collect(_library()["sets"]); by=first)
        push!(sets, _category_set_summary(name, config))
    end
    sets
end

# Persist the whole library, validating the document before it lands on disk.
# This is the single gate through which every filter and set write passes, so
# a dangling filter reference is caught here rather than reaching the SQL
# generator, where it would silently classify every row as 'Unassigned'.
# Returns the library, or an HTTP.Response error.
function _save_library(lib::Dict)
    errors = CompositionLibrary.validate(lib)
    isempty(errors) || return json_error(400, "invalid_library", join(errors, "; "))
    _write_library(lib)
    lib
end

# Save a filter config verbatim into the library. Returns the saved filter
# config, or an HTTP.Response error.
function _save_filter_unlocked(name::String, config::AbstractDict)
    lib = _library()
    lib["filters"][name] = config
    result = _save_library(lib)
    result isa HTTP.Response && return result
    @info "Saved filter: $name (library $(_library_path()))"
    result["filters"][name]
end

# Remove a filter from the library; refuses when a set still references it,
# naming the referencing sets. Returns the deleted name, or an HTTP.Response.
function _delete_filter_unlocked(name::String)
    lib = _library()
    users = CompositionLibrary.filter_in_use(lib, name)
    isempty(users) || return json_error(409, "filter_in_use",
        "Filter '$name' is used by: " * join(users, ", "))
    haskey(lib["filters"], name) || return json_error(404, "filter_not_found",
        "Filter '$name' not found")
    delete!(lib["filters"], name)
    _write_library(lib)
    @info "Deleted filter: $name (library $(_library_path()))"
    name
end

# Save a whole set config verbatim into the library. This is the composition
# builder's direct-save path, distinct from `_save_category_set`'s
# save-as-a-recolour path used by the existing category-sets routes. Returns
# the saved set config, or an HTTP.Response error.
function _save_composition_set_unlocked(name::String, config::AbstractDict)
    lib = _library()
    lib["sets"][name] = config
    result = _save_library(lib)
    result isa HTTP.Response && return result
    @info "Saved category set: $name (library $(_library_path()))"
    result["sets"][name]
end

## Filter to SQL translation
# These helpers now delegate to the Categories module, which is self-contained
# and does not depend on FuncDBAnnotation or server-state globals.
_col_translate_map(source::String, present=nothing) =
    Categories.col_translate_map(source, present)

_filter_to_sql_conditions(filter_config::Dict, col_set::Set{String},
                          table_alias::String;
                          col_map::Dict{String,String}=Dict{String,String}()) =
    Categories.filter_to_sql_conditions(filter_config, col_set, table_alias; col_map)

## Live composition summary

"""
    _composition_summary(study, run, category_set_name, subgroup;
                         group, params, table="merged", tag="category", value) -> HTTP.Response

Per-label row and read counts for one results table. `tag` follows the chart
routes: `"category"` groups by the `Category__<value>` column (backfilled when
absent), `"rank"` groups by the taxonomy rank `value`, blank ranks counting as
Unclassified as in `aggregate_by_taxon`.

`subgroup` is either `nothing` (all sample columns) or a prefix string
(columns matching `"<prefix>_"`). Returns a 400 when the subgroup prefix
matches no sample columns.

The returned JSON matches the `CompositionBuildResult` frontend type:
`{ table, source, category_set, tag, value, total_rows, total_reads, categories }`.
"""
function _composition_summary(study::String, run::String, category_set_name::String,
                               subgroup::Union{String,Nothing};
                               group::Union{String,Nothing}=nothing,
                               params::Dict{String,String}=Dict{String,String}(),
                               table::String="merged",
                               tag::String="category",
                               value::String=category_set_name)
    tag in ("category", "rank") || return json_error(400, "bad_tag",
        "tag must be 'category' or 'rank'")
    Validation.is_safe_name(table) || return json_error(400, "invalid_table",
        "Table name must contain only letters, numbers, dots, hyphens, and underscores")

    lib = _library()
    if tag == "category"
        # The name is interpolated into a quoted SQL identifier below, via
        # `Categories.column_name`. Library membership holds only because every
        # API write charset-checks the key, so the guard is repeated here.
        Validation.is_safe_name(value) ||
            return json_error(400, "invalid_name",
                "Category set name must contain only letters, numbers, dots, hyphens, and underscores")
        haskey(lib["sets"], value) || return json_error(404, "category_set_not_found",
            "Category set '$value' not found")
    end

    source = _tagging_source(study, run; group)

    _with_analysis_results_table(study, run, table; group,
                                 readonly = tag != "category") do con, columns
        label_expr = if tag == "category"
            Categories.ensure_columns!(con, table, source, [value]; library=lib, suffixed=_suffixed())
            "\"$(_category_column(value))\""
        else
            rank_col = _rank_column(columns, value)
            isnothing(rank_col) && return json_error(400, "bad_rank",
                "Unknown rank '$value' in table '$table'")
            "COALESCE(NULLIF(TRIM(\"$rank_col\"), ''), 'Unclassified')"
        end

        all_sample_cols = Analysis.sample_columns(con, table)

        # Resolve the retained sample columns for the requested subgroup.
        retained = _filter_by_prefix(all_sample_cols, subgroup)
        isempty(retained) && return json_error(400, "no_subgroup_samples",
            "Sub-group selection '$(something(subgroup, ""))' matches no sample columns in '$table'")

        # The summary applies no row filters, so both read-count bases measure
        # the same thing here: each sample's whole library.
        retained = _retain_sample_columns(con, table, retained, params, "", Any[])
        isempty(retained) && return json_error(400, "no_samples",
            "No samples pass the sample read-count filter")

        # Per-row read sum (for the WHERE clause) and per-group read sum (for SELECT).
        row_sum  = join(["COALESCE(\"$c\", 0)" for c in retained], " + ")
        grp_sum  = join(["COALESCE(SUM(\"$c\"), 0)" for c in retained], " + ")

        # Rows with no reads in the retained columns are excluded; HAVING drops
        # any label whose aggregate is zero.
        sql = """
            SELECT $label_expr AS cat, COUNT(*) AS rows, ($grp_sum) AS reads
            FROM \"$table\"
            WHERE ($row_sum) > 0
            GROUP BY cat
            HAVING ($grp_sum) > 0
            ORDER BY reads DESC, cat ASC
        """
        df = DataFrame(DBInterface.execute(con, sql))

        total_rows  = sum(df.rows;  init=0)
        total_reads = sum(df.reads; init=0)

        # Reads descending, but pinned as the figure legends pin them, so a
        # category holds one place across the whole view.
        cat_stats = OrderedDict{String,Any}()
        for row in eachrow(df[pinned_segment_order(string.(df.cat)), :])
            cat_stats[string(row.cat)] = Dict(
                "rows"          => Int(row.rows),
                "reads"         => Int(row.reads),
                "reads_percent" => total_reads > 0 ?
                    round(100.0 * row.reads / total_reads; digits=2) : 0.0,
            )
        end

        json(Dict(
            "table"        => table,
            "source"       => source,
            "category_set" => category_set_name,
            "tag"          => tag,
            "value"        => value,
            "total_rows"   => total_rows,
            "total_reads"  => total_reads,
            "samples"      => retained,
            "categories"   => cat_stats,
        ))
    end
end

## Routes

# The whole composition library: filters plus sets, for the frontend's
# composition builder (Tasks 7-8).
@get "/api/v1/composition" function(req)
    json(_library())
end

# Create or update a named filter. The whole library is validated before the
# write lands on disk, so a name violating the charset or a value that is not
# a Dict is rejected with 400 "invalid_library".
@post "/api/v1/composition/filters/{name}" function(req, name::String)
    body = _to_plain(JSON3.read(String(req.body)))
    result = _save_filter(name, body)
    result isa HTTP.Response && return result
    json(result)
end

# Delete a filter. Refused with 409 "filter_in_use" (naming the referencing
# sets) while any set still references it; 404 "filter_not_found" when absent.
@delete "/api/v1/composition/filters/{name}" function(req, name::String)
    result = _delete_filter(name)
    result isa HTTP.Response && return result
    json((; deleted=result))
end

# Create or update a set from its whole config (label, description,
# unassigned_colour, categories). A category naming a filter absent from the
# library fails validation and is rejected with 400 "invalid_library".
@post "/api/v1/composition/sets/{name}" function(req, name::String)
    body = _to_plain(JSON3.read(String(req.body)))
    result = _save_composition_set(name, body)
    result isa HTTP.Response && return result
    json(result)
end

# Delete a set; `default` is protected. Shares `_delete_category_set` since
# the library is the single store regardless of which route deletes a set.
@delete "/api/v1/composition/sets/{name}" function(req, name::String)
    result = _delete_category_set(name)
    result isa HTTP.Response && return result
    json((; deleted=result))
end

# List category sets
@get "/api/v1/category-sets" function(req)
    json(_list_category_sets())
end

# Emit one saved category set in the shape `_list_category_sets` produces.
function _category_set_summary(name::String, config)
    cs = get(config, "categories", [])
    cats = [Dict("name" => get(c, "name", ""), "colour" => get(c, "colour", nothing))
            for c in cs]
    Dict(
        "name"             => name,
        "label"            => get(config, "label", get(config, "name", name)),
        "description"      => get(config, "description", ""),
        "categories"       => cats,
        "unassigned_colour" => string(get(config, "unassigned_colour", _UNASSIGNED_COLOUR)),
    )
end

# Clone the base set's structure verbatim, overriding only the colours by
# category name, and write the updated set into the library at
# `config/composition.yml`. Returns the saved set summary, or an
# `HTTP.Response` error mirroring the route's failure modes.
function _save_category_set_unlocked(name::String, base::String,
                            colours::Dict{String,String};
                            label::Union{String,Nothing}=nothing,
                            description::Union{String,Nothing}=nothing)
    Validation.is_safe_name(name) || return json_error(400, "invalid_name",
        "Name must contain only letters, numbers, dots, hyphens, and underscores")
    isempty(base) && return json_error(400, "missing_base", "Body must include 'base'")

    lib = _library()
    base_config = get(lib["sets"], base, nothing)
    isnothing(base_config) && return json_error(404, "base_not_found",
        "Base category set '$base' not found")

    for clr in values(colours)
        occursin(CompositionLibrary.COLOUR_RE, clr) || return json_error(400, "invalid_colour",
            "Colour '$clr' must be a hex triplet like '#3498db'")
    end

    # Rebuild the category list preserving order, `filter`, and `funcdb_require`
    # byte-for-byte; only `colour` is replaced where an override is supplied.
    base_cats = get(base_config, "categories", [])
    out_cats = OrderedDict{String,Any}[]
    for cat in base_cats
        cat_name = string(get(cat, "name", ""))
        entry = OrderedDict{String,Any}("name" => cat_name)
        new_colour = get(colours, cat_name, get(cat, "colour", nothing))
        isnothing(new_colour) || (entry["colour"] = string(new_colour))
        haskey(cat, "filter") && (entry["filter"] = cat["filter"])
        haskey(cat, "funcdb_require") && (entry["funcdb_require"] = cat["funcdb_require"])
        push!(out_cats, entry)
    end

    out = OrderedDict{String,Any}(
        "label"       => isnothing(label) ? get(base_config, "label", name) : label,
        "description" => isnothing(description) ? get(base_config, "description", "") : description,
        "categories"  => out_cats,
    )

    # Route the "Unassigned" colour override to the top-level field rather than a category entry.
    unassigned = get(colours, _UNASSIGNED_CATEGORY,
                     get(base_config, "unassigned_colour", nothing))
    isnothing(unassigned) || (out["unassigned_colour"] = string(unassigned))

    lib["sets"][name] = out
    result = _save_library(lib)
    result isa HTTP.Response && return result
    @info "Saved category set: $name (library $(_library_path()))"
    _category_set_summary(name, out)
end

# Remove a set from the library; `default` is protected. Returns the deleted
# name, or an `HTTP.Response` error.
function _delete_category_set_unlocked(name::String)
    Validation.is_safe_name(name) || return json_error(400, "invalid_name",
        "Name must contain only letters, numbers, dots, hyphens, and underscores")
    name == "default" && return json_error(400, "protected_set",
        "The 'default' category set cannot be deleted")

    lib = _library()
    haskey(lib["sets"], name) || return json_error(404, "set_not_found",
        "Category set '$name' not found")
    delete!(lib["sets"], name)
    _write_library(lib)
    @info "Deleted category set: $name (library $(_library_path()))"
    name
end

# Save a category set, mirroring the filter-preset save apparatus.
@post "/api/v1/category-sets/{name}" function(req, name::String)
    body = JSON3.read(String(req.body))
    base = string(get(body, :base, ""))
    colours = Dict{String,String}()
    let c = get(body, :colours, nothing)
        isnothing(c) || for (cat, clr) in pairs(c)
            colours[string(cat)] = string(clr)
        end
    end
    label = let l = get(body, :label, nothing); isnothing(l) ? nothing : string(l) end
    description = let d = get(body, :description, nothing); isnothing(d) ? nothing : string(d) end

    result = _save_category_set(name, base, colours; label, description)
    result isa HTTP.Response && return result
    json(result)
end

# Delete a category set; `default` is protected and cannot be removed.
@delete "/api/v1/category-sets/{name}" function(req, name::String)
    result = _delete_category_set(name)
    result isa HTTP.Response && return result
    json((; deleted=result))
end

# Live composition summary: per-label row and read counts for a results table
# (default merged), by category set or by rank.
@post "/api/v1/studies/{study}/runs/{run}/composition/summary" function(req,
                                                                         study::String,
                                                                         run::String)
    err = _validate_run_request(study, run)
    isnothing(err) || return err

    body = JSON3.read(String(req.body))
    category_set = string(get(body, :category_set, "default"))
    subgroup = let s = get(body, :subgroup, nothing)
        isnothing(s) ? nothing : string(s)
    end
    table = string(get(body, :table, "merged"))
    tag   = string(get(body, :tag, "category"))
    value = string(get(body, :value, category_set))

    try
        _composition_summary(study, run, category_set, subgroup;
                             group=_req_group(req), params=_body_filter_params(body),
                             table, tag, value)
    catch e
        @error "Composition summary failed" study run category_set subgroup exception=(e, catch_backtrace())
        json_error(500, "composition_summary_failed",
            "Failed to compute composition summary: $(sprint(showerror, e))")
    end
end

# Paginated query of a results table (default merged), ensuring the requested
# category set column exists.
@post "/api/v1/studies/{study}/runs/{run}/composition/{source}/query" function(req,
                                                                                study::String,
                                                                                run::String,
                                                                                source::String)
    err = _require_study_run_source(study, run, source)
    !isnothing(err) && return err
    group = _req_group(req)

    body = JSON3.read(String(req.body))
    category_set = string(get(body, :category_set, "default"))

    table = string(get(body, :table, "merged"))
    Validation.is_safe_name(table) || return json_error(400, "invalid_table",
        "Table name must contain only letters, numbers, dots, hyphens, and underscores")

    tag_src = _tagging_source(study, run; group)
    _with_analysis_results_table(study, run, table; group, readonly=false) do con, _
        Categories.ensure_columns!(con, table, tag_src, [category_set];
                                   library=_library(), suffixed=_suffixed())
        _duckdb_paginated_query(con, table, body)
    end
end

# Distinct values for a results-table column (default merged); ensures the
# category column when requested.
@post "/api/v1/studies/{study}/runs/{run}/composition/{source}/distinct/{column}" function(req,
                                                                                            study::String,
                                                                                            run::String,
                                                                                            source::String,
                                                                                            column::String)
    err = _require_study_run_source(study, run, source)
    !isnothing(err) && return err
    group = _req_group(req)

    body = JSON3.read(String(req.body))
    category_set = string(get(body, :category_set, "default"))

    table = string(get(body, :table, "merged"))
    Validation.is_safe_name(table) || return json_error(400, "invalid_table",
        "Table name must contain only letters, numbers, dots, hyphens, and underscores")

    tag_src = _tagging_source(study, run; group)
    _with_analysis_results_table(study, run, table; group, readonly=false) do con, _
        Categories.ensure_columns!(con, table, tag_src, [category_set];
                                   library=_library(), suffixed=_suffixed())
        _duckdb_distinct(con, table, column, body)
    end
end

## Library edits hold the config lock so a load-modify-write cannot interleave with another save.
_save_filter(args...; kwargs...) = lock(() -> _save_filter_unlocked(args...; kwargs...), _config_file_lock)
_delete_filter(args...; kwargs...) = lock(() -> _delete_filter_unlocked(args...; kwargs...), _config_file_lock)
_save_composition_set(args...; kwargs...) = lock(() -> _save_composition_set_unlocked(args...; kwargs...), _config_file_lock)
_save_category_set(args...; kwargs...) = lock(() -> _save_category_set_unlocked(args...; kwargs...), _config_file_lock)
_delete_category_set(args...; kwargs...) = lock(() -> _delete_category_set_unlocked(args...; kwargs...), _config_file_lock)
