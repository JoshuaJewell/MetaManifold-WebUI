# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: on-demand analysis (alpha diversity, taxa bars, pipeline stats, NMDS, PERMANOVA)

using CSV, DataFrames

using ..Analysis: alpha_chart,
                             sample_columns, taxonomy_levels, taxon_column, filtered_counts,
                             aggregate_by_taxon, combined_counts_across_runs,
                             alpha_boxplot, nmds_chart, run_nmds, run_permanova, r_available,
                             sequence_column_name, normalise_counts, venn_taxa_present,
                             pool_columns, auto_min_depth, bar_chart, faceted_bar_chart,
                             natural_sort
using ..Categories, ..CompositionLibrary

function _validate_run_request(study::String, run::String)
    study in _study_names() || return json_error(404, "study_not_found", "Study '$study' not found")
    run in _all_run_names(study) || return json_error(404, "run_not_found", "Run '$run' not found")
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
    cfg_cache = Dict{Tuple{String,String,String}, Any}()
    norm = _alpha_normalisation(study; run, group, config_cache=cfg_cache)

    prefix = let p = get(body, :prefix, nothing); isnothing(p) ? nothing : string(p) end

    _with_analysis_results_table(study, run, table; group) do con, columns
        scols = _filter_by_prefix(sample_columns(con, table), prefix)
        isempty(scols) && return json_error(400, "no_samples", "No sample columns found in table '$table'")
        exclude_conditions = _exclusion_conditions(study, run, group, columns; surface="diversity", config_cache=cfg_cache)
        where_clause, where_params = _analysis_where_clause(params, columns; exclude_conditions)
        scols = _retain_sample_columns(con, table, scols, params, where_clause, where_params)
        isempty(scols) && return json_error(400, "no_samples",
                                            "No samples pass the sample read-count filter")
        mat = filtered_counts(con, table, scols, where_clause, where_params)
        alpha = alpha_diversity(mat; norm.method, norm.depth,
                                seed=_study_seed(study), norm.iterations)
        scols = scols[alpha.kept]
        isempty(scols) && return json_error(400, "no_samples",
                                            "All samples fall below the normalisation depth")

        fig = alpha_chart(scols, alpha.richness, alpha.shannon, alpha.simpson; short_labels=true)
        HTTP.Response(200, ["Content-Type" => "application/json"],
                      body=JSON3.write(fig))
    end
end

## Unified chart helpers

# Resolve the composition library without depending on composition.jl (loaded
# after this file).
#
# A single chart request reads the library from three points (row exclusions,
# the Category__ backfill, and the legend colours), and the cross-run chart does
# so once per run, so parsing the YAML afresh at each was 2N+1 parses of one
# small file per figure. The parse is memoised on the file's exact bytes: the
# read still happens on every call, but the read is not the cost, and the key
# cannot go stale. An (mtime, size) key would be cheaper still and is what one
# reaches for first, but it is unsound here: `mtime` is quantised coarsely
# enough that two writes in quick succession share a timestamp, and a same-size
# edit would then be served from the cache. A figure rendered against a
# superseded category set is wrong in a way nothing downstream would catch, so
# the key is the content itself.
#
# Readers must treat the returned Dict as immutable, since they now share one;
# every mutating path goes through composition.jl's `_library`, which parses
# its own copy.
const _CHART_LIB_LOCK  = ReentrantLock()
const _CHART_LIB_CACHE = Ref{Tuple{String,Vector{UInt8},Dict}}(("", UInt8[], Dict{String,Any}()))

function _chart_library()
    path  = joinpath(dirname(ServerState.projects_dir()), "config", "composition.yml")
    bytes = isfile(path) ? read(path) : UInt8[]
    lock(_CHART_LIB_LOCK) do
        cached_path, cached_bytes, cached_lib = _CHART_LIB_CACHE[]
        (cached_path == path && cached_bytes == bytes) && return cached_lib
        lib = CompositionLibrary.load(path)
        _CHART_LIB_CACHE[] = (path, bytes, lib)
        lib
    end
end

## Colour map for a category set: category name -> hex colour.
# Falls back to _UNASSIGNED_COLOUR for "Unassigned" and "#95a5a6" for anything else.
function _category_colour_for(set_name::String)
    cfg = get(_chart_library()["sets"], set_name, nothing)
    colour_map = Dict{String,String}()
    if !isnothing(cfg)
        for cat in get(cfg, "categories", [])
            name = get(cat, "name", nothing)
            col  = get(cat, "colour", nothing)
            (!isnothing(name) && !isnothing(col)) && (colour_map[string(name)] = string(col))
        end
        unassigned_col = get(cfg, "unassigned_colour", "#95a5a6")
        colour_map["Unassigned"] = string(unassigned_col)
    end
    lbl -> get(colour_map, lbl, "#95a5a6")
end

## Core chart-data helper (pure: no HTTP, testable directly).
#
# Resolves label column, applies scope/pooling, drops zero-read rows, and returns
# (segment_labels, sample_names, counts_matrix) or an error dict.
#
# Counts are raw read sums per (label, sample), and the chart routes plot them
# as they come. Depth normalisation (rarefy/SRS) is deliberately NOT applied
# here, unlike the alpha/NMDS/PERMANOVA surfaces: those estimate richness and
# distance, which are biased by uneven sampling effort, so equalising library
# size earns its cost there. A composition bar answers "how many reads, and in
# what proportions", and rarefying breaks both halves of that. It rescales every
# sample to the smallest positive library in the set, so the absolute bars stop
# being read counts and all stand at the same height; and it quantises the
# relative bars to multiples of 1/depth. A single shallow sample is enough to
# do it -- a 4-read library drags every other sample down to a 4-read total,
# which renders as bars that look like a count of four ASVs.
#
# Parameters:
#   con              - DuckDB connection (already open on the annotation table)
#   table            - table name
#   columns          - all column names present in the table
#   all_scols        - all sample columns discovered for this run/table
#   where_clause     - SQL WHERE fragment (may be "")
#   where_params     - positional parameters for where_clause
#   tag              - "rank" or "category"
#   value            - rank name or category-set name
#   subgroup         - nothing | a prefix string | "__pool__"
#   pool_groups      - sub-group labels to pool into when subgroup == "__pool__"
#   pool_fallback    - label used when pool_groups is empty (e.g. the run name)
#   top_n            - top-N limit; ignored for categories (uses typemax(Int))
#   source           - annotation source ("VSEARCH"|"DADA2"), for ensure_columns!
function _chart_data(con, table::String, columns::Vector{String},
                     all_scols::Vector{String},
                     where_clause::String, where_params::Vector,
                     tag::String, value::String,
                     subgroup::Union{String,Nothing},
                     pool_groups::Vector{String},
                     pool_fallback::String,
                     top_n::Int,
                     source::String)
    # Scope sample columns by subgroup.
    scols = if isnothing(subgroup) || subgroup == "__pool__"
        all_scols
    else
        _filter_by_prefix(all_scols, subgroup)
    end
    isempty(scols) && return json_error(400, "no_samples",
                                        "No sample columns for the requested scope")

    label_col = _label_column(con, table, columns, tag, value, source)
    label_col isa HTTP.Response && return label_col

    # For categories, never collapse into "Other".
    effective_top_n = (tag == "category") ? typemax(Int) : top_n

    agg = aggregate_by_taxon(con, table, scols, label_col, where_clause, where_params)
    nrow(agg) == 0 && return json_error(400, "no_data", "No data after filtering")

    segment_labels = String.(agg.taxon)
    counts = Matrix{Float64}(agg[:, scols])  # segments x samples

    # Zero-read guard: drop label rows whose sum across retained samples is zero.
    row_totals = vec(sum(counts, dims=2))
    keep_rows  = row_totals .> 0
    if !all(keep_rows)
        segment_labels = segment_labels[keep_rows]
        counts         = counts[keep_rows, :]
    end
    isempty(segment_labels) && return json_error(400, "no_data",
                                                  "All labels have zero reads")

    # Pool sample columns when requested.
    if subgroup == "__pool__"
        scols, counts = pool_columns(scols, counts, pool_groups;
                                     fallback_label=pool_fallback)
    end

    (; segment_labels, sample_names=scols, counts, effective_top_n)
end

# The column a chart groups rows by: a rank of the request's classifier, or a
# category set's tags (backfilled for older runs). Returns the column name or an
# HTTP error response.
function _label_column(con, table::String, columns::Vector{String},
                       tag::String, value::String, tag_source::String)
    if tag == "rank"
        rank = if isempty(value)
            levels = taxonomy_levels(con, table)
            isempty(levels) && return json_error(400, "no_taxonomy", "No taxonomy columns found")
            last(levels)
        else
            value
        end
        col = _rank_column(columns, rank)
        isnothing(col) || return col
        src = _analysis_source()
        return json_error(400, "bad_rank", isnothing(src) ? "Unknown rank '$rank'" :
                                           "Table '$table' has no $src taxonomy at rank '$rank'")
    end
    Categories.ensure_columns!(con, table, tag_source, [value];
                               library=_chart_library(), suffixed=_suffixed())
    col = _category_column(value)
    # Still absent when the set is missing from the library (ensure_columns! skips it).
    col in _duckdb_columns(con, table) ||
        return json_error(404, "category_set_not_found", "Category set '$value' not found")
    col
end

## Per-run unified chart
@post "/api/v1/studies/{study}/runs/{run}/analysis/chart" function(req,
                                                                    study::String,
                                                                    run::String)
    group = _req_group(req)
    dir, err = _require_run_results(study, run; group)
    !isnothing(err) && return err

    body = JSON3.read(String(req.body))
    table    = string(get(body, :table, "merged"))
    tag      = string(get(body, :tag, "rank"))
    value    = string(get(body, :value, ""))
    top_n    = Int(get(body, :top_n, 15))
    relative = Bool(get(body, :relative, true))
    mode     = string(get(body, :mode, "stacked"))
    keep_empty = Bool(get(body, :keep_empty, false))
    subgroup = let v = get(body, :subgroup, nothing)
        (isnothing(v) || v === nothing) ? nothing : string(v)
    end
    params = _body_filter_params(body)
    cfg_cache = Dict{Tuple{String,String,String}, Any}()

    # Resolve sub-group prefixes for pooling (used when subgroup == "__pool__").
    run_data_dir = joinpath(ServerState.data_dir(), study,
                            isnothing(group) ? run : joinpath(group, run))
    pool_groups = _is_pooled(run_data_dir) ? _subgroup_names(run_data_dir) : String[]

    # Only the category branch backfills Category__ columns; rank charts read-only.
    _with_analysis_results_table(study, run, table; group,
                                 readonly = tag != "category") do con, columns
        all_scols = sample_columns(con, table)
        isempty(all_scols) && return json_error(400, "no_samples",
                                                "No sample columns found in table '$table'")
        # Category-tagged rendering is the composition surface; rank is taxa.
        surface = tag == "category" ? "composition" : "taxa"
        exclude_conditions = _exclusion_conditions(study, run, group, columns; surface, config_cache=cfg_cache)
        where_clause, where_params = _analysis_where_clause(params, columns;
                                                             exclude_conditions)
        all_scols = _retain_sample_columns(con, table, all_scols, params,
                                           where_clause, where_params)
        isempty(all_scols) && return json_error(400, "no_samples",
                                                "No samples pass the sample read-count filter")

        # Source resolves the taxonomy side for category backfill; rank charts ignore it.
        tag_source = tag == "category" ? _tagging_source(study, run; group) : ""
        result = _chart_data(con, table, columns, all_scols,
                             where_clause, where_params,
                             tag, value, subgroup,
                             pool_groups, run, top_n, tag_source)
        result isa HTTP.Response && return result

        # Composition bars plot raw read sums: see the note on _chart_data.
        counts = result.counts
        scols  = result.sample_names

        colour_for = (tag == "category") ? _category_colour_for(value) : nothing
        # Pooled bars are sub-group names, not samples: nothing to trim.
        fig = bar_chart(result.segment_labels, scols, counts;
                        top_n=result.effective_top_n, relative, mode, colour_for,
                        keep_empty, short_labels = subgroup != "__pool__")
        HTTP.Response(200, ["Content-Type" => "application/json"],
                      body=JSON3.write(fig))
    end
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

    table = get(HTTP.queryparams(req), "table", "merged")
    _with_analysis_results_table(study, run, table; group) do con, columns
        # Only the ranks the requested classifier actually has.
        json(filter(r -> !isnothing(_rank_column(columns, r)), taxonomy_levels(con, table)))
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
    source = _opt_string(run_spec, :source)
    # Names come from the request body, so they must name a run of this study.
    isnothing(group) || (_valid_name(group) && group in _group_names(study)) || return nothing
    run_name in _all_run_names(study) || return nothing
    dir = _require_duckdb(study, run_name; group)
    isnothing(dir) && return nothing
    label = if !isnothing(prefix)
        prefix
    elseif !isnothing(group)
        "$group/$run_name"
    else
        run_name
    end
    (; run=run_name, group, dir, label, prefix, source)
end

"""Filter sample columns to those matching a prefix (e.g. "SubgroupName_")."""
function _filter_by_prefix(scols::Vector{String}, prefix::Union{String,Nothing})
    isnothing(prefix) && return scols
    pfx = prefix * "_"
    filter(c -> startswith(c, pfx), scols)
end

# Keep columns matching any of several sub-group prefixes, preserving input order.
# An empty prefix list means no restriction, mirroring the `nothing` case above.
function _filter_by_prefix(scols::Vector{String}, prefixes::Vector{String})
    isempty(prefixes) && return scols
    pfxs = [p * "_" for p in prefixes]
    filter(c -> any(p -> startswith(c, p), pfxs), scols)
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

# Resolve a run's config cascade, memoised in `config_cache` by (study, run, group).
function _resolved_run_config(study::String, run, group, config_cache)
    key = (study, something(run, ""), something(group, ""))
    if !isnothing(config_cache) && haskey(config_cache, key)
        return config_cache[key]
    end
    cfg = _resolve_config(study, run, group)
    !isnothing(config_cache) && (config_cache[key] = cfg)
    cfg
end

## Category-based figure exclusion
# The chart surfaces an exclusion may apply to. "diversity" covers alpha, NMDS
# and PERMANOVA; "taxa" is the rank-tagged chart; "composition" is the
# category-tagged chart; "venn" the presence/absence sets; "differential" the
# per-taxon differential abundance test.
const _EXCLUSION_SURFACES = Set(["diversity", "taxa", "composition", "venn", "differential"])

# Resolve the configured `analysis.exclude_categories` for a run into a vector of
# (set, category, surfaces) specs. Each names a composition category set and the
# single category within it whose taxa should be dropped. `surfaces` is the set
# of chart surfaces the spec applies to, parsed from an optional `apply_to` list;
# when `apply_to` is absent the spec applies to every surface (`surfaces` is
# `nothing`). Unknown surface names are dropped with a warning.
function _analysis_exclude_categories(study::String;
                                      run::Union{String,Nothing}=nothing,
                                      group::Union{String,Nothing}=nothing,
                                      config_cache::Union{Dict,Nothing}=nothing)
    resolved = _resolved_run_config(study, run, group, config_cache)
    entry = get(resolved, "analysis.exclude_categories", (; value=[], source="default"))
    raw = get(entry, :value, [])
    specs = @NamedTuple{set::String, category::String, surfaces::Union{Nothing,Set{String}}}[]
    raw isa AbstractVector || return specs
    for item in raw
        item isa AbstractDict || continue
        set = string(get(item, "set", ""))
        category = string(get(item, "category", ""))
        (isempty(set) || isempty(category)) && continue
        surfaces = _parse_apply_to(get(item, "apply_to", nothing); set, category)
        push!(specs, (; set, category, surfaces))
    end
    specs
end

# Parse an `apply_to` value into a set of known surfaces, or `nothing` (all).
# A missing or malformed value means "all surfaces"; an explicit list is filtered
# to the known surfaces, and unknown entries are warned about and dropped.
function _parse_apply_to(raw; set::String, category::String)
    raw isa AbstractVector || return nothing
    surfaces = Set{String}()
    for s in raw
        name = lowercase(strip(string(s)))
        if name in _EXCLUSION_SURFACES
            push!(surfaces, name)
        else
            @warn "exclude_categories: unknown apply_to surface; ignoring" set category surface=name
        end
    end
    surfaces
end

# Build SQL WHERE fragments that drop rows classified into an excluded category,
# for the given chart `surface` (see `_EXCLUSION_SURFACES`). A spec applies only
# when its `surfaces` is `nothing` (all) or contains `surface`, so an exclusion
# can be kept out of, say, the composition chart while still cleaning diversity.
# The composition category CASE expression is evaluated inline and compared to
# the excluded category name, so no Category__ column need be materialised and
# the query stays read-only. Column references are bare (no table alias), which
# matches the context in which these fragments are appended. `strict=true` gates
# on exact realisability so a degraded filter drops nothing rather than the wrong
# rows; a warning is logged for each skipped spec.
function _exclusion_conditions(study::String, run, group, columns::Vector{String};
                               surface::String, config_cache::Union{Dict,Nothing}=nothing)
    isnothing(run) && return String[]
    specs = filter(s -> isnothing(s.surfaces) || surface in s.surfaces,
                   _analysis_exclude_categories(study; run, group, config_cache))
    isempty(specs) && return String[]
    source = _tagging_source(study, run; group)
    lib = _chart_library()
    filters = lib["filters"]
    col_set = Set(columns)
    conds = String[]
    for spec in specs
        cfg = get(lib["sets"], spec.set, nothing)
        isnothing(cfg) && continue
        cats = get(cfg, "categories", [])
        case = Categories.category_case_when(cats, col_set, source;
                                             filters, table_alias="", strict=true)
        # Drop rows only when the set is exactly realisable AND the excluded
        # category is reachable in the CASE, either as its own branch (THEN) or
        # as the catch-all (ELSE); otherwise the condition would silently exclude
        # nothing or the wrong rows.
        cat_lit = Categories._sql_str(spec.category)
        if isnothing(case) || !(occursin("THEN $cat_lit", case) || occursin("ELSE $cat_lit", case))
            @warn "exclude_categories: category not classifiable for this table; skipping exclusion" set=spec.set category=spec.category
            continue
        end
        push!(conds, "($case) != $cat_lit")
    end
    conds
end

# How alpha diversity normalises counts before richness, Shannon and Simpson.
function _alpha_normalisation(study::String;
                              run::Union{String,Nothing}=nothing,
                              group::Union{String,Nothing}=nothing,
                              config_cache::Union{Dict,Nothing}=nothing)
    resolved = _resolved_run_config(study, run, group, config_cache)
    val(k, d) = something(get(get(resolved, k, (; value=d)), :value, d), d)
    (; method = String(val("analysis.alpha.normalisation", "none")),
       depth = Int(val("analysis.alpha.depth", 0)),
       iterations = Int(val("analysis.alpha.rarefaction_iterations", 1)))
end

# How NMDS and PERMANOVA prepare counts for Bray-Curtis. hellinger is a
# transform; rarefy and srs are depth normalisation.
function _beta_normalisation(study::String;
                             run::Union{String,Nothing}=nothing,
                             group::Union{String,Nothing}=nothing,
                             config_cache::Union{Dict,Nothing}=nothing)
    resolved = _resolved_run_config(study, run, group, config_cache)
    val(k, d) = something(get(get(resolved, k, (; value=d)), :value, d), d)
    method = String(val("analysis.beta.normalisation", "hellinger"))
    method == "hellinger" ? (; method="none", depth=0, transform="hellinger") :
                            (; method, depth=Int(val("analysis.beta.depth", 0)), transform="none")
end

function _analysis_where_clause(params, columns::Vector{String};
                                exclude_conditions::Vector{String}=String[])
    where_clause, where_params = _build_where(params, columns)
    for cond in exclude_conditions
        where_clause = isempty(where_clause) ? "WHERE $cond" : "$where_clause AND $cond"
    end
    where_clause, where_params
end

# Open a run's merged results DB for analysis and invoke `f(con, columns)`.
# The merged table carries taxonomy ranks, Category__<set> columns and per-sample
# counts, so every analysis reads it directly; the per-source annotation DBs of
# the quarantined FuncDB workflow are not consulted. Pass `readonly=false` when
# the body backfills Category__ columns via `Categories.ensure_columns!`.
function _with_analysis_results_table(f::Function, study::String, run::String, table::String;
                                      group::Union{String,Nothing}=nothing,
                                      readonly::Bool=true)
    merge_dir = _require_duckdb(study, run; group)
    isnothing(merge_dir) && return json_error(404, "no_results",
        "No results for run '$run' - run the pipeline first")
    opener = readonly ? with_results_db : with_results_db_write
    opener(merge_dir) do con
        columns = _duckdb_columns(con, table)
        isempty(columns) && return json_error(404, "table_not_found",
                                              "Table '$table' not found in results database")
        f(con, columns)
    end
end

# Cross-run analogue: open the merged results DB of a run already resolved by
# `_resolve_run_duckdb` (which carries its `dir`), invoking `f(con, columns)`.
# Returns `nothing` when the table is absent so accumulating loops simply skip it.
function _with_resolved_results_table(f::Function, resolved, table::String;
                                      readonly::Bool=true)
    opener = readonly ? with_results_db : with_results_db_write
    opener(resolved.dir) do con
        columns = _duckdb_columns(con, table)
        isempty(columns) && return nothing
        f(con, columns)
    end
end

function _study_alpha_plot_config(study::String)
    resolved = _resolve_config(study)
    get_value(key, default) = get(get(resolved, key, (; value=default, source="default")), :value, default)
    (;
        show_points = Bool(get_value("analysis.alpha.show_points", true)),
        annotate_significance = Bool(get_value("analysis.alpha.annotate_significance", false)),
        pairwise_brackets = Bool(get_value("analysis.alpha.pairwise_brackets", false)),
        paired_samples = Bool(get_value("analysis.alpha.paired_samples", false)),
        paired_lines = Bool(get_value("analysis.alpha.paired_lines", false)),
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

function _expand_comparison_run_specs(study::String, runs_spec; aggregate::Bool=false)
    expanded = NamedTuple[]
    for spec in runs_spec
        run_name = string(get(spec, :run, ""))
        isempty(run_name) && continue
        group  = _opt_string(spec, :group)
        prefix = _opt_string(spec, :prefix)
        source = _opt_string(spec, :source)
        if aggregate
            push!(expanded, (; run=run_name, group, prefix=nothing, source))
            continue
        end

        if !isnothing(prefix)
            push!(expanded, (; run=run_name, group, prefix, source))
            continue
        end

        data_rel = isnothing(group) ? run_name : joinpath(group, run_name)
        data_dir = joinpath(ServerState.data_dir(), study, data_rel)
        if _is_pooled(data_dir)
            for subgroup in _subgroup_names(data_dir)
                push!(expanded, (; run=run_name, group, prefix=subgroup, source))
            end
        else
            push!(expanded, (; run=run_name, group, prefix=nothing, source))
        end
    end
    expanded
end

function _collect_cross_run_asvs(study::String, runs_spec, table, params; aggregate::Bool=false)
    run_data = Tuple{String, Vector{String}, DataFrame}[]
    group_labels = String[]
    run_labels = String[]
    individual_ids = String[]
    missing_seq = false
    cfg_cache = Dict{Tuple{String,String,String}, Any}()

    for spec in _expand_comparison_run_specs(study, runs_spec; aggregate)
        resolved = _resolve_run_duckdb(study, spec)
        isnothing(resolved) && continue
        _with_resolved_results_table(resolved, table) do con, columns
            scols = _filter_by_prefix(sample_columns(con, table), resolved.prefix)
            isempty(scols) && return

            seq_col = sequence_column_name(columns)
            if isnothing(seq_col)
                missing_seq = true
                return
            end

            exclude_conditions = _exclusion_conditions(study, resolved.run, resolved.group, columns; surface="diversity", config_cache=cfg_cache)
            where_clause, where_params = _analysis_where_clause(params, columns; exclude_conditions)
            scols = _retain_sample_columns(con, table, scols, params, where_clause, where_params)
            isempty(scols) && return
            seq_filter = isempty(where_clause) ?
                "WHERE \"$seq_col\" IS NOT NULL AND TRIM(\"$seq_col\") != ''" :
                "$where_clause AND \"$seq_col\" IS NOT NULL AND TRIM(\"$seq_col\") != ''"

            sum_exprs = join(["SUM(COALESCE(\"$c\", 0)) AS \"$c\"" for c in scols], ", ")
            sql = """
                SELECT "$seq_col" AS "sequence", $sum_exprs
                FROM "$table" $seq_filter
                GROUP BY "$seq_col"
                ORDER BY "$seq_col"
            """
            df = DataFrame(DBInterface.execute(con, sql, where_params))
            nrow(df) == 0 && return

            push!(run_data, (resolved.label, scols, df))
            group_label = _resolved_group_label(resolved)
            run_label = _resolved_run_label(resolved)
            append!(group_labels, fill(group_label, length(scols)))
            append!(run_labels, fill(run_label, length(scols)))
            append!(individual_ids, [_comparison_sample_id(s, resolved) for s in scols])
        end
    end

    run_data, group_labels, run_labels, individual_ids, missing_seq
end

function _drop_empty_samples(mat::Matrix{Float64}, aligned::AbstractVector...)
    keep = vec(sum(mat; dims=2) .> 0)
    filtered = Any[mat[keep, :]]
    for values in aligned
        push!(filtered, values[keep])
    end
    Tuple(filtered), count(keep), length(keep)
end

## Resolve the rarefaction depth shared across runs: a configured value wins;
## auto (0) resolves to the pooled minimum positive library size. Skipped for
## method "none", where normalise_counts ignores the depth entirely.
function _resolve_cross_run_depth(method::String, configured_depth::Int, matrices)
    (method == "none" || configured_depth != 0) && return configured_depth
    auto_min_depth(Iterators.flatten(vec(sum(m; dims=2)) for m in matrices))
end

function _apply_normalisation(mat::Matrix{Float64}, method::String, depth::Int,
                              seed::Int, aligned::AbstractVector...)
    norm = normalise_counts(mat; method, depth, seed)
    filtered = Any[norm.mat]
    for values in aligned
        push!(filtered, values[norm.kept])
    end
    Tuple(filtered), length(norm.kept)
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

    aggregate = Bool(get(body, :aggregate, false))
    expanded_specs = _expand_comparison_run_specs(study, runs_spec; aggregate)
    resolved_specs = filter(!isnothing, [_resolve_run_duckdb(study, spec) for spec in expanded_specs])
    isempty(resolved_specs) && return json_error(400, "no_data", "No data found for any run")

    compare_mode = length(unique(_resolved_run_label(r) for r in resolved_specs)) > 1 ? :run : :subgroup
    cfg_cache = Dict{Tuple{String,String,String}, Any}()
    norm = _alpha_normalisation(study)

    # First pass: collect raw matrices (fetched inside annotation DB closures)
    collected = Tuple{Matrix{Float64}, Vector{String}, Any}[]
    for resolved in resolved_specs
        _with_resolved_results_table(resolved, table) do con, columns
            scols = _filter_by_prefix(sample_columns(con, table), resolved.prefix)
            isempty(scols) && return
            exclude_conditions = _exclusion_conditions(study, resolved.run, resolved.group, columns; surface="diversity", config_cache=cfg_cache)
            where_clause, where_params = _analysis_where_clause(params, columns; exclude_conditions)
            scols = _retain_sample_columns(con, table, scols, params, where_clause, where_params)
            isempty(scols) && return
            mat = filtered_counts(con, table, scols, where_clause, where_params)
            push!(collected, (mat, scols, resolved))
        end
    end

    isempty(collected) && return json_error(400, "no_data", "No data found for any run")

    norm_depth = _resolve_cross_run_depth(norm.method, norm.depth,
                                          (mat for (mat, _, _) in collected))

    seed = _study_seed(study)

    # Second pass: normalise each run's matrix and compute diversity metrics
    groups = Dict{String, NamedTuple{(:sample_ids, :richness, :shannon, :simpson),
                                     Tuple{Vector{String}, Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()
    for (raw_mat, raw_scols, resolved) in collected
        alpha = alpha_diversity(raw_mat; norm.method, depth=norm_depth, seed,
                                norm.iterations)
        scols = raw_scols[alpha.kept]
        isempty(scols) && continue

        r, sh, si = alpha.richness, alpha.shannon, alpha.simpson
        sample_ids = [_comparison_pair_key(s, resolved; mode=compare_mode) for s in scols]
        label = compare_mode == :run ? _resolved_run_label(resolved) : _resolved_group_label(resolved)
        bucket = get!(groups, label, (; sample_ids=String[], richness=Float64[], shannon=Float64[], simpson=Float64[]))
        append!(bucket.sample_ids, sample_ids)
        append!(bucket.richness, r)
        append!(bucket.shannon, sh)
        append!(bucket.simpson, si)
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
                        paired_lines=alpha_cfg.paired_lines,
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
    aggregate = Bool(get(body, :aggregate, false))
    (norm_method, norm_depth, transform) = _beta_normalisation(study)

    run_data, group_labels, run_name_labels, _, missing_seq = _collect_cross_run_asvs(study, runs_spec, table, params; aggregate)

    missing_seq && return json_error(400, "no_sequence_data",
                                         "Table '$table' does not contain a DNA sequence column required for cross-run NMDS/PERMANOVA (expected `sequence` or `Sequence`, as in `merged` or `asv_counts`)")
    length(run_data) < 2 && return json_error(400, "too_few_runs",
                                                   "Need data from at least 2 runs for NMDS")
    mat, all_samples, _, _ = combined_asv_counts_across_runs(run_data)
    (filtered, kept, total) = _drop_empty_samples(mat, all_samples, group_labels, run_name_labels)
    mat, all_samples, group_labels, run_name_labels = filtered
    kept < total && @info "NMDS: dropped $(total - kept) empty samples before ordination"
    (filtered, nkept) = _apply_normalisation(mat, norm_method, norm_depth,
                                             _study_seed(study),
                                             all_samples, group_labels, run_name_labels)
    mat, all_samples, group_labels, run_name_labels = filtered
    nkept < kept &&
        @info "NMDS: dropped $(kept - nkept) samples below normalisation depth"
    size(mat, 1) < 3 && return json_error(400, "too_few_samples",
                                               "Need at least 3 non-empty samples for NMDS")

    coords, stress = run_nmds(mat; seed=_study_seed(study), transform)
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
    aggregate = Bool(get(body, :aggregate, false))
    (norm_method, norm_depth, transform) = _beta_normalisation(study)

    run_data, group_labels, run_labels, individual_ids, missing_seq = _collect_cross_run_asvs(study, runs_spec, table, params; aggregate)

    missing_seq && return json_error(400, "no_sequence_data",
                                         "Table '$table' does not contain a DNA sequence column required for cross-run NMDS/PERMANOVA (expected `sequence` or `Sequence`, as in `merged` or `asv_counts`)")
    length(run_data) < 2 && return json_error(400, "too_few_runs",
                                                   "Need data from at least 2 runs")
    mat, all_samples, _, _ = combined_asv_counts_across_runs(run_data)
    (filtered, kept, total) = _drop_empty_samples(mat, all_samples, group_labels, run_labels, individual_ids)
    mat, all_samples, group_labels, run_labels, individual_ids = filtered
    kept < total && @info "PERMANOVA: dropped $(total - kept) empty samples before analysis"
    (filtered, nkept) = _apply_normalisation(mat, norm_method, norm_depth,
                                             _study_seed(study),
                                             all_samples, group_labels, run_labels, individual_ids)
    mat, all_samples, group_labels, run_labels, individual_ids = filtered
    nkept < kept &&
        @info "PERMANOVA: dropped $(kept - nkept) samples below normalisation depth"
    size(mat, 1) < 3 && return json_error(400, "too_few_samples",
                                               "Need at least 3 non-empty samples for PERMANOVA")

    # Samples from the same individual are not exchangeable, so permutations
    # stay within each individual's block.
    metadata = DataFrame(sample=all_samples, group=group_labels, run=run_labels)
    result = run_permanova(mat, metadata; seed=_study_seed(study), transform,
                           blocks=individual_ids)
    isnothing(result) && return json_error(500, "permanova_failed", "PERMANOVA computation failed")
    hasproperty(result, :message) && return json_error(500, "permanova_failed",
                                                           "PERMANOVA computation failed: $(result.message)")
    json(result)
end

## Cross-run unified chart
@post "/api/v1/studies/{study}/analysis/chart" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    body = JSON3.read(String(req.body))
    runs_spec = get(body, :runs, [])
    isempty(runs_spec) && return json_error(400, "no_runs", "Provide at least one run")
    table    = string(get(body, :table, "merged"))
    tag      = string(get(body, :tag, "rank"))
    value    = string(get(body, :value, ""))
    top_n    = Int(get(body, :top_n, 15))
    relative = Bool(get(body, :relative, true))
    mode     = string(get(body, :mode, "stacked"))
    keep_empty = Bool(get(body, :keep_empty, false))
    subgroup = let v = get(body, :subgroup, nothing)
        (isnothing(v) || v === nothing) ? nothing : string(v)
    end
    params   = _body_filter_params(body)
    cfg_cache = Dict{Tuple{String,String,String}, Any}()

    # Collect per-run aggregated counts.
    collected = Tuple{String, Vector{String}, DataFrame}[]
    first_error = Ref{Union{HTTP.Response,Nothing}}(nothing)
    effective_top_n = (tag == "category") ? typemax(Int) : top_n

    for spec in runs_spec
        resolved = _resolve_run_duckdb(study, spec)
        isnothing(resolved) && continue
        # Only the category branch backfills Category__ columns; rank charts read-only.
        _with_resolved_results_table(resolved, table; readonly = tag != "category") do con, columns

            # Scope sample columns to the requested subgroup prefix (not "__pool__" here).
            prefix_scope = (isnothing(subgroup) || subgroup == "__pool__") ?
                            resolved.prefix : subgroup
            scols = _filter_by_prefix(sample_columns(con, table), prefix_scope)
            isempty(scols) && return nothing
            surface = tag == "category" ? "composition" : "taxa"
            exclude_conditions = _exclusion_conditions(study, resolved.run, resolved.group, columns; surface, config_cache=cfg_cache)
            where_clause, where_params = _analysis_where_clause(params, columns;
                                                                exclude_conditions)
            scols = _retain_sample_columns(con, table, scols, params,
                                           where_clause, where_params)
            isempty(scols) && return nothing

            label_col = _label_column(con, table, columns, tag, value,
                                      _tagging_source(study, resolved.run; group=resolved.group))
            if label_col isa HTTP.Response
                isnothing(first_error[]) && (first_error[] = label_col)
                return nothing
            end

            agg = aggregate_by_taxon(con, table, scols, label_col,
                                     where_clause, where_params)
            nrow(agg) == 0 && return nothing

            # Zero-read guard: drop label rows that are all-zero across scols.
            row_totals = vec(sum(Matrix{Float64}(agg[:, scols]), dims=2))
            keep_rows  = row_totals .> 0
            if !all(keep_rows)
                agg = agg[keep_rows, :]
            end
            nrow(agg) == 0 && return nothing

            label = _resolved_group_label(resolved)
            push!(collected, (label, scols, agg))
            nothing
        end
    end

    isempty(collected) && return something(first_error[], json_error(400, "no_data", "No data found for any run"))

    # Combine across runs via the existing helper.
    mat, all_samples, segment_labels, run_labels = combined_counts_across_runs(collected)

    # Transpose to segments x samples. Raw read sums: see the note on _chart_data.
    counts = permutedims(mat)

    # Zero-read guard on the combined matrix.
    row_totals = vec(sum(counts, dims=2))
    keep_rows  = row_totals .> 0
    if !all(keep_rows)
        segment_labels = segment_labels[keep_rows]
        counts         = counts[keep_rows, :]
    end

    # Pool: cross-run always pools by run label when subgroup == "__pool__" or by default.
    # Pool on each sample's own run label rather than `pool_columns`' name-prefix
    # match: only a pooled run's sub-group columns carry a "<label>_" prefix, so
    # the prefix match sent every sample of an ordinary run to "Other" and
    # folded all the runs being compared into that one indistinguishable bar.
    unique_labels = unique(run_labels)
    pooled = isnothing(subgroup) || subgroup == "__pool__"
    if pooled
        scols  = unique_labels
        counts = reduce(hcat, [vec(sum(counts[:, run_labels .== l], dims=2))
                               for l in unique_labels])
    else
        scols = all_samples
    end

    # Pooled bars keep the order the runs were selected in and their full labels.
    colour_for = (tag == "category") ? _category_colour_for(value) : nothing
    fig = bar_chart(segment_labels, scols, counts;
                    top_n=effective_top_n, relative, mode, colour_for, keep_empty,
                    natural_order = !pooled, short_labels = !pooled)
    HTTP.Response(200, ["Content-Type" => "application/json"],
                  body=JSON3.write(fig))
end

## Faceted composition grid
# The dimensions a facet axis can be keyed on. "run" and "group" are the study's
# own structure; "subgroup" is the per-prefix split of a pooled run.
const _FACET_DIMENSIONS = Set(["run", "group", "subgroup"])

const _FACET_UNSET = Dict("group" => "(ungrouped)", "subgroup" => "(all)")

# The value a resolved run spec takes on one facet dimension. A spec that carries
# no group or prefix still needs a cell to live in, so it gets an explicit
# placeholder rather than being dropped from the grid.
function _facet_value(resolved, dimension::String)
    if dimension == "run"
        resolved.run
    elseif dimension == "group"
        something(resolved.group, _FACET_UNSET["group"])
    else
        something(resolved.prefix, _FACET_UNSET["subgroup"])
    end
end

# Fold the run specs that share a facet cell into one panel. Sample column names
# are unique only within their own run, so once a cell draws on more than one
# spec the names are qualified with the spec label; left bare, two runs' samples
# that happen to share a name would collapse onto one categorical slot.
#
# `counts` comes back samples-by-segments, the orientation
# `combined_counts_across_runs` produces; the caller transposes it to the
# segments-by-samples that `faceted_bar_chart` plots.
function _facet_panel(row::String, col::String,
                      entries::Vector{Tuple{String, Vector{String}, DataFrame}})
    mat, sample_names, segment_labels, entry_labels = combined_counts_across_runs(entries)
    qualified = length(entries) > 1 ?
        ["$(entry_labels[i])/$(sample_names[i])" for i in eachindex(sample_names)] :
        sample_names
    (; row, col, segment_labels, sample_names=qualified, counts=mat)
end

@post "/api/v1/studies/{study}/analysis/chart-facet" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    body = JSON3.read(String(req.body))
    runs_spec = get(body, :runs, [])
    isempty(runs_spec) && return json_error(400, "no_runs", "Provide at least one run")

    table      = string(get(body, :table, "merged"))
    tag        = string(get(body, :tag, "rank"))
    value      = string(get(body, :value, ""))
    top_n      = Int(get(body, :top_n, 15))
    relative   = Bool(get(body, :relative, true))
    mode       = string(get(body, :mode, "stacked"))
    keep_empty = Bool(get(body, :keep_empty, false))
    row_dim    = string(get(body, :rows, "run"))
    col_dim    = string(get(body, :cols, "subgroup"))

    row_dim in _FACET_DIMENSIONS || return json_error(400, "bad_facet",
        "Unknown row facet '$row_dim' (expected run, group or subgroup)")
    col_dim in _FACET_DIMENSIONS || return json_error(400, "bad_facet",
        "Unknown column facet '$col_dim' (expected run, group or subgroup)")
    row_dim == col_dim && return json_error(400, "bad_facet",
        "Row and column facets must differ; both are '$row_dim'")

    params = _body_filter_params(body)
    cfg_cache = Dict{Tuple{String,String,String}, Any}()
    effective_top_n = (tag == "category") ? typemax(Int) : top_n

    # Collect each expanded spec's aggregated counts into its facet cell.
    cells = Dict{Tuple{String,String}, Vector{Tuple{String, Vector{String}, DataFrame}}}()
    first_error = Ref{Union{HTTP.Response,Nothing}}(nothing)
    for spec in _expand_comparison_run_specs(study, runs_spec)
        resolved = _resolve_run_duckdb(study, spec)
        isnothing(resolved) && continue
        _with_resolved_results_table(resolved, table; readonly = tag != "category") do con, columns
            scols = _filter_by_prefix(sample_columns(con, table), resolved.prefix)
            isempty(scols) && return nothing

            surface = tag == "category" ? "composition" : "taxa"
            exclude_conditions = _exclusion_conditions(study, resolved.run, resolved.group,
                                                       columns; surface, config_cache=cfg_cache)
            where_clause, where_params = _analysis_where_clause(params, columns;
                                                                exclude_conditions)
            scols = _retain_sample_columns(con, table, scols, params,
                                           where_clause, where_params)
            isempty(scols) && return nothing

            label_col = _label_column(con, table, columns, tag, value,
                                      _tagging_source(study, resolved.run; group=resolved.group))
            if label_col isa HTTP.Response
                isnothing(first_error[]) && (first_error[] = label_col)
                return nothing
            end

            agg = aggregate_by_taxon(con, table, scols, label_col, where_clause, where_params)
            nrow(agg) == 0 && return nothing

            cell = (_facet_value(resolved, row_dim), _facet_value(resolved, col_dim))
            push!(get!(cells, cell, Tuple{String, Vector{String}, DataFrame}[]),
                  (resolved.label, scols, agg))
            nothing
        end
    end

    isempty(cells) && return something(first_error[], json_error(400, "no_data", "No data found for any run"))

    panels = [_facet_panel(row, col, entries) for ((row, col), entries) in cells]

    # Raw read sums, panel by panel: see the note on _chart_data. Every panel is
    # on the one read-count scale already, which is what makes a grid comparable.
    prepared = NamedTuple[]
    for p in panels
        isempty(p.sample_names) && continue
        counts = permutedims(p.counts)               # segments x samples
        push!(prepared, (; p.row, p.col, p.segment_labels, p.sample_names, counts))
    end
    isempty(prepared) && return json_error(400, "no_samples",
                                           "No samples in any facet panel")

    row_labels = natural_sort(unique(String[p.row for p in prepared]))
    col_labels = natural_sort(unique(String[p.col for p in prepared]))

    colour_for = (tag == "category") ? _category_colour_for(value) : nothing
    fig = faceted_bar_chart(prepared, row_labels, col_labels;
                            top_n=effective_top_n, relative, mode, colour_for,
                            keep_empty, row_title=row_dim, col_title=col_dim,
                            short_labels=true)
    HTTP.Response(200, ["Content-Type" => "application/json"],
                  body=JSON3.write(fig))
end

## Cross-run taxon Venn sets
"""
    _venn_condition_labels(conditions) -> Vector{String}

Build a user-facing label for each condition in the Venn diagram, including
only the identity components (group / run / prefix) that actually differ
across the selected conditions. This avoids ambiguity when e.g. multiple
groups share the same pooled-run subgroup prefixes.
"""
function _venn_condition_labels(conditions)
    groups_differ   = length(unique(c.group   for c in conditions)) > 1
    runs_differ     = length(unique(c.run     for c in conditions)) > 1
    prefixes_differ = length(unique(c.prefix  for c in conditions)) > 1
    labels = String[]
    for c in conditions
        parts = String[]
        groups_differ   && !isnothing(c.group)  && push!(parts, c.group)
        runs_differ                              && push!(parts, c.run)
        prefixes_differ && !isnothing(c.prefix) && push!(parts, c.prefix)
        if isempty(parts)
            # Nothing varies (only one condition); fall back to the most specific identifier.
            push!(parts, something(c.prefix, c.run))
        end
        push!(labels, join(parts, "/"))
    end
    labels
end

@post "/api/v1/studies/{study}/analysis/venn" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    body   = JSON3.read(String(req.body))
    runs_spec = get(body, :runs, [])
    isempty(runs_spec) && return json_error(400, "no_runs", "Provide at least one run")
    table  = string(get(body, :table, "merged"))
    rank   = string(get(body, :rank,  "Genus"))
    params = _body_filter_params(body)

    # Collect each condition with its identity components so we can compute
    # disambiguated labels after seeing the full set of selected conditions.
    aggregate = Bool(get(body, :aggregate, false))
    conditions = NamedTuple{
        (:group, :run, :prefix, :taxa),
        Tuple{Union{String,Nothing}, String, Union{String,Nothing}, Vector{String}}
    }[]
    cfg_cache  = Dict{Tuple{String,String,String}, Any}()

    for spec in _expand_comparison_run_specs(study, runs_spec; aggregate)
        resolved = _resolve_run_duckdb(study, spec)
        isnothing(resolved) && continue

        _with_resolved_results_table(resolved, table) do con, columns
            levels = taxonomy_levels(con, table)
            isempty(levels) && return
            effective_rank = rank in levels ? rank : last(levels)
            rank_col = _rank_column(columns, effective_rank)
            isnothing(rank_col) && return

            scols = _filter_by_prefix(sample_columns(con, table), resolved.prefix)
            isempty(scols) && return

            exclude_conditions = _exclusion_conditions(study, resolved.run, resolved.group, columns; surface="venn", config_cache=cfg_cache)
            where_clause, where_params = _analysis_where_clause(params, columns;
                                                                exclude_conditions)
            scols = _retain_sample_columns(con, table, scols, params,
                                           where_clause, where_params)
            isempty(scols) && return nothing
            taxa = venn_taxa_present(con, table, scols, rank_col,
                                     where_clause, where_params)
            push!(conditions, (; group=resolved.group, run=resolved.run,
                                  prefix=resolved.prefix, taxa))
        end
    end

    length(conditions) < 2 && return json_error(400, "too_few_conditions",
        "Need data from at least 2 conditions for the taxon overlap diagram")

    labels = _venn_condition_labels(conditions)
    sets = [(; name=labels[i], taxa=conditions[i].taxa) for i in eachindex(conditions)]

    json(Dict("sets" => sets, "rank" => rank))
end

## Differential abundance between two conditions
"""
    _differential_config(study) -> DifferentialConfig

The study's `analysis.differential` settings, resolved through the config
cascade.
"""
function _differential_config(study::String)
    resolved = _resolved_run_config(study, nothing, nothing, nothing)
    val(k, d) = something(get(get(resolved, k, (; value=d)), :value, d), d)
    DifferentialConfig(String(val("analysis.differential.offset", "tss")),
                       Float64(val("analysis.differential.min_prevalence", 0.0));
                       method=String(val("analysis.differential.method", "nb_glm")),
                       replacement_delta=Float64(val("analysis.differential.replacement_delta", 0.65)))
end

"""
    _differential(study, body) -> HTTP.Response

Run a differential abundance test between exactly two conditions (runs or
pooled-run sub-groups) and return the per-taxon table with a volcano plot. The
first condition is the reference, so a positive estimate means more abundant in
the second. Every refusal is an explicit JSON error.
"""
function _differential(study::String, body)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    runs_spec = get(body, :runs, [])
    table     = string(get(body, :table, "merged"))
    rank      = string(get(body, :rank,  "Genus"))
    params    = _body_filter_params(body)
    aggregate = Bool(get(body, :aggregate, false))

    specs = _expand_comparison_run_specs(study, runs_spec; aggregate)
    length(specs) == 2 || return json_error(400, "two_conditions_required",
        "Differential abundance compares exactly 2 conditions; the selection gives $(length(specs))")

    config = try
        _differential_config(study)
    catch e
        e isa ArgumentError || rethrow()
        return json_error(400, "invalid_config", sprint(showerror, e))
    end

    conditions = NamedTuple[]
    effective_ranks = String[]
    for spec in specs
        label = something(spec.prefix, spec.run)
        resolved = _resolve_run_duckdb(study, spec)
        isnothing(resolved) && return json_error(404, "results_not_found",
            "No results database for '$label'")
        problem = _with_resolved_results_table(resolved, table) do con, columns
            levels = taxonomy_levels(con, table)
            isempty(levels) && return "table '$table' of '$label' has no taxonomy columns"
            effective_rank = rank in levels ? rank : last(levels)
            rank_col = _rank_column(columns, effective_rank)
            isnothing(rank_col) && return "rank '$effective_rank' is not a column of '$label'"
            scols = _filter_by_prefix(sample_columns(con, table), resolved.prefix)
            isempty(scols) && return "'$label' has no sample columns"
            exclude_conditions = _exclusion_conditions(study, resolved.run, resolved.group, columns;
                                                       surface="differential")
            where_clause, where_params = _analysis_where_clause(params, columns; exclude_conditions)
            scols = _retain_sample_columns(con, table, scols, params, where_clause, where_params)
            isempty(scols) && return "every sample of '$label' was removed by the read-count bounds"
            df = aggregate_by_taxon(con, table, scols, rank_col, where_clause, where_params)
            push!(conditions, (; group=resolved.group, run=resolved.run,
                                  prefix=resolved.prefix, scols, df))
            push!(effective_ranks, effective_rank)
            nothing
        end
        isnothing(problem) || return json_error(400, "condition_unusable", problem)
    end
    length(conditions) == 2 || return json_error(404, "results_not_found",
        "Results table '$table' is missing for at least one condition")
    length(unique(effective_ranks)) == 1 || return json_error(400, "rank_mismatch",
        "The two conditions resolved to different ranks: $(join(effective_ranks, ", "))")

    labels = _venn_condition_labels(conditions)
    labels[1] == labels[2] && return json_error(400, "identical_conditions",
        "Both conditions are '$(labels[1])'; choose two different runs or sub-groups")
    run_data = Tuple{String, Vector{String}, DataFrame}[
        (labels[i], conditions[i].scols, conditions[i].df) for i in 1:2]
    mat, samples, taxa, groups = combined_counts_across_runs(run_data)

    result = try
        differential_abundance(mat, samples, taxa, groups;
                               reference=labels[1], contrast=labels[2], config)
    catch e
        e isa ArgumentError    && return json_error(400, "invalid_input", sprint(showerror, e))
        e isa ScalingRefusal   && return json_error(422, "scaling_refused", sprint(showerror, e))
        e isa MASSUnavailable  && return json_error(503, "r_unavailable", sprint(showerror, e))
        e isa ErrorException   && return json_error(422, "no_taxon_fitted", e.msg)
        rethrow()
    end
    result["rank"] = first(effective_ranks)
    result["table"] = table
    result["figure"] = volcano_chart(result)
    HTTP.Response(200, ["Content-Type" => "application/json"],
                  body=JSON3.write(result))
end

## Cross-run differential abundance (negative binomial + BH)
@post "/api/v1/studies/{study}/analysis/differential" function(req, study::String)
    _differential(study, JSON3.read(String(req.body)))
end
