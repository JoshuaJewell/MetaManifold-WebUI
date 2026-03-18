module Analysis

# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

using DataFrames, JSON3, DuckDB, DBInterface
using ..DiversityMetrics: richness, shannon, simpson

export sample_columns, filtered_counts, filtered_df, taxonomy_levels,
       aggregate_by_taxon, combined_counts_across_runs,
       alpha_chart, taxa_bar_chart, pipeline_stats_chart,
       alpha_boxplot, nmds_chart,
       run_nmds, run_permanova, r_available

## DuckDB query helpers for analysis
"""
    sample_columns(con, table) -> Vector{String}

Identify per-sample count columns in a DuckDB table.
Returns all numeric columns except known non-count names.
"""
function sample_columns(con, table::String)
    result = DataFrame(DBInterface.execute(con,
        "SELECT column_name, data_type FROM information_schema.columns WHERE table_name = ?", [table]))
    non_count = Set([
        "SeqName", "Pident", "Accession", "rRNA", "Organellum", "specimen",
        "sequence", "OTU", "ASV",
        "Domain", "Supergroup", "Division", "Subdivision",
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    ])
    numeric_types = Set(["BIGINT", "INTEGER", "DOUBLE", "FLOAT", "HUGEINT", "SMALLINT", "TINYINT"])
    cols = String[]
    for row in eachrow(result)
        col = string(row.column_name)
        dtype = string(row.data_type)
        col in non_count && continue
        endswith(col, "_dada2") && continue
        endswith(col, "_boot") && continue
        endswith(col, "_vsearch") && continue
        dtype in numeric_types || continue
        push!(cols, col)
    end
    cols
end

"""
    taxonomy_levels(con, table) -> Vector{String}

Detect which taxonomy rank columns exist in the table.
"""
function taxonomy_levels(con, table::String)
    result = DataFrame(DBInterface.execute(con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = ?", [table]))
    existing = Set(string.(result.column_name))
    known_ranks = [
        "Domain", "Supergroup", "Division", "Subdivision",
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    ]
    filter(r -> r in existing, known_ranks)
end

"""
    filtered_counts(con, table, sample_cols, where_clause, where_params) -> Matrix{Float64}

Return a samples-by-features count matrix from DuckDB with filters applied.
Rows = samples, columns = ASV/OTU rows matching the filter.
"""
function filtered_counts(con, table::String, sample_cols::Vector{String},
                         where_clause::String, where_params::Vector)
    col_sql = join(["\"$c\"" for c in sample_cols], ", ")
    sql = "SELECT $col_sql FROM \"$table\" $where_clause"
    result = DataFrame(DBInterface.execute(con, sql, where_params))
    nrow(result) == 0 && return zeros(Float64, length(sample_cols), 0)
    n_features = nrow(result)
    n_samples = length(sample_cols)
    mat = zeros(Float64, n_samples, n_features)
    for (j, row) in enumerate(eachrow(result))
        for (i, col) in enumerate(sample_cols)
            v = row[Symbol(col)]
            mat[i, j] = ismissing(v) ? 0.0 : Float64(v)
        end
    end
    mat
end

"""
    filtered_df(con, table, where_clause, where_params) -> DataFrame

Return the full filtered table as a DataFrame.
"""
function filtered_df(con, table::String,
                     where_clause::String, where_params::Vector)
    sql = "SELECT * FROM \"$table\" $where_clause"
    DataFrame(DBInterface.execute(con, sql, where_params))
end

"""
    aggregate_by_taxon(con, table, sample_cols, rank, where_clause, where_params) -> DataFrame

Group rows by a taxonomy rank, summing sample counts.
"""
function aggregate_by_taxon(con, table::String, sample_cols::Vector{String},
                            rank::String, where_clause::String, where_params::Vector)
    sum_exprs = join(["SUM(COALESCE(\"$c\", 0)) AS \"$c\"" for c in sample_cols], ", ")
    sql = """
        SELECT COALESCE(NULLIF(TRIM("$rank"), ''), 'Unclassified') AS taxon, $sum_exprs
        FROM "$table" $where_clause
        GROUP BY taxon
        ORDER BY ($(join(["SUM(COALESCE(\"$c\", 0))" for c in sample_cols], " + "))) DESC
    """
    DataFrame(DBInterface.execute(con, sql, where_params))
end

"""
    combined_counts_across_runs(run_data) -> (Matrix{Float64}, Vector{String}, Vector{String}, Vector{String})

Build a combined taxonomy-aggregated count matrix across multiple runs for NMDS/PERMANOVA.
`run_data` is a vector of `(label, sample_cols, taxon_counts_df)` tuples.

Returns (matrix, all_sample_names, taxon_labels, run_labels_per_sample).
"""
function combined_counts_across_runs(
    run_data::Vector{Tuple{String, Vector{String}, DataFrame}}
)
    taxa_counts = Dict{String, Dict{String, Float64}}()
    all_samples = String[]
    run_labels = String[]

    for (label, sample_cols, df) in run_data
        append!(all_samples, sample_cols)
        append!(run_labels, fill(label, length(sample_cols)))
        for row in eachrow(df)
            taxon = string(row.taxon)
            td = get!(taxa_counts, taxon, Dict{String, Float64}())
            for col in sample_cols
                v = row[Symbol(col)]
                td[col] = get(td, col, 0.0) + (ismissing(v) ? 0.0 : Float64(v))
            end
        end
    end

    taxa_labels = sort(collect(keys(taxa_counts)))
    n_samples = length(all_samples)
    n_taxa = length(taxa_labels)
    mat = zeros(Float64, n_samples, n_taxa)

    for (j, taxon) in enumerate(taxa_labels)
        td = taxa_counts[taxon]
        for (i, sample) in enumerate(all_samples)
            mat[i, j] = get(td, sample, 0.0)
        end
    end

    mat, all_samples, taxa_labels, run_labels
end

## Plotly chart builders
# All functions return a Dict suitable for JSON3.write -> Plotly JSON.

## Colour palette
function _palette_hex(n::Int)::Vector{String}
    base = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7"]
    n <= length(base) && return base[1:n]
    colours = copy(base)
    for i in (length(base) + 1):n
        hue = mod(i * 137.508, 360.0)
        c = 0.85 * 0.65
        x = c * (1.0 - abs(mod(hue / 60.0, 2.0) - 1.0))
        m = 0.85 - c
        r, g, b = hue < 60  ? (c + m, x + m, m) :
                hue < 120 ? (x + m, c + m, m) :
                hue < 180 ? (m, c + m, x + m) :
                hue < 240 ? (m, x + m, c + m) :
                hue < 300 ? (x + m, m, c + m) :
                            (c + m, m, x + m)
        hex(v) = lpad(string(round(Int, v * 255), base=16), 2, '0')
        push!(colours, "#" * hex(r) * hex(g) * hex(b))
    end
    colours[1:n]
end

const _GREY_HEX = "#B3B3B3"

function _apply_grey!(colours::Vector{String}, labels::Vector{String})
    for (i, l) in enumerate(labels)
        l in ("Unclassified", "Other") && (colours[i] = _GREY_HEX)
    end
end

## Alpha diversity chart
"""
    alpha_chart(sample_names, richness, shannon, simpson) -> Dict

Three-panel bar chart of alpha diversity metrics per sample.
"""
function alpha_chart(sample_names::Vector{String},
                     richness_values::Vector{Int},
                     shannon_values::Vector{Float64},
                     simpson_values::Vector{Float64})
    colour = _palette_hex(1)[1]
    metrics = [
        ("y",  "x",  "Richness (observed ASVs)", Float64.(richness_values)),
        ("y2", "x2", "Shannon index",            shannon_values),
        ("y3", "x3", "Simpson index",            simpson_values),
    ]
    traces = [Dict{String,Any}(
        "type" => "bar", "x" => sample_names, "y" => vals,
        "xaxis" => xax, "yaxis" => yax,
        "showlegend" => false, "marker" => Dict("colour" => colour),
    ) for (yax, xax, _, vals) in metrics]

    layout = Dict{String,Any}(
        "title" => Dict("text" => "Alpha diversity"),
        "grid" => Dict("rows" => 3, "columns" => 1, "pattern" => "independent"),
        "yaxis"  => Dict("title" => "Richness (observed ASVs)"),
        "yaxis2" => Dict("title" => "Shannon index"),
        "yaxis3" => Dict("title" => "Simpson index"),
        "xaxis"  => Dict("tickangle" => -45),
        "xaxis2" => Dict("tickangle" => -45),
        "xaxis3" => Dict("title" => "Sample", "tickangle" => -45),
    )
    Dict("data" => traces, "layout" => layout)
end

## Taxa bar chart
"""
    taxa_bar_chart(taxon_labels, sample_names, counts; top_n, relative) -> Dict

Stacked bar chart of taxonomic composition.
`counts` is a taxa-by-samples matrix (from aggregate_by_taxon output).
"""
function taxa_bar_chart(taxon_labels::Vector{String},
                        sample_names::Vector{String},
                        counts::Matrix{Float64};
                        top_n::Int=15, relative::Bool=true)
    n_taxa = length(taxon_labels)
    n_samples = length(sample_names)

    totals = vec(sum(counts, dims=2))
    order = sortperm(totals, rev=true)

    if n_taxa > top_n
        keep = order[1:top_n]
        other = vec(sum(counts[order[(top_n+1):end], :], dims=1))
        counts = vcat(counts[keep, :], other')
        final_labels = vcat(taxon_labels[keep], ["Other"])
    else
        counts = counts[order, :]
        final_labels = taxon_labels[order]
    end

    if relative
        col_sums = vec(sum(counts, dims=1))
        for j in 1:n_samples
            col_sums[j] > 0 && (counts[:, j] ./= col_sums[j])
        end
    end

    n_final = length(final_labels)
    colours = _palette_hex(n_final)
    _apply_grey!(colours, final_labels)
    ylabel = relative ? "Relative abundance" : "Read count"
    title = relative ? "Taxonomic composition" : "Taxonomic composition (absolute)"

    traces = [Dict{String,Any}(
        "type" => "bar", "name" => final_labels[i],
        "x" => sample_names, "y" => collect(counts[i, :]),
        "marker" => Dict("colour" => colours[i]),
    ) for i in 1:n_final]

    layout = Dict{String,Any}(
        "barmode" => "stack",
        "title" => Dict("text" => title),
        "xaxis" => Dict("title" => "Sample", "tickangle" => -45),
        "yaxis" => Dict{String,Any}("title" => ylabel),
        "legend" => Dict("traceorder" => "normal"),
    )
    relative && (layout["yaxis"]["range"] = [0, 1])
    Dict("data" => traces, "layout" => layout)
end

## Pipeline stats chart
const _STAGE_ORDER = ["input", "filtered", "denoisedF", "denoisedR", "merged", "nochim"]
const _STAGE_LABELS = Dict(
    "input" => "Input", "filtered" => "Filtered",
    "denoisedF" => "Denoised (F)", "denoisedR" => "Denoised (R)",
    "merged" => "Merged", "nochim" => "Chimera-free",
)

"""
    pipeline_stats_chart(stats_df) -> Dict

Grouped bar chart of read counts at each DADA2 pipeline stage.
`stats_df` must have a `sample` column plus stage count columns.
"""
function pipeline_stats_chart(stats_df::DataFrame)
    all_cols = names(stats_df)
    stage_cols = [c for c in _STAGE_ORDER if c in all_cols]
    reads_cols = sort([c for c in all_cols if startswith(c, "reads_")])
    append!(stage_cols, reads_cols)
    isempty(stage_cols) && return nothing

    sample_names = String.(stats_df.sample)
    stage_labels = [get(_STAGE_LABELS, c, c) for c in stage_cols]
    colours = _palette_hex(length(sample_names))

    traces = [Dict{String,Any}(
        "type" => "bar", "name" => sample_names[si],
        "x" => stage_labels,
        "y" => [let v = stats_df[si, Symbol(c)]; ismissing(v) ? 0.0 : Float64(v) end
                for c in stage_cols],
        "marker" => Dict("colour" => colours[si]),
    ) for si in eachindex(sample_names)]

    layout = Dict{String,Any}(
        "barmode" => "group",
        "title" => Dict("text" => "Reads through pipeline stages"),
        "xaxis" => Dict("title" => "Pipeline stage", "tickangle" => -45),
        "yaxis" => Dict("title" => "Read count"),
    )
    Dict("data" => traces, "layout" => layout)
end

## Alpha diversity boxplot (cross-run)
"""
    alpha_boxplot(groups) -> Dict

Three-panel boxplot comparing alpha diversity across groups.
`groups` is a vector of `(label, richness, shannon, simpson)` tuples.
"""
function alpha_boxplot(groups::Vector{Tuple{String, Vector{Int}, Vector{Float64}, Vector{Float64}}})
    colours = _palette_hex(length(groups))
    panels = [
        ("y",  "Richness (observed ASVs)", 1),
        ("y2", "Shannon index",            2),
        ("y3", "Simpson index",            3),
    ]
    traces = Any[]
    for (pi, (yax, _, panel_idx)) in enumerate(panels)
        for (gi, (label, richness_values, shannon_values, simpson_values)) in enumerate(groups)
            vals = pi == 1 ? Float64.(richness_values) : pi == 2 ? shannon_values : simpson_values
            push!(traces, Dict{String,Any}(
                "type" => "box", "name" => label,
                "y" => vals, "yaxis" => yax,
                "xaxis" => panel_idx == 1 ? "x" : "x$panel_idx",
                "marker" => Dict("colour" => colours[gi]),
                "showlegend" => pi == 1,
                "legendgroup" => label,
            ))
        end
    end
    layout = Dict{String,Any}(
        "title" => Dict("text" => "Alpha diversity comparison"),
        "grid" => Dict("rows" => 3, "columns" => 1, "pattern" => "independent"),
        "yaxis"  => Dict("title" => "Richness (observed ASVs)"),
        "yaxis2" => Dict("title" => "Shannon index"),
        "yaxis3" => Dict("title" => "Simpson index"),
    )
    Dict("data" => traces, "layout" => layout)
end

## NMDS
"""
    nmds_chart(coords, labels; colour_by, stress) -> Dict

NMDS ordination scatter plot. `coords` is an Nx2 matrix.
"""
const _MARKER_SYMBOLS = [
    "circle", "square", "diamond", "triangle-up", "triangle-down",
    "pentagon", "hexagon", "star", "cross", "x",
]

function nmds_chart(coords::Matrix{Float64}, labels::Vector{String};
                    colour_by::Union{Vector{String},Nothing}=nothing,
                    shape_by::Union{Vector{String},Nothing}=nothing,
                    stress::Union{Float64,Nothing}=nothing)
    n = size(coords, 1)
    cgroups = isnothing(colour_by) ? ["all"] : unique(colour_by)
    sgroups = isnothing(shape_by) ? ["all"] : unique(shape_by)
    colours = _palette_hex(length(cgroups))
    cb = isnothing(colour_by) ? fill("all", n) : colour_by
    sb = isnothing(shape_by) ? fill("all", n) : shape_by

    traces = Any[]
    for (ci, cg) in enumerate(cgroups)
        for (si, sg) in enumerate(sgroups)
            mask = (cb .== cg) .& (sb .== sg)
            any(mask) || continue
            name = length(sgroups) <= 1 ? cg : "$cg ($sg)"
            sym = _MARKER_SYMBOLS[mod1(si, length(_MARKER_SYMBOLS))]
            push!(traces, Dict{String,Any}(
                "type" => "scatter", "mode" => "markers", "name" => name,
                "x" => coords[mask, 1], "y" => coords[mask, 2],
                "text" => labels[mask],
                "legendgroup" => cg,
                "marker" => Dict("colour" => colours[ci], "size" => 10,
                                 "symbol" => sym),
            ))
        end
    end

    annotations = Any[]
    if !isnothing(stress)
        push!(annotations, Dict{String,Any}(
            "text" => "stress = $(round(stress; digits=3))",
            "xref" => "paper", "yref" => "paper",
            "x" => 0.02, "y" => 0.98,
            "xanchor" => "left", "yanchor" => "top",
            "showarrow" => false, "font" => Dict("size" => 12),
        ))
    end

    layout = Dict{String,Any}(
        "title" => Dict("text" => "NMDS ordination"),
        "xaxis" => Dict("title" => "NMDS1"),
        "yaxis" => Dict("title" => "NMDS2"),
        "annotations" => annotations,
    )
    Dict("data" => traces, "layout" => layout)
end

## NMDS and PERMANOVA
const _r_lock = ReentrantLock()
const _r_loaded = Ref(false)

function _ensure_r()
    _r_loaded[] && return true
    lock(_r_lock) do
        _r_loaded[] && return true
        try
            @eval using RCall
            getfield(@__MODULE__, :RCall).reval("suppressPackageStartupMessages(library(vegan))")
            _r_loaded[] = true
            return true
        catch e
            @warn "R/vegan not available - NMDS and PERMANOVA disabled" exception=e
            return false
        end
    end
end

r_available() = _ensure_r()

"""
    run_nmds(mat) -> (coords::Matrix{Float64}, stress::Float64)

NMDS via vegan::metaMDS with Bray-Curtis distance.
`mat` is samples-by-features. Returns NaN-filled results on failure.
"""
function run_nmds(mat::Matrix{Float64})
    _ensure_r() || return (fill(NaN, size(mat, 1), 2), NaN)
    lock(_r_lock) do
        rcall = getfield(@__MODULE__, :RCall)
        rcall.globalEnv[:mat] = mat
        rcall.reval("""
            set.seed(42)
            nmds_res <- tryCatch(
                metaMDS(mat, distance = "bray", k = 2, trymax = 200,
                        autotransform = FALSE, trace = 0),
                error = function(e) NULL
            )
            if (!is.null(nmds_res)) {
                nmds_coords <- nmds_res\$points
                nmds_stress <- nmds_res\$stress
            } else {
                nmds_coords <- matrix(NA_real_, nrow = nrow(mat), ncol = 2)
                nmds_stress <- NA_real_
            }
        """)
        coords = rcall.rcopy(rcall.reval("nmds_coords"))::Matrix{Float64}
        stress = rcall.rcopy(rcall.reval("nmds_stress"))::Float64
        rcall.reval("rm(mat, nmds_res, nmds_coords, nmds_stress); gc()")
        coords, stress
    end
end

"""
    run_permanova(mat, metadata) -> Union{NamedTuple, Nothing}

PERMANOVA via vegan::adonis2.
`metadata` is a DataFrame with one row per sample and covariate columns.
"""
function run_permanova(mat::Matrix{Float64}, metadata::DataFrame)
    _ensure_r() || return nothing
    covariates = [c for c in names(metadata) if lowercase(c) != "sample"]
    isempty(covariates) && return nothing
    formula_rhs = join(covariates, " + ")

    lock(_r_lock) do
        rcall = getfield(@__MODULE__, :RCall)
        meta_r = copy(metadata)
        rcall.globalEnv[:mat] = mat
        rcall.globalEnv[:meta_r] = meta_r
        rcall.globalEnv[:formula_rhs] = formula_rhs
        rcall.reval("""
            set.seed(42)
            dist_mat <- vegdist(mat, method = "bray")
            form <- as.formula(paste("dist_mat ~", formula_rhs))
            perm_res <- tryCatch(
                adonis2(form, data = meta_r, permutations = 999),
                error = function(e) NULL
            )
            if (!is.null(perm_res)) {
                perm_text <- paste(capture.output(print(perm_res)), collapse = "\\n")
                perm_r2 <- perm_res\$R2[1]
                perm_f <- perm_res\$F[1]
                perm_p <- perm_res[["Pr(>F)"]][1]
            } else {
                perm_text <- NA_character_
                perm_r2 <- NA_real_
                perm_f <- NA_real_
                perm_p <- NA_real_
            }
        """)
        txt = rcall.rcopy(rcall.reval("perm_text"))
        r2 = rcall.rcopy(rcall.reval("perm_r2"))
        f_stat = rcall.rcopy(rcall.reval("perm_f"))
        p_val = rcall.rcopy(rcall.reval("perm_p"))
        rcall.reval("rm(mat, meta_r, formula_rhs, dist_mat, form, perm_res, perm_text, perm_r2, perm_f, perm_p); gc()")

        (ismissing(txt) || txt == "NA") && return nothing
        (; text=txt,
           r2=ismissing(r2) ? nothing : r2,
           f_statistic=ismissing(f_stat) ? nothing : f_stat,
           p_value=ismissing(p_val) ? nothing : p_val)
    end
end

end
