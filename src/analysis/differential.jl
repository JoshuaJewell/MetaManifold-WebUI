# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

## Differential abundance
# Per-taxon negative-binomial GLMs comparing two groups of samples, fitted with
# MASS::glm.nb, with a library-size offset and Benjamini-Hochberg correction
# across the taxa that produced a p-value.
#
# Count matrices here are samples x taxa, as everywhere else in `Analysis`.
# The model is fitted on the observed integer read counts: there are no
# pseudocounts, and anything that is not a non-negative integer is refused
# rather than rounded. A taxon whose fit fails is reported with its reason and
# kept out of the multiple-testing family; it never receives a stand-in p-value.
module Differential

using DataFrames, RCall
using ..RRuntime: with_r_lock
using ..Analysis: _palette_hex, R_WAIT_SECONDS

export DifferentialConfig, ScalingRefusal, MASSUnavailable, OFFSET_METHODS,
       tss_factors, rle_factors, size_factors, bh_adjust, validate_counts,
       fit_nb, differential_abundance, volcano_chart

# Offset (size factor) methods accepted under `analysis.differential.offset`.
const OFFSET_METHODS = ("tss", "rle")

# Above this the NB dispersion parameter has run to its upper bound: the counts
# show no overdispersion and the model has collapsed to a Poisson. Below the
# lower bound the dispersion is degenerate the other way.
const THETA_UPPER = 1e7
const THETA_LOWER = 1e-8

"""
    DifferentialConfig(offset::AbstractString="tss", min_prevalence::Real=0.0)

Settings for one differential abundance analysis, read from
`analysis.differential` in the pipeline config. `offset` names the size factor
method (`"tss"` or `"rle"`). `min_prevalence` is the fraction of samples, in
[0, 1], in which a taxon must have at least one read to be tested.
"""
struct DifferentialConfig
    offset         :: String
    min_prevalence :: Float64
    function DifferentialConfig(offset::AbstractString="tss", min_prevalence::Real=0.0)
        o = lowercase(strip(String(offset)))
        o in OFFSET_METHODS || throw(ArgumentError(
            "analysis.differential.offset must be one of $(join(OFFSET_METHODS, ", ")), not '$offset'"))
        (isfinite(min_prevalence) && 0 <= min_prevalence <= 1) || throw(ArgumentError(
            "analysis.differential.min_prevalence must be a fraction in [0, 1], not $min_prevalence"))
        new(o, Float64(min_prevalence))
    end
end

"""
    ScalingRefusal(method, reason)

Raised when a size factor method cannot be computed on the data it was given,
for instance total-sum scaling on an empty sample. The reason names the sample
or the property at fault.
"""
struct ScalingRefusal <: Exception
    method :: String
    reason :: String
end

Base.showerror(io::IO, e::ScalingRefusal) =
    print(io, "$(uppercase(e.method)) scaling refused: $(e.reason)")

"""
    MASSUnavailable()

Raised when the R package MASS, which provides `glm.nb`, cannot be loaded.
"""
struct MASSUnavailable <: Exception end

Base.showerror(io::IO, ::MASSUnavailable) = print(io,
    "the R package MASS is not available, so negative-binomial models cannot be " *
    "fitted; install R with its recommended packages (MASS ships with them)")

"""
    _geomean(x) -> Float64

Geometric mean of a vector of positive numbers.
"""
_geomean(x) = exp(sum(log.(x)) / length(x))

"""
    _median(x) -> Float64

Median of a non-empty vector, the mean of the two middle values when its length
is even.
"""
function _median(x)
    s = sort(collect(Float64, x))
    m = length(s) ÷ 2
    isodd(length(s)) ? s[m + 1] : (s[m] + s[m + 1]) / 2
end

"""
    tss_factors(counts, samples) -> Vector{Float64}

Total-sum scaling factors, one per row of the samples x taxa matrix `counts`:
each sample's library size divided by the geometric mean of all library sizes,
so the factors have geometric mean 1. Throws `ScalingRefusal` naming the first
sample with no reads.
"""
function tss_factors(counts::AbstractMatrix{<:Real}, samples::AbstractVector{<:AbstractString})
    lib = vec(sum(counts; dims=2))
    for (i, l) in enumerate(lib)
        l > 0 || throw(ScalingRefusal("tss",
            "sample '$(samples[i])' has no reads, so it has no library size to scale by"))
    end
    lib ./ _geomean(lib)
end

"""
    rle_factors(counts, samples) -> Vector{Float64}

Relative log expression (median-of-ratios) scaling factors, one per row of the
samples x taxa matrix `counts`. Each sample's factor is the median, over the
taxa with at least one read in every sample, of its count divided by that
taxon's geometric mean across samples; the factors are then centred to
geometric mean 1. Throws `ScalingRefusal` when no taxon is present in every
sample, since every geometric mean would then be zero.
"""
function rle_factors(counts::AbstractMatrix{<:Real}, samples::AbstractVector{<:AbstractString})
    size(counts, 1) == length(samples) ||
        throw(ArgumentError("$(length(samples)) sample names for $(size(counts, 1)) rows"))
    shared = [j for j in axes(counts, 2) if all(>(0), view(counts, :, j))]
    isempty(shared) && throw(ScalingRefusal("rle",
        "no taxon has reads in every sample, so the per-taxon geometric mean the " *
        "median-of-ratios needs is zero for all of them; use offset: tss"))
    sub = Float64.(counts[:, shared])
    gm = [_geomean(view(sub, :, j)) for j in axes(sub, 2)]
    raw = [_median(view(sub, i, :) ./ gm) for i in axes(sub, 1)]
    raw ./ _geomean(raw)
end

"""
    size_factors(counts, samples, method) -> Vector{Float64}

Size factors by the named `method` (`"tss"` or `"rle"`).
"""
function size_factors(counts::AbstractMatrix{<:Real}, samples::AbstractVector{<:AbstractString},
                      method::AbstractString)
    method == "tss" && return tss_factors(counts, samples)
    method == "rle" && return rle_factors(counts, samples)
    throw(ArgumentError("unknown offset method '$method'; expected one of $(join(OFFSET_METHODS, ", "))"))
end

"""
    bh_adjust(p) -> Vector{Float64}

Benjamini-Hochberg adjusted p-values, matching R's `p.adjust(p, "BH")`. Every
input must be a finite probability in [0, 1]; anything else throws an
`ArgumentError` rather than being dropped, because dropping it would silently
shrink the family. An empty input gives an empty result.
"""
function bh_adjust(p::AbstractVector{<:Real})
    n = length(p)
    n == 0 && return Float64[]
    for (i, x) in enumerate(p)
        (isfinite(x) && 0 <= x <= 1) || throw(ArgumentError(
            "p-value $i is $x: Benjamini-Hochberg needs finite probabilities in [0, 1]"))
    end
    order = sortperm(p; rev=true)
    out = Vector{Float64}(undef, n)
    running = 1.0
    for (k, idx) in enumerate(order)
        rank = n - k + 1
        running = min(running, n / rank * Float64(p[idx]))
        out[idx] = min(running, 1.0)
    end
    out
end

"""
    validate_counts(counts, samples, taxa) -> Nothing

Check that every entry of the samples x taxa matrix `counts` is a finite,
non-negative integer. Throws an `ArgumentError` naming the first offending
sample and taxon: the model is defined on read counts, and a normalised,
rarefied-and-averaged or pseudocounted value would be silently misfitted.
"""
function validate_counts(counts::AbstractMatrix{<:Real}, samples::AbstractVector{<:AbstractString},
                         taxa::AbstractVector{<:AbstractString})
    size(counts) == (length(samples), length(taxa)) || throw(ArgumentError(
        "count matrix is $(size(counts)) but there are $(length(samples)) samples and $(length(taxa)) taxa"))
    for j in axes(counts, 2), i in axes(counts, 1)
        x = counts[i, j]
        (isfinite(x) && x >= 0 && isinteger(x)) || throw(ArgumentError(
            "the count for taxon '$(taxa[j])' in sample '$(samples[i])' is $x: the " *
            "negative-binomial model needs non-negative integer read counts"))
    end
    nothing
end

# The per-taxon fit loop. Variables are prefixed `da_` so they cannot collide
# with another analysis's globals, and are removed afterwards.
const _FIT_R = raw"""
da_g <- factor(da_group, levels = da_levels)
da_term <- paste0("da_g", da_levels[2])
da_n <- ncol(da_counts)
da_result <- data.frame(status = rep("ok", da_n), note = rep("", da_n),
                        estimate = rep(NA_real_, da_n), se = rep(NA_real_, da_n),
                        statistic = rep(NA_real_, da_n), pvalue = rep(NA_real_, da_n),
                        theta = rep(NA_real_, da_n), stringsAsFactors = FALSE)
for (da_j in seq_len(da_n)) {
  da_y <- da_counts[, da_j]
  if (length(unique(da_y)) < 2L) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- "the counts are constant across samples: there is nothing to estimate"
    next
  }
  # Complete separation: with no reads in every sample of one group the group
  # coefficient runs off to infinity, its standard error with it, and the Wald
  # p-value collapses towards 1 although the fit reports convergence. That p is
  # not evidence of no difference, so the taxon is refused with its reason.
  da_absent <- levels(da_g)[vapply(levels(da_g), function(l) all(da_y[da_g == l] == 0), logical(1))]
  if (length(da_absent) > 0L) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- paste0("no reads in any '", da_absent[[1]], "' sample (complete ",
                                      "separation): the fold change is unbounded and the Wald test ",
                                      "is undefined, so no p-value is reported")
    next
  }
  da_warns <- character(0)
  da_fit <- tryCatch(
    withCallingHandlers(
      MASS::glm.nb(da_y ~ da_g + offset(da_offset), control = glm.control(maxit = 100)),
      warning = function(w) {
        da_warns <<- c(da_warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  if (inherits(da_fit, "error")) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- paste("glm.nb stopped with an error:", conditionMessage(da_fit))
    next
  }
  da_warns <- unique(da_warns)
  if (!isTRUE(da_fit[["converged"]]) || !is.null(da_fit[["th.warn"]]) ||
      any(grepl("iteration limit|alternation limit|did not converge|NaNs produced", da_warns))) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- paste("the fit did not converge:",
                                     paste(c(da_fit[["th.warn"]], da_warns), collapse = "; "))
    next
  }
  da_co <- summary(da_fit)[["coefficients"]]
  if (!(da_term %in% rownames(da_co))) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- "the group coefficient is not estimable (aliased)"
    next
  }
  da_row <- da_co[da_term, ]
  if (!all(is.finite(da_row)) || da_row[[4]] < 0 || da_row[[4]] > 1) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- "the fit returned a non-finite estimate, standard error or p-value"
    next
  }
  da_theta <- da_fit[["theta"]]
  da_result[da_j, c("estimate", "se", "statistic", "pvalue")] <- unname(da_row[1:4])
  da_result[da_j, "theta"] <- da_theta
  if (!is.finite(da_theta) || da_theta >= da_theta_upper) {
    da_result[da_j, "status"] <- "boundary"
    da_result[da_j, "note"] <- "the dispersion parameter theta reached its upper bound: these counts show no overdispersion, so the fit is effectively Poisson"
  } else if (da_theta <= da_theta_lower) {
    da_result[da_j, "status"] <- "boundary"
    da_result[da_j, "note"] <- "the dispersion parameter theta reached its lower bound: the variance is extreme relative to the mean"
  } else if (length(da_warns) > 0L) {
    da_result[da_j, "note"] <- paste("R warned:", paste(da_warns, collapse = "; "))
  }
}
"""

"""
    _num(x) -> Union{Float64,Nothing}

A finite number as `Float64`, or `nothing` for R's NA, NaN or an infinity.
"""
_num(x) = (ismissing(x) || isnothing(x) || !isfinite(x)) ? nothing : Float64(x)

"""
    fit_nb(counts, groups, offset; levels) -> DataFrame

Fit `MASS::glm.nb(y ~ group + offset(offset))` to each column of the samples x
taxa integer matrix `counts`. `levels` is `(reference, contrast)`, so the
estimate is the natural-log fold change of `contrast` over `reference`.
Returns one row per taxon with `status` (`ok`, `boundary` or `failed`), `note`,
`estimate`, `se`, `statistic`, `pvalue` and `theta`; a failed taxon has
`nothing` for every statistic. Throws `MASSUnavailable` when MASS cannot be
loaded and `RBusyError` when the R runtime stays busy.
"""
function fit_nb(counts::AbstractMatrix{<:Integer}, groups::AbstractVector{<:AbstractString},
                offset::AbstractVector{<:Real}; levels::NTuple{2,String})
    size(counts, 1) == length(groups) == length(offset) || throw(ArgumentError(
        "counts have $(size(counts, 1)) samples, groups $(length(groups)), offset $(length(offset))"))
    raw = with_r_lock(; timeout=R_WAIT_SECONDS[]) do
        RCall.rcopy(RCall.reval("requireNamespace('MASS', quietly = TRUE)")) || throw(MASSUnavailable())
        RCall.globalEnv[:da_counts] = Matrix{Int}(counts)
        RCall.globalEnv[:da_group] = String.(groups)
        RCall.globalEnv[:da_levels] = collect(levels)
        RCall.globalEnv[:da_offset] = Float64.(offset)
        RCall.globalEnv[:da_theta_upper] = THETA_UPPER
        RCall.globalEnv[:da_theta_lower] = THETA_LOWER
        try
            RCall.reval(_FIT_R)
            DataFrame(RCall.rcopy(RCall.reval("da_result")))
        finally
            RCall.reval("rm(list = intersect(ls(), c('da_counts', 'da_group', 'da_levels', " *
                        "'da_offset', 'da_theta_upper', 'da_theta_lower', 'da_g', 'da_term', " *
                        "'da_n', 'da_result', 'da_j', 'da_y', 'da_warns', 'da_fit', 'da_co', " *
                        "'da_row', 'da_theta', 'da_absent'))); invisible(gc())")
        end
    end
    DataFrame(status    = String.(raw.status),
              note      = String.(raw.note),
              estimate  = _num.(raw.estimate),
              se        = _num.(raw.se),
              statistic = _num.(raw.statistic),
              pvalue    = _num.(raw.pvalue),
              theta     = _num.(raw.theta))
end

"""
    differential_abundance(counts, samples, taxa, groups; reference, contrast,
                           config=DifferentialConfig()) -> Dict{String,Any}

Test every taxon (column) of the samples x taxa read-count matrix for a
difference in abundance between the samples labelled `contrast` and those
labelled `reference` in `groups`.

Size factors are computed on the full matrix by `config.offset` and enter each
model as `log(factor)`. Taxa present in fewer than `config.min_prevalence` of
the samples are reported as `filtered` and not fitted. The Benjamini-Hochberg
family is the fitted taxa that produced a p-value; failed and filtered taxa
carry `padj = nothing`.

Throws `ArgumentError` for unusable input, `ScalingRefusal` when the offset
cannot be computed, `MASSUnavailable` when MASS is missing, and an
`ErrorException` when no taxon could be fitted at all.
"""
function differential_abundance(counts::AbstractMatrix{<:Real},
                                samples::AbstractVector{<:AbstractString},
                                taxa::AbstractVector{<:AbstractString},
                                groups::AbstractVector{<:AbstractString};
                                reference::AbstractString, contrast::AbstractString,
                                config::DifferentialConfig=DifferentialConfig())
    reference == contrast && throw(ArgumentError(
        "the two groups must differ, but both are '$reference'"))
    length(groups) == length(samples) || throw(ArgumentError(
        "$(length(groups)) group labels for $(length(samples)) samples"))
    isempty(taxa) && throw(ArgumentError("there are no taxa to test"))
    validate_counts(counts, samples, taxa)
    unknown = setdiff(unique(groups), (reference, contrast))
    isempty(unknown) || throw(ArgumentError(
        "samples belong to groups other than '$reference' and '$contrast': $(join(unknown, ", "))"))
    n_ref = count(==(reference), groups)
    n_con = count(==(contrast), groups)
    (n_ref >= 1 && n_con >= 1) || throw(ArgumentError(
        "each group needs at least one sample; '$reference' has $n_ref and '$contrast' has $n_con"))
    n_ref + n_con >= 3 || throw(ArgumentError(
        "a group effect and a dispersion cannot be estimated from $(n_ref + n_con) samples; at least 3 are needed"))

    factors = size_factors(counts, samples, config.offset)
    offset = log.(factors)

    n = length(samples)
    prevalence = [count(>(0), view(counts, :, j)) / n for j in axes(counts, 2)]
    tested = findall(>=(config.min_prevalence), prevalence)

    fits = isempty(tested) ? nothing :
        fit_nb(Int.(counts[:, tested]), groups, offset; levels=(String(reference), String(contrast)))

    rows = [Dict{String,Any}(
                "taxon" => String(taxa[j]), "status" => "filtered",
                "note" => "present in $(count(>(0), view(counts, :, j))) of $n samples, " *
                          "below analysis.differential.min_prevalence = $(config.min_prevalence)",
                "estimate" => nothing, "log2_fold_change" => nothing,
                "standard_error" => nothing, "statistic" => nothing,
                "pvalue" => nothing, "padj" => nothing, "dispersion_theta" => nothing,
                "prevalence" => prevalence[j])
            for j in axes(counts, 2)]
    if !isnothing(fits)
        for (k, j) in enumerate(tested)
            f = fits[k, :]
            est = f.estimate
            merge!(rows[j], Dict{String,Any}(
                "status" => f.status, "note" => f.note,
                "estimate" => est,
                "log2_fold_change" => isnothing(est) ? nothing : est / log(2),
                "standard_error" => f.se, "statistic" => f.statistic,
                "pvalue" => f.pvalue, "dispersion_theta" => f.theta))
        end
    end

    family = findall(r -> r["status"] in ("ok", "boundary") && !isnothing(r["pvalue"]), rows)
    if isempty(family)
        notes = unique(r["note"] for r in rows)
        error("no taxon could be fitted: " * join(first(notes, 5), " | "))
    end
    padj = bh_adjust([rows[j]["pvalue"] for j in family])
    for (k, j) in enumerate(family)
        rows[j]["padj"] = padj[k]
    end

    sort!(rows; by = r -> (isnothing(r["padj"]) ? 2.0 : r["padj"],
                           isnothing(r["pvalue"]) ? 2.0 : r["pvalue"], r["taxon"]))

    n_failed = count(r -> r["status"] == "failed", rows)
    n_boundary = count(r -> r["status"] == "boundary", rows)
    n_filtered = count(r -> r["status"] == "filtered", rows)
    Dict{String,Any}(
        "status" => n_failed == 0 ? "ok" : "partial",
        "method" => "Negative-binomial GLM per taxon (MASS::glm.nb, log link) with a " *
                    "$(uppercase(config.offset)) size-factor offset; Wald test of the group " *
                    "coefficient; Benjamini-Hochberg adjustment over the fitted taxa",
        "groups" => Dict("reference" => String(reference), "contrast" => String(contrast)),
        "n_samples" => Dict(String(reference) => n_ref, String(contrast) => n_con),
        "size_factors" => [Dict("sample" => String(samples[i]), "group" => String(groups[i]),
                                "factor" => factors[i]) for i in eachindex(samples)],
        "config" => Dict("offset" => config.offset, "min_prevalence" => config.min_prevalence),
        "diagnostics" => Dict("n_taxa" => length(taxa), "n_tested" => length(family),
                              "n_failed" => n_failed, "n_boundary" => n_boundary,
                              "n_filtered" => n_filtered),
        "rows" => rows,
    )
end

"""
    volcano_chart(result; alpha=0.05) -> Dict

A Plotly volcano plot of a `differential_abundance` result: log2 fold change
against -log10 of the raw p-value, coloured by whether the BH-adjusted p-value
is below `alpha`. Taxa without a p-value are omitted here; they remain in the
results table with their reason.
"""
function volcano_chart(result::AbstractDict; alpha::Real=0.05)
    ref = result["groups"]["reference"]
    con = result["groups"]["contrast"]
    fitted = filter(r -> !isnothing(r["padj"]), result["rows"])
    colours = _palette_hex(2)
    traces = Any[]
    for (sig, name, colour) in ((false, "padj ≥ $alpha", colours[2]),
                                (true,  "padj < $alpha", colours[1]))
        sel = filter(r -> (r["padj"] < alpha) == sig, fitted)
        isempty(sel) && continue
        # A p-value that underflowed to 0 is drawn at the smallest positive
        # double and says so in its hover text.
        y = [-log10(max(r["pvalue"], floatmin(Float64))) for r in sel]
        text = [string(r["taxon"], "<br>padj = ", round(r["padj"]; sigdigits=3),
                       r["pvalue"] == 0 ? "<br>p underflowed to 0" : "") for r in sel]
        push!(traces, Dict{String,Any}(
            "type" => "scatter", "mode" => "markers", "name" => name,
            "x" => [r["log2_fold_change"] for r in sel], "y" => y,
            "text" => text, "hoverinfo" => "text+x+y",
            "marker" => Dict("color" => colour, "size" => 8),
        ))
    end
    layout = Dict{String,Any}(
        "title" => Dict("text" => "Differential abundance: $con vs $ref"),
        "xaxis" => Dict("title" => "log2 fold change ($con / $ref)", "zeroline" => true),
        "yaxis" => Dict("title" => "-log10 p"),
    )
    Dict("data" => traces, "layout" => layout)
end

end # module Differential
