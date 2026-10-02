# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

## Differential abundance
# Per-taxon models comparing two groups of samples, with Benjamini-Hochberg
# correction across the taxa that produced a p-value. Two methods:
#
# - nb_glm: a negative-binomial GLM per taxon (MASS::glm.nb) on the observed
#   integer read counts, with a library-size offset. There are no pseudocounts.
# - clr_lm: a Gaussian linear model per taxon (stats::lm) on centred log-ratios,
#   after multiplicative replacement of the zeros (`ZeroReplacement`). The
#   estimate is a difference in CLR units, not a fold change.
#
# Count matrices here are samples x taxa, as everywhere else in `Analysis`.
# Anything that is not a non-negative integer count is refused rather than
# rounded. A taxon whose fit fails is reported with its reason and kept out of
# the multiple-testing family; it never receives a stand-in p-value.
module Differential

using DataFrames, RCall
using ..RRuntime: with_r_lock
using ..Analysis: _palette_hex, R_WAIT_SECONDS
using ..ZeroReplacement: DEFAULT_DELTA, multiplicative_replacement

export DifferentialConfig, ScalingRefusal, MASSUnavailable, OFFSET_METHODS, METHODS,
       tss_factors, rle_factors, size_factors, bh_adjust, validate_counts,
       clr_transform, fit_nb, fit_lm, differential_abundance, volcano_chart

# Methods accepted under `analysis.differential.method`.
const METHODS = ("nb_glm", "clr_lm")

# Offset (size factor) methods accepted under `analysis.differential.offset`.
const OFFSET_METHODS = ("tss", "rle")

# Above this the NB dispersion parameter has run to its upper bound: the counts
# show no overdispersion and the model has collapsed to a Poisson. Below the
# lower bound the dispersion is degenerate the other way.
const THETA_UPPER = 1e7
const THETA_LOWER = 1e-8

"""
    DifferentialConfig(offset="tss", min_prevalence=0.0; method="nb_glm",
                       replacement_delta=0.65)

Settings for one differential abundance analysis, read from
`analysis.differential` in the pipeline config. `method` is `"nb_glm"` or
`"clr_lm"`. `offset` names the size factor method (`"tss"` or `"rle"`) and is
used by `nb_glm` only: a library-size offset has no meaning for a Gaussian
model of log-ratios, which are already free of library size.
`min_prevalence` is the fraction of samples, in [0, 1], in which a taxon must
have at least one read to be tested. `replacement_delta`, in (0, 1), is the
fraction of a taxon's detection limit given to a replaced zero by `clr_lm`.
"""
struct DifferentialConfig
    method            :: String
    offset            :: String
    min_prevalence    :: Float64
    replacement_delta :: Float64
    function DifferentialConfig(offset::AbstractString="tss", min_prevalence::Real=0.0;
                                method::AbstractString="nb_glm",
                                replacement_delta::Real=DEFAULT_DELTA)
        m = lowercase(strip(String(method)))
        m in METHODS || throw(ArgumentError(
            "analysis.differential.method must be one of $(join(METHODS, ", ")), not '$method'"))
        o = lowercase(strip(String(offset)))
        o in OFFSET_METHODS || throw(ArgumentError(
            "analysis.differential.offset must be one of $(join(OFFSET_METHODS, ", ")), not '$offset'"))
        (isfinite(min_prevalence) && 0 <= min_prevalence <= 1) || throw(ArgumentError(
            "analysis.differential.min_prevalence must be a fraction in [0, 1], not $min_prevalence"))
        (isfinite(replacement_delta) && 0 < replacement_delta < 1) || throw(ArgumentError(
            "analysis.differential.replacement_delta must be in (0, 1), not $replacement_delta"))
        new(m, o, Float64(min_prevalence), Float64(replacement_delta))
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
sample and taxon: the models are defined on read counts, and a normalised,
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
            "models are defined on non-negative integer read counts"))
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
                        "'da_row', 'da_theta'))); invisible(gc())")
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
    clr_transform(x) -> Matrix{Float64}

Centred log-ratio of each row (sample) of the samples x taxa matrix `x`:
`log(x_ij)` minus the mean of `log(x_i·)`, so every row sums to zero. Every
value must be positive; replace zeros first.
"""
function clr_transform(x::AbstractMatrix{<:Real})
    all(v -> isfinite(v) && v > 0, x) || throw(ArgumentError(
        "the centred log-ratio needs positive finite values; replace zeros first"))
    l = log.(Float64.(x))
    l .- sum(l; dims=2) ./ size(l, 2)
end

# The per-taxon Gaussian fit loop, with the same `da_` naming and statuses as
# `_FIT_R` but no `boundary`: an lm has no dispersion parameter to run to a bound.
const _FIT_LM_R = raw"""
da_g <- factor(da_group, levels = da_levels)
da_term <- paste0("da_g", da_levels[2])
da_n <- ncol(da_values)
da_result <- data.frame(status = rep("ok", da_n), note = rep("", da_n),
                        estimate = rep(NA_real_, da_n), se = rep(NA_real_, da_n),
                        statistic = rep(NA_real_, da_n), pvalue = rep(NA_real_, da_n),
                        stringsAsFactors = FALSE)
for (da_j in seq_len(da_n)) {
  da_y <- da_values[, da_j]
  if (length(unique(da_y)) < 2L) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- "the values are constant across samples: there is nothing to estimate"
    next
  }
  da_warns <- character(0)
  da_co <- tryCatch(
    withCallingHandlers(
      summary(stats::lm(da_y ~ da_g))[["coefficients"]],
      warning = function(w) {
        da_warns <<- c(da_warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  if (inherits(da_co, "error")) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- paste("lm stopped with an error:", conditionMessage(da_co))
    next
  }
  if (!(da_term %in% rownames(da_co))) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- "the group coefficient is not estimable (aliased)"
    next
  }
  da_row <- da_co[da_term, ]
  if (!all(is.finite(da_row)) || da_row[[2]] <= 0 || da_row[[4]] < 0 || da_row[[4]] > 1) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- paste(c("the fit returned a non-finite or zero standard error, so there is no test (the residual variance is zero)",
                                       unique(da_warns)), collapse = "; ")
    next
  }
  da_result[da_j, c("estimate", "se", "statistic", "pvalue")] <- unname(da_row[1:4])
  if (length(da_warns) > 0L) {
    da_result[da_j, "note"] <- paste("R warned:", paste(unique(da_warns), collapse = "; "))
  }
}
"""

"""
    fit_lm(values, groups; levels) -> DataFrame

Fit `stats::lm(y ~ group)` to each column of the samples x taxa matrix
`values` (for `clr_lm`, centred log-ratios). `levels` is
`(reference, contrast)`, so the estimate is the mean of `contrast` minus the
mean of `reference`, and the test is the two-sample t-test with pooled
variance. Returns one row per taxon with `status` (`ok` or `failed`), `note`,
`estimate`, `se`, `statistic` (t) and `pvalue`; a failed taxon has `nothing`
for every statistic. Throws `RBusyError` when the R runtime stays busy.
"""
function fit_lm(values::AbstractMatrix{<:Real}, groups::AbstractVector{<:AbstractString};
                levels::NTuple{2,String})
    size(values, 1) == length(groups) || throw(ArgumentError(
        "values have $(size(values, 1)) samples and groups $(length(groups))"))
    raw = with_r_lock(; timeout=R_WAIT_SECONDS[]) do
        RCall.globalEnv[:da_values] = Matrix{Float64}(values)
        RCall.globalEnv[:da_group] = String.(groups)
        RCall.globalEnv[:da_levels] = collect(levels)
        try
            RCall.reval(_FIT_LM_R)
            DataFrame(RCall.rcopy(RCall.reval("da_result")))
        finally
            RCall.reval("rm(list = intersect(ls(), c('da_values', 'da_group', 'da_levels', " *
                        "'da_g', 'da_term', 'da_n', 'da_result', 'da_j', 'da_y', 'da_warns', " *
                        "'da_co', 'da_row'))); invisible(gc())")
        end
    end
    DataFrame(status    = String.(raw.status),
              note      = String.(raw.note),
              estimate  = _num.(raw.estimate),
              se        = _num.(raw.se),
              statistic = _num.(raw.statistic),
              pvalue    = _num.(raw.pvalue))
end

"""
    _row(taxon, prevalence, note) -> Dict{String,Any}

A results row with no statistics, as for a taxon that was not fitted.
"""
_row(taxon, prevalence, note) = Dict{String,Any}(
    "taxon" => String(taxon), "status" => "filtered", "note" => note,
    "estimate" => nothing, "log2_fold_change" => nothing,
    "standard_error" => nothing, "statistic" => nothing,
    "pvalue" => nothing, "padj" => nothing, "dispersion_theta" => nothing,
    "prevalence" => prevalence)

"""
    differential_abundance(counts, samples, taxa, groups; reference, contrast,
                           config=DifferentialConfig()) -> Dict{String,Any}

Test every taxon (column) of the samples x taxa read-count matrix for a
difference in abundance between the samples labelled `contrast` and those
labelled `reference` in `groups`, by `config.method`.

`nb_glm`: size factors are computed on the full matrix by `config.offset` and
enter each negative-binomial model as `log(factor)`; the estimate is a
natural-log fold change.

`clr_lm`, in this order: (1) taxa with no reads in any of the samples are set
aside, because a taxon never observed has no detection limit and is not part of
the composition; (2) the zeros of the remaining taxa are replaced
multiplicatively with `config.replacement_delta`; (3) each sample is
transformed to centred log-ratios over all of those taxa; (4) only then is
`min_prevalence` applied, choosing which taxa are tested. Filtering before the
transform would change every sample's geometric mean and so every other
taxon's value. The estimate is the difference in mean CLR, contrast minus
reference: a difference in centred-log-ratio units, not a fold change, and it
is relative to the geometric mean of the observed taxa. Replaced zeros are not
measurements, so a taxon with many zeros carries an estimate that depends on
`replacement_delta`.

For both methods, taxa present in fewer than `config.min_prevalence` of the
samples are reported as `filtered` and not fitted, and the Benjamini-Hochberg
family is the fitted taxa that produced a p-value; failed and filtered taxa
carry `padj = nothing`. The result's `effect` names the row key that holds the
effect size and its label.

Throws `ArgumentError` for unusable input (including a zero replacement that
cannot be made), `ScalingRefusal` when the offset cannot be computed,
`MASSUnavailable` when MASS is missing, and an `ErrorException` when no taxon
could be fitted at all.
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
        "a group effect and a variance cannot be estimated from $(n_ref + n_con) samples; at least 3 are needed"))

    n = length(samples)
    present = [count(>(0), view(counts, :, j)) for j in axes(counts, 2)]
    prevalence = present ./ n
    rows = [_row(taxa[j], prevalence[j],
                 "present in $(present[j]) of $n samples, below " *
                 "analysis.differential.min_prevalence = $(config.min_prevalence)")
            for j in axes(counts, 2)]
    levels = (String(reference), String(contrast))
    diagnostics = Dict{String,Any}()

    if config.method == "nb_glm"
        factors = size_factors(counts, samples, config.offset)
        tested = findall(>=(config.min_prevalence), prevalence)
        fits = isempty(tested) ? nothing :
            fit_nb(Int.(counts[:, tested]), groups, log.(factors); levels)
        if !isnothing(fits)
            for (k, j) in enumerate(tested)
                f = fits[k, :]
                merge!(rows[j], Dict{String,Any}(
                    "status" => f.status, "note" => f.note, "estimate" => f.estimate,
                    "log2_fold_change" => isnothing(f.estimate) ? nothing : f.estimate / log(2),
                    "standard_error" => f.se, "statistic" => f.statistic,
                    "pvalue" => f.pvalue, "dispersion_theta" => f.theta))
            end
        end
        method_text = "Negative-binomial GLM per taxon (MASS::glm.nb, log link) with a " *
                      "$(uppercase(config.offset)) size-factor offset; Wald test of the group " *
                      "coefficient; Benjamini-Hochberg adjustment over the fitted taxa"
        effect = Dict("key" => "log2_fold_change", "label" => "log2 fold change")
        size_factor_rows = [Dict("sample" => String(samples[i]), "group" => String(groups[i]),
                                 "factor" => factors[i]) for i in eachindex(samples)]
        config_echo = Dict{String,Any}("method" => config.method, "offset" => config.offset,
                                       "min_prevalence" => config.min_prevalence)
    else
        observed = findall(>(0), present)
        for j in setdiff(axes(counts, 2), observed)
            rows[j]["note"] = "no reads in any sample of either group, so it has no " *
                              "detection limit and is not part of the composition"
        end
        length(observed) >= 2 || throw(ArgumentError(
            "the centred log-ratio needs at least 2 taxa with reads; there are $(length(observed))"))
        replaced = multiplicative_replacement(counts[:, observed];
                                              delta=config.replacement_delta,
                                              samples, taxa=taxa[observed])
        clr = clr_transform(replaced.counts)
        tested = findall(k -> prevalence[observed[k]] >= config.min_prevalence, eachindex(observed))
        fits = isempty(tested) ? nothing : fit_lm(clr[:, tested], groups; levels)
        if !isnothing(fits)
            for (k, c) in enumerate(tested)
                f = fits[k, :]
                merge!(rows[observed[c]], Dict{String,Any}(
                    "status" => f.status, "note" => f.note, "estimate" => f.estimate,
                    "standard_error" => f.se, "statistic" => f.statistic,
                    "pvalue" => f.pvalue))
            end
        end
        method_text = "Gaussian linear model per taxon (stats::lm) on centred log-ratios, after " *
                      "multiplicative replacement of zeros (delta = $(config.replacement_delta) of " *
                      "each taxon's smallest observed count); t test of the group coefficient; " *
                      "Benjamini-Hochberg adjustment over the fitted taxa. Estimates are " *
                      "differences in CLR units, not fold changes, and replaced zeros are not measurements"
        effect = Dict("key" => "estimate", "label" => "CLR difference")
        size_factor_rows = Dict{String,Any}[]
        config_echo = Dict{String,Any}("method" => config.method,
                                       "min_prevalence" => config.min_prevalence,
                                       "replacement_delta" => config.replacement_delta)
        diagnostics["zero_replacement"] = Dict{String,Any}(
            "delta" => replaced.delta, "zeros_replaced" => replaced.zeros_replaced,
            "n_taxa_in_composition" => length(observed),
            "n_taxa_unobserved" => length(taxa) - length(observed),
            "max_imputed_fraction" => maximum(replaced.imputed_mass))
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
    merge!(diagnostics, Dict{String,Any}(
        "n_taxa" => length(taxa), "n_tested" => length(family), "n_failed" => n_failed,
        "n_boundary" => count(r -> r["status"] == "boundary", rows),
        "n_filtered" => count(r -> r["status"] == "filtered", rows)))
    Dict{String,Any}(
        "status" => n_failed == 0 ? "ok" : "partial",
        "method" => method_text,
        "effect" => effect,
        "groups" => Dict("reference" => String(reference), "contrast" => String(contrast)),
        "n_samples" => Dict(String(reference) => n_ref, String(contrast) => n_con),
        "size_factors" => size_factor_rows,
        "config" => config_echo,
        "diagnostics" => diagnostics,
        "rows" => rows,
    )
end

"""
    volcano_chart(result; alpha=0.05) -> Dict

A Plotly volcano plot of a `differential_abundance` result: the effect named by
`result["effect"]` (log2 fold change, or CLR difference) against -log10 of the
raw p-value, coloured by whether the BH-adjusted p-value is below `alpha`. Taxa
without a p-value are omitted here; they remain in the results table with their
reason.
"""
function volcano_chart(result::AbstractDict; alpha::Real=0.05)
    ref = result["groups"]["reference"]
    con = result["groups"]["contrast"]
    effect = get(result, "effect", Dict("key" => "log2_fold_change", "label" => "log2 fold change"))
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
            "x" => [r[effect["key"]] for r in sel], "y" => y,
            "text" => text, "hoverinfo" => "text+x+y",
            "marker" => Dict("color" => colour, "size" => 8),
        ))
    end
    layout = Dict{String,Any}(
        "title" => Dict("text" => "Differential abundance: $con vs $ref"),
        "xaxis" => Dict("title" => "$(effect["label"]) ($con vs $ref)", "zeroline" => true),
        "yaxis" => Dict("title" => "-log10 p"),
    )
    Dict("data" => traces, "layout" => layout)
end

end # module Differential
