# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

## Differential abundance
# Per-taxon models comparing two groups of samples, with Benjamini-Hochberg
# correction across the taxa that produced a p-value. Two methods:
#
# - nb_glm: a negative-binomial GLM per taxon (MASS::glm.nb) on the observed
#   integer read counts, with a library-size offset. There are no pseudocounts.
# - clr_welch: Welch's two-sample t-test per taxon (stats::t.test) on centred
#   log-ratios, over the taxa that pass `min_prevalence`, after their zeros are
#   replaced by zCompositions::cmultRepl (Bayesian-multiplicative, for counts).
#   The estimate is a difference in CLR units, not a fold change.
#
# Count matrices here are samples x taxa, as everywhere else in `Analysis`.
# Anything that is not a non-negative integer count is refused rather than
# rounded. A taxon whose fit fails is reported with its reason and kept out of
# the multiple-testing family; it never receives a stand-in p-value.
module Differential

using DataFrames, RCall
using ..RRuntime: with_r_lock
using ..Analysis: _palette_hex, R_WAIT_SECONDS
export DifferentialConfig, ScalingRefusal, MASSUnavailable, ZCompositionsUnavailable,
       OFFSET_METHODS, METHODS, tss_factors, rle_factors, size_factors, bh_adjust,
       validate_counts, replace_zeros, clr_transform, fit_nb, fit_welch,
       differential_abundance, volcano_chart

# Methods accepted under `analysis.differential.method`.
const METHODS = ("nb_glm", "clr_welch")

# Offset (size factor) methods accepted under `analysis.differential.offset`.
const OFFSET_METHODS = ("tss", "rle")

# Above this the NB dispersion parameter has run to its upper bound: the counts
# show no overdispersion and the model has collapsed to a Poisson. Below the
# lower bound the dispersion is degenerate the other way.
const THETA_UPPER = 1e7
const THETA_LOWER = 1e-8

# zCompositions::cmultRepl's default `frac`, as in Martín-Fernández et al. (2003).
const DEFAULT_DELTA = 0.65

"""
    DifferentialConfig(offset="tss", min_prevalence=0.0; method="nb_glm",
                       replacement_delta=0.65)

Settings for one differential abundance analysis, read from
`analysis.differential` in the pipeline config. `method` is `"nb_glm"` or
`"clr_welch"`. `offset` names the size factor method (`"tss"` or `"rle"`) and is
used by `nb_glm` only: a library-size offset has no meaning for a test on
log-ratios, which are already free of library size.
`min_prevalence` is the fraction of samples, in [0, 1], in which a taxon must
have at least one read to be tested. `replacement_delta`, in (0, 1), is the
`frac` passed to `zCompositions::cmultRepl` by `clr_welch`: an imputed proportion
above its taxon's smallest observed proportion is lowered to this fraction of it.
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
    ZCompositionsUnavailable()

Raised when the R package zCompositions, which provides `cmultRepl`, cannot be
loaded.
"""
struct ZCompositionsUnavailable <: Exception end

Base.showerror(io::IO, ::ZCompositionsUnavailable) = print(io,
    "the R package zCompositions is not available, so zeros cannot be replaced for " *
    "clr_welch; restore the R library from renv.lock (Rscript -e 'renv::restore()')")

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

# zCompositions::cmultRepl on the taxa that will be tested. z.warning = 1 turns
# off its own sparsity screen, and z.delete = FALSE stops it from dropping a
# row or column: either would silently change the taxa or misalign the samples
# with their groups. Which taxa are sparse enough to leave out is decided by
# `min_prevalence` before this runs. Errors and warnings are captured, not raised.
const _REPLACE_R = raw"""
da_warns <- character(0)
da_out <- tryCatch(
  withCallingHandlers(
    as.matrix(zCompositions::cmultRepl(da_counts, label = 0, method = "GBM", output = "prop",
                                       frac = da_frac, z.warning = 1, z.delete = FALSE,
                                       suppress.print = TRUE)),
    warning = function(w) {
      da_warns <<- c(da_warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
  error = function(e) conditionMessage(e))
"""

"""
    _check_compositions(out, zero_mask, samples) -> out

Check that cmultRepl's output `out` is one valid composition per sample: the
same shape as the input whose zeros `zero_mask` marks, and every row finite,
positive and summing to 1. Throws `ArgumentError` naming the first sample that
fails, pointing at `min_prevalence`.
"""
function _check_compositions(out::AbstractMatrix{<:Real}, zero_mask::AbstractMatrix{Bool},
                             samples::AbstractVector{<:AbstractString})
    size(out) == size(zero_mask) || error(
        "cmultRepl returned a $(size(out)) table for $(size(zero_mask)) input; a row or column was dropped")
    for i in axes(out, 1)
        row = view(out, i, :)
        if !all(v -> isfinite(v) && v > 0, row) || abs(sum(row) - 1) > 1e-8
            throw(ArgumentError(
                "zero replacement left sample '$(samples[i])' without a valid composition: " *
                "its $(count(view(zero_mask, i, :))) replaced zeros would take the whole sample. " *
                "Raise analysis.differential.min_prevalence so fewer rare taxa are tested"))
        end
    end
    out
end

"""
    replace_zeros(counts; frac=DEFAULT_DELTA, samples, taxa) -> NamedTuple

Replace the zeros of the samples x taxa read-count matrix `counts` with
`zCompositions::cmultRepl(counts, label = 0, method = "GBM", output = "prop",
frac = frac)`, the geometric Bayesian-multiplicative replacement for count
data. Returns `(proportions, zeros_replaced, imputed_mass, warnings)`:
`proportions` has the same shape as `counts`, each row summing to 1 with no
zero left; `imputed_mass` is the share of each sample's total held by replaced
values; `warnings` is whatever R warned.

The caller chooses the taxa: everything passed in is replaced and kept, none is
dropped. Throws `ArgumentError` when there are fewer than 2 taxa, when a taxon
has reads in fewer than 2 samples (the GBM prior is then undefined), when a
sample has no reads, or when cmultRepl refuses or returns a value that is not a
positive proportion; the messages name the taxon or sample and point at
`min_prevalence`. Throws `ZCompositionsUnavailable` when zCompositions cannot
be loaded and `RBusyError` when the R runtime stays busy.
"""
function replace_zeros(counts::AbstractMatrix{<:Real};
                       frac::Real=DEFAULT_DELTA,
                       samples::AbstractVector{<:AbstractString}=["sample $i" for i in axes(counts, 1)],
                       taxa::AbstractVector{<:AbstractString}=["taxon $j" for j in axes(counts, 2)])
    size(counts) == (length(samples), length(taxa)) || throw(ArgumentError(
        "the table is $(size(counts)) but there are $(length(samples)) samples and $(length(taxa)) taxa"))
    size(counts, 2) >= 2 || throw(ArgumentError(
        "the centred log-ratio needs at least 2 taxa; there are $(size(counts, 2))"))
    n = size(counts, 1)
    sparse = [j for j in axes(counts, 2) if count(>(0), view(counts, :, j)) < 2]
    isempty(sparse) || throw(ArgumentError(
        "$(length(sparse)) taxa have reads in fewer than 2 samples, so cmultRepl's " *
        "Bayesian-multiplicative prior cannot be estimated for them (" *
        join(("'$(taxa[j])'" for j in first(sparse, 5)), ", ") *
        (length(sparse) > 5 ? ", ..." : "") * "); raise analysis.differential.min_prevalence " *
        "to at least $(round(2 / n; digits=3)) so that only taxa seen in 2 or more samples are tested"))
    for i in axes(counts, 1)
        sum(view(counts, i, :)) > 0 || throw(ArgumentError(
            "sample '$(samples[i])' has no reads in the $(size(counts, 2)) taxa that pass " *
            "min_prevalence, so it has no composition to transform"))
    end

    zero_mask = counts .== 0
    n_zero = count(zero_mask)
    out, warns = with_r_lock(; timeout=R_WAIT_SECONDS[]) do
        RCall.rcopy(RCall.reval("requireNamespace('zCompositions', quietly = TRUE)")) ||
            throw(ZCompositionsUnavailable())
        # cmultRepl stops when there is no zero to replace; the closure of the
        # data is then the whole answer.
        n_zero == 0 && return (Float64.(counts) ./ sum(counts; dims=2), String[])
        RCall.globalEnv[:da_counts] = Matrix{Float64}(counts)
        RCall.globalEnv[:da_frac] = Float64(frac)
        try
            RCall.reval(_REPLACE_R)
            res = RCall.rcopy(RCall.reval("da_out"))
            res isa AbstractString && throw(ArgumentError(
                "zCompositions::cmultRepl refused the table: $res. Raising " *
                "analysis.differential.min_prevalence leaves fewer, better-observed taxa"))
            (Matrix{Float64}(res), String.(RCall.rcopy(Vector{String}, RCall.reval("da_warns"))))
        finally
            RCall.reval("rm(list = intersect(ls(), c('da_counts', 'da_frac', 'da_out', 'da_warns'))); " *
                        "invisible(gc())")
        end
    end
    _check_compositions(out, zero_mask, samples)
    mass = [sum(out[i, j] for j in axes(out, 2) if zero_mask[i, j]; init=0.0) for i in axes(out, 1)]
    (proportions=out, zeros_replaced=n_zero, imputed_mass=mass, warnings=unique(warns))
end

# The per-taxon Welch t-test loop, with the same `da_` naming and statuses as
# `_FIT_R` but no `boundary`: a t-test has no dispersion parameter to run to a bound.
const _FIT_WELCH_R = raw"""
da_g <- factor(da_group, levels = da_levels)
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
  da_tt <- tryCatch(
    withCallingHandlers(
      stats::t.test(da_y[da_g == da_levels[2]], da_y[da_g == da_levels[1]], var.equal = FALSE),
      warning = function(w) {
        da_warns <<- c(da_warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  if (inherits(da_tt, "error")) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- paste("t.test stopped with an error:", conditionMessage(da_tt))
    next
  }
  da_row <- c(unname(da_tt[["estimate"]][1] - da_tt[["estimate"]][2]), unname(da_tt[["stderr"]]),
              unname(da_tt[["statistic"]]), da_tt[["p.value"]])
  if (!all(is.finite(da_row)) || da_row[[2]] <= 0 || da_row[[4]] < 0 || da_row[[4]] > 1) {
    da_result[da_j, "status"] <- "failed"
    da_result[da_j, "note"] <- paste(c("the test returned a non-finite or zero standard error, so there is no test (both groups have zero variance)",
                                       unique(da_warns)), collapse = "; ")
    next
  }
  da_result[da_j, c("estimate", "se", "statistic", "pvalue")] <- da_row
  if (length(da_warns) > 0L) {
    da_result[da_j, "note"] <- paste("R warned:", paste(unique(da_warns), collapse = "; "))
  }
}
"""

"""
    fit_welch(values, groups; levels) -> DataFrame

Welch's two-sample t-test, `stats::t.test(contrast, reference, var.equal =
FALSE)`, on each column of the samples x taxa matrix `values` (for `clr_welch`,
centred log-ratios). `levels` is `(reference, contrast)`, so the estimate is
the mean of `contrast` minus the mean of `reference`; the variances of the two
groups are not assumed equal, since the groups often differ in size. Returns
one row per taxon with `status` (`ok` or `failed`), `note`, `estimate`, `se`,
`statistic` (Welch's t) and `pvalue`; a failed taxon has `nothing` for every
statistic. Throws `RBusyError` when the R runtime stays busy.
"""
function fit_welch(values::AbstractMatrix{<:Real}, groups::AbstractVector{<:AbstractString};
                   levels::NTuple{2,String})
    size(values, 1) == length(groups) || throw(ArgumentError(
        "values have $(size(values, 1)) samples and groups $(length(groups))"))
    raw = with_r_lock(; timeout=R_WAIT_SECONDS[]) do
        RCall.globalEnv[:da_values] = Matrix{Float64}(values)
        RCall.globalEnv[:da_group] = String.(groups)
        RCall.globalEnv[:da_levels] = collect(levels)
        try
            RCall.reval(_FIT_WELCH_R)
            DataFrame(RCall.rcopy(RCall.reval("da_result")))
        finally
            RCall.reval("rm(list = intersect(ls(), c('da_values', 'da_group', 'da_levels', " *
                        "'da_g', 'da_n', 'da_result', 'da_j', 'da_y', 'da_warns', 'da_tt', " *
                        "'da_row'))); invisible(gc())")
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

`clr_welch`, in this order: (1) `min_prevalence` chooses the taxa, and taxa with
no reads in any sample are always left out; (2) the zeros of those taxa only
are replaced by `zCompositions::cmultRepl` (GBM, `frac =
config.replacement_delta`), so the rare taxa that are not tested contribute no
imputed values; (3) each sample is transformed to centred log-ratios over the
same taxa; (4) each taxon is tested by Welch's t-test. Replacing zeros across
every observed taxon first would let imputed values for rare taxa take over
the samples of a sparse table. The estimate is the difference in mean CLR,
contrast minus reference: a difference in centred-log-ratio units, not a fold
change, relative to the geometric mean of the tested taxa. Replaced zeros are
not measurements; `diagnostics.zero_replacement` reports their share of each
sample. Each group needs at least 2 samples.

For both methods, taxa present in fewer than `config.min_prevalence` of the
samples are reported as `filtered` and not fitted, and the Benjamini-Hochberg
family is the fitted taxa that produced a p-value; failed and filtered taxa
carry `padj = nothing`. The result's `effect` names the row key that holds the
effect size and its label.

Throws `ArgumentError` for unusable input (including a zero replacement that
cannot be made), `ScalingRefusal` when the offset cannot be computed,
`MASSUnavailable` or `ZCompositionsUnavailable` when the R package a method
needs is missing, and an `ErrorException` when no taxon could be fitted at all.
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
        (n_ref >= 2 && n_con >= 2) || throw(ArgumentError(
            "Welch's t-test needs at least 2 samples in each group to estimate its variance; " *
            "'$reference' has $n_ref and '$contrast' has $n_con"))
        for j in findall(==(0), present)
            rows[j]["note"] = "no reads in any sample of either group, so it is not part of the composition"
        end
        tested = findall(j -> present[j] > 0 && prevalence[j] >= config.min_prevalence, axes(counts, 2))
        length(tested) >= 2 || throw(ArgumentError(
            "the centred log-ratio needs at least 2 taxa with reads at or above " *
            "analysis.differential.min_prevalence = $(config.min_prevalence); there are $(length(tested))"))
        replaced = replace_zeros(counts[:, tested]; frac=config.replacement_delta,
                                 samples, taxa=taxa[tested])
        clr = clr_transform(replaced.proportions)
        fits = fit_welch(clr, groups; levels)
        for (k, j) in enumerate(tested)
            f = fits[k, :]
            merge!(rows[j], Dict{String,Any}(
                "status" => f.status, "note" => f.note, "estimate" => f.estimate,
                "standard_error" => f.se, "statistic" => f.statistic,
                "pvalue" => f.pvalue))
        end
        method_text = "Welch's two-sample t-test per taxon (stats::t.test, var.equal = FALSE) on " *
                      "centred log-ratios over the taxa passing min_prevalence, after their zeros " *
                      "were replaced by zCompositions::cmultRepl (GBM, frac = $(config.replacement_delta)); " *
                      "Benjamini-Hochberg adjustment over the fitted taxa. Estimates are " *
                      "differences in CLR units, not fold changes, and replaced zeros are not measurements"
        effect = Dict("key" => "estimate", "label" => "CLR difference")
        size_factor_rows = Dict{String,Any}[]
        config_echo = Dict{String,Any}("method" => config.method,
                                       "min_prevalence" => config.min_prevalence,
                                       "replacement_delta" => config.replacement_delta)
        sorted_mass = sort(replaced.imputed_mass)
        diagnostics["zero_replacement"] = Dict{String,Any}(
            "method" => "zCompositions::cmultRepl (GBM)",
            "delta" => config.replacement_delta, "zeros_replaced" => replaced.zeros_replaced,
            "n_taxa_in_composition" => length(tested),
            "n_taxa_unobserved" => count(==(0), present),
            "median_imputed_fraction" => _median(sorted_mass),
            "max_imputed_fraction" => last(sorted_mass),
            "warnings" => replaced.warnings)
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
