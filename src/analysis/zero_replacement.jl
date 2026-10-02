# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

## Zero replacement for compositional transforms
# Multiplicative replacement (Martín-Fernández, Barceló-Vidal & Pawlowsky-Glahn
# 2003, Mathematical Geology 35(3):253-278), the operator of
# `zCompositions::multRepl` with `frac = delta` and a per-taxon detection limit.
# A log-ratio transform is undefined at zero, so a zero has to be given some
# positive value first; this operator does that while keeping each sample's
# total and the ratios among its observed taxa exactly as they were, which a
# pseudocount does not.
#
# Matrices here are samples x taxa, as everywhere else in `Analysis` (the
# reference and the original paper write compositions as rows too).
#
# A replaced value is not a measurement. The replacement is biased, and no rule
# that sees only the observed data can avoid that; the result records it.
module ZeroReplacement

export DEFAULT_DELTA, ReplacementOutcome, multiplicative_replacement,
       replacement_invariants_hold

"""
    DEFAULT_DELTA

The fraction of a taxon's detection limit that a replaced zero receives by
default, 0.65, as in `zCompositions::multRepl` and the original paper.
"""
const DEFAULT_DELTA = 0.65

# Relative tolerance of the runtime invariant checks.
const TOLERANCE = 1e-9

"""
    ReplacementOutcome

The result of `multiplicative_replacement`: the replaced samples x taxa table,
the `delta` used, the per-taxon `detection_limits`, the number of
`zeros_replaced`, and `imputed_mass`, the fraction of each sample's total now
held by replaced values.
"""
struct ReplacementOutcome
    counts           :: Matrix{Float64}
    delta            :: Float64
    detection_limits :: Vector{Float64}
    zeros_replaced   :: Int
    imputed_mass     :: Vector{Float64}
end

"""
    _invariant_violation(original, replaced) -> Union{String,Nothing}

The first way in which `replaced` fails to be a multiplicative replacement of
`original`, as a sentence, or `nothing` when it is one: same shape, every
sample total preserved, every observed taxon in a sample scaled by the same
factor, and no zero left.
"""
function _invariant_violation(original::AbstractMatrix{<:Real}, replaced::AbstractMatrix{<:Real})
    size(original) == size(replaced) ||
        return "the replaced table is $(size(replaced)) for an input of $(size(original))"
    for i in axes(original, 1)
        before = sum(view(original, i, :))
        before > 0 || return "sample $i has no reads"
        after = sum(view(replaced, i, :))
        abs(after - before) <= TOLERANCE * before ||
            return "sample $i changed total from $before to $after"
        scale = nothing
        for j in axes(original, 2)
            if original[i, j] > 0
                r = replaced[i, j] / original[i, j]
                scale = something(scale, r)
                abs(r - scale) <= TOLERANCE * max(1.0, abs(scale)) ||
                    return "sample $i scaled taxon $j by $r and an earlier taxon by $scale"
            elseif !(replaced[i, j] > 0)
                return "sample $i still holds $(replaced[i, j]) for taxon $j"
            end
        end
    end
    nothing
end

"""
    replacement_invariants_hold(original, replaced) -> Bool

Whether `replaced` keeps every sample total of `original`, scales the observed
taxa of each sample by one common factor, and leaves no zero.
"""
replacement_invariants_hold(original::AbstractMatrix{<:Real}, replaced::AbstractMatrix{<:Real}) =
    isnothing(_invariant_violation(original, replaced))

"""
    multiplicative_replacement(counts; delta=DEFAULT_DELTA, samples, taxa) -> ReplacementOutcome

Replace the zeros of the samples x taxa table `counts`. With `DL_j` the
detection limit of taxon `j` (its smallest observed value anywhere in the
table), `S_i` the total of sample `i`, `Z_i` its zero taxa and
`Δ_i = Σ_{j ∈ Z_i} delta·DL_j / S_i`:

    x̃_ij = delta · DL_j        for j ∈ Z_i
    x̃_ij = (1 - Δ_i) · x_ij    otherwise

Each sample total and the ratios among each sample's observed taxa are kept
exactly; every replaced value is positive and below its taxon's detection
limit. On a table of proportions this equals `zCompositions::multRepl(X,
label = 0, dl, frac = delta)`; on counts with unequal totals, dividing each row
by its total gives the reference's output.

`samples` and `taxa` name rows and columns in error messages. Throws
`ArgumentError` when `delta` is not in (0, 1), when a value is negative or not
finite, when a taxon is never observed (its detection limit cannot be
estimated), when a sample has no reads, or when a sample's replaced values
would take its whole total (`Δ_i ≥ 1`); that message gives the largest
`delta` the sample admits.
"""
function multiplicative_replacement(counts::AbstractMatrix{<:Real};
                                    delta::Real=DEFAULT_DELTA,
                                    samples::AbstractVector{<:AbstractString}=["sample $i" for i in axes(counts, 1)],
                                    taxa::AbstractVector{<:AbstractString}=["taxon $j" for j in axes(counts, 2)])
    δ = Float64(delta)
    0 < δ < 1 || throw(ArgumentError(
        "the replacement delta must be in (0, 1), not $delta: it is the fraction of a " *
        "taxon's detection limit that a replaced zero receives (0.65 is the published default)"))
    size(counts) == (length(samples), length(taxa)) || throw(ArgumentError(
        "the table is $(size(counts)) but there are $(length(samples)) samples and $(length(taxa)) taxa"))
    X = Matrix{Float64}(counts)
    for j in axes(X, 2), i in axes(X, 1)
        (isfinite(X[i, j]) && X[i, j] >= 0) || throw(ArgumentError(
            "the value for taxon '$(taxa[j])' in sample '$(samples[i])' is $(X[i, j]); " *
            "zero replacement needs non-negative finite values"))
    end

    dl = Vector{Float64}(undef, size(X, 2))
    for j in axes(X, 2)
        observed = filter(>(0), view(X, :, j))
        isempty(observed) && throw(ArgumentError(
            "taxon '$(taxa[j])' has no reads in any sample, so its detection limit " *
            "cannot be estimated from this table"))
        dl[j] = minimum(observed)
    end

    out = copy(X)
    mass = zeros(Float64, size(X, 1))
    n_zero = 0
    for i in axes(X, 1)
        total = sum(view(X, i, :))
        total > 0 || throw(ArgumentError("sample '$(samples[i])' has no reads"))
        zeros_i = findall(==(0), view(X, i, :))
        isempty(zeros_i) && continue
        limit = sum(dl[zeros_i])
        Δ = δ * limit / total
        Δ < 1 || throw(ArgumentError(
            "sample '$(samples[i])' has $(length(zeros_i)) zeros whose detection limits sum " *
            "to $limit against a total of $total, so delta = $δ would give replaced values " *
            "$(round(100Δ; digits=1))% of the sample. The largest delta this sample admits " *
            "is below $(total / limit); lower delta or raise min_prevalence"))
        for j in axes(X, 2)
            out[i, j] = X[i, j] == 0 ? δ * dl[j] : (1 - Δ) * X[i, j]
        end
        mass[i] = Δ
        n_zero += length(zeros_i)
    end

    problem = _invariant_violation(X, out)
    isnothing(problem) || error("multiplicative replacement broke its own invariant: $problem")
    ReplacementOutcome(out, δ, dl, n_zero, mass)
end

end # module ZeroReplacement
