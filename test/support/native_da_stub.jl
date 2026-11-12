# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: AGPL-3.0-or-later

## A stand-in for CompositionalDA
#
# `CompositionalDA` is not in the General registry, so it is not a declared
# dependency and `Pkg.instantiate()` cannot install it: in CI,
# `compositional_da_available()` is false, and the adapter in
# `src/analysis/differential.jl` that turns a `DAResult` into result rows is
# unreachable through the package.
#
# `differential_abundance(...; native_backend=…)` is the seam: the adapter calls
# a module, so a test can supply this one. The module below offers the two entry
# points the adapter uses (`aldex2`, `ancombc`) and the error type it translates
# (`CompositionError`), with the column names and the verdict strings of the
# real package, so the adapter code under test is the code that runs in
# production. Only the estimator behind the boundary differs, and the result
# records which module produced it under `diagnostics.native_method.package`,
# so a record made with this stand-in cannot be read as one made with the
# package.
#
# The values are deliberately simple arithmetic on the feature's position, so
# every assertion in `test/unit/test_differential_native.jl` can name the number
# it expects rather than compare two implementations of the same model.
#
# Loaded from `test/runtests.jl` at top level: `module` is not valid inside the
# `@testset` that includes the unit files.

module NativeDAStub

using DataFrames

struct CompositionError <: Exception
    msg::String
end

Base.showerror(io::IO, e::CompositionError) = print(io, "CompositionError: ", e.msg)

# Every call the adapter made, in order, with the arguments it used. The tests
# read these to check the wiring (which taxa were sent, in what order, with
# which labels and settings) rather than inferring it from the result.
const CALLS = Any[]

# What the next call should do: `:none` answers normally; `:composition` raises
# the package's own error; `:other` raises something else; `:nonfinite` answers
# with a NaN statistic for the first feature.
const FAILURE = Ref{Symbol}(:none)

# Append a feature the adapter never asked about, to check that a row the caller
# has no place for is skipped instead of invented.
const EXTRA_FEATURE = Ref{Bool}(false)

function reset!()
    empty!(CALLS)
    FAILURE[] = :none
    EXTRA_FEATURE[] = false
    nothing
end

last_call() = last(CALLS)

function _record(method::Symbol, counts, features, groups, settings::Dict{Symbol,Any})
    push!(CALLS, (method=method, counts=copy(counts),
                  features=collect(String.(features)), groups=collect(String.(groups)),
                  settings=settings))
    nothing
end

function _dispatch_failure()
    FAILURE[] === :composition && throw(CompositionError("stub refused the composition"))
    FAILURE[] === :other && error("stub backend failed for a reason of its own")
    nothing
end

"""
    aldex2(counts, features, groups; mc_samples=128, seed=1)

The same signature and column names as `CompositionalDA.aldex2`. `diff_btw` is
the estimate the adapter reports, `effect` the within-instance effect it keeps
beside it, and `qvalue` is deliberately not `bh_adjust(pvalue)`: the adapter must
publish ALDEx2's expected within-instance BH values rather than adjust a second
time, and a stand-in whose two agreed could not tell the two apart.
"""
function aldex2(counts::AbstractMatrix{<:Integer}, features::AbstractVector,
                groups::AbstractVector; mc_samples::Integer=128, seed::Integer=1)
    _record(:aldex2, counts, features, groups,
            Dict{Symbol,Any}(:mc_samples => mc_samples, :seed => seed))
    _dispatch_failure()
    n = length(features)
    p = [0.01 * j for j in 1:n]
    table = DataFrame(
        feature   = String.(features),
        effect    = [0.25 * j for j in 1:n],
        pvalue    = p,
        qvalue    = p .* 0.5,
        wi_pvalue = p .* 0.75,
        wi_qvalue = p .* 0.375,
        diff_btw  = [0.5 * j for j in 1:n],
        diff_win  = [0.125 * j for j in 1:n],
        rab_all   = [1.0 * j for j in 1:n])
    FAILURE[] === :nonfinite && (table.pvalue[1] = NaN)
    EXTRA_FEATURE[] && push!(table,
        ("not_requested", 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0))
    (table = table,
     provenance = (package_version = v"0.0.0-stub",
                   parameters = Dict{Symbol,Any}(:mc_samples => mc_samples, :seed => seed,
                                                 :prior => 0.5)))
end

"""
    ancombc(counts, features, groups; prevalence=0.0, structural_zeros=true, conserve=true)

The same signature, columns and verdict strings as `CompositionalDA.ancombc`. A
feature that is zero in every sample of one condition is `"structural_zero"`,
one that is zero in both is `"structural_zero_both"`, and the rest are
`"tested"`.

A feature named `structural_both` is reported as structurally absent in both
conditions even when it has reads. Real data cannot produce that verdict through
`differential_abundance` (an all-zero column is filtered out before the backend
is called), but the adapter still has to map the verdict, and a mapping nobody
exercises is where a wrong message hides.
"""
function ancombc(counts::AbstractMatrix{<:Integer}, features::AbstractVector,
                 groups::AbstractVector; prevalence::Real=0.0,
                 structural_zeros::Bool=true, conserve::Bool=true)
    _record(:ancombc, counts, features, groups,
            Dict{Symbol,Any}(:prevalence => prevalence,
                             :structural_zeros => structural_zeros,
                             :conserve => conserve))
    _dispatch_failure()
    labels = String.(groups)
    n = length(features)
    status = map(1:n) do j
        String(features[j]) == "structural_both" && return "structural_zero_both"
        zero0 = sum(counts[labels .== "0", j]) == 0
        zero1 = sum(counts[labels .== "1", j]) == 0
        zero0 && zero1 && return "structural_zero_both"
        (zero0 || zero1) && return "structural_zero"
        "tested"
    end
    # `effect / log(2)` is the adapter's log2 fold change, so the effect is a
    # half-integer multiple of log(2) and the expected value is exactly j/2.
    effect = [log(2) * j / 2 for j in 1:n]
    se = [0.1 * j for j in 1:n]
    # Structural decisions carry no numerical test: NaN, exactly as the package
    # leaves them. The adapter must branch on the verdict *before* it looks at
    # finiteness, or these rows would raise instead of being reported.
    p = [status[j] == "tested" ? 0.02 * j : NaN for j in 1:n]
    table = DataFrame(
        feature   = String.(features),
        effect    = effect,
        se        = se,
        statistic = effect ./ se,
        pvalue    = p,
        qvalue    = p .* 0.25,
        status    = status)
    if FAILURE[] === :nonfinite
        first_tested = findfirst(==("tested"), status)
        first_tested === nothing || (table.pvalue[first_tested] = NaN)
    end
    EXTRA_FEATURE[] && push!(table, ("not_requested", 0.0, 1.0, 0.0, 1.0, 1.0, "tested"))
    (table = table,
     provenance = (package_version = v"0.0.0-stub",
                   parameters = Dict{Symbol,Any}(:prevalence => prevalence,
                                                 :structural_zeros => structural_zeros,
                                                 :conserve => conserve)))
end

end # module NativeDAStub
