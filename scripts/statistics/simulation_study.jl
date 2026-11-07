#!/usr/bin/env julia
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Seeded simulation study for the validated statistics layer.
#
# This is a study of the implementation, not of an estimator. It asks two
# questions with pre-set tolerances:
#
#   1. Exact arithmetic: over random integer matrices, does every proportion
#      reconstruct its count exactly and sum to exactly one? (tolerance 0)
#   2. High precision: how large is the absolute error of Shannon entropy and
#      Gini-Simpson diversity at 128/256/512 bits, measured against an
#      separately evaluated 4096-bit computation from the same exact proportions?
#
# Tolerances are declared before the run: `tolerance(bits) = 2^-(bits-16)`,
# This is an absolute-error budget, not a claim about units in the last place.
# Both calculations use MPFR; this study checks precision consistency, while
# test/references/check_biological_reference.py uses an independent implementation.
# The study fails if any replicate
# exceeds its tolerance; it never widens a tolerance after seeing a result.
#
# Usage:
#   julia --project=. scripts/statistics/simulation_study.jl \
#       [--replicates 40] [--seed 20261023] [--output FILE]
#
# Exit status: 0 = every check within tolerance, 1 = a tolerance was exceeded
# or an exact identity failed, 2 = bad invocation.

using Random
using JSON3
using Dates

using MetaManifold
using MetaManifold.NumericPolicy
using MetaManifold.ExactSummaries
using MetaManifold.HighPrecision

const REFERENCE_BITS = 4096
const PRECISIONS = (128, 256, 512)
const EXACT_REPLICATES = 200

"""Allowed absolute error at `bits`, declared before the run."""
tolerance(bits::Integer) = big(2.0)^(-(Int(bits) - 16))

"""Draw a random counts matrix. Zeros are included deliberately: a composition
with an absent taxon is a different shape from a composition with a tiny one,
and the zero-total case must stay undefined rather than become zero."""
function draw_counts(rng, features::Int, samples::Int)
    counts = Matrix{BigInt}(undef, features, samples)
    for j in 1:samples, i in 1:features
        counts[i, j] = rand(rng) < 0.12 ? big(0) : BigInt(rand(rng, 1:1_000_000))
    end
    # Every sample has at least one read; a zero-total sample is covered by the
    # dedicated zero-state tests, not by a division by zero here.
    counts[1, :] .= BigInt(rand(rng, 1:1_000_000))
    counts
end

"""Higher-precision MPFR reference from exact proportions."""
function reference_metrics(column::AbstractVector{BigInt})
    # The do-block binds to the call it is written against, so `first(f() do … end)`
    # would pass the function to `first`. Bind the result first, then unwrap it.
    computed = with_precision(; bits=REFERENCE_BITS) do
        total = sum(column)
        shannon = zero(BigFloat)
        squares = big(0)
        for count in column
            squares += count * count
            iszero(count) && continue
            p = BigFloat(count) / BigFloat(total)
            shannon -= p * log(p)
        end
        simpson = BigFloat((total * total - squares) // (total * total))
        (shannon=shannon, simpson=simpson)
    end
    first(computed)
end

"""Measure errors at reference precision, independent of the caller's MPFR state."""
function metric_errors(observed::BigFloat, reference::BigFloat)
    computed = with_precision(; bits=REFERENCE_BITS) do
        absolute = abs(observed - reference)
        (absolute=absolute, relative=absolute / max(one(BigFloat), abs(reference)))
    end
    first(computed)
end

function main(args::Vector{String})
    replicates = 40
    seed = 20261023
    output = nothing
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument == "--replicates" && index < length(args)
            replicates = parse(Int, args[index + 1]); index += 2
        elseif argument == "--seed" && index < length(args)
            seed = parse(Int, args[index + 1]); index += 2
        elseif argument == "--output" && index < length(args)
            output = args[index + 1]; index += 2
        else
            println(stderr, "simulation_study: unknown or incomplete argument '$argument'")
            return 2
        end
    end
    replicates >= 1 || (println(stderr, "simulation_study: --replicates must be positive"); return 2)

    rng = MersenneTwister(seed)
    exact_failures = String[]
    metric_failures = String[]
    maxima = Dict{Int,Dict{String,BigFloat}}(
        bits => Dict("shannon_abs" => big(0.0), "shannon_rel" => big(0.0),
                     "simpson_abs" => big(0.0), "simpson_rel" => big(0.0))
        for bits in PRECISIONS)

    # 1. Exact identity study: tolerance is zero, not a small number.
    exact_policy = numeric_policy(; mode=:exact_counts, max_denominator_bits=4096)
    for replicate in 1:EXACT_REPLICATES
        features = rand(rng, 1:8)
        samples = rand(rng, 1:5)
        counts = draw_counts(rng, features, samples)
        summary = exact_summary(counts;
            sample_labels=["s$j" for j in 1:samples],
            feature_labels=["f$i" for i in 1:features], policy=exact_policy)
        for (j, sample) in pairs(summary.samples)
            for i in 1:features
                if sample.proportions[i] * sample.total != counts[i, j]
                    push!(exact_failures, "replicate $replicate sample $j feature $i: proportion does not reconstruct the count")
                end
            end
            if sum(sample.proportions) != big(1) // big(1)
                push!(exact_failures, "replicate $replicate sample $j: proportions do not sum to one")
            end
        end
    end

    # 2. Precision study against a higher-precision MPFR computation.
    for replicate in 1:replicates
        features = rand(rng, 2:8)
        counts = draw_counts(rng, features, 1)
        column = counts[:, 1]
        reference = reference_metrics(column)
        for bits in PRECISIONS
            policy = numeric_policy(; mode=:high_precision, precision_bits=bits)
            result = high_precision_summary(counts;
                sample_labels=["s"], feature_labels=["f$i" for i in 1:features],
                policy)
            sample = result.samples[1]
            allowed = tolerance(bits)
            shannon_errors = metric_errors(sample.shannon, reference.shannon)
            simpson_errors = metric_errors(sample.gini_simpson, reference.simpson)
            shannon_error, shannon_relative = shannon_errors.absolute, shannon_errors.relative
            simpson_error, simpson_relative = simpson_errors.absolute, simpson_errors.relative
            maxima[bits]["shannon_abs"] = max(maxima[bits]["shannon_abs"], shannon_error)
            maxima[bits]["simpson_abs"] = max(maxima[bits]["simpson_abs"], simpson_error)
            maxima[bits]["shannon_rel"] = max(maxima[bits]["shannon_rel"], shannon_relative)
            maxima[bits]["simpson_rel"] = max(maxima[bits]["simpson_rel"], simpson_relative)
            shannon_error <= allowed || push!(metric_failures,
                "replicate $replicate at $bits bits: Shannon absolute error $shannon_error exceeds $allowed")
            simpson_error <= allowed || push!(metric_failures,
                "replicate $replicate at $bits bits: Gini-Simpson absolute error $simpson_error exceeds $allowed")
            sample.shannon >= 0 || push!(metric_failures, "replicate $replicate at $bits bits: negative Shannon entropy")
            sample.gini_simpson >= 0 || push!(metric_failures, "replicate $replicate at $bits bits: negative Gini-Simpson")
        end
    end

    complete = isempty(exact_failures) && isempty(metric_failures)
    report = Dict{String,Any}(
        "schema" => "metamanifold.statistics-simulation-study.v1",
        "seed" => seed,
        "exact_replicates" => EXACT_REPLICATES,
        "precision_replicates" => replicates,
        "reference_bits" => REFERENCE_BITS,
        "comparison_bits" => REFERENCE_BITS,
        "precisions" => collect(PRECISIONS),
        "tolerance_rule" => "2^-(bits-16) absolute, declared before the run",
        "tolerances" => Dict(string(bits) => string(tolerance(bits)) for bits in PRECISIONS),
        "max_errors" => Dict(string(bits) => Dict(key => string(value) for (key, value) in maxima[bits])
                             for bits in PRECISIONS),
        "exact_failures" => exact_failures,
        "metric_failures" => metric_failures,
        "complete" => complete,
        "software" => Dict(
            "julia" => string(VERSION),
            "mpfr" => string(Base.MPFR.version()),
            "generated" => string(Dates.now(UTC)),
        ),
    )
    text = JSON3.write(report)
    if isnothing(output)
        println(text)
    else
        open(output, "w") do io
            JSON3.pretty(io, report)
            println(io)
        end
        println("simulation_study: $(EXACT_REPLICATES) exact and $(replicates) precision replicates; ",
                complete ? "all within tolerance" : "TOLERANCE EXCEEDED", "; report ", output)
    end
    if !complete
        for failure in Iterators.flatten((exact_failures, metric_failures))
            println(stderr, "simulation_study: ", failure)
        end
    end
    return complete ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
