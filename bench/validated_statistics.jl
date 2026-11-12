# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

#!/usr/bin/env julia
# Lightweight, reproducible benchmark for issue #1. This is not an estimator or
# correctness validation; the unit/reference suites provide those checks.
using JSON3
using MetaManifold
using MetaManifold.NumericPolicy
using MetaManifold.ExactSummaries
using MetaManifold.HighPrecision

function counts_fixture(n_features::Int, n_samples::Int)
    labels = ["sample_$j" for j in 1:n_samples]
    features = ["feature_$i" for i in 1:n_features]
    counts = [BigInt((37i + 101j + 17i*j) % 1000)
              for i in 1:n_features, j in 1:n_samples]
    counts, labels, features
end

function measure(f::Function; repetitions::Int=5)
    f() # warm up compilation outside all reported timings
    GC.gc()
    trials = [@timed f() for _ in 1:repetitions]
    times = sort([t.time for t in trials])
    bytes = sort([t.bytes for t in trials])
    middle = (repetitions + 1) ÷ 2
    Dict{String,Any}(
        "repetitions" => repetitions,
        "seconds" => [t.time for t in trials],
        "allocated_bytes" => [t.bytes for t in trials],
        "median_seconds" => times[middle],
        "median_allocated_bytes" => bytes[middle],
    )
end

function benchmark()
    measurements = Dict{String,Any}[]
    for (features, samples) in ((8, 2), (100, 8), (1000, 8))
        counts, labels, feature_labels = counts_fixture(features, samples)
        groups = ["group_$(j % 2)" for j in 1:samples]
        policy = numeric_policy(:exact_counts)
        result = measure(; repetitions=features == 1000 ? 3 : 5) do
            exact_summary(counts; sample_labels=labels,
                feature_labels=feature_labels, groups, policy)
        end
        result["operation"] = "exact_summary"
        result["precision_bits"] = nothing
        result["features"] = features
        result["samples"] = samples
        push!(measurements, result)
    end

    counts, labels, feature_labels = counts_fixture(100, 8)
    for bits in (64, 256, 1024)
        policy = numeric_policy(mode=:high_precision, precision_bits=bits)
        result = measure() do
            high_precision_summary(counts; sample_labels=labels,
                feature_labels=feature_labels, policy)
        end
        result["operation"] = "high_precision_summary"
        result["precision_bits"] = bits
        result["features"] = 100
        result["samples"] = 8
        push!(measurements, result)
    end

    # Concurrency evidence on the same input and precision, not a claim that a
    # single MPFR operation is parallel. Run the benchmark with -t2 or more.
    policy = numeric_policy(mode=:high_precision, precision_bits=256)
    work = () -> high_precision_summary(counts; sample_labels=labels,
        feature_labels=feature_labels, policy)
    sequential = measure(; repetitions=3) do
        foreach(_ -> work(), 1:4)
    end
    sequential["operation"] = "four_high_precision_summaries_sequential"
    concurrent = measure(; repetitions=3) do
        @sync for _ in 1:4
            Threads.@spawn work()
        end
    end
    concurrent["operation"] = "four_high_precision_summaries_concurrent"
    push!(measurements, sequential, concurrent)

    Dict{String,Any}(
        "schema" => "metamanifold.validated-statistics-benchmark.v1",
        "commit" => get(ENV, "GITHUB_SHA", "unknown"),
        "software" => Dict{String,Any}(
            "julia" => string(VERSION), "mpfr" => string(Base.MPFR.version())),
        "platform" => Dict{String,Any}(
            "os" => string(Sys.KERNEL), "arch" => string(Sys.ARCH),
            "cpu" => string(Sys.CPU_NAME), "julia_threads" => Threads.nthreads()),
        "interpretation" => "descriptive runtime/allocation baseline only; no statistical validation or universal speed claim",
        "measurements" => measurements,
    )
end

output = nothing
if !isempty(ARGS)
    if length(ARGS) == 2 && ARGS[1] == "--output"
        output = ARGS[2]
    else
        error("usage: julia --project=. -t2 bench/validated_statistics.jl [--output PATH]")
    end
end
result = benchmark()
encoded = JSON3.write(result)
if isnothing(output)
    println(encoded)
else
    mkpath(dirname(abspath(output)))
    write(output, encoded)
    println("Wrote validated-statistics benchmark to $output")
end
