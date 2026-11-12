# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Property-based checks for the exact and high-precision layers. The random
# matrices are seeded, so a failure is reproducible from the printed trial
# number; shrinking is done by re-running a single trial with its parameters.
using Test
using Random
using MetaManifold
using MetaManifold.NumericPolicy
using MetaManifold.ExactSummaries
using MetaManifold.HighPrecision

const _PROP_BITS = 512
const _PROP_TOLERANCE = big"1e-70"

# A do-block binds to the call it is written against, so `first(f() do … end)`
# would hand the function to `first`. Each helper binds and then unwraps.
function _reference_shannon(column)
    computed = with_precision(; bits=_PROP_BITS) do
        total = sum(big.(column))
        iszero(total) && return nothing
        h = zero(BigFloat)
        for count in column
            iszero(count) && continue
            p = BigFloat(count) / BigFloat(total)
            h -= p * log(p)
        end
        h
    end
    first(computed)
end

function _reference_simpson(column)
    computed = with_precision(; bits=_PROP_BITS) do
        total = sum(big.(column))
        iszero(total) && return nothing
        squares = sum(big(c)^2 for c in column)
        BigFloat((total * total - squares) // (total * total))
    end
    first(computed)
end

_close(reference, candidate) =
    abs(candidate - reference) <= _PROP_TOLERANCE * max(big(1), abs(reference))

function _prop_matrix(rng, features, samples)
    matrix = Matrix{BigInt}(undef, features, samples)
    for i in 1:features, j in 1:samples
        matrix[i, j] = rand(rng) < 0.12 ? big(2)^rand(rng, 70:90) :
                                          BigInt(rand(rng, 0:10^6))
    end
    matrix
end

@testset "Property: exact proportions reconstruct counts" begin
    rng = MersenneTwister(20261023)
    policy = numeric_policy(; mode=:exact_counts, round_digits=8)
    for trial in 1:25
        features = rand(rng, 1:5)
        samples = rand(rng, 1:4)
        counts = _prop_matrix(rng, features, samples)
        labels = ["s$j" for j in 1:samples]
        names = ["f$i" for i in 1:features]
        summary = exact_summary(counts; sample_labels=labels,
            feature_labels=names, policy)

        for (j, sample) in pairs(summary.samples)
            column = counts[:, j]
            @test sample.total == sum(column)
            @test sample.counts == column
            if iszero(sample.total)
                @test all(isnothing, sample.proportions)
                @test sample.warnings ==
                      ["'$(labels[j])' has zero reads: relative abundances are undefined, not zero"]
                continue
            end
            @test all(p -> !isnothing(p), sample.proportions)
            for i in 1:features
                @test sample.proportions[i] == counts[i, j] // sample.total
                @test sample.proportions[i] * sample.total == counts[i, j]
            end
            # Exact rationals sum to exactly one; no floating point is involved.
            @test sum(sample.proportions) == big(1) // big(1)
        end

        # Scaling every count by k leaves proportions and ordering untouched.
        k = rand(rng, 2:9)
        scaled = exact_summary(k .* counts; sample_labels=labels,
            feature_labels=names, policy)
        for j in eachindex(labels)
            @test scaled.samples[j].total == k * summary.samples[j].total
            @test scaled.samples[j].proportions == summary.samples[j].proportions
        end
        # Storage strings round-trip without a lossy numeric conversion.
        stored = summary_to_storage(summary)
        for (j, entry) in pairs(stored["samples"])
            total_text = entry["total"]
            @test (total_text isa AbstractString ? parse_exact_count(total_text) :
                   BigInt(total_text)) == summary.samples[j].total
            for (i, text) in pairs(entry["counts"])
                # Small exact integers remain JSON numbers; big ones are strings.
                parsed = text isa AbstractString ? parse_exact_count(text) : BigInt(text)
                @test parsed == summary.samples[j].counts[i]
            end
            for (i, text) in pairs(entry["proportions"])
                isnothing(text) && continue
                @test parse_exact_rational(String(text)) == summary.samples[j].proportions[i]
            end
        end
    end
end

@testset "Property: Shannon and Gini-Simpson satisfy independent references" begin
    rng = MersenneTwister(424242)
    policy = numeric_policy(; mode=:high_precision, precision_bits=256)
    for trial in 1:20
        features = rand(rng, 1:6)
        samples = rand(rng, 1:4)
        counts = _prop_matrix(rng, features, samples)
        summary = high_precision_summary(counts;
            sample_labels=["s$j" for j in 1:samples],
            feature_labels=["f$i" for i in 1:features], policy)

        for (j, sample) in pairs(summary.samples)
            column = counts[:, j]
            @test sample.total == sum(column)
            @test sample.richness == count(!iszero, column)
            if iszero(sample.total)
                @test isnothing(sample.shannon)
                @test isnothing(sample.gini_simpson)
                continue
            end
            shannon = _reference_shannon(column)
            simpson = _reference_simpson(column)
            @test shannon isa BigFloat
            @test _close(shannon, sample.shannon)
            @test _close(simpson, sample.gini_simpson)
            @test sample.shannon >= 0
            @test sample.shannon <= log(BigFloat(features)) + big"1e-70"
            @test sample.gini_simpson >= 0
            @test sample.gini_simpson <= big(1)
            @test sample.gini_simpson <= (features - 1) / features + big"1e-70"
        end

        # A zero-count feature adds no diversity and does not change any value.
        padded = vcat(counts, zeros(BigInt, 1, samples))
        padded_summary = high_precision_summary(padded;
            sample_labels=["s$j" for j in 1:samples],
            feature_labels=vcat(["f$i" for i in 1:features], ["absent"]), policy)
        for j in eachindex(summary.samples)
            @test padded_summary.samples[j].richness == summary.samples[j].richness
            if !isnothing(summary.samples[j].shannon)
                @test padded_summary.samples[j].shannon == summary.samples[j].shannon
                @test padded_summary.samples[j].gini_simpson ==
                      summary.samples[j].gini_simpson
            end
        end

        # Permuting features permutes the values but not the totals or metrics.
        order = randperm(rng, features)
        permuted = high_precision_summary(counts[order, :];
            sample_labels=["s$j" for j in 1:samples],
            feature_labels=["p$i" for i in 1:features], policy)
        for j in eachindex(summary.samples)
            @test permuted.samples[j].total == summary.samples[j].total
            @test permuted.samples[j].richness == summary.samples[j].richness
            if !isnothing(summary.samples[j].shannon)
                # Summation order changes with the permutation, so the last bit
                # of a 256-bit sum may differ. The metric must agree to within
                # the declared tolerance, not bit-for-bit: asserting equality
                # here would test the summation order, not the statistic.
                @test _close(summary.samples[j].shannon, permuted.samples[j].shannon)
                @test _close(summary.samples[j].gini_simpson,
                             permuted.samples[j].gini_simpson)
            end
        end
    end
end

@testset "Property: uniform and single-species cases have closed forms" begin
    policy = numeric_policy(; mode=:high_precision, precision_bits=256)
    for features in (1, 2, 3, 8, 33)
        summary = high_precision_summary(fill(BigInt(7), features, 2);
            sample_labels=["u1", "u2"], feature_labels=["f$i" for i in 1:features],
            policy)
        for sample in summary.samples
            @test sample.richness == features
            @test _close(log(BigFloat(features)), sample.shannon)
            @test _close(BigFloat(features - 1) / features, sample.gini_simpson)
        end
    end
    single = high_precision_summary(BigInt[5 0; 0 9];
        sample_labels=["only", "other"], feature_labels=["f1", "f2"], policy)
    @test single.samples[1].richness == 1
    @test single.samples[1].shannon == 0
    @test single.samples[1].gini_simpson == 0
    @test single.samples[2].shannon == 0
end

@testset "Property: rational budgets bound exact representations" begin
    counts = BigInt[1 big(2)^200; big(2)^200 1]
    policy = numeric_policy(; mode=:exact_counts, max_denominator_bits=64)
    @test_throws ResourceLimitError exact_summary(counts;
        sample_labels=["a", "b"], feature_labels=["f1", "f2"], policy)
    generous = numeric_policy(; mode=:exact_counts,
        max_denominator_bits=DEFAULT_MAX_DENOMINATOR_BITS)
    summary = exact_summary(counts; sample_labels=["a", "b"],
        feature_labels=["f1", "f2"], policy=generous)
    @test summary.samples[1].proportions[1] == 1 // (big(2)^200 + 1)
    @test summary.samples[1].proportions[2] == big(2)^200 // (big(2)^200 + 1)
    @test sum(summary.samples[1].proportions) == big(1) // big(1)
end

@testset "Property: high-precision storage round-trips at achieved precision" begin
    rng = MersenneTwister(99)
    for bits in (64, 128, 256, 1024)
        counts = _prop_matrix(rng, 3, 3)
        policy = numeric_policy(; mode=:high_precision, precision_bits=bits)
        summary = high_precision_summary(counts;
            sample_labels=["a", "b", "c"], feature_labels=["f1", "f2", "f3"], policy)
        @test summary.achieved_bits == bits
        stored = high_precision_to_storage(summary)
        @test stored["achieved_precision_bits"] == bits
        for (j, entry) in pairs(stored["samples"])
            iszero(summary.samples[j].total) && continue
            for metric in ("shannon", "gini_simpson")
                value = entry[metric]
                @test value["kind"] == "approximate"
                @test value["precision_bits"] == bits
                parsed = parse_approximate_decimal(String(value["value"]); bits=bits)
                original = metric == "shannon" ? summary.samples[j].shannon :
                                                 summary.samples[j].gini_simpson
                @test parsed.value == original
                # The display is explicitly a rounded, non-authoritative view.
                @test value["display_kind"] == "rounded"
                display = String(value["display"])
                @test occursin("approximate", display)
                @test occursin("$(bits)-bit", display)
                @test occursin("rounded", display)
            end
        end
    end
end
