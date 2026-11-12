# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Property-based checks for the exact layer. The random
# matrices are seeded, so a failure is reproducible from the printed trial
# number; shrinking is done by re-running a single trial with its parameters.
using Test
using Random
using MetaManifold
using MetaManifold.NumericPolicy
using MetaManifold.ExactSummaries

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
