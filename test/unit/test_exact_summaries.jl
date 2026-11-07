# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

using Test
using JSON3
using MetaManifold.NumericPolicy
using MetaManifold.ExactSummaries
using MetaManifold.Canonical

const EXACT_FIXTURE = [1 2; 1 0; 1 0]
const EXACT_LABELS = ["A", "B"]
const EXACT_FEATURES = ["f1", "f2", "f3"]

@testset "Exact sample and group proportions" begin
    s = exact_summary(EXACT_FIXTURE; sample_labels=EXACT_LABELS,
        feature_labels=EXACT_FEATURES, groups=["g1", "g2"])
    @test !s.approximate
    @test s.samples[1].total == 3
    @test s.samples[1].proportions == fill(big(1)//3, 3)
    # Sample B has [2, 0, 0] reads: a zero count in a sample that has reads is a
    # defined proportion of zero. Only a zero *total* leaves proportions
    # undefined, and this sample's total is 2.
    @test s.samples[2].proportions == Union{Rational{BigInt},Nothing}[big(1), big(0), big(0)]
    @test length(s.groups) == 2
    @test s.groups[1].label == "g1"
    @test s.groups[2].counts == BigInt[2, 0, 0]
    @test all(has_defined_proportions, s.samples)
    @test length(s.input_sha256) == 64
    @test s.input_sha256 == counts_fingerprint(EXACT_FIXTURE;
        sample_labels=EXACT_LABELS, feature_labels=EXACT_FEATURES, groups=["g1", "g2"])

    storage = summary_to_storage(s)
    @test storage["status"] == "success"
    @test storage["input_status"] == "exact"
    @test storage["samples"][1]["proportions"] == ["1/3", "1/3", "1/3"]
    @test storage["samples"][1]["total"] == 3
    @test storage["samples"][2]["proportions"] == ["1/1", "0/1", "0/1"]
    @test length(jcs_sha256(storage)) == 64
    display = summary_to_display(s; digits=8)
    @test occursin("exact", display)
    @test occursin("0.33333333", display)
    @test !occursin("significant", lowercase(display))
end

@testset "Exact arithmetic beyond Float64 safe integers" begin
    large = big(2)^80 + 17
    counts = reshape(BigInt[large, 1, large, 3], 2, 2)
    s = exact_summary(counts; sample_labels=EXACT_LABELS, feature_labels=["x", "y"])
    @test s.samples[1].total == large + 1
    @test s.samples[1].counts[1] == large
    @test s.samples[1].proportions[1] == large//(large+1)
    storage = summary_to_storage(s)
    @test storage["samples"][1]["counts"][1] == string(large)
    @test storage["samples"][1]["total"] == string(large+1)
    @test jcs_bytes(storage) == jcs_bytes(storage)

    # Metamorphic property: multiplying every count leaves all proportions fixed.
    scaled = exact_summary(counts .* big(1000); sample_labels=EXACT_LABELS,
                           feature_labels=["x", "y"])
    @test [x.proportions for x in scaled.samples] == [x.proportions for x in s.samples]
    @test scaled.input_sha256 != s.input_sha256
end

@testset "Zero totals are undefined, never zero" begin
    s = exact_summary([0 0; 0 0]; sample_labels=EXACT_LABELS, feature_labels=["x", "y"],
                      groups=["g", "g"])
    for x in vcat(s.samples, s.groups)
        @test x.total == 0
        @test all(isnothing, x.proportions)
        @test !has_defined_proportions(x)
        @test any(w -> occursin("undefined", w), x.warnings)
    end
    storage = summary_to_storage(s)
    @test storage["status"] == "partial"
    @test storage["diagnostics"]["n_undefined_proportions"] == 3
    @test storage["samples"][1]["status"] == "undefined_zero_total"
    @test storage["samples"][1]["proportions"] == [nothing, nothing]
    @test occursin("undefined (zero total) — not zero", summary_to_display(s))
end

@testset "Floating-point provenance is never erased" begin
    s = exact_summary(Float64.(EXACT_FIXTURE); sample_labels=EXACT_LABELS,
        feature_labels=EXACT_FEATURES)
    @test s.approximate
    @test all(x -> x.approximate, s.samples)
    @test all(x -> any(w -> occursin("approximations", w), x.warnings), s.samples)
    @test summary_to_storage(s)["input_status"] == "approximate"
    @test occursin("remain approximate", summary_to_display(s))
    @test_throws UnsupportedRepresentationError exact_summary(Float64.(EXACT_FIXTURE);
        sample_labels=EXACT_LABELS, feature_labels=EXACT_FEATURES, allow_approximate=false)
end

@testset "Malformed tables and resource budgets fail closed" begin
    @test_throws ArgumentError exact_summary(EXACT_FIXTURE;
        sample_labels=["A"], feature_labels=EXACT_FEATURES)
    @test_throws ArgumentError exact_summary(EXACT_FIXTURE;
        sample_labels=EXACT_LABELS, feature_labels=["x"])
    @test_throws ArgumentError exact_summary(EXACT_FIXTURE;
        sample_labels=["A", "A"], feature_labels=EXACT_FEATURES)
    @test_throws ArgumentError exact_summary(EXACT_FIXTURE;
        sample_labels=EXACT_LABELS, feature_labels=["x", "x", "z"])
    @test_throws ArgumentError exact_summary([1 -1; 2 3];
        sample_labels=EXACT_LABELS, feature_labels=["x", "y"])
    @test_throws UnsupportedRepresentationError exact_summary([1.5 2; 2 3];
        sample_labels=EXACT_LABELS, feature_labels=["x", "y"])
    @test_throws UnsupportedRepresentationError exact_summary([Inf 2; 2 3];
        sample_labels=EXACT_LABELS, feature_labels=["x", "y"])
    @test_throws ArgumentError exact_summary(EXACT_FIXTURE;
        sample_labels=EXACT_LABELS, feature_labels=EXACT_FEATURES, groups=["g"])
    @test_throws ResourceLimitError exact_summary(zeros(Int, 1, MAX_CELLS + 1);
        sample_labels=["only-one"])
    tight = numeric_policy(mode=:exact_counts, max_denominator_bits=32)
    @test_throws ResourceLimitError exact_summary(reshape(BigInt[1, big(2)^40-1], 2, 1);
        sample_labels=["huge"], feature_labels=["x", "y"], policy=tight)
end

@testset "Generative invariance and exact proportion identities" begin
    seed = 0x5eed
    state = UInt(seed)
    next_uint() = (state = state * 6364136223846793005 + 1442695040888963407)
    for n_features in (1, 2, 7), n_samples in (1, 4)
        counts = [Int(next_uint() % 20) for _ in 1:n_features, _ in 1:n_samples]
        labels = ["s$j" for j in 1:n_samples]
        features = ["f$i" for i in 1:n_features]
        s = exact_summary(counts; sample_labels=labels, feature_labels=features)
        for sample in s.samples
            if iszero(sample.total)
                @test all(isnothing, sample.proportions)
            else
                @test exact_rational_sum(sample.proportions) == 1
                @test all(p -> numerator(p) >= 0 && denominator(p) > 0, sample.proportions)
            end
        end
        # Feature permutation changes coordinates, not each feature's fraction.
        order = reverse(1:n_features)
        permuted = exact_summary(counts[order, :]; sample_labels=labels,
                                 feature_labels=features[order])
        for j in 1:n_samples
            @test permuted.samples[j].proportions == s.samples[j].proportions[order]
            @test permuted.samples[j].total == s.samples[j].total
        end
    end
end

@testset "Independent Python fractions reference" begin
    python = Sys.which("python3")
    if isnothing(python)
        @info "Skipping Independent Python fractions reference: python3 is unavailable"
    else
        payloads = [
            Dict{String,Any}(
                "sample_labels" => ["A", "B", "empty"],
                "feature_labels" => ["x", "y", "z"],
                "groups" => ["g1", "g1", "g2"],
                "counts" => [["4", "2", "0"], ["6", "4", "0"], ["0", "0", "0"]],
            ),
            Dict{String,Any}(
                "sample_labels" => ["large", "small"],
                "feature_labels" => ["x", "y"],
                "groups" => nothing,
                "counts" => [[string(big(2)^80 + 17), "1"], ["1", "3"]],
            ),
        ]
        script = joinpath(@__DIR__, "..", "references", "exact_summary_reference.py")
        strings(values) = [isnothing(value) ? nothing : string(value) for value in values]
        for payload in payloads
            labels = payload["sample_labels"]
            features = payload["feature_labels"]
            counts = [parse(BigInt, payload["counts"][i][j])
                      for i in eachindex(features), j in eachindex(labels)]
            summary = exact_summary(counts; sample_labels=labels,
                feature_labels=features, groups=payload["groups"])
            stored = summary_to_storage(summary)
            mktempdir() do directory
                input = joinpath(directory, "input.json")
                write(input, JSON3.write(payload))
                reference = JSON3.read(read(`$python $script $input`, String))
                for kind in ("samples", "groups")
                    for (i, actual) in enumerate(stored[kind])
                        expected = reference[kind][i]
                        @test actual["label"] == expected["label"]
                        @test string(actual["total"]) == expected["total"]
                        @test strings(actual["counts"]) == strings(expected["counts"])
                        @test strings(actual["proportions"]) == strings(expected["proportions"])
                    end
                end
            end
        end
    end
end

@testset "Legacy analysis remains detached from opt-in statistics" begin
    for file in ("diversity.jl", "analysis.jl", "differential.jl")
        source = lowercase(read(joinpath(@__DIR__, "..", "..", "src", "analysis", file), String))
        @test !occursin("exactsummaries", source)
        @test !occursin("highprecision", source)
        @test !occursin("numericpolicy", source)
    end
end
