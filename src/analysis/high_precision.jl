# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Opt-in high-precision descriptive diversity. MPFR approximations are not exact
# and no sampling-confidence or significance claim is made.
module HighPrecision

import ..ExactSummaries
import ..NumericPolicy
using ..ExactSummaries: validate_coordinates, validated_count_columns
using ..NumericPolicy: NumericPolicySpec, numeric_policy, assert_mode, exact_fraction,
    approximate_value, to_storage, to_display, policy_to_storage, policy_fingerprint,
    check_cancelled!, check_cancelled_periodically!

export HighPrecisionSummary, HighPrecisionSample, high_precision_summary,
       high_precision_to_storage, high_precision_to_display

struct HighPrecisionSample
    label::String
    total::BigInt
    richness::Int
    shannon::Union{BigFloat,Nothing}
    gini_simpson::Union{BigFloat,Nothing}
    warnings::Vector{String}
end

struct HighPrecisionSummary
    samples::Vector{HighPrecisionSample}
    policy::NumericPolicySpec
    achieved_bits::Int
    warnings::Vector{String}
    input_sha256::String
end

function _metric(column::Vector{BigInt}, policy::NumericPolicySpec; cancellation=nothing)
    check_cancelled!(cancellation)
    total = sum(column; init=big(0))
    richness = count(!iszero, column)
    if iszero(total)
        return (total=total, richness=richness, shannon=nothing, gini_simpson=nothing)
    end

    # Exact rational conversion happens once per proportion. Shannon's logarithm
    # is necessarily approximate; Gini-Simpson has an exact rational intermediate,
    # avoiding a final 1 - x cancellation in floating point.
    shannon = zero(BigFloat)
    sum_squares = big(0)
    for (i, count) in enumerate(column)
        check_cancelled_periodically!(cancellation, i)
        sum_squares += count * count
        iszero(count) && continue
        proportion = exact_fraction(count, total; policy)
        p = BigFloat(proportion)
        shannon -= p * log(p)
    end
    simpson_rational = exact_fraction(total * total - sum_squares, total * total; policy)
    (total=total, richness=richness, shannon=shannon,
     gini_simpson=BigFloat(simpson_rational))
end

"""Compute observed richness, Shannon H = -Σp log(p), and Gini-Simpson
1 - Σp² per sample at an explicit MPFR precision. Inputs must be exact integer
(or unit-rational) counts; floating-point counts are refused, never upgraded.
Zero-total samples have undefined Shannon/Simpson values (`nothing`)."""
function high_precision_summary(counts::AbstractMatrix{<:Real};
                                sample_labels::AbstractVector{<:AbstractString},
                                feature_labels::AbstractVector{<:AbstractString}=
                                    ["feature $i" for i in 1:size(counts, 1)],
                                policy::NumericPolicySpec=numeric_policy(:high_precision),
                                cancellation=nothing)
    assert_mode(policy, :high_precision)
    labels, features = validate_coordinates(counts, sample_labels, feature_labels)
    columns, _ = validated_count_columns(
        counts, labels, features; allow_approximate=false, cancellation)

    metrics, achieved = NumericPolicy.with_precision(; bits=policy.precision_bits) do
        results = Vector{Any}(undef, length(columns))
        for j in eachindex(columns)
            check_cancelled!(cancellation)
            results[j] = _metric(columns[j], policy; cancellation)
        end
        results
    end
    samples = Vector{HighPrecisionSample}(undef, length(labels))
    for j in eachindex(labels)
        check_cancelled!(cancellation)
        warnings = String[]
        if iszero(metrics[j].total)
            push!(warnings, "'$(labels[j])' has zero reads: Shannon and Gini-Simpson are undefined, not zero")
        end
        samples[j] = HighPrecisionSample(
            labels[j], metrics[j].total, metrics[j].richness,
            metrics[j].shannon, metrics[j].gini_simpson, warnings)
    end
    fingerprint = ExactSummaries._fingerprint(columns, labels, features, nothing)
    HighPrecisionSummary(samples, policy, achieved, String[], fingerprint)
end

function _metric_storage(x::BigFloat, policy::NumericPolicySpec)
    wrapped = approximate_value(x; bits=precision(x))
    Dict{String,Any}(
        "kind" => "approximate",
        "value" => to_storage(wrapped),
        "precision_bits" => precision(x),
        "display" => to_display(wrapped; digits=policy.round_digits),
        "display_kind" => "rounded",
    )
end

function high_precision_to_storage(s::HighPrecisionSummary)
    samples = [Dict{String,Any}(
        "label" => x.label,
        "status" => iszero(x.total) ? "undefined_zero_total" : "success",
        "total" => to_storage(NumericPolicy.exact_value(x.total)),
        "richness" => x.richness,
        "shannon" => isnothing(x.shannon) ? nothing : _metric_storage(x.shannon, s.policy),
        "gini_simpson" => isnothing(x.gini_simpson) ? nothing : _metric_storage(x.gini_simpson, s.policy),
        "warnings" => x.warnings,
    ) for x in s.samples]
    undefined = count(x -> iszero(x.total), s.samples)
    Dict{String,Any}(
        "schema" => "metamanifold.high-precision-diversity.v1",
        "status" => undefined == 0 ? "success" : "partial",
        "input_status" => "exact",
        "claim" => "MPFR approximations from exact counts; higher precision is not exactness; no inference is claimed",
        "numeric_policy" => policy_to_storage(s.policy),
        "policy_fingerprint" => policy_fingerprint(s.policy),
        "mode" => String(s.policy.mode),
        "achieved_precision_bits" => s.achieved_bits,
        "input_sha256" => s.input_sha256,
        "software" => Dict{String,String}(
            "julia" => string(VERSION), "mpfr" => string(Base.MPFR.version())),
        "samples" => samples,
        "warnings" => s.warnings,
        "diagnostics" => Dict{String,Int}(
            "n_samples" => length(s.samples), "n_undefined_metrics" => undefined),
    )
end

function high_precision_to_display(s::HighPrecisionSummary; digits::Integer=s.policy.round_digits)
    io = IOBuffer()
    println(io, "High-precision descriptive diversity — MPFR approximations at ",
            s.achieved_bits, " bits, not exact values. No inference is claimed.")
    for x in s.samples
        println(io, "  sample ", x.label, " — total ", x.total,
                "; observed richness ", x.richness)
        if iszero(x.total)
            println(io, "    Shannon and Gini-Simpson: undefined (zero total) — not zero")
        else
            println(io, "    Shannon: ", to_display(
                approximate_value(x.shannon; bits=precision(x.shannon)); digits))
            println(io, "    Gini-Simpson: ", to_display(
                approximate_value(x.gini_simpson; bits=precision(x.gini_simpson)); digits))
        end
        foreach(warning -> println(io, "  note: ", warning), x.warnings)
    end
    foreach(warning -> println(io, "  note: ", warning), s.warnings)
    String(take!(io))
end

end # module HighPrecision
