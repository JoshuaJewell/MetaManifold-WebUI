# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Catalogue item 1: exact descriptive counts and proportions. No inference.
module ExactSummaries

using ..Canonical: jcs_sha256
using ..NumericPolicy: NumericPolicySpec, numeric_policy, assert_mode, require_count,
    checked_count_sum, exact_relative_abundance, exact_value, to_storage, to_display,
    policy_fingerprint, policy_to_storage, UnsupportedRepresentationError,
    ResourceLimitError, MAX_CELLS, MAX_COUNT_BITS, check_cancelled!,
    check_cancelled_periodically!

export ExactSummary, SampleSummary, GroupSummary, exact_summary,
       has_defined_proportions, summary_to_storage, summary_to_display,
       counts_fingerprint

struct SampleSummary
    label::String
    features::Vector{String}
    counts::Vector{BigInt}
    total::BigInt
    proportions::Vector{Union{Rational{BigInt},Nothing}}
    approximate::Bool
    warnings::Vector{String}
end

struct GroupSummary
    label::String
    members::Vector{String}
    features::Vector{String}
    counts::Vector{BigInt}
    total::BigInt
    proportions::Vector{Union{Rational{BigInt},Nothing}}
    approximate::Bool
    warnings::Vector{String}
end

struct ExactSummary
    samples::Vector{SampleSummary}
    groups::Vector{GroupSummary}
    policy::NumericPolicySpec
    approximate::Bool
    warnings::Vector{String}
    input_sha256::String
end

has_defined_proportions(s::Union{SampleSummary,GroupSummary}) = !iszero(s.total)

function validate_coordinates(counts::AbstractMatrix, sample_labels, feature_labels)
    n_features, n_samples = size(counts)
    n_features * n_samples <= MAX_CELLS || throw(ResourceLimitError(
        "table cells", MAX_CELLS, n_features * n_samples))
    length(sample_labels) == n_samples || throw(ArgumentError(
        "expected $n_samples sample labels, got $(length(sample_labels))"))
    length(feature_labels) == n_features || throw(ArgumentError(
        "expected $n_features feature labels, got $(length(feature_labels))"))
    labels = String[String(s) for s in sample_labels]
    features = String[String(f) for f in feature_labels]
    all(s -> !isempty(strip(s)), labels) || throw(ArgumentError("sample labels must not be blank"))
    all(f -> !isempty(strip(f)), features) || throw(ArgumentError("feature labels must not be blank"))
    length(unique(labels)) == length(labels) || throw(ArgumentError("sample labels must be unique"))
    length(unique(features)) == length(features) || throw(ArgumentError("feature labels must be unique"))
    labels, features
end

function _exact_count(value, feature::AbstractString, sample::AbstractString;
                      allow_approximate::Bool)
    allow_approximate || !(value isa AbstractFloat) || throw(UnsupportedRepresentationError(
        "an integer count from an exact source", string(value),
        "feature '$feature', sample '$sample'; floating-point input is not upgraded"))
    count = try
        require_count(value)
    catch err
        (err isa UnsupportedRepresentationError || err isa ArgumentError) || rethrow()
        detail = err isa UnsupportedRepresentationError ? err.value : sprint(showerror, err)
        throw(UnsupportedRepresentationError("an exactly-known count", string(value),
            "feature '$feature', sample '$sample' ($detail)"))
    end
    count < 0 && throw(ArgumentError(
        "negative count $count at feature '$feature', sample '$sample'"))
    bits = iszero(count) ? 1 : ndigits(abs(count), base=2)
    bits <= MAX_COUNT_BITS || throw(ResourceLimitError("count bits", MAX_COUNT_BITS, bits))
    BigInt(count)
end

"""Shared, lossless count-matrix validation used by both opt-in summary modes.
Returns columns of BigInt counts and whether floating-point input was accepted."""
function validated_count_columns(counts::AbstractMatrix{<:Real}, labels::Vector{String},
                                 features::Vector{String}; allow_approximate::Bool,
                                 cancellation=nothing)
    approximate = false
    cell_index = 0
    columns = Vector{Vector{BigInt}}(undef, length(labels))
    for j in eachindex(labels)
        check_cancelled!(cancellation)
        columns[j] = Vector{BigInt}(undef, length(features))
        for i in eachindex(features)
            cell_index += 1
            check_cancelled_periodically!(cancellation, cell_index)
            value = counts[i, j]
            approximate |= value isa AbstractFloat
            columns[j][i] = _exact_count(value, features[i], labels[j]; allow_approximate)
        end
    end
    columns, approximate
end

function _summarise(label::String, counts::Vector{BigInt}, policy::NumericPolicySpec,
                    approximate::Bool, warnings::Vector{String};
                    cancellation=nothing)
    check_cancelled!(cancellation)
    total = checked_count_sum(counts; on_overflow=:widen)
    proportions = Vector{Union{Rational{BigInt},Nothing}}(undef, length(counts))
    if iszero(total)
        fill!(proportions, nothing)
        push!(warnings, "'$label' has zero reads: relative abundances are undefined, not zero")
    else
        for (i, count) in enumerate(counts)
            check_cancelled_periodically!(cancellation, i)
            proportions[i] = exact_relative_abundance(count, total; policy)
        end
    end
    total, proportions, warnings
end

function _sample_summary(label, features, counts, policy, approximate;
                         cancellation=nothing)
    warnings = String[]
    approximate && push!(warnings,
        "'$label' came from floating-point input: source values remain approximations")
    total, proportions, warnings = _summarise(
        label, counts, policy, approximate, warnings; cancellation)
    SampleSummary(label, features, counts, total, proportions, approximate, warnings)
end

function _group_summary(label, members, features, row_counts, policy, approximate;
                        cancellation=nothing)
    warnings = String[]
    approximate && push!(warnings,
        "'$label' aggregates floating-point input: source values remain approximations")
    counts = Vector{BigInt}(undef, length(row_counts))
    for (i, row) in enumerate(row_counts)
        check_cancelled_periodically!(cancellation, i)
        counts[i] = checked_count_sum(row; on_overflow=:widen)
    end
    total, proportions, warnings = _summarise(
        label, counts, policy, approximate, warnings; cancellation)
    GroupSummary(label, members, features, counts, total, proportions, approximate, warnings)
end

function _fingerprint(columns::Vector{Vector{BigInt}}, labels, features, groups)
    rows = [[string(columns[j][i]) for j in eachindex(labels)] for i in eachindex(features)]
    jcs_sha256(Dict{String,Any}(
        "schema" => "metamanifold.exact-count-matrix.v1",
        "orientation" => "features_by_samples",
        "feature_labels" => features,
        "sample_labels" => labels,
        "groups" => groups,
        "counts" => rows,
    ))
end

function counts_fingerprint(counts::AbstractMatrix{<:Integer};
                            sample_labels, feature_labels, groups=nothing)
    labels, features = validate_coordinates(counts, sample_labels, feature_labels)
    columns, _ = validated_count_columns(counts, labels, features; allow_approximate=false)
    _fingerprint(columns, labels, features,
                 isnothing(groups) ? nothing : String[String(g) for g in groups])
end

"""Exactly summarise a features-by-samples count table.

Integer and unit-denominator rational inputs are exact. Integral floats at most
2^53 are accepted for compatibility only when `allow_approximate=true`; the
result then labels the source as approximate. HTTP boundaries use `false`, so a
requested exact analysis never silently upgrades a JSON/CSV float.
"""
function exact_summary(counts::AbstractMatrix{<:Real};
                       sample_labels::AbstractVector{<:AbstractString},
                       feature_labels::AbstractVector{<:AbstractString}=
                           ["feature $i" for i in 1:size(counts, 1)],
                       groups::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
                       policy::NumericPolicySpec=numeric_policy(:exact_counts),
                       allow_approximate::Bool=true, cancellation=nothing)
    assert_mode(policy, :exact_counts)
    labels, features = validate_coordinates(counts, sample_labels, feature_labels)
    columns, approximate = validated_count_columns(
        counts, labels, features; allow_approximate, cancellation)

    group_labels = nothing
    if !isnothing(groups)
        length(groups) == length(labels) || throw(ArgumentError(
            "expected $(length(labels)) group labels, got $(length(groups))"))
        group_labels = String[String(g) for g in groups]
        all(g -> !isempty(strip(g)), group_labels) ||
            throw(ArgumentError("group labels must not be blank"))
    end

    samples = Vector{SampleSummary}(undef, length(labels))
    for j in eachindex(labels)
        check_cancelled!(cancellation)
        samples[j] = _sample_summary(labels[j], features, columns[j], policy, approximate;
                                     cancellation)
    end
    group_summaries = GroupSummary[]
    if !isnothing(group_labels)
        for label in unique(group_labels)
            check_cancelled!(cancellation)
            selected = [j for j in eachindex(labels) if group_labels[j] == label]
            members = labels[selected]
            row_counts = [BigInt[columns[j][i] for j in selected] for i in eachindex(features)]
            push!(group_summaries, _group_summary(
                label, members, features, row_counts, policy, approximate; cancellation))
        end
    end
    fingerprint = _fingerprint(columns, labels, features, group_labels)
    ExactSummary(samples, group_summaries, policy, approximate, String[], fingerprint)
end

function _counts_to_storage(x::Union{SampleSummary,GroupSummary})
    Dict{String,Any}(
        "label" => x.label,
        "status" => has_defined_proportions(x) ? "success" : "undefined_zero_total",
        "total" => to_storage(exact_value(x.total)),
        "richness" => count(!iszero, x.counts),
        "counts" => [to_storage(exact_value(c)) for c in x.counts],
        "proportions" => [isnothing(p) ? nothing : to_storage(exact_value(p))
                          for p in x.proportions],
        "proportions_defined" => has_defined_proportions(x),
        "value_kind" => x.approximate ? "exact arithmetic from approximate input" : "exact",
        "warnings" => x.warnings,
    )
end

function summary_to_storage(s::ExactSummary)
    samples = [_counts_to_storage(x) for x in s.samples]
    groups = [_counts_to_storage(x) for x in s.groups]
    for (i, x) in enumerate(s.samples)
        samples[i]["features"] = x.features
    end
    for (i, x) in enumerate(s.groups)
        groups[i]["features"] = x.features
        groups[i]["members"] = x.members
    end
    undefined = count(x -> !has_defined_proportions(x), vcat(s.samples, s.groups))
    Dict{String,Any}(
        "schema" => "metamanifold.exact-descriptive-summary.v1",
        "status" => undefined == 0 ? "success" : "partial",
        "input_status" => s.approximate ? "approximate" : "exact",
        "claim" => s.approximate ?
            "exact integer/rational arithmetic after conversion; floating-point sources remain approximate; no inference" :
            "counts and proportions are exact; no comparison, test or significance is claimed",
        "numeric_policy" => policy_to_storage(s.policy),
        "policy_fingerprint" => policy_fingerprint(s.policy), # compatibility/readability
        "mode" => String(s.policy.mode),
        "input_sha256" => s.input_sha256,
        "software" => Dict{String,String}("julia" => string(VERSION)),
        "samples" => samples,
        "groups" => groups,
        "warnings" => s.warnings,
        "diagnostics" => Dict{String,Int}(
            "n_samples" => length(s.samples),
            "n_groups" => length(s.groups),
            "n_undefined_proportions" => undefined,
        ),
    )
end

function summary_to_display(s::ExactSummary; digits::Integer=s.policy.round_digits)
    digits >= 0 || throw(ArgumentError("digits must not be negative"))
    io = IOBuffer()
    if s.approximate
        println(io, "Descriptive summary — exact arithmetic after conversion, but floating-point sources remain approximate. No inference is claimed.")
    else
        println(io, "Descriptive summary — counts and proportions are exact. No comparison, no test and no significance is claimed.")
    end
    for (heading, entries) in (("sample", s.samples), ("group", s.groups))
        for x in entries
            println(io, "  $heading ", x.label, " — total ", to_display(exact_value(x.total)))
            if !has_defined_proportions(x)
                println(io, "    proportions: undefined (zero total) — not zero")
                continue
            end
            for (feature, proportion) in zip(x.features, x.proportions)
                rendered = isnothing(proportion) ? "undefined (not observed)" :
                    to_display(exact_value(proportion); digits)
                println(io, "    ", feature, ": ", rendered)
            end
        end
    end
    notes = copy(s.warnings)
    for entry in vcat(s.samples, s.groups), warning in entry.warnings
        push!(notes, warning)
    end
    foreach(warning -> println(io, "  note: ", warning), notes)
    String(take!(io))
end

end # module ExactSummaries
