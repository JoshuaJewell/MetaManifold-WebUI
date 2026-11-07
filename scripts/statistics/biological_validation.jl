#!/usr/bin/env julia
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

# Produce descriptive results from a counts-only CSV. The independent checker
# reads that same CSV separately; this producer does not validate itself.
using CSV, JSON3, SHA
using MetaManifold.NumericPolicy
using MetaManifold.ExactSummaries
using MetaManifold.HighPrecision
using MetaManifold.Canonical

function biological_results(input::String, provenance::String, output::String)
    isfile(input) || error("counts CSV does not exist: $input")
    isfile(provenance) || error("pipeline provenance is required")
    filesize(input) <= MAX_REQUEST_BYTES || error("counts CSV exceeds the input byte limit")
    # header=false prevents CSV.jl from silently renaming duplicate samples.
    table = collect(CSV.File(input; header=false, types=String, strict=true))
    length(table) >= 2 || error("CSV needs a header and at least one feature")
    header = collect(first(table))
    length(header) >= 2 || error("CSV needs at least one sample")
    # R's write.csv gives the row-name/feature column an empty heading. It is
    # not a sample; only the sample headings must be present and nonblank.
    all(x -> x isa AbstractString, header[2:end]) || error("missing sample CSV header")
    labels = String.(header[2:end])
    (length(table) - 1) * length(labels) <= MAX_CELLS || error("counts CSV exceeds the cell limit")
    features = String[]
    counts = Matrix{BigInt}(undef, length(table) - 1, length(labels))
    for (i, row) in enumerate(table[2:end])
        row[1] isa AbstractString || error("missing feature identifier at row $(i + 1)")
        push!(features, row[1])
        for j in eachindex(labels)
            cell = row[j + 1]
            cell isa AbstractString || error("missing count at row $(i + 1), sample $(labels[j])")
            counts[i, j] = parse_exact_count(cell)
        end
    end
    exact = exact_summary(counts; sample_labels=labels, feature_labels=features,
                          allow_approximate=false)
    high = high_precision_summary(counts; sample_labels=labels, feature_labels=features,
        policy=numeric_policy(; mode=:high_precision, precision_bits=256))
    report = Dict{String,Any}(
        "schema" => "metamanifold.descriptive-validation-input.v1",
        "csv_sha256" => bytes2hex(sha256(read(input))),
        "provenance_sha256" => bytes2hex(sha256(read(provenance))),
        "provenance" => ijson_read(read(provenance)),
        "exact" => summary_to_storage(exact),
        "high_precision" => high_precision_to_storage(high),
        "independent_check" => "not performed by this producer; run check_biological_reference.py",
        "limitations" => ["descriptive counts only", "no taxonomy validation", "no population representativeness claim",
                           "scientific review, data-use approval and release approval are separate gates"],
    )
    open(output, "w") do io
        write(io, jcs_bytes(report))
        write(io, '\n')
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 3 || error("usage: biological_validation.jl COUNTS.csv PROVENANCE.json RESULT.json")
    biological_results(ARGS...)
end
