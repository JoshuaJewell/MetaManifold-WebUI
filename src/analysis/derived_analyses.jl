# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Content-addressed, immutable derived-analysis records. Source result databases
# are never modified. Publication is atomic and an existing conflicting record is
# a data-integrity failure, not an overwrite opportunity.
module DerivedAnalyses

using JSON3
using ..Canonical: jcs_bytes, jcs_sha256, ijson_read

export DerivedAnalysisError, derived_analysis_id, persist_derived_analysis,
       load_derived_analysis, MAX_DERIVED_BYTES

const SCHEMA = "metamanifold.derived-analysis.v1"
const MAX_DERIVED_BYTES = 64 * 1024 * 1024

struct DerivedAnalysisError <: Exception
    code::String
    reason::String
end
Base.showerror(io::IO, e::DerivedAnalysisError) = print(io, e.reason)

function derived_analysis_id(source, config, result)
    jcs_sha256(Dict{String,Any}(
        "schema" => SCHEMA, "source" => source, "config" => config, "result" => result))
end

function _record(id::String, source, config, result)
    Dict{String,Any}(
        "schema" => SCHEMA,
        "analysis_id" => id,
        "source" => source,
        "config" => config,
        "result" => result,
    )
end

_get(x, key::String, default=nothing) = haskey(x, key) ? x[key] :
    haskey(x, Symbol(key)) ? x[Symbol(key)] : default

function _validate_bytes(bytes::Vector{UInt8})
    length(bytes) <= MAX_DERIVED_BYTES || throw(DerivedAnalysisError(
        "invalid_record", "derived analysis exceeds the $MAX_DERIVED_BYTES-byte read limit"))
    record = try
        ijson_read(bytes)
    catch err
        throw(DerivedAnalysisError("invalid_record",
            "derived analysis is not valid JSON: $(sprint(showerror, err))"))
    end
    (record isa JSON3.Object || record isa AbstractDict) || throw(DerivedAnalysisError(
        "invalid_record", "derived analysis root must be a JSON object"))
    jcs_bytes(record) == bytes || throw(DerivedAnalysisError(
        "invalid_record", "derived analysis is not canonical I-JSON/JCS"))
    _get(record, "schema") == SCHEMA || throw(DerivedAnalysisError(
        "invalid_record", "derived analysis has an unsupported schema"))
    id = _get(record, "analysis_id")
    id isa AbstractString && occursin(r"^[0-9a-f]{64}$", id) ||
        throw(DerivedAnalysisError("invalid_record", "derived analysis has no valid content ID"))
    source = _get(record, "source")
    config = _get(record, "config")
    result = _get(record, "result")
    (isnothing(source) || isnothing(config) || isnothing(result)) &&
        throw(DerivedAnalysisError("invalid_record", "derived analysis is incomplete"))
    all(x -> x isa JSON3.Object || x isa AbstractDict, (source, config, result)) ||
        throw(DerivedAnalysisError("invalid_record", "derived source, config and result must be objects"))
    expected = derived_analysis_id(source, config, result)
    expected == id || throw(DerivedAnalysisError("conflict",
        "derived analysis content does not match its ID; refusing a tampered or conflicting record"))
    record
end

function _existing(path::String, expected_id::String)
    filesize(path) <= MAX_DERIVED_BYTES || throw(DerivedAnalysisError(
        "conflict", "existing derived analysis exceeds the read limit; refusing to overwrite it"))
    bytes = try
        read(path)
    catch err
        throw(DerivedAnalysisError("conflict",
            "cannot read existing derived analysis: $(sprint(showerror, err))"))
    end
    record = try
        _validate_bytes(bytes)
    catch err
        err isa DerivedAnalysisError || rethrow()
        throw(DerivedAnalysisError("conflict",
            "existing derived analysis is invalid; refusing to overwrite it: $(err.reason)"))
    end
    _get(record, "analysis_id") == expected_id || throw(DerivedAnalysisError(
        "conflict", "existing derived analysis has a different content ID"))
    record
end

"""Persist one canonical record and return `(id, path, created)`.

The temp file is linked into place, and hardlink requires the destination not to
exist. Thus two concurrent callers can publish the same ID but neither can
truncate an existing record. Repeating an analysis verifies and reuses its file.
"""
function persist_derived_analysis(directory::AbstractString; source, config, result)
    id = derived_analysis_id(source, config, result)
    record = _record(id, source, config, result)
    bytes = jcs_bytes(record)
    length(bytes) <= MAX_DERIVED_BYTES || throw(DerivedAnalysisError(
        "resource_limit", "derived analysis exceeds the $MAX_DERIVED_BYTES-byte storage limit"))
    mkpath(directory)
    final = joinpath(directory, "$id.json")
    if isfile(final)
        _existing(final, id)
        return (id=id, path=final, created=false)
    end

    temp, io = mktemp(directory; cleanup=false)
    created = false
    try
        write(io, bytes)
        flush(io)
        close(io)
        try
            Base.hardlink(temp, final)
            created = true
        catch err
            err isa InterruptException && rethrow()
            # A concurrent publisher may have won the race. It is acceptable only
            # if its complete record verifies to the same ID; otherwise fail closed.
            if isfile(final)
                _existing(final, id)
            else
                rethrow()
            end
        end
    finally
        isopen(io) && close(io)
        isfile(temp) && rm(temp; force=true)
    end
    (id=id, path=final, created=created)
end

function load_derived_analysis(directory::AbstractString, id::AbstractString)
    occursin(r"^[0-9a-f]{64}$", id) || throw(DerivedAnalysisError(
        "not_found", "derived analysis ID must be 64 lowercase hexadecimal characters"))
    path = joinpath(directory, "$id.json")
    isfile(path) || throw(DerivedAnalysisError(
        "not_found", "derived analysis '$id' does not exist"))
    filesize(path) <= MAX_DERIVED_BYTES || throw(DerivedAnalysisError(
        "invalid_record", "derived analysis exceeds the read limit"))
    _existing(path, String(id))
end

end # module DerivedAnalyses
