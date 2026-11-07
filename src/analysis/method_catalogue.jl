# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The bounded implementation catalogue of validated statistics methods. The
# review and evidence record is docs/statistics/validation-evidence.md; the
# machine-readable source of truth is config/defaults/validated_statistics.yml
# Scientific review and upstream release approval remain explicit gates.
#
# This module does not implement methods. It makes the published scope
# checkable, so a supported entry without evidence, or a limit that drifts from
# NumericPolicy, fails a test instead of silently becoming a claim.
module MethodCatalogue

using SHA, YAML
using ..Canonical: jcs_sha256
using ..NumericPolicy: MAX_CELLS, MAX_PRECISION_BITS, MAX_DENOMINATOR_BITS,
    MAX_COUNT_BITS, MAX_ROUND_DIGITS, MAX_REQUEST_BYTES

export CATALOGUE_SCHEMA, CATALOGUE_STATUSES, catalogue_path,
       load_method_catalogue, validate_catalogue, catalogue_fingerprint,
       catalogue_to_storage, supported_method_ids, method_by_id

const CATALOGUE_SCHEMA = "metamanifold.statistics-method-catalogue.v1"
const CATALOGUE_STATUSES = ("supported", "unsupported", "deferred", "out_of_scope")

# The limits the catalogue must state exactly. The catalogue is not allowed to
# advertise a larger matrix, precision or rational budget than the code enforces.
const _REQUIRED_LIMITS = Dict(
    "max_cells" => MAX_CELLS,
    "max_precision_bits" => MAX_PRECISION_BITS,
    "max_denominator_bits" => MAX_DENOMINATOR_BITS,
    "max_count_bits" => MAX_COUNT_BITS,
    "max_round_digits" => MAX_ROUND_DIGITS,
    "max_request_bytes" => MAX_REQUEST_BYTES,
)

"""Path of the shipped YAML catalogue. This is an application repository, so
the configuration file is part of the installed tree rather than a generated
artifact."""
catalogue_path() = normpath(joinpath(@__DIR__, "..", "..", "config", "defaults",
                                    "validated_statistics.yml"))

function _string_keys(value)
    if value isa AbstractDict
        result = Dict{String,Any}()
        for (key, entry) in pairs(value)
            name = String(key)
            haskey(result, name) && throw(ArgumentError(
                "catalogue repeats the key '$name'; duplicate keys are not a valid catalogue"))
            result[name] = _string_keys(entry)
        end
        return result
    elseif value isa AbstractVector
        return Any[_string_keys(entry) for entry in value]
    end
    value
end

function _require_string(mapping::AbstractDict, key::String, context::String)
    value = get(mapping, key, nothing)
    (value isa AbstractString && !isempty(strip(value))) || throw(ArgumentError(
        "$context must contain a non-empty string '$key'"))
    String(value)
end

function _require_strings(mapping::AbstractDict, key::String, context::String)
    value = get(mapping, key, nothing)
    value isa AbstractVector || throw(ArgumentError(
        "$context must contain an array '$key'"))
    all(entry -> entry isa AbstractString && !isempty(strip(entry)), value) ||
        throw(ArgumentError("$context.$key must contain only non-empty strings"))
    String[String(entry) for entry in value]
end

function _validate_method(method::Any, index::Int)
    method isa AbstractDict || throw(ArgumentError(
        "catalogue method $index must be an object"))
    context = "catalogue method $index"
    id = _require_string(method, "id", context)
    context = "catalogue method '$id'"
    status = _require_string(method, "status", context)
    status in CATALOGUE_STATUSES || throw(ArgumentError(
        "$context has unknown status '$status'; expected one of $(join(CATALOGUE_STATUSES, ", "))"))
    _require_string(method, "summary", context)
    _require_strings(method, "response_types", context)
    _require_strings(method, "designs", context)
    _require_strings(method, "supported_conditions", context)
    _require_strings(method, "unsupported_conditions", context)
    _require_string(method, "uncertainty", context)
    _require_strings(method, "diagnostics", context)
    _require_string(method, "effect_sizes", context)
    _require_string(method, "multiple_testing", context)
    _require_strings(method, "warnings", context)
    _require_strings(method, "references", context)
    evidence = _require_strings(method, "evidence", context)
    if status == "supported"
        _require_string(method, "implementation", context)
        isempty(evidence) && throw(ArgumentError(
            "$context is supported but names no evidence"))
    else
        _require_string(method, "reason", context)
    end
    method
end

"""Validate the parsed catalogue and return it unchanged.

Validation is deliberately strict. A catalogue entry with an unknown status, a
supported method without evidence, or a limit that differs from the enforcement
constants is a published claim that the implementation cannot back.
"""
function validate_catalogue(document::AbstractDict)
    schema = _require_string(document, "schema", "catalogue")
    schema == CATALOGUE_SCHEMA || throw(ArgumentError(
        "catalogue schema is '$schema'; expected '$CATALOGUE_SCHEMA'"))

    catalogue = get(document, "catalogue", nothing)
    catalogue isa AbstractDict || throw(ArgumentError("catalogue.catalogue must be an object"))
    _require_string(catalogue, "id", "catalogue.catalogue")
    status = _require_string(catalogue, "status", "catalogue.catalogue")
    status == "pending_review" || throw(ArgumentError(
        "catalogue must retain pending_review until independent review and release approval"))
    _require_string(catalogue, "review_status", "catalogue.catalogue")
    _require_string(catalogue, "review_reference", "catalogue.catalogue")

    methods = get(document, "methods", nothing)
    methods isa AbstractVector && !isempty(methods) || throw(ArgumentError(
        "catalogue.methods must be a non-empty array"))
    ids = String[]
    for (index, method) in enumerate(methods)
        _validate_method(method, index)
        push!(ids, String(method["id"]))
    end
    length(unique(ids)) == length(ids) || throw(ArgumentError(
        "catalogue method ids must be unique; got duplicates in $(join(ids, ", "))"))

    limits = get(document, "limits", nothing)
    limits isa AbstractDict || throw(ArgumentError("catalogue.limits must be an object"))
    for (key, expected) in _REQUIRED_LIMITS
        value = get(limits, key, nothing)
        (value isa Integer && !(value isa Bool)) || throw(ArgumentError(
            "catalogue.limits.$key must be an integer"))
        Int(value) == expected || throw(ArgumentError(
            "catalogue.limits.$key is $value but the code enforces $expected"))
    end
    document
end

"""Load and validate the catalogue, returning its path, parsed document and
canonical fingerprint."""
function load_method_catalogue(; path::AbstractString=catalogue_path())
    isfile(path) || throw(ArgumentError(
        "statistics method catalogue not found at '$path'; the installed tree is incomplete"))
    raw = YAML.load_file(path)
    raw isa AbstractDict || throw(ArgumentError("catalogue root must be a mapping"))
    document = _string_keys(raw)
    validate_catalogue(document)
    (path=String(path), document=document, fingerprint=catalogue_fingerprint(document))
end

function _payload(document::AbstractDict)
    Dict{String,Any}(
        "schema" => document["schema"],
        "catalogue" => document["catalogue"],
        "methods" => document["methods"],
        "limits" => document["limits"],
    )
end

"""Content fingerprint of the catalogue. It covers parsed data, not comments or
whitespace, so a reformat that preserves meaning does not change the ID."""
catalogue_fingerprint(document::AbstractDict) = jcs_sha256(_payload(document))

function catalogue_to_storage(document::AbstractDict)
    result = _payload(document)
    result["catalogue_sha256"] = catalogue_fingerprint(document)
    result
end

supported_method_ids(document::AbstractDict) =
    String[String(method["id"]) for method in document["methods"]
           if method["status"] == "supported"]

function method_by_id(document::AbstractDict, id::AbstractString)
    for method in document["methods"]
        method["id"] == id && return method
    end
    nothing
end

end # module MethodCatalogue
