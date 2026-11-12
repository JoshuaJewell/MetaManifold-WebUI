# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Opt-in validated descriptive statistics. These routes deliberately do not use
# Analysis.filtered_counts (Float64) or SQL SUM (fixed-width overflow).
using JSON3, DataFrames, DuckDB, DBInterface
using ..NumericPolicy
using ..ExactSummaries
using ..HighPrecision
using ..DerivedAnalyses
using ..MethodCatalogue
using ..Canonical: jcs_bytes, ijson_read

# The limit is owned by NumericPolicy so the catalogue, route and tests
# cannot drift apart.
const _STATISTICS_REQUEST_LIMIT = MAX_REQUEST_BYTES
const _EXACT_JSON_INTEGER_MAX = 2^53 - 1
const _DUCKDB_INTEGER_TYPES = Set([
    "TINYINT", "SMALLINT", "INTEGER", "BIGINT", "HUGEINT",
])
# DuckDB integer types that `Analysis.sample_columns` does not admit as counts.
# A column of one of these is ambiguous rather than absent: it must be rejected
# (`_validated_integer_sample_columns`) instead of silently omitted, whether the
# omission comes from its unsignedness or from a reserved metadata name. BIGNUM
# belongs here because its width is unbounded, so it cannot be read as an exact
# 64-bit count.
const _UNSUPPORTED_INTEGER_TYPES = Set([
    "UTINYINT", "USMALLINT", "UINTEGER", "UBIGINT", "UHUGEINT", "BIGNUM",
])

# Responses are I-JSON values in RFC 8785 canonical bytes. `JSON3.write` emits
# neither key order nor ECMAScript number formatting, so it cannot be used here:
# the content hashes for this surface depend on deterministic canonical bytes.
function _statistics_json(payload; status::Int=200)
    HTTP.Response(status, ["Content-Type" => "application/json"], jcs_bytes(payload))
end

function _statistics_guard(f::Function)
    try
        f()
    catch err
        err isa InterruptException && rethrow()
        err isa OperationCancelledError && return json_error(
            408, "statistics_cancelled", sprint(showerror, err))
        err isa ResourceLimitError && return json_error(
            413, "numeric_resource_limit", sprint(showerror, err))
        err isa UnsupportedRepresentationError && return json_error(
            422, "invalid_numeric_input", sprint(showerror, err))
        err isa DerivedAnalysisError && return json_error(
            err.code == "not_found" ? 404 : err.code == "conflict" ? 409 :
            err.code == "resource_limit" ? 413 : 500,
            "derived_analysis_$(err.code)", err.reason)
        err isa ArgumentError && return json_error(
            400, "invalid_statistics_request", sprint(showerror, err))
        rethrow()
    end
end

function _statistics_fields(body, allowed)
    for key in keys(body)
        string(key) in allowed || throw(ArgumentError("unsupported statistics field '$key'; no filter or transformation was applied"))
    end
end

function _statistics_body(req; reject_noninteger_tokens::Bool=false)
    body_bytes = req.body isa AbstractString ? ncodeunits(req.body) : length(req.body)
    body_bytes <= _STATISTICS_REQUEST_LIMIT || throw(ResourceLimitError(
        "request bytes", _STATISTICS_REQUEST_LIMIT, body_bytes))
    source = req.body isa AbstractString ? String(req.body) : String(copy(req.body))
    parsed = ijson_read(source)
    parsed isa JSON3.Object || throw(ArgumentError("request body must be a JSON object"))
    reject_noninteger_tokens && _reject_noninteger_json_tokens(source)
    parsed
end

"""Reject decimal/exponent tokens whose representation JSON3 would erase.
JSON3 1.14 normalises integral `1.0` and `1e0` tokens to Int64. Every numeric
field in a direct validated-statistics request is contractually integral, so a
lexical pass over the already-validated JSON preserves that boundary."""
function _reject_noninteger_json_tokens(source::AbstractString)
    bytes = codeunits(source)
    position = 1
    while position <= length(bytes)
        byte = bytes[position]
        if byte == UInt8('"')
            position += 1
            closed = false
            while position <= length(bytes)
                byte = bytes[position]
                if byte == UInt8('\\')
                    position += 2 # JSON3 has already validated the escape.
                elseif byte == UInt8('"')
                    closed = true
                    position += 1
                    break
                else
                    position += 1
                end
            end
            closed || throw(ArgumentError("unterminated JSON string"))
        elseif 0x30 <= byte <= 0x39 || byte == UInt8('-') || byte == UInt8('+')
            token_start = position
            noninteger = false
            while position <= length(bytes)
                byte = bytes[position]
                if 0x30 <= byte <= 0x39 || byte == UInt8('+') || byte == UInt8('-') ||
                   byte == UInt8('.') || byte == UInt8('e') || byte == UInt8('E')
                    noninteger |= byte == UInt8('.') || byte == UInt8('e') || byte == UInt8('E')
                    position += 1
                else
                    break
                end
            end
            if noninteger
                token = String(bytes[token_start:position - 1])
                throw(UnsupportedRepresentationError(
                    "an integer JSON token or whole-number string", token,
                    "validated direct statistics; decimal and exponent JSON numbers are not exact input"))
            end
        else
            position += 1
        end
    end
    nothing
end

function _policy_integer(raw, key::Symbol, default)
    value = get(raw, key, default)
    (value isa Integer && !(value isa Bool)) || throw(ArgumentError(
        "numeric_policy.$key must be an integer, not $(repr(value))"))
    typemin(Int) <= value <= typemax(Int) || throw(ArgumentError(
        "numeric_policy.$key is outside the configuration integer range"))
    Int(value)
end

function _statistics_policy(body)
    raw = get(body, :numeric_policy, nothing)
    isnothing(raw) && throw(ArgumentError("numeric_policy is required; ordinary mode is not implicit"))
    (raw isa JSON3.Object || raw isa AbstractDict) || throw(ArgumentError(
        "numeric_policy must be an object"))
    mode_value = get(raw, :mode, nothing)
    mode_value isa AbstractString || throw(ArgumentError("numeric_policy.mode must be a string"))
    mode_value in ("exact_counts", "high_precision") || throw(ArgumentError(
        "validated statistics supports only exact_counts or high_precision; use the existing ordinary endpoints otherwise"))
    allowed = Set(["mode", "precision_bits", "max_denominator_bits", "round_digits"])
    for key in keys(raw)
        string(key) in allowed || throw(ArgumentError("unknown numeric_policy field '$key'"))
    end
    precision = haskey(raw, :precision_bits) ?
                _policy_integer(raw, :precision_bits, 0) : nothing
    numeric_policy(; mode=Symbol(mode_value), precision_bits=precision,
        max_denominator_bits=_policy_integer(raw, :max_denominator_bits, DEFAULT_MAX_DENOMINATOR_BITS),
        round_digits=_policy_integer(raw, :round_digits, 6))
end

function _statistics_labels(raw, name::String)
    raw isa AbstractVector || throw(ArgumentError("$name must be an array of strings"))
    all(x -> x isa AbstractString, raw) || throw(ArgumentError("$name must contain only strings"))
    String[String(x) for x in raw]
end

function _json_exact_count(value, feature::String, sample::String)
    if value isa AbstractString
        try
            return BigInt(parse_exact_count(value))
        catch err
            err isa UnsupportedRepresentationError || rethrow()
            throw(UnsupportedRepresentationError("an integer count", string(value),
                "feature '$feature', sample '$sample'; use a canonical whole-number string"))
        end
    end
    # Non-integral decimal/exponent tokens remain floats after JSON3 parsing.
    # Integral ones are rejected by the raw-token pass before this function.
    (value isa Integer && !(value isa Bool)) || throw(UnsupportedRepresentationError(
        "an integer JSON token or whole-number string", repr(value),
        "feature '$feature', sample '$sample'; floating-point JSON is not exact input"))
    0 <= value <= _EXACT_JSON_INTEGER_MAX || throw(UnsupportedRepresentationError(
        "a JSON-safe integer", repr(value),
        "feature '$feature', sample '$sample'; counts above 2^53-1 must be strings"))
    BigInt(value)
end

function _validated_filter_params(body)
    filter_text = get(body, :filter, nothing)
    (isnothing(filter_text) || filter_text isa AbstractString) ||
        throw(ArgumentError("filter must be a string"))
    filters = get(body, :colFilters, nothing)
    normalized_filters = Dict{String,Any}()
    if !isnothing(filters)
        (filters isa JSON3.Object || filters isa AbstractDict) ||
            throw(ArgumentError("colFilters must be an object"))
        for (column, specification) in pairs(filters)
            string(column) == SAMPLE_READS_FILTER_KEY && throw(ArgumentError(
                "sample read filters are not supported for validated statistics"))
            (specification isa JSON3.Object || specification isa AbstractDict) ||
                throw(ArgumentError("each column filter must be an object"))
            normalized = Dict{Symbol,Any}()
            for (field, value) in pairs(specification)
                name = string(field)
                name in ("include", "exclude", "text", "min", "max") ||
                    throw(ArgumentError("unsupported column filter field '$name'"))
                isnothing(value) && continue # cleared UI fields are inactive
                if name in ("include", "exclude")
                    (value isa AbstractVector && all(v -> v isa AbstractString, value)) ||
                        throw(ArgumentError("$name must be an array of strings"))
                    isempty(value) && continue
                    # The shared query helper transports these through a pipe-
                    # separated string. Refuse values it cannot represent intact.
                    any(v -> isempty(v) || occursin('|', v), value) &&
                        throw(UnsupportedRepresentationError("a lossless column filter", repr(value),
                            "validated statistics; empty or pipe-containing filter values are unsupported"))
                elseif name == "text"
                    value isa AbstractString || throw(ArgumentError("text filter must be a string"))
                else
                    (value isa Real && !(value isa Bool) && isfinite(value)) ||
                        throw(ArgumentError("$name filter must be a finite number"))
                end
                normalized[Symbol(name)] = value
            end
            lower, upper = get(specification, :min, nothing), get(specification, :max, nothing)
            (!isnothing(lower) && !isnothing(upper) && lower > upper) &&
                throw(ArgumentError("column filter min must not exceed max"))
            normalized_filters[string(column)] = normalized
        end
    end
    try
        _body_filter_params(Dict(:filter => filter_text, :colFilters => normalized_filters))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("invalid filter parameters: $(sprint(showerror, err))"))
    end
end

function _direct_count_input(body; cancellation=nothing)
    sample_raw = get(body, :sample_labels, nothing)
    feature_raw = get(body, :feature_labels, nothing)
    rows = get(body, :counts, nothing)
    isnothing(sample_raw) && throw(ArgumentError("sample_labels is required"))
    isnothing(feature_raw) && throw(ArgumentError("feature_labels is required"))
    rows isa AbstractVector || throw(ArgumentError("counts must be a features-by-samples array"))
    labels = _statistics_labels(sample_raw, "sample_labels")
    features = _statistics_labels(feature_raw, "feature_labels")
    isempty(labels) && throw(ArgumentError("at least one sample is required"))
    isempty(features) && throw(ArgumentError("at least one feature is required"))
    length(rows) == length(features) || throw(ArgumentError(
        "counts has $(length(rows)) rows but feature_labels has $(length(features))"))
    length(rows) * length(labels) <= MAX_CELLS || throw(ResourceLimitError(
        "table cells", MAX_CELLS, length(rows) * length(labels)))
    counts = Matrix{BigInt}(undef, length(features), length(labels))
    cell_index = 0
    for i in eachindex(features)
        row = rows[i]
        row isa AbstractVector || throw(ArgumentError("counts row $i must be an array"))
        length(row) == length(labels) || throw(ArgumentError(
            "counts row $i has $(length(row)) cells; expected $(length(labels))"))
        for j in eachindex(labels)
            cell_index += 1
            check_cancelled_periodically!(cancellation, cell_index)
            counts[i, j] = _json_exact_count(row[j], features[i], labels[j])
        end
    end
    groups_raw = get(body, :groups, nothing)
    groups = isnothing(groups_raw) ? nothing : _statistics_labels(groups_raw, "groups")
    !isnothing(groups) && length(groups) != length(labels) && throw(ArgumentError(
        "groups has $(length(groups)) labels; expected $(length(labels))"))
    # Shared checks cover blank/duplicate labels before any analysis is claimed.
    ExactSummaries.validate_coordinates(counts, labels, features)
    counts, labels, features, groups
end

function _run_statistics(counts, labels, features, groups, policy; cancellation=nothing)
    if policy.mode === :exact_counts
        summary = exact_summary(counts; sample_labels=labels, feature_labels=features,
            groups, policy, allow_approximate=false, cancellation)
        return summary_to_storage(summary)
    end
    isnothing(groups) || throw(ArgumentError(
        "groups are supported for exact_counts only; high_precision is per-sample"))
    summary = high_precision_summary(counts; sample_labels=labels,
        feature_labels=features, policy, cancellation)
    high_precision_to_storage(summary)
end

"""Optional cooperative deadline for one request.

The deadline is operational, not part of the numeric policy: it does not change
a value or a precision, and it is not persisted in a derived record. It is
bounded to one hour so a client cannot pin a worker indefinitely.
"""
function _statistics_cancellation(body)
    raw = get(body, :time_limit_seconds, nothing)
    isnothing(raw) && return nothing
    value = if raw isa Integer && !(raw isa Bool)
        Float64(raw)
    elseif raw isa AbstractFloat
        Float64(raw)
    else
        throw(ArgumentError("time_limit_seconds must be a number, not $(repr(raw))"))
    end
    isfinite(value) && value > 0 && value <= 3600 || throw(ArgumentError(
        "time_limit_seconds must be greater than 0 and at most 3600"))
    CancellationToken(; deadline_seconds=value)
end

function _validated_integer_sample_columns(con, table::String, columns::Vector{String})
    schema = DataFrame(DBInterface.execute(con,
        "SELECT column_name, data_type FROM information_schema.columns WHERE table_name = ?",
        [table]))
    types = Dict(string(row.column_name) => uppercase(string(row.data_type)) for row in eachrow(schema))
    # DECIMAL and BOOLEAN are refused outright: the first is exact but not an
    # integer count DuckDB will sum as one, and the second is not a count at
    # all. Neither may be silently dropped while claiming a complete summary.
    for (column, datatype) in types
        if startswith(datatype, "DECIMAL") || datatype == "BOOLEAN"
            throw(UnsupportedRepresentationError("an unambiguous supported count schema",
                "$column has type $datatype", "validated statistics; this column cannot be silently omitted"))
        end
    end
    samples = Analysis.sample_columns(con, table)
    # The shared detector reserves metadata names and suffixes, but FASTQ
    # sample stems can use them too. Without an authoritative column-role map
    # an excluded integer column is ambiguous: neither silently omit it nor
    # reinterpret possible metadata as reads. Direct labelled counts remain
    # available when the source schema cannot establish the sample set.
    selected = Set(samples)
    for column in columns
        datatype = get(types, column, "")
        if (datatype in _DUCKDB_INTEGER_TYPES || datatype in _UNSUPPORTED_INTEGER_TYPES) &&
                !(column in selected)
            throw(UnsupportedRepresentationError("an unambiguous supported count schema",
                "integer column '$column' ($datatype) is excluded by the sample detector",
                "validated statistics; use an explicit labelled count matrix to distinguish counts from metadata"))
        end
    end
    isempty(samples) && throw(ArgumentError("no sample count columns were found"))
    for sample in samples
        type = get(types, sample, "")
        type in _DUCKDB_INTEGER_TYPES || throw(UnsupportedRepresentationError(
            "an integer DuckDB count column", "$sample has type $type",
            "validated run statistics; FLOAT/DOUBLE/DECIMAL storage is never upgraded"))
    end
    samples
end

function _reject_float_count_filters(sample_columns::Vector{String}, params::Dict{String,String})
    for sample in sample_columns, prefix in ("col_min.", "col_max.")
        key = prefix * sample
        haskey(params, key) && throw(UnsupportedRepresentationError(
            "an exact integer count-filter predicate", "$key=$(params[key])",
            "validated run statistics; legacy col_min/col_max predicates use Float64"))
    end
    nothing
end

"""Read and aggregate an integer-safe run table. SQL performs only projection
and filtering; feature grouping and all count addition use Julia BigInt. No SQL
SUM, arbitrary grouped value, or Float64 matrix is used."""
function _validated_duckdb_counts(con, table::String, columns::Vector{String};
                                  prefix=nothing, taxonomic_level=nothing,
                                  where_clause::String="", where_params=Any[],
                                  cancellation=nothing)
    check_cancelled!(cancellation)
    samples = _validated_integer_sample_columns(con, table, columns)
    samples = _filter_by_prefix(samples, prefix)
    isempty(samples) && throw(ArgumentError("no sample columns remain after the prefix filter"))
    sample_symbols = Symbol.(samples)

    feature_expression = if isnothing(taxonomic_level)
        identity = Analysis.sequence_column_name(columns)
        if isnothing(identity) && "SeqName" in columns
            identity = "SeqName"
        end
        isnothing(identity) && throw(ArgumentError(
            "taxonomic_level is required because this table has no sequence identifier column"))
        _sql_identifier(identity)
    else
        taxon = Analysis.taxon_column(columns, taxonomic_level)
        isnothing(taxon) && throw(ArgumentError(
            "taxonomic level '$taxonomic_level' is not present in table '$table'"))
        # Missing/blank taxonomy is an explicit biological category, not a count
        # missing value. Count NULLs below are always rejected.
        "COALESCE(NULLIF(CAST($(_sql_identifier(taxon)) AS VARCHAR), ''), 'Unassigned')"
    end
    count_columns = join([_sql_identifier(c) for c in samples], ", ")
    # DuckDB identifiers are case-insensitive. Reserve a label alias outside the
    # source schema so a sample named `feature` (or our internal name) keeps its
    # own count column and cannot be mistaken for the projected feature label.
    source_names = Set(lowercase.(columns))
    feature_alias = "__validated_feature"
    while lowercase(feature_alias) in source_names
        feature_alias *= "_"
    end
    feature_symbol = Symbol(feature_alias)
    sql = "SELECT $feature_expression AS $(_sql_identifier(feature_alias)), $count_columns FROM $(_sql_identifier(table)) " *
          where_clause

    # A SQL GROUP BY would either overflow fixed-width SUM or require an
    # arbitrary per-group value. Stream the raw rows instead, and bound the
    # aggregated matrix as features arrive rather than materialising the raw
    # table in Julia. DuckDB.Tables is used through the direct DuckDB dependency.
    features = String[]
    rows_by_feature = Vector{Vector{BigInt}}()
    feature_rows = Dict{String,Int}()
    raw_rows = 0
    statement = DBInterface.prepare(con, sql, DuckDB.StreamResult)
    try
        result = DBInterface.execute(statement, where_params)
        try
            for chunk in DuckDB.Tables.partitions(result)
                chunk_columns = DuckDB.Tables.columns(chunk)
                for row in DuckDB.Tables.rows(chunk_columns)
                    raw_rows += 1
                    check_cancelled_periodically!(cancellation, raw_rows)
                    feature_value = row[feature_symbol]
                    ismissing(feature_value) && throw(UnsupportedRepresentationError(
                        "a non-missing feature label", "missing",
                        "validated run statistics row $raw_rows"))
                    feature = string(feature_value)
                    i = get(feature_rows, feature, 0)
                    if iszero(i)
                        n_features = length(features) + 1
                        n_features <= MAX_CELLS || throw(ResourceLimitError(
                            "features", MAX_CELLS, n_features))
                        n_features * length(samples) <= MAX_CELLS || throw(ResourceLimitError(
                            "table cells", MAX_CELLS, n_features * length(samples)))
                        push!(features, feature)
                        push!(rows_by_feature, zeros(BigInt, length(samples)))
                        i = length(features)
                        feature_rows[feature] = i
                    end
                    for (j, sample) in enumerate(samples)
                        value = row[sample_symbols[j]]
                        ismissing(value) && throw(UnsupportedRepresentationError(
                            "a non-missing integer count", "NULL",
                            "feature '$feature', sample '$sample'"))
                        (value isa Integer && !(value isa Bool)) || throw(UnsupportedRepresentationError(
                            "an integer count", repr(value), "feature '$feature', sample '$sample'"))
                        value >= 0 || throw(ArgumentError(
                            "negative count $value at feature '$feature', sample '$sample'"))
                        rows_by_feature[i][j] += BigInt(value)
                    end
                end
            end
        finally
            DBInterface.close!(result)
        end
    finally
        DBInterface.close!(statement)
    end

    isempty(features) && throw(ArgumentError("no rows remain after filtering"))
    order = sortperm(features)
    features = features[order]
    rows_by_feature = rows_by_feature[order]
    counts = Matrix{BigInt}(undef, length(features), length(samples))
    for i in eachindex(features), j in eachindex(samples)
        counts[i, j] = rows_by_feature[i][j]
    end
    ExactSummaries.validate_coordinates(counts, samples, features)
    counts, samples, features
end

function _direct_source(result, labels, features, groups, counts)
    # Persisted direct analyses are self-contained: the original exact input can
    # be reconstructed without trusting a browser number or an external file.
    Dict{String,Any}(
        "schema" => "metamanifold.direct-count-input.v1",
        "source_format" => "direct-json",
        "input_sha256" => result["input_sha256"],
        "input" => Dict{String,Any}(
            "orientation" => "features_by_samples",
            "sample_labels" => labels,
            "feature_labels" => features,
            "groups" => groups,
            "counts" => [[string(counts[i, j]) for j in eachindex(labels)]
                         for i in eachindex(features)],
        ),
    )
end

# A run inside a group is addressed either as `/runs/{run}` (the group is
# resolved from the directory layout) or with an explicit `?group=`. Both reach
# the same project directory, so the resolved group -- never the request's
# optional parameter -- is what a derived record is filed under and checked
# against. Recording the unresolved value would misidentify the record's source
# and make the explicit form fail to retrieve a record the implicit form
# created.
_resolved_run_group(study::String, run::String, group::Union{String,Nothing}) =
    isnothing(group) ? _run_group(study, run) : group

function _derived_run_directory(study::String, run::String; group=nothing)
    joinpath(_run_project_dir(study, run; group), "derived", "validated-statistics")
end

function _persisted_statistics_response(directory, source, config, result)
    saved = persist_derived_analysis(directory; source, config, result)
    response = copy(result)
    response["review_status"] = "independent scientific review and upstream release approval pending"
    response["analysis_id"] = saved.id
    response["derived_analysis"] = Dict{String,Any}(
        "id" => saved.id, "created" => saved.created, "immutable" => true)
    _statistics_json(response)
end

@post "/api/v1/statistics/summaries" function(req)
    _statistics_guard() do
        body = _statistics_body(req; reject_noninteger_tokens=true)
        _statistics_fields(body, ("numeric_policy", "sample_labels", "feature_labels",
                                  "counts", "groups", "persist", "time_limit_seconds"))
        policy = _statistics_policy(body)
        cancellation = _statistics_cancellation(body)
        counts, labels, features, groups = _direct_count_input(body; cancellation)
        result = _run_statistics(counts, labels, features, groups, policy; cancellation)
        persist = get(body, :persist, false)
        persist isa Bool || throw(ArgumentError("persist must be a Boolean"))
        if persist
            source = _direct_source(result, labels, features, groups, counts)
            config = Dict{String,Any}("numeric_policy" => policy_to_storage(policy))
            directory = joinpath(ServerState.projects_dir(), "_derived", "direct")
            return _persisted_statistics_response(directory, source, config, result)
        end
        _statistics_json(result)
    end
end

@post "/api/v1/studies/{study}/runs/{run}/statistics/summaries" function(req, study::String, run::String)
    _statistics_guard() do
        err = _validate_run_request(study, run)
        !isnothing(err) && return err
        group = _resolved_run_group(study, run, _req_group(req))
        body = _statistics_body(req)
        _statistics_fields(body, ("numeric_policy", "table", "prefix", "groups",
                                  "time_limit_seconds", "taxonomic_level", "rank", "colFilters", "filter"))
        policy = _statistics_policy(body)
        cancellation = _statistics_cancellation(body)
        haskey(body, :sampleReads) && throw(ArgumentError(
            "sampleReads filtering uses the legacy Float64/SQL-SUM path and is not supported for validated statistics"))
        table = get(body, :table, "merged")
        (table isa AbstractString && Validation.is_safe_name(table)) ||
            throw(ArgumentError("invalid table name"))
        table = String(table)
        prefix_raw = get(body, :prefix, nothing)
        (isnothing(prefix_raw) || prefix_raw isa AbstractString) ||
            throw(ArgumentError("prefix must be a string"))
        prefix = isnothing(prefix_raw) ? nothing : String(prefix_raw)
        level_raw = get(body, :taxonomic_level, get(body, :rank, nothing))
        (isnothing(level_raw) || level_raw isa AbstractString) ||
            throw(ArgumentError("taxonomic_level must be a string"))
        level = isnothing(level_raw) ? nothing : String(level_raw)
        params = _validated_filter_params(body)
        any(startswith(key, "sample_" ) for key in keys(params)) && throw(ArgumentError(
            "sample read filters use floating-point aggregation and are not supported"))
        groups_raw = get(body, :groups, nothing)
        groups = isnothing(groups_raw) ? nothing : _statistics_labels(groups_raw, "groups")

        outcome = _with_analysis_results_table(study, run, table; group, readonly=true) do con, columns
            for key in keys(params)
                (startswith(key, "col_") || startswith(key, "col.")) || continue
                column = split(key, '.'; limit=2)[end]
                column in columns || throw(ArgumentError("unknown filter column '$column'"))
            end
            candidate_samples = _validated_integer_sample_columns(con, table, columns)
            _reject_float_count_filters(candidate_samples, params)
            exclusions = String[]
            where_clause, where_params = _analysis_where_clause(params, columns;
                exclude_conditions=exclusions)
            counts, labels, features = _validated_duckdb_counts(con, table, columns;
                prefix, taxonomic_level=level, where_clause, where_params, cancellation)
            !isnothing(groups) && length(groups) != length(labels) && throw(ArgumentError(
                "groups has $(length(groups)) labels; expected $(length(labels)) after sample selection"))
            result = _run_statistics(counts, labels, features, groups, policy; cancellation)
            source = Dict{String,Any}(
                "schema" => "metamanifold.run-count-input.v1",
                "source_format" => "duckdb-read-only-bigint-aggregation",
                "study_id" => study, "run_id" => run, "run_group" => group,
                "table" => table, "input_sha256" => result["input_sha256"],
            )
            config = Dict{String,Any}(
                "numeric_policy" => policy_to_storage(policy),
                "taxonomic_level" => level, "prefix" => prefix,
                "groups" => groups, "filters" => params,
                "inherited_exclusion_conditions" => exclusions,
                "sample_columns" => labels,
            )
            (; source, config, result)
        end
        outcome isa HTTP.Response && return outcome
        directory = _derived_run_directory(study, run; group)
        _persisted_statistics_response(directory, outcome.source, outcome.config, outcome.result)
    end
end

function _check_derived_source(record, key::String, expected)
    source = get(record, :source, nothing)
    isnothing(source) && throw(DerivedAnalysisError("invalid_record", "derived record has no source"))
    get(source, Symbol(key), nothing) == expected || throw(DerivedAnalysisError(
        "not_found", "derived analysis does not belong to the requested source"))
end

@get "/api/v1/statistics/catalogue" function(req)
    _statistics_guard() do
        loaded = load_method_catalogue()
        payload = catalogue_to_storage(loaded.document)
        payload["supported_method_ids"] = supported_method_ids(loaded.document)
        payload["method_count"] = length(loaded.document["methods"])
        payload["file"] = "config/defaults/validated_statistics.yml"
        _statistics_json(payload)
    end
end

@get "/api/v1/statistics/derived/{id}" function(req, id::String)
    _statistics_guard() do
        directory = joinpath(ServerState.projects_dir(), "_derived", "direct")
        record = load_derived_analysis(directory, id)
        _check_derived_source(record, "source_format", "direct-json")
        _statistics_json(record)
    end
end

@get "/api/v1/studies/{study}/runs/{run}/statistics/derived/{id}" function(req, study::String, run::String, id::String)
    _statistics_guard() do
        err = _validate_run_request(study, run)
        !isnothing(err) && return err
        group = _resolved_run_group(study, run, _req_group(req))
        directory = _derived_run_directory(study, run; group)
        record = load_derived_analysis(directory, id)
        _check_derived_source(record, "study_id", study)
        _check_derived_source(record, "run_id", run)
        _check_derived_source(record, "run_group", group)
        _statistics_json(record)
    end
end
