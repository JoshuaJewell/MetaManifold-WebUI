# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Negative controls for the validated statistics surface. Every case here is
# expected to be refused with a specific status. A case that starts returning
# 200 means a claim in VALIDATED_STATISTICS.md stopped being true.
using Test
using JSON3
using MetaManifold
using MetaManifold.Canonical
using MetaManifold.NumericPolicy
using MetaManifold.ExactSummaries
using MetaManifold.HighPrecision

const _ADV_SV = MetaManifold.Server

_adv_text(response) = response.body isa AbstractString ? String(response.body) :
                      String(copy(response.body))

function _adv_post(path::String, body; raw::Bool=false)
    payload = raw ? body : JSON3.write(body)
    req = _ADV_SV.HTTP.Request("POST", path, ["Content-Type" => "application/json"], payload)
    _ADV_SV.Oxygen.internalrequest(req)
end

function _adv_get(path::String)
    req = _ADV_SV.HTTP.Request("GET", path, String[], "")
    _ADV_SV.Oxygen.internalrequest(req)
end

function _adv_with_root(f::Function)
    root = mktempdir()
    old = _ADV_SV.ServerState._root[]
    try
        _ADV_SV.ServerState.set_root!(root)
        mkpath(joinpath(root, "projects"))
        f(root)
    finally
        _ADV_SV.ServerState._root[] = old
        rm(root; recursive=true, force=true)
    end
end

const _ADV_BASE = Dict{String,Any}(
    "numeric_policy" => Dict{String,Any}("mode" => "exact_counts"),
    "sample_labels" => ["A", "B"],
    "feature_labels" => ["f1", "f2"],
    "counts" => [[1, 2], [3, 4]],
)

_adv_body(overrides::Dict) = merge(_ADV_BASE, overrides)
# The call sites pass pairs, not keyword arguments: `_adv_body("k" => v)` is
# one positional Pair, and without this method it is a MethodError rather than
# a request body.
_adv_body(pairs::Pair...) = merge(_ADV_BASE, Dict{String,Any}(pairs...))
_adv_body(key::String, value) = merge(_ADV_BASE, Dict(key => value))

"""Assert a refused request: documented status, documented code, canonical body."""
function _adv_refused(response, status::Int, code::AbstractString)
    @test response.status == status
    text = _adv_text(response)
    parsed = JSON3.read(text)
    @test String(parsed.error) == code
    @test haskey(parsed, :message)
    # The error body is the same I-JSON/JCS contract as a success body.
    @test String(jcs_bytes(parsed)) == text
    # A refusal is a refusal, not a stack trace or an internal detail.
    @test !occursin("Stacktrace", text)
    @test !occursin("MethodError", text)
    parsed
end

@testset "Time limits are validated and never persisted" begin
    _adv_with_root() do root
        # 3600.5 and 1e9 are float tokens; the raw lexical pass rejects every
        # non-integral JSON number in a validated-statistics request (422), so
        # they are refused before the deadline is even read.
        for bad in (0, -1, 3601, "5", true)
            response = _adv_post("/api/v1/statistics/summaries",
                _adv_body("time_limit_seconds" => bad))
            @test response.status == 400
            @test response.status != 500
        end
        # NaN and infinity are not representable in I-JSON and are refused at
        # parse time rather than becoming a permanent deadline.
        for token in ("NaN", "Infinity", "-Infinity")
            raw = replace(JSON3.write(_ADV_BASE), "\"counts\"" =>
                "\"time_limit_seconds\":$token,\"counts\"")
            @test _adv_post("/api/v1/statistics/summaries", raw; raw=true).status in (400, 422)
        end
        for token in ("3600.5", "1e2")
            raw = replace(JSON3.write(_ADV_BASE), "\"counts\"" =>
                "\"time_limit_seconds\":$token,\"counts\"")
            @test _adv_post("/api/v1/statistics/summaries", raw; raw=true).status == 422
        end

        # An explicit null is the absence of a deadline, not a zero-second one.
        @test _adv_post("/api/v1/statistics/summaries",
            _adv_body("time_limit_seconds" => nothing)).status == 200

        accepted = _adv_post("/api/v1/statistics/summaries",
            _adv_body("time_limit_seconds" => 3600))
        @test accepted.status == 200
        # The deadline is operational: it must not appear in the response or
        # in a derived record, so it cannot change a stored value.
        payload = JSON3.read(_adv_text(accepted))
        @test !haskey(payload, :time_limit_seconds)
        @test !ispath(joinpath(root, "projects", "_derived", "direct"))
    end
end

@testset "Cancellation is observable and specific" begin
    counts = BigInt[1 2; 3 4]
    labels = ["A", "B"]
    features = ["f1", "f2"]
    policy = numeric_policy(; mode=:exact_counts)

    token = CancellationToken()
    cancel!(token)
    @test_throws OperationCancelledError check_cancelled!(token)
    @test_throws OperationCancelledError check_cancelled!(token)
    @test_throws OperationCancelledError exact_summary(counts;
        sample_labels=labels, feature_labels=features, policy, cancellation=token)

    # The deadline is armed, not merely sampled: a token created with a 1 ns
    # budget is still expired after sleeping.
    expired = CancellationToken(; deadline_seconds=1e-9)
    sleep(0.01)
    @test_throws OperationCancelledError check_cancelled!(expired)
    @test_throws OperationCancelledError check_cancelled!(expired) # remains armed

    pre = CancellationToken()
    cancel!(pre)
    @test_throws OperationCancelledError exact_summary(counts;
        sample_labels=labels, feature_labels=features, policy, cancellation=pre)
    pre_high = CancellationToken()
    cancel!(pre_high)
    @test_throws OperationCancelledError high_precision_summary(counts;
        sample_labels=labels, feature_labels=features,
        policy=numeric_policy(; mode=:high_precision, precision_bits=256),
        cancellation=pre_high)

    # The route guard turns the same exception into a documented 408, not a 500.
    response = _ADV_SV._statistics_guard() do
        throw(OperationCancelledError("time limit exceeded"))
    end
    body = _adv_refused(response, 408, "statistics_cancelled")
    @test occursin("no partial result was returned", lowercase(String(body.message)))

    # Incremental checking means a token cancelled mid-loop stops the work.
    # The interval contract is independent of scheduling and is asserted here;
    # the end-to-end deadline path is proven by the hosted CI job.
    @test MetaManifold.NumericPolicy._CANCELLATION_CHECK_INTERVAL == 256
    half = CancellationToken()
    cancel!(half)
    @test_throws OperationCancelledError check_cancelled_periodically!(half, 256)
    @test isnothing(check_cancelled_periodically!(CancellationToken(), 255))
end

@testset "Numeric policy limits are enforced before work starts" begin
    _adv_with_root() do _
        cases = [
            (Dict{String,Any}("mode" => "high_precision", "precision_bits" => 1), 400),
            (Dict{String,Any}("mode" => "high_precision",
                              "precision_bits" => MAX_PRECISION_BITS + 1), 400),
            (Dict{String,Any}("mode" => "high_precision", "precision_bits" => 53), 400),
            (Dict{String,Any}("mode" => "exact_counts", "precision_bits" => 256), 400),
            (Dict{String,Any}("mode" => "exact_counts", "max_denominator_bits" => 31), 400),
            (Dict{String,Any}("mode" => "exact_counts",
                              "max_denominator_bits" => MAX_DENOMINATOR_BITS + 1), 400),
            (Dict{String,Any}("mode" => "exact_counts", "round_digits" => -1), 400),
            (Dict{String,Any}("mode" => "exact_counts",
                              "round_digits" => MAX_ROUND_DIGITS + 1), 400),
            (Dict{String,Any}("mode" => "ordinary"), 400),
            (Dict{String,Any}("mode" => "Exact_Counts"), 400),
            (Dict{String,Any}("mode" => "exact_counts", "unknown" => 1), 400),
            (Dict{String,Any}("mode" => "exact_counts", "round_digits" => 8.0), 422),
        ]
        for (policy, status) in cases
            response = _adv_post("/api/v1/statistics/summaries",
                _adv_body("numeric_policy" => policy))
            @test response.status == status
            @test response.status != 500
        end
        # Missing numeric_policy is refused: ordinary mode is never implicit.
        without = Dict("sample_labels" => ["A"], "feature_labels" => ["f"],
                       "counts" => [[1]])
        @test _adv_post("/api/v1/statistics/summaries", without).status == 400
    end
end

@testset "Matrix shape and label hygiene" begin
    _adv_with_root() do _
        cases = [
            ("sample_labels", String[]),
            ("sample_labels", ["A", "A"]),
            ("sample_labels", ["A", " "]),
            ("sample_labels", ["A", 2]),
            ("feature_labels", String[]),
            ("feature_labels", ["f1", "f1"]),
            ("feature_labels", ["f1", ""]),
            ("sample_labels", "A"),
            ("counts", [[1]]),                    # rows != features
            ("counts", [[1, 2, 3], [3, 4]]),      # ragged
            ("counts", [[1], [3, 4]]),
            ("counts", "[[1,2],[3,4]]"),
            ("groups", ["g1"]),
            ("groups", ["g1", "g2", "g3"]),
            ("groups", ["g1", 2]),
        ]
        for (key, value) in cases
            response = _adv_post("/api/v1/statistics/summaries", _adv_body(key, value))
            @test response.status == 400
            @test response.status != 500
        end

        # Booleans are not counts. `Any[Any[true, 2], Any[3, 4]]` keeps the
        # Boolean in the document: `[[true, 2], [3, 4]]` is an Int array in
        # Julia, so it would serialise as 1 and test nothing. A JSON float
        # token and a negative integer are refused for the same reason.
        for value in (true, false)
            body = _adv_body("counts" => Any[Any[value, 2], Any[3, 4]])
            status = _adv_post("/api/v1/statistics/summaries", body).status
            @test (repr(value), status) == (repr(value), 422)
        end
        for value in (1.5, -1)
            body = _adv_body("counts" => [[value, 2], [3, 4]])
            status = _adv_post("/api/v1/statistics/summaries", body).status
            @test (repr(value), status) == (repr(value), 422)
        end

        # A JSON number beyond 2^53-1 is not I-JSON-exact, so the canonical
        # reader refuses the document itself (400) before the statistics layer
        # sees a count. The same magnitude is accepted when it arrives as a
        # whole-number string, which is the documented route for large counts.
        huge_number = _adv_post("/api/v1/statistics/summaries",
            _adv_body("counts" => [[9007199254740992, 2], [3, 4]])).status
        @test huge_number == 400
        huge_string = _adv_post("/api/v1/statistics/summaries",
            _adv_body("counts" => [["9007199254740993", 2], [3, 4]])).status
        @test huge_string == 200
        # Non-integral or negative strings are refused by the statistics layer
        # (422), not by the JSON reader: they are well-formed JSON.
        for bad in ("-1", "1.5", "1e3")
            body = _adv_body("counts" => [[bad, 2], [3, 4]])
            status = _adv_post("/api/v1/statistics/summaries", body).status
            @test (repr(bad), status) == (repr(bad), 422)
        end
    end
end

@testset "Resource limits bound the work" begin
    _adv_with_root() do _
        # Exactly MAX_CELLS is accepted; one feature more is refused with 413.
        samples = ["s$i" for i in 1:250]
        features = ["f$i" for i in 1:400]
        exact = _adv_body("feature_labels" => features, "sample_labels" => samples,
            "counts" => [zeros(Int, 250) for _ in 1:400])
        @test length(features) * length(samples) == MAX_CELLS
        @test _adv_post("/api/v1/statistics/summaries", exact).status == 200

        over = _adv_body("feature_labels" => vcat(features, "f401"),
            "sample_labels" => samples,
            "counts" => [zeros(Int, 250) for _ in 1:401])
        _adv_refused(_adv_post("/api/v1/statistics/summaries", over),
            413, "numeric_resource_limit")

        # A count whose bit length exceeds the reviewed budget is refused, and
        # the textual check happens before allocating the integer.
        huge = "9" ^ (MAX_COUNT_BITS + 17)
        @test_throws ResourceLimitError parse_exact_count(huge)
        body = _adv_body("counts" => [[huge, 2], [3, 4]])
        @test _adv_post("/api/v1/statistics/summaries", body).status == 413

        # The request body cap is on bytes, not characters.
        @test_throws ResourceLimitError _ADV_SV._statistics_body(
            (; body="{\"pad\":\"" * repeat("x", MAX_REQUEST_BYTES) * "\"}"))
        multibyte = "{\"pad\":\"" * repeat("é", div(MAX_REQUEST_BYTES, 2) + 8) * "\"}"
        @test ncodeunits(multibyte) > MAX_REQUEST_BYTES
        @test length(multibyte) < MAX_REQUEST_BYTES
        @test_throws ResourceLimitError _ADV_SV._statistics_body((; body=multibyte))
    end
end

@testset "Derived identifiers cannot escape their directory" begin
    _adv_with_root() do _
        for id in ("../../etc/passwd", "..%2fetc%2fpasswd", "ABC", "0" ^ 63,
                   repeat("A", 64), repeat("0", 65), "../" * repeat("a", 64))
            response = _adv_get("/api/v1/statistics/derived/$id")
            @test response.status == 404
            @test response.status != 500
        end
        # Unknown but well-formed IDs are a 404, not an empty record.
        _adv_refused(_adv_get("/api/v1/statistics/derived/" * repeat("0", 64)),
            404, "derived_analysis_not_found")
        for study in ("..%2f..%2fetc", "studyV/../..", ".")
            @test _adv_get("/api/v1/studies/$study/runs/run1/statistics/derived/" *
                           repeat("0", 64)).status == 404
        end
    end
end

@testset "Error bodies are canonical I-JSON" begin
    _adv_with_root() do _
        for response in (
            _adv_get("/api/v1/statistics/derived/nope"),
            _adv_post("/api/v1/statistics/summaries", "{bad"; raw=true),
            _adv_post("/api/v1/statistics/summaries",
                _adv_body("numeric_policy" => Dict("mode" => "ordinary"))),
        )
            text = _adv_text(response)
            parsed = JSON3.read(text)
            @test haskey(parsed, :error)
            @test String(jcs_bytes(parsed)) == text
            @test !occursin("\n", text)
        end
        # The catalogue is a GET and is safe/idempotent; it must not mutate.
        first_response = _adv_get("/api/v1/statistics/catalogue")
        second_response = _adv_get("/api/v1/statistics/catalogue")
        @test first_response.status == 200
        @test _adv_text(first_response) == _adv_text(second_response)
    end
end
