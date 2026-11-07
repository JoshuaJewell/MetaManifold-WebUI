# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: AGPL-3.0-or-later

using Test
using JSON3
using MetaManifold
const VALIDATED_SV = MetaManifold.Server

function _vs_request(method::String, path::String, body; raw::Bool=false)
    payload = raw ? body : JSON3.write(body)
    req = VALIDATED_SV.HTTP.Request(method, path, ["Content-Type" => "application/json"], payload)
    VALIDATED_SV.Oxygen.internalrequest(req)
end
_vs_post(path, body; raw=false) = _vs_request("POST", path, body; raw)
_vs_get(path) = _vs_request("GET", path, "")

function _vs_with_root(f::Function)
    root = mktempdir()
    old = VALIDATED_SV.ServerState._root[]
    try
        VALIDATED_SV.ServerState.set_root!(root)
        mkpath(joinpath(root, "projects"))
        f(root)
    finally
        VALIDATED_SV.ServerState._root[] = old
        rm(root; recursive=true, force=true)
    end
end

const _VS_EXACT_BODY = Dict(
    "numeric_policy" => Dict("mode" => "exact_counts", "round_digits" => 8),
    "sample_labels" => ["A", "B"],
    "feature_labels" => ["f1", "f2", "f3"],
    "counts" => [[1, "9007199254740992"], [1, 1], [1, 1]],
    "groups" => ["g1", "g2"],
    "persist" => true,
)

@testset "Validated direct JSON statistics routes" begin
    _vs_with_root() do root
        response = _vs_post("/api/v1/statistics/summaries", _VS_EXACT_BODY)
        @test response.status == 200
        body = JSON3.read(String(response.body))
        @test body.status == "success"
        @test body.input_status == "exact"
        @test collect(body.samples[1].proportions) == ["1/3", "1/3", "1/3"]
        @test body.samples[2].total == "9007199254740994"
        @test body.samples[2].counts[1] == "9007199254740992"
        @test length(body.analysis_id) == 64
        @test body.derived_analysis.created
        id = String(body.analysis_id)

        repeat_response = _vs_post("/api/v1/statistics/summaries", _VS_EXACT_BODY)
        @test repeat_response.status == 200
        repeated = JSON3.read(String(repeat_response.body))
        @test repeated.analysis_id == id
        @test !repeated.derived_analysis.created
        loaded_response = _vs_get("/api/v1/statistics/derived/$id")
        @test loaded_response.status == 200
        loaded = JSON3.read(String(loaded_response.body))
        @test loaded.source.source_format == "direct-json"
        @test loaded.source.input.counts[1][2] == "9007199254740992"
        @test loaded.result.input_sha256 == body.input_sha256

        no_store = merge(_VS_EXACT_BODY, Dict("persist" => false))
        response = _vs_post("/api/v1/statistics/summaries", no_store)
        @test response.status == 200
        @test !haskey(JSON3.read(String(response.body)), :analysis_id)
    end
end

@testset "High-precision JSON response and zero-total state" begin
    _vs_with_root() do _
        body = Dict(
            "numeric_policy" => Dict("mode" => "high_precision", "precision_bits" => 512),
            "sample_labels" => ["uniform", "empty"],
            "feature_labels" => ["x", "y"],
            "counts" => [[1, 0], [1, 0]],
        )
        response = _vs_post("/api/v1/statistics/summaries", body)
        @test response.status == 200
        result = JSON3.read(String(response.body))
        @test result.status == "partial"
        @test result.achieved_precision_bits == 512
        @test result.samples[1].richness == 2
        @test result.samples[1].shannon.kind == "approximate"
        @test result.samples[1].shannon.precision_bits == 512
        @test result.samples[1].shannon.value isa AbstractString
        @test result.samples[2].shannon === nothing
        @test result.samples[2].status == "undefined_zero_total"
        grouped = merge(body, Dict("groups" => ["a", "b"]))
        @test _vs_post("/api/v1/statistics/summaries", grouped).status == 400
    end
end

@testset "Numeric JSON representations fail closed" begin
    _vs_with_root() do _
        for bad in (1.5, -1)
            payload = merge(_VS_EXACT_BODY, Dict("persist" => false,
                "counts" => [[bad, 1], [1, 0], [1, 0]]))
            response = _vs_post("/api/v1/statistics/summaries", payload)
            @test response.status == 422
            @test JSON3.read(String(response.body)).error == "invalid_numeric_input"
        end
        # Decimal and exponent notation are float tokens even though their
        # values are integral. These go through raw bodies: a Dict round-trip
        # cannot represent the token distinction (JSON3.write(1.0) is "1").
        for token in ("1.0", "1e0")
            raw = """{"numeric_policy":{"mode":"exact_counts"},"sample_labels":["A"],"feature_labels":["f"],"counts":[[$token]]}"""
            @test _vs_post("/api/v1/statistics/summaries", raw; raw=true).status == 422
            # String(::Vector{UInt8}) can consume the request buffer. Parsing and
            # lexical validation must use the same string for byte-backed bodies.
            @test _vs_post("/api/v1/statistics/summaries",
                Vector{UInt8}(codeunits(raw)); raw=true).status == 422
            integer_raw = replace(raw, token => "1")
            @test _vs_post("/api/v1/statistics/summaries",
                Vector{UInt8}(codeunits(integer_raw)); raw=true).status == 200
        end

        # Fractional input stays invalid as well, even when it is nearly 1.
        fractional = """{"numeric_policy":{"mode":"exact_counts"},"sample_labels":["A"],"feature_labels":["f"],"counts":[[1.0000000001]]}"""
        @test _vs_post("/api/v1/statistics/summaries", fractional; raw=true).status == 422
        # Strict I-JSON is checked on the raw body, before JSON3 can collapse
        # duplicate keys or normalise an oversized integer token. Unlike a
        # negative count, an oversized JSON number is invalid at parse time.
        for raw in ("""{"numeric_policy":{"mode":"exact_counts","mode":"high_precision"},"sample_labels":["A"],"feature_labels":["f"],"counts":[[1]]}""",
                    """{"numeric_policy":{"mode":"exact_counts"},"sample_labels":["A"],"feature_labels":["f"],"counts":[[9007199254740992]]}""")
            @test _vs_post("/api/v1/statistics/summaries", raw; raw=true).status == 400
        end
        negative_string = merge(_VS_EXACT_BODY, Dict("persist" => false,
            "counts" => [["-1", 1], [1, 0], [1, 0]]))
        @test _vs_post("/api/v1/statistics/summaries", negative_string).status == 422
        float_policy = merge(_VS_EXACT_BODY, Dict("persist" => false,
            "numeric_policy" => Dict("mode" => "exact_counts", "round_digits" => 8.0)))
        @test _vs_post("/api/v1/statistics/summaries", float_policy).status == 422
        ordinary = merge(_VS_EXACT_BODY, Dict("numeric_policy" => Dict("mode" => "ordinary")))
        @test _vs_post("/api/v1/statistics/summaries", ordinary).status == 400
        ragged = merge(_VS_EXACT_BODY, Dict("counts" => [[1], [1, 2], [1, 0]]))
        @test _vs_post("/api/v1/statistics/summaries", ragged).status == 400
    end
end

@testset "Statistics JSON parser fuzz and malformed objects" begin
    _vs_with_root() do _
        malformed = ["", "{", "[1", "{\"counts\":", "null", "[]", "1", "\"text\"",
                     "{bad}", "{\"numeric_policy\":[]}",
                     "{\"numeric_policy\":{\"mode\":\"exact_counts\"},\"counts\":{}}"]
        alphabet = collect("{}[],:\"0123456789eE+-abc ")
        state = UInt(0x1234)
        for _ in 1:80
            n = Int(state % 50)
            chars = Char[]
            for _ in 1:n
                state = state * 6364136223846793005 + 1
                push!(chars, alphabet[Int(state % length(alphabet)) + 1])
            end
            push!(malformed, String(chars))
        end
        for raw in malformed
            response = _vs_post("/api/v1/statistics/summaries", raw; raw=true)
            @test response.status in (400, 413, 422)
            @test response.status != 500
        end
    end
end

function _vs_make_run(root::String)
    data_run = joinpath(root, "data", "studyV", "run1")
    mkpath(data_run)
    touch(joinpath(data_run, "s1_R1.fastq.gz"))
    touch(joinpath(data_run, "s1_R2.fastq.gz"))
    write(joinpath(data_run, "pipeline.yml"), "analysis:\n  exclude_categories: []\n")
    merge_dir = joinpath(root, "projects", "studyV", "run1", "merged")
    mkpath(merge_dir)
    db = DuckDB.DB(joinpath(merge_dir, "results.duckdb"))
    con = DBInterface.connect(db)
    try
        DBInterface.execute(con, """CREATE TABLE merged (
            SeqName VARCHAR, Genus VARCHAR, Pident DOUBLE,
            s1 BIGINT, s2 HUGEINT, Empty BIGINT)""")
        huge = big(2)^126
        DBInterface.execute(con, """INSERT INTO merged VALUES
            ('seq2', 'G', 99.0, 2, $huge, 0),
            ('seq1', 'G', 98.0, 1, $huge, 0),
            ('seq3', 'H', 97.0, 4, 1, 0)""")
        DBInterface.execute(con, """CREATE TABLE merged_float (
            SeqName VARCHAR, Genus VARCHAR, s1 DOUBLE, s2 BIGINT)""")
        DBInterface.execute(con, "INSERT INTO merged_float VALUES ('x', 'G', 1.0, 2)")
        DBInterface.execute(con, """CREATE TABLE merged_negative (
            SeqName VARCHAR, Genus VARCHAR, s1 BIGINT)""")
        DBInterface.execute(con, "INSERT INTO merged_negative VALUES ('x', 'G', -1)")
        DBInterface.execute(con, """CREATE TABLE merged_null (
            SeqName VARCHAR, Genus VARCHAR, s1 BIGINT)""")
        DBInterface.execute(con, "INSERT INTO merged_null VALUES ('x', 'G', NULL)")
        DBInterface.execute(con, "CREATE TABLE merged_decimal (SeqName VARCHAR, Genus VARCHAR, s1 BIGINT, s2 DECIMAL(18,0))")
        DBInterface.execute(con, "INSERT INTO merged_decimal VALUES ('x', 'G', 1, 2)")
        for saved_name in ("filtered-v1", "2026.01")
            DBInterface.execute(con, "CREATE TABLE \"$(saved_name)\" AS SELECT * FROM merged")
        end
        DBInterface.execute(con, """CREATE TABLE merged_aliases (
            SeqName VARCHAR, Genus VARCHAR, feature BIGINT,
            __VALIDATED_FEATURE BIGINT, __validated_feature_ BIGINT, "odd""sample" HUGEINT)""")
        DBInterface.execute(con, """INSERT INTO merged_aliases VALUES
            ('seq2', 'G', 2, 20, 200, 9007199254740993),
            ('seq1', 'G', 1, 10, 100, 1), ('seq3', 'H', 4, 40, 400, 3)""")
        # An unbounded integer is an integer the detector cannot name, so it
        # must be refused rather than silently dropped.
        DBInterface.execute(con, "CREATE TABLE merged_bignum (SeqName VARCHAR, Genus VARCHAR, s1 BIGINT, counted BIGNUM)")
        DBInterface.execute(con, "INSERT INTO merged_bignum VALUES ('x', 'G', 1, 2)")
        DBInterface.execute(con, """CREATE TABLE merged_subgroups (
            SeqName VARCHAR, Genus VARCHAR, A_one BIGINT, A_two BIGINT, AB_one BIGINT)""")
        DBInterface.execute(con, """INSERT INTO merged_subgroups VALUES
            ('seq1', 'G', 1, 4, 100), ('seq2', 'H', 2, 5, 200)""")
        for (index, name) in enumerate(("total", "total_A", "sample_boot", "sample_dada2",
                                       "sample_vsearch", "Pident", "consensus_score", "Genus"))
            DBInterface.execute(con, "CREATE TABLE ambiguous_$index (SeqName VARCHAR, s1 BIGINT, \"$name\" BIGINT)")
            DBInterface.execute(con, "INSERT INTO ambiguous_$index VALUES ('seq1', 1, 5)")
        end
    finally
        DBInterface.close!(con)
        close(db)
    end
end

@testset "Read-only run statistics and BigInt aggregation" begin
    _vs_with_root() do root
        _vs_make_run(root)
        path = "/api/v1/studies/studyV/runs/run1/statistics/summaries"
        request = Dict("numeric_policy" => Dict("mode" => "exact_counts"),
                       "table" => "merged", "taxonomic_level" => "Genus")
        response = _vs_post(path, request)
        @test response.status == 200
        result = JSON3.read(String(response.body))
        @test result.status == "partial" # Empty has undefined proportions.
        @test collect(result.samples[1].features) == ["G", "H"]
        s1 = only(filter(x -> x.label == "s1", result.samples))
        s2 = only(filter(x -> x.label == "s2", result.samples))
        empty = only(filter(x -> x.label == "Empty", result.samples))
        @test collect(s1.counts) == [3, 4]
        @test s1.total == 7
        @test collect(s1.proportions) == ["3/7", "4/7"]
        @test s2.counts[1] == string(big(2)^127) # would overflow SQL SUM(HUGEINT)
        @test s2.total == string(big(2)^127 + 1)
        @test empty.total == 0
        @test collect(empty.proportions) == [nothing, nothing]
        id = String(result.analysis_id)
        @test isfile(joinpath(root, "projects", "studyV", "run1", "derived",
                              "validated-statistics", "$id.json"))
        @test !isfile(joinpath(root, "projects", "studyV", "run1", "merged", "$id.json"))

        # Omitted table selects the pipeline output, while saved composition
        # tables share the same safe-name contract as their creation routes.
        default_request = filter(pair -> first(pair) != "table", request)
        default_response = _vs_post(path, default_request)
        @test default_response.status == 200
        @test JSON3.read(String(default_response.body)).analysis_id == id
        for saved_name in ("filtered-v1", "2026.01")
            saved_response = _vs_post(path, merge(request, Dict("table" => saved_name)))
            @test saved_response.status == 200
            saved = JSON3.read(String(saved_response.body))
            @test saved.samples == result.samples
            @test saved.analysis_id != id # provenance identifies the source table
        end
        for invalid in ("", "merged\n", "merged;DROP TABLE merged", "../merged", 123, true, nothing, ["merged"])
            invalid_response = _vs_post(path, merge(request, Dict("table" => invalid)))
            @test invalid_response.status == 400
        end
        @test _vs_post(path, request).status == 200 # invalid names did not mutate the source
        sequence_request = filter(pair -> first(pair) != "taxonomic_level", request)
        sequence_response = _vs_post(path, merge(sequence_request, Dict("time_limit_seconds" => 30)))
        @test sequence_response.status == 200
        sequence_result = JSON3.read(String(sequence_response.body))
        @test Set(sequence_result.samples[1].features) == Set(["seq1", "seq2", "seq3"])
        @test only(filter(x -> x.label == "s1", sequence_result.samples)).total == 7

        repeated = JSON3.read(String(_vs_post(path, request).body))
        @test repeated.analysis_id == id
        @test !repeated.derived_analysis.created
        loaded_response = _vs_get("/api/v1/studies/studyV/runs/run1/statistics/derived/$id")
        @test loaded_response.status == 200
        loaded = JSON3.read(String(loaded_response.body))
        @test loaded.source.source_format == "duckdb-read-only-bigint-aggregation"
        @test loaded.config.numeric_policy.mode == "exact_counts"
        @test !haskey(loaded.result, :analysis_id) # response-only field is not in the stored result

        # Run filters retain decimal support for numeric metadata columns.
        filtered = _vs_post(path, merge(request, Dict(
            "colFilters" => Dict("Pident" => Dict("min" => 99.0)))))
        @test filtered.status == 200
        filtered_result = JSON3.read(String(filtered.body))
        @test length(filtered_result.samples) == 3
        filtered_s2 = only(filter(x -> x.label == "s2", filtered_result.samples))
        @test collect(filtered_s2.counts) == [string(big(2)^126)]

        selected = _vs_post(path, merge(request, Dict("colFilters" => Dict(
            "Genus" => Dict("include" => ["G"], "exclude" => ["H"], "text" => "G")))))
        @test selected.status == 200
        @test collect(JSON3.read(String(selected.body)).samples[1].features) == ["G"]
        for specification in (Dict("include" => "G"), Dict("exclude" => 1),
            Dict("include" => [1]), Dict("text" => true), Dict("min" => "1"),
            Dict("max" => true), Dict("min" => 2, "max" => 1), Dict("incldue" => ["G"]), 1)
            rejected = _vs_post(path, merge(request, Dict("colFilters" => Dict("Genus" => specification))))
            @test rejected.status == 400
        end
        @test _vs_post(path, merge(request, Dict("colFilters" => ["G"]))).status == 400
        @test _vs_post(path, merge(request, Dict("filter" => ["G"]))).status == 400
        @test _vs_post(path, merge(request, Dict("prefix" => true))).status == 400
        @test _vs_post(path, merge(request, Dict("taxonomic_level" => 1))).status == 400
        for values in (["G|H"], [""])
            @test _vs_post(path, merge(request, Dict("colFilters" => Dict(
                "Genus" => Dict("include" => values))))).status == 422
        end
        cleared = Dict("include" => String[], "exclude" => nothing, "text" => "", "min" => nothing)
        @test _vs_post(path, merge(request, Dict("colFilters" => Dict("Genus" => cleared)))).status == 200

        @test _vs_post(path, merge(request, Dict("table" => "merged_float"))).status == 422
        @test _vs_post(path, merge(request, Dict("table" => "merged_negative"))).status == 400
        @test _vs_post(path, merge(request, Dict("table" => "merged_null"))).status == 422
        @test _vs_post(path, merge(request, Dict("table" => "merged_decimal"))).status == 422
        @test _vs_post(path, merge(request, Dict("colFilters" => Dict("unknown" => Dict("min" => 1))))).status == 400
        @test _vs_post(path, merge(request, Dict("taxonomic_level" => "Species"))).status == 400
        @test _vs_post(path, merge(request, Dict("sampleReads" => Dict("min" => 1)))).status == 400
        @test _vs_post(path, merge(request, Dict(
            "colFilters" => Dict("s1" => Dict("min" => 1))))).status == 422
        @test _vs_get("/api/v1/statistics/derived/$id").status == 404
    end
end

@testset "Run summaries preserve colliding sample names and subgroup boundaries" begin
    _vs_with_root() do root
        _vs_make_run(root)
        path = "/api/v1/studies/studyV/runs/run1/statistics/summaries"
        for level in (nothing, "Genus")
            response = _vs_post(path, Dict("numeric_policy" => Dict("mode" => "exact_counts"),
                "table" => "merged_aliases", "taxonomic_level" => level))
            @test response.status == 200
            response.status == 200 || continue
            result = JSON3.read(String(response.body))
            @test Set(s.label for s in result.samples) ==
                Set(["feature", "__VALIDATED_FEATURE", "__validated_feature_", "odd\"sample"])
            for (label, factor) in (("feature", 1), ("__VALIDATED_FEATURE", 10), ("__validated_feature_", 100))
                sample = only(filter(s -> s.label == label, result.samples))
                @test sample.total == 7 * factor
                @test collect(sample.features) == (isnothing(level) ? ["seq1", "seq2", "seq3"] : ["G", "H"])
                @test collect(sample.counts) == factor .* (isnothing(level) ? [1, 2, 4] : [3, 4])
                @test collect(sample.proportions) == (isnothing(level) ? ["1/7", "2/7", "4/7"] : ["3/7", "4/7"])
            end
            quoted = only(filter(s -> s.label == "odd\"sample", result.samples))
            @test quoted.total == "9007199254740997"
            @test collect(quoted.counts) == (isnothing(level) ? Any[1, "9007199254740993", 3] : Any["9007199254740994", 3])
            loaded = JSON3.read(String(_vs_get("/api/v1/studies/studyV/runs/run1/statistics/derived/$(result.analysis_id)").body))
            @test loaded.result.samples == result.samples
        end

        for mode in ("exact_counts", "high_precision"), prefix in ("A", "AB", nothing)
            response = _vs_post(path, Dict("numeric_policy" => Dict("mode" => mode),
                "table" => "merged_subgroups", "prefix" => prefix))
            @test response.status == 200
            response.status == 200 || continue
            result = JSON3.read(String(response.body))
            expected = prefix == "A" ? Dict("A_one" => 3, "A_two" => 9) :
                       prefix == "AB" ? Dict("AB_one" => 300) :
                       Dict("A_one" => 3, "A_two" => 9, "AB_one" => 300)
            @test Dict(s.label => s.total for s in result.samples) == expected
            loaded = JSON3.read(String(_vs_get("/api/v1/studies/studyV/runs/run1/statistics/derived/$(result.analysis_id)").body))
            @test loaded.config.prefix == prefix
            @test Set(loaded.config.sample_columns) == Set(keys(expected))
            @test loaded.result.samples == result.samples
        end
    end
end

@testset "Reserved integer names cannot silently remove samples" begin
    _vs_with_root() do root
        _vs_make_run(root)
        path = "/api/v1/studies/studyV/runs/run1/statistics/summaries"
        for index in 1:8, mode in ("exact_counts", "high_precision")
            response = _vs_post(path, Dict("numeric_policy" => Dict("mode" => mode),
                "table" => "ambiguous_$index"))
            @test response.status == 422
            @test occursin("excluded by the sample detector", String(response.body))
        end
        directory = joinpath(root, "projects", "studyV", "run1", "derived", "validated-statistics")
        @test !isdir(directory) || isempty(readdir(directory))
        # Float metadata remains usable; rejecting ambiguous integers must not
        # prevent the existing ordinary metadata/count schema from working.
        @test _vs_post(path, Dict("numeric_policy" => Dict("mode" => "exact_counts"))).status == 200
    end
end

@testset "Unknown filter columns fail before any derived record is written" begin
    _vs_with_root() do root
        _vs_make_run(root)
        path = "/api/v1/studies/studyV/runs/run1/statistics/summaries"
        request = Dict("numeric_policy" => Dict("mode" => "exact_counts"), "table" => "merged")
        directory = joinpath(root, "projects", "studyV", "run1", "derived", "validated-statistics")
        for criterion in (Dict("text" => "G"), Dict("include" => ["G"]),
                          Dict("exclude" => ["G"]), Dict("min" => 1), Dict("max" => 1))
            response = _vs_post(path, merge(request, Dict("colFilters" => Dict("unknown" => criterion))))
            @test response.status == 400
            @test !isdir(directory) || isempty(readdir(directory))
        end
        response = _vs_post(path, merge(request, Dict("colFilters" => Dict("Genus" => Dict("text" => "H")))))
        @test response.status == 200
        result = JSON3.read(String(response.body))
        @test only(filter(s -> s.label == "s1", result.samples)).total == 4
        @test collect(result.samples[1].features) == ["seq3"]
        loaded = JSON3.read(String(_vs_get("/api/v1/studies/studyV/runs/run1/statistics/derived/$(result.analysis_id)").body))
        @test loaded.config.filters["col.Genus"] == "H"
        @test loaded.result.samples == result.samples
    end
end

@testset "Concurrent identical direct requests publish one immutable record" begin
    _vs_with_root() do _
        responses = Vector{Any}(undef, 8)
        @sync for i in eachindex(responses)
            Threads.@spawn responses[i] = _vs_post("/api/v1/statistics/summaries", _VS_EXACT_BODY)
        end
        @test all(r -> r.status == 200, responses)
        bodies = [JSON3.read(String(r.body)) for r in responses]
        @test length(unique([String(b.analysis_id) for b in bodies])) == 1
        @test count(b -> b.derived_analysis.created, bodies) == 1
    end
end

@testset "Statistics request body preserves byte buffers and counts UTF-8 bytes" begin
    @test_throws ErrorException VALIDATED_SV._statistics_guard(() -> error("unexpected failure"))
    @test_throws InterruptException VALIDATED_SV._statistics_guard(() -> throw(InterruptException()))
    _vs_with_root() do _
        escaped = merge(_VS_EXACT_BODY, Dict("sample_labels" => ["sample\"1", "sample\\2"]))
        # Escapes and digits inside labels must not be read as number tokens.
        escaped_response = _vs_post("/api/v1/statistics/summaries", escaped)
        @test escaped_response.status == 200
        @test [s.label for s in JSON3.read(String(escaped_response.body)).samples] == escaped["sample_labels"]
    end
    raw_bytes = Vector{UInt8}(codeunits("{\"counts\":[[1,2]]}"))
    parsed_bytes = VALIDATED_SV._statistics_body((; body = raw_bytes); reject_noninteger_tokens = true)
    @test collect(parsed_bytes.counts[1]) == [1, 2]
    # Julia 1.12 `String(::Vector{UInt8})` truncates in place; `_statistics_body` must copy first.
    @test length(raw_bytes) == ncodeunits("{\"counts\":[[1,2]]}")

    sub_body = SubString("  {\"label\":\"ok\"}  ", 3, 16)
    parsed_sub = VALIDATED_SV._statistics_body((; body = sub_body))
    @test String(parsed_sub.label) == "ok"

    # Multi-byte UTF-8 input whose codepoint count is under the limit while its
    # UTF-8 byte count exceeds `_STATISTICS_REQUEST_LIMIT` must fail on bytes.
    limit = VALIDATED_SV._STATISTICS_REQUEST_LIMIT
    multibyte_payload = SubString("{\"v\":\"" * repeat("é", div(limit, 2) + 8) * "\"}", 1)
    @test length(multibyte_payload) < limit
    @test ncodeunits(multibyte_payload) > limit
    @test_throws MetaManifold.NumericPolicy.ResourceLimitError VALIDATED_SV._statistics_body((; body = multibyte_payload))
end

@testset "Saved-table filter columns are quoted, not interpolated" begin
    _vs_with_root() do root
        _vs_make_run(root)
        path = "/api/v1/studies/studyV/runs/run1/statistics/summaries"
        # A column name the safe-name contract for *tables* still admits as a
        # saved table's column: `odd"sample` already exists in merged_aliases.
        request = Dict("numeric_policy" => Dict("mode" => "exact_counts"),
                       "table" => "merged_aliases", "taxonomic_level" => "Genus")
        for filters in (Dict("filter" => "G"),
                        Dict("colFilters" => Dict("Genus" => Dict("text" => "G"))),
                        Dict("colFilters" => Dict("Genus" => Dict("include" => ["G"]))),
                        Dict("colFilters" => Dict("Genus" => Dict("exclude" => ["H"]))),
                        Dict("colFilters" => Dict("Pident" => Dict("min" => 1, "max" => 99))))
            response = _vs_post(path, merge(request, filters))
            @test response.status == 200
        end
        # The quoted identifier itself also builds valid SQL. `min` on a sample
        # column is refused by the exact-integer count-filter guard, which runs
        # after the predicate is built, so reaching that refusal -- and not a
        # 500 from malformed SQL -- is the regression under test.
        quoted = _vs_post(path, merge(request,
            Dict("colFilters" => Dict("odd\"sample" => Dict("min" => 1)))))
        @test quoted.status == 422
        @test occursin("count-filter", String(quoted.body))
    end
end

@testset "An unbounded integer is never silently omitted" begin
    _vs_with_root() do root
        _vs_make_run(root)
        path = "/api/v1/studies/studyV/runs/run1/statistics/summaries"
        # `merged_bignum` has one recognised count column (s1) and one BIGNUM
        # column. Two outcomes are defensible -- refuse the ambiguous integer,
        # or include it under the existing count bounds -- but a summary that
        # reports success while omitting it is not. Assert the invariant rather
        # than one mechanism, so a later decision to support BIGNUM does not
        # turn this into a test of the wrong thing.
        complete = false
        for mode in ("exact_counts", "high_precision")
            response = _vs_post(path, Dict("numeric_policy" => Dict("mode" => mode),
                "table" => "merged_bignum"))
            body = String(response.body)
            if response.status == 200
                @test occursin("counted", body) # the column is summarised
                complete = true
            else
                @test response.status == 422
                @test occursin("counted", body) || occursin("BIGNUM", body)
            end
        end
        if !complete
            directory = joinpath(root, "projects", "studyV", "run1", "derived", "validated-statistics")
            @test !isdir(directory) || isempty(readdir(directory))
        end
    end
end

@testset "A grouped run files and retrieves its record under the resolved group" begin
    _vs_with_root() do root
        _vs_make_run(root)
        # A second run that lives inside a group directory.
        grouped = joinpath(root, "data", "studyV", "grpA", "runG")
        mkpath(grouped)
        touch(joinpath(grouped, "g1_R1.fastq.gz"))
        write(joinpath(grouped, "pipeline.yml"), "analysis:\n  exclude_categories: []\n")
        merge_dir = joinpath(root, "projects", "studyV", "grpA", "runG", "merged")
        mkpath(merge_dir)
        db = DuckDB.DB(joinpath(merge_dir, "results.duckdb"))
        con = DBInterface.connect(db)
        try
            DBInterface.execute(con, """CREATE TABLE merged (
                SeqName VARCHAR, Genus VARCHAR, g1 BIGINT, g2 BIGINT)""")
            DBInterface.execute(con, "INSERT INTO merged VALUES ('seq1', 'G', 1, 2), ('seq2', 'H', 3, 4)")
        finally
            DBInterface.close!(con)
            close(db)
        end
        request = Dict("numeric_policy" => Dict("mode" => "exact_counts"),
                       "table" => "merged", "taxonomic_level" => "Genus")
        implicit = _vs_post("/api/v1/studies/studyV/runs/runG/statistics/summaries", request)
        @test implicit.status == 200
        id = String(JSON3.read(String(implicit.body)).analysis_id)
        @test isfile(joinpath(root, "projects", "studyV", "grpA", "runG", "derived",
                              "validated-statistics", "$id.json"))
        # The record names the group it was actually read from, so the explicit
        # query form retrieves it instead of failing the source check.
        loaded = _vs_get("/api/v1/studies/studyV/runs/runG/statistics/derived/$id?group=grpA")
        @test loaded.status == 200
        @test JSON3.read(String(loaded.body)).source.run_group == "grpA"
        @test _vs_get("/api/v1/studies/studyV/runs/runG/statistics/derived/$id").status == 200
    end
end
