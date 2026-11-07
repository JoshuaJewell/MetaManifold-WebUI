# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Boundary audit for the descriptive-statistics port (upstream #65). Each boundary is crossed with a
# real partner where one is available: DuckDB here, Node and Rscript when the
# runner provides them. A missing partner is skipped loudly and never counted
# as evidence that the boundary is safe.
using Test
using JSON3
using DuckDB, DBInterface
using MetaManifold
using MetaManifold.Canonical
using MetaManifold.NumericPolicy
using MetaManifold.ExactSummaries

const _BOUNDARY_SV = MetaManifold.Server

function _boundary_post(path::String, body; raw::Bool=false)
    payload = raw ? body : JSON3.write(body)
    req = _BOUNDARY_SV.HTTP.Request("POST", path, ["Content-Type" => "application/json"], payload)
    _BOUNDARY_SV.Oxygen.internalrequest(req)
end

function _boundary_root(f::Function)
    root = mktempdir()
    old = _BOUNDARY_SV.ServerState._root[]
    try
        _BOUNDARY_SV.ServerState.set_root!(root)
        mkpath(joinpath(root, "projects"))
        f(root)
    finally
        _BOUNDARY_SV.ServerState._root[] = old
        rm(root; recursive=true, force=true)
    end
end

@testset "JSON integer and rational round trips" begin
    @test to_storage(exact_value(big(2)^80 + 17)) == string(big(2)^80 + 17)
    @test parse_exact_count(to_storage(exact_value(big(2)^80 + 17))) == big(2)^80 + 17
    @test parse_exact_rational(to_storage(exact_value(big(-7)//9))) == big(-7)//9
    @test parse_exact_rational(to_storage(exact_value(big(2)^80//(big(2)^80 + 1)))) ==
          big(2)^80 // (big(2)^80 + 1)
    @test parse_exact_count("9007199254740991") == 9007199254740991
    @test parse_exact_count("9007199254740993") == big(9007199254740993)
    @test to_storage(exact_value(9007199254740991)) == 9007199254740991
    @test to_storage(exact_value(big(9007199254740993))) == "9007199254740993"

    for text in ("1.0", "1e0", "9.007199254740992e15", "9007199254740993.0",
                 "-1", "+-1", " 1 2 ", "NaN", "Inf", "1/2")
        @test_throws UnsupportedRepresentationError parse_exact_count(text)
    end
    for text in ("1/0", "0/0", "1.0", "1e-3", "1/2/3", "NaN")
        @test_throws UnsupportedRepresentationError parse_exact_rational(text)
    end

    value, bits = with_precision(() -> BigFloat("0.5"); bits=512)
    @test bits == 512
    encoded = to_storage(approximate_value(value))
    parsed = parse_approximate_decimal(encoded; bits=512)
    @test parsed.value == value
    @test parsed.precision_bits == 512
    @test encoded == "0.5"
end

@testset "Direct JSON number tokens are exact or refused" begin
    _boundary_root() do _
        base = Dict(
            "numeric_policy" => Dict("mode" => "exact_counts"),
            "sample_labels" => ["A"],
            "feature_labels" => ["f"],
            "counts" => [[1]],
        )
        safe = merge(base, Dict("counts" => [[9007199254740991]]))
        response = _boundary_post("/api/v1/statistics/summaries", safe)
        @test response.status == 200
        body = JSON3.read(String(copy(response.body)))
        @test body.samples[1].counts[1] == 9007199254740991
        @test String(jcs_bytes(body)) == String(copy(response.body))

        # 2^53 is outside the exact I-JSON integer range and is rejected on the
        # raw text before JSON3 can normalise it.
        raw = """{"numeric_policy":{"mode":"exact_counts"},"sample_labels":["A"],"feature_labels":["f"],"counts":[[9007199254740992]]}"""
        @test _boundary_post("/api/v1/statistics/summaries", raw; raw=true).status == 400

        # A whole-number string is the only representation allowed above 2^53-1.
        string_body = merge(base, Dict("counts" => [["9007199254740993"]]))
        response = _boundary_post("/api/v1/statistics/summaries", string_body)
        @test response.status == 200
        body = JSON3.read(String(copy(response.body)))
        @test body.samples[1].counts[1] == "9007199254740993"
        @test body.samples[1].total == "9007199254740993"
        @test body.samples[1].proportions == ["1/1"]

        for token in ("1.0", "1e0")
            raw = """{"numeric_policy":{"mode":"exact_counts"},"sample_labels":["A"],"feature_labels":["f"],"counts":[[$token]]}"""
            @test _boundary_post("/api/v1/statistics/summaries", raw; raw=true).status == 422
        end
    end
end

@testset "DuckDB integer storage boundary" begin
    db = DuckDB.DB()
    con = DBInterface.connect(db)
    try
        DBInterface.execute(con, """CREATE TABLE integer_counts (
            SeqName VARCHAR, s1 TINYINT, s2 SMALLINT, s3 INTEGER,
            s4 BIGINT, s5 HUGEINT)""")
        DBInterface.execute(con, """INSERT INTO integer_counts VALUES
            ('seq1', 1, 2, 3, 4, $(big(2)^100))""")
        columns = String["SeqName", "s1", "s2", "s3", "s4", "s5"]
        @test _BOUNDARY_SV._validated_integer_sample_columns(con, "integer_counts", columns) ==
              ["s1", "s2", "s3", "s4", "s5"]

        DBInterface.execute(con, """CREATE TABLE double_counts (
            SeqName VARCHAR, s1 DOUBLE, s2 BIGINT)""")
        @test_throws UnsupportedRepresentationError _BOUNDARY_SV._validated_integer_sample_columns(
            con, "double_counts", String["SeqName", "s1", "s2"])

        pre_cancelled = CancellationToken()
        cancel!(pre_cancelled)
        @test_throws OperationCancelledError _BOUNDARY_SV._validated_duckdb_counts(
            con, "integer_counts", columns; cancellation=pre_cancelled)
    finally
        DBInterface.close!(con)
        close(db)
    end
end

@testset "Browser JSON.parse boundary (Node)" begin
    node = Sys.which("node")
    if isnothing(node)
        @info "Skipping browser JSON boundary: node is not on PATH"
    else
        script = """
        const lossy = JSON.parse('{"n":9007199254740993}').n;
        const exact = JSON.parse('{"n":"9007199254740993"}').n;
        const above = JSON.parse('{"n":9007199254740992}').n;
        console.log(String(lossy));
        console.log(String(exact));
        console.log(String(above));
        """
        output = split(strip(read(`$node -e $script`, String)), '\n')
        @test output[1] == "9007199254740992"   # one short, unmarked
        @test output[2] == "9007199254740993"   # string crosses byte for byte
        @test output[3] == "9007199254740992"   # 2^53 is representable, 2^53+1 is not
        @test parse_exact_count(output[2]) == big(9007199254740993)
    end
end

@testset "R double/integer boundary (Rscript)" begin
    rscript = Sys.which("Rscript")
    if isnothing(rscript)
        @info "Skipping R numeric boundary: Rscript is not on PATH"
    else
        r_code = """
        cat(sprintf("%.0f\\n", 9007199254740993))
        cat("9007199254740993\\n")
        cat(is.na(as.integer(2^31)), "\\n")
        """
        output = split(strip(read(`$rscript --vanilla -e $r_code`, String)), '\n')
        # R has no integer type wider than 32 bits and stores every numeric as a
        # double, so the literal arrives one short. The string is exact.
        @test output[1] == "9007199254740992"
        @test output[2] == "9007199254740993"
        @test output[3] == "TRUE"
    end
end

@testset "Canonical JSON round trips" begin
    payload = Dict{String,Any}(
        "schema" => "metamanifold.boundary-test.v1",
        "counts" => ["9007199254740993", 1, 2],
        "proportions" => ["1/3", nothing],
        "flags" => [true, false],
        "nested" => Dict("b" => 2, "a" => 1),
    )
    bytes = jcs_bytes(payload)
    reparsed = ijson_read(bytes)
    @test jcs_bytes(reparsed) == bytes
    @test !endswith(String(copy(bytes)), "\n")
    @test String(copy(bytes)) == """{"counts":["9007199254740993",1,2],"flags":[true,false],"nested":{"a":1,"b":2},"proportions":["1/3",null],"schema":"metamanifold.boundary-test.v1"}"""
end
