# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Tests against stored canonical strings independently checked using
# ECMAScript JSON.stringify and UTF-16 key order. Rust fixture generation
# provenance has not been reproduced in this audit.
using MetaManifold
using MetaManifold.Canonical: jcs_bytes, jcs_sha256, ijson_read
using JSON3, SHA, Test

@testset "Untrusted JSON nesting is bounded before parsing" begin
    @test_throws ArgumentError ijson_read(repeat("[", 65) * "0" * repeat("]", 65))
    @test ijson_read(JSON3.write(Dict("label" => repeat("[", 1000)))).label == repeat("[", 1000)
    @test_throws ArgumentError ijson_read("]")
end

@testset "canonical JSON (JCS)" begin
    # Byte-for-byte against the stored expected strings, over the fixture
    # corpus: key sorting (including the astral-before-high-BMP UTF-16
    # case), ES number formatting, string escaping, nesting, literals.
    fixtures = joinpath(@__DIR__, "..", "fixtures", "jcs_oracle.jsonl")
    for (i, line) in enumerate(eachline(fixtures))
        case = JSON3.read(line)
        @testset "oracle case $i" begin
            value = case.input
            # The input arrives parsed; the expected is the stored canonical text.
            @test String(jcs_bytes(value)) == case.expected
        end
    end

    @testset "key order is irrelevant, bytes are not" begin
        a = Dict("b" => 1, "a" => Dict("y" => [1.0, 2], "x" => "s"))
        b = Dict("a" => Dict("x" => "s", "y" => [1.0, 2]), "b" => 1)
        @test jcs_bytes(a) == jcs_bytes(b)
        # ...but array order is preserved, as JCS requires.
        c = Dict("a" => [2, 1])
        @test jcs_bytes(a) != jcs_bytes(c)
    end

    @testset "the astral key sorts before the high BMP key" begin
        # U+1F600 arrives as the surrogate pair D83D DE00; UTF-16 code-unit
        # order puts it before U+FF21. Naive code-point and UTF-8 byte order would put the BMP key
        # first instead -- this case distinguishes the UTF-16 rule.
        d = Dict("\uff21" => "high bmp", "\U0001f600" => "astral")
        @test String(jcs_bytes(d)) == "{\"😀\":\"astral\",\"Ａ\":\"high bmp\"}"
    end

    @testset "float notation follows ECMAScript" begin
        @test String(jcs_bytes([4.50, 2.4499999999999998e16, 1e30, 1e-6, 1e-7,
                                1e21, 1e20, -0.0, 5e-324])) ==
              "[4.5,24500000000000000,1e+30,0.000001,1e-7,1e+21,100000000000000000000,0,5e-324]"
    end

    @testset "refusals" begin
        @test_throws ArgumentError jcs_bytes(2^53)
        @test_throws ArgumentError jcs_bytes(-(2^53))
        @test_throws ArgumentError jcs_bytes(Inf)
        @test_throws ArgumentError jcs_bytes(NaN)
        @test jcs_bytes(2^53 - 1) == b"9007199254740991"
    end

    @testset "a content address is stable" begin
        # The use this exists for: a chart-like spec built in any key
        # order hashes identically, and the hash can be recorded where the
        # bytes cannot.
        spec = Dict(
            "layout" => Dict("title" => Dict("text" => "CO₂")),
            "data" => [Dict("x" => [1, 2, 3], "y" => [4.5, 5.5, 6.5],
                             "type" => "scatter")])
        same = Dict(
            "data" => [Dict("type" => "scatter", "y" => [4.5, 5.5, 6.5],
                            "x" => [1, 2, 3])],
            "layout" => Dict("title" => Dict("text" => "CO₂")))
        @test jcs_sha256(spec) == jcs_sha256(same)
        @test jcs_sha256(spec) == bytes2hex(sha256(jcs_bytes(spec)))
    end
end

@testset "canonical input validation regressions" begin
    @test String(jcs_bytes(Dict(:b => 2, :a => 1))) == "{\"a\":1,\"b\":2}"
    @test String(jcs_bytes(copy(JSON3.read("{\"a\":1}")))) == "{\"a\":1}"
    @test_throws ArgumentError jcs_bytes(Dict{Any,Int}(:a => 1, "a" => 2))
    @test_throws ArgumentError jcs_bytes(JSON3.read("{\"a\":1,\"a\":2}"))
    @test_throws ArgumentError jcs_bytes(Dict(1 => "not a JSON key"))
    invalid = String(UInt8[0xff])
    @test_throws ArgumentError jcs_bytes(invalid)
    @test_throws ArgumentError jcs_bytes(Dict(invalid => 1))
    @test_throws ArgumentError jcs_bytes(BigFloat(1) + eps(BigFloat))
    # Records in the statistics layer are built as named tuples; a named tuple is
    # an object with unique keys, so it must hash exactly like the equivalent Dict.
    @test jcs_bytes((b = 2, a = 1)) == jcs_bytes(Dict("a" => 1, "b" => 2))
    @test String(jcs_bytes((; a = [1, nothing], b = "s"))) == "{\"a\":[1,null],\"b\":\"s\"}"
    @test_throws ArgumentError jcs_bytes((a = BigFloat(1) + eps(BigFloat),))
    @test String(jcs_bytes([Float16(1.5), Float32(1.5), Float64(1.5)])) == "[1.5,1.5,1.5]"
    for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
        @test String(jcs_bytes(T(1))) == "1"
    end
    @test String(jcs_bytes(UInt8[0, 1, 255])) == "[0,1,255]"
end

@testset "rounding fixture is explicitly binary64 input" begin
    cases = readlines(joinpath(@__DIR__, "..", "fixtures", "jcs_oracle.jsonl"))
    @test JSON3.read(cases[5]).input[2] isa AbstractFloat
end

@testset "strict I-JSON input and JCS output (ijson-jcs)" begin
    for input in ("{\"a\":1,\"a\":2}",
                  "{\"a\":1,\"\\u0061\":2}",
                  "{\"nested\":[{\"b\":0,\"b\":1}]}",
                  "{\"n\":9007199254740992}",
                  "{\"n\":9.007199254740992e15}",
                  "{\"n\":-9.007199254740992e15}",
                  "{\"n\":2.4499999999999998e16}",
                  "[-9007199254740992]", "[1e400]", "{\"key\":\"\\ud800\"}",
                  "{\"key\":\"\\udc00\"}",
                  "{\"\\ufdd0\":0}", "[\"\\uFFFF\"]",
                  "01", "1.", "-", "[1] x", "{\"a\":1}{\"a\":2}")
        @test_throws ArgumentError ijson_read(input)
    end
    @test_throws ArgumentError ijson_read(UInt8[0xff])
    @test_throws ArgumentError jcs_bytes(Dict("bad" => "\ufdd0"))
    @test_throws ArgumentError jcs_bytes(Dict("\U0001ffff" => 1))
    for input in ("{\"a\":1,\"b\":[true,null]}",
                  "{\"emoji\":\"\\ud83d\\ude00\"}",
                  "{\"exact\":\"9007199254740992\"}",
                  "{\"float\":1e-7}")
        @test ijson_read(String(jcs_bytes(ijson_read(input)))) isa JSON3.Object
    end
    @test ijson_read("{\"n\":1e21}").n == 1e21 # JCS retains exponent notation
    @test ijson_read("{\"n\":9e15}").n == 9e15 # within exact range
    @test String(jcs_bytes(ijson_read("{\"b\":2,\"a\":1}"))) == "{\"a\":1,\"b\":2}"
    @test ijson_read(Vector{UInt8}(codeunits("{\"a\":1}"))) isa JSON3.Object
    @test ijson_read("[1]  \n") isa JSON3.Array
end
