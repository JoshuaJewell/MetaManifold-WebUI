# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Unit tests for the JCS (RFC 8785) canonical JSON writer. The oracle is the
# estate's ijson-jcs crate (Rust): every fixture in
# test/fixtures/jcs_oracle.jsonl was PRODUCED by that crate, and this file
# asserts MetaManifold.Canonical reproduces it byte-for-byte. Any
# disagreement is a bug in one of the two implementations, and the fixture
# names the case.
using MetaManifold
using MetaManifold.Canonical: jcs_bytes, jcs_sha256
using JSON3, SHA, Test

@testset "canonical JSON (JCS)" begin
    # Byte-for-byte against the Rust crate's output, over the fixture
    # corpus: key sorting (including the astral-before-high-BMP UTF-16
    # case), ES number formatting, string escaping, nesting, literals.
    fixtures = joinpath(@__DIR__, "..", "fixtures", "jcs_oracle.jsonl")
    for (i, line) in enumerate(eachline(fixtures))
        case = JSON3.read(line)
        @testset "oracle case $i" begin
            value = case.input
            # The input arrives already parsed; the expected is the oracle's text.
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
        # order puts it before U+FF21. Naive code-point order would agree
        # here; naive byte order would too -- this case pins the rule.
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
