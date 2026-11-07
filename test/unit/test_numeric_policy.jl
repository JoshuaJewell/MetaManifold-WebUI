# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

using Test
using MetaManifold.NumericPolicy
using MetaManifold.Canonical

@testset "Numeric storage and precision metadata across float widths" begin
    for (value, bits) in ((Float16(0.5), 11), (Float32(0.5), 24), (0.5, 53))
        wrapped = approximate_value(value)
        @test wrapped.precision_bits == bits
        @test to_storage(wrapped) === 0.5
        @test encode_json_number(value) === 0.5
    end
    @test encode_json_number(42) == 42
    @test encode_json_number(big(2)^53) == "9007199254740992"
    @test encode_json_number(BigFloat(0.5)) == "0.5"
    for value in (BigFloat(Inf), BigFloat(NaN), true, "1", nothing)
        @test_throws UnsupportedRepresentationError encode_json_number(value)
    end
    @test_throws UnsupportedRepresentationError require_count(1//2)
    @test occursin("count overflow", sprint(showerror, CountOverflowError(big(1), big(2))))
end

@testset "Policy validation and ambient rounding cannot change the claim" begin
    invalid = NumericPolicySpec(:high_precision, MAX_PRECISION_BITS + 1, 4096, 6)
    @test_throws ArgumentError assert_mode(invalid, :high_precision)
    before_bits = precision(BigFloat)
    before_rounding = rounding(BigFloat)
    setrounding(BigFloat, RoundDown) do
        value, bits = with_precision(; bits=256) do
            @test rounding(BigFloat) == RoundNearest
            BigFloat(1) / BigFloat(3)
        end
        @test bits == 256
        @test precision(value) == 256
        @test rounding(BigFloat) == RoundDown
        @test_throws ErrorException with_precision(() -> error("test interruption"); bits=512)
        @test rounding(BigFloat) == RoundDown
    end
    @test precision(BigFloat) == before_bits
    @test rounding(BigFloat) == before_rounding
end

@testset "Numeric policy and claim separation" begin
    @test numeric_policy().mode == :ordinary
    @test numeric_policy(:exact_counts).precision_bits == 53
    @test numeric_policy(:high_precision).precision_bits == 256
    @test numeric_policy(mode=:high_precision, precision_bits=512).precision_bits == 512
    @test_throws ArgumentError numeric_policy(mode=:exact_counts, precision_bits=128)
    @test_throws ArgumentError numeric_policy(mode=:high_precision, precision_bits=53)
    @test_throws ArgumentError numeric_policy(mode=:ordinary, precision_bits=128)
    @test_throws ArgumentError numeric_policy(mode=:unknown)
    @test_throws ArgumentError numeric_policy(mode=:high_precision, precision_bits=1)
    @test_throws ArgumentError numeric_policy(max_denominator_bits=31)
    @test_throws ArgumentError numeric_policy(round_digits=-1)
    @test_throws ArgumentError numeric_policy(precision_bits=big(typemax(UInt128)))
    @test length(policy_fingerprint(numeric_policy(:exact_counts))) == 64
    @test policy_fingerprint(numeric_policy(:exact_counts)) !=
          policy_fingerprint(numeric_policy(:high_precision))
    @test policy_to_storage(numeric_policy(:high_precision))["precision_bits"] == 256
end

@testset "Exact counts, overflow and rational budgets" begin
    @test require_count(7) == 7
    @test require_count(UInt128(9)) == UInt128(9)
    @test require_count(7.0) == 7
    @test require_count(big(10)^30) == big(10)^30
    @test require_count(6//3) == 2
    @test_throws UnsupportedRepresentationError require_count(7.5)
    @test_throws UnsupportedRepresentationError require_count(Inf)
    @test_throws UnsupportedRepresentationError require_count(NaN)
    @test_throws UnsupportedRepresentationError require_count(true)
    @test_throws UnsupportedRepresentationError require_count(2.0^54)
    @test_throws UnsupportedRepresentationError require_count("7")

    @test checked_count_sum([1, 2, 3]) == 6
    @test checked_count_sum(Int[]) == 0
    @test checked_count_sum([typemax(Int) - 1, 1]) == typemax(Int)
    @test_throws CountOverflowError checked_count_sum([typemax(Int), 1])
    @test checked_count_sum([typemax(Int), 1]; on_overflow=:widen) == big(typemax(Int)) + 1
    @test checked_count_sum([typemax(Int), big(1)]) == big(typemax(Int)) + 1

    @test exact_fraction(1, 3) == big(1)//big(3)
    @test exact_relative_abundance(1, 3) == big(1)//big(3)
    @test isnothing(exact_relative_abundance(0, 0))
    @test_throws ArgumentError exact_relative_abundance(4, 3)
    @test_throws UnsupportedRepresentationError exact_fraction(1, 0)
    @test exact_rational_sum([big(1)//3, big(1)//6]) == big(1)//2
    @test isnothing(exact_rational_sum(Rational{BigInt}[]))
    tight = numeric_policy(mode=:exact_counts, max_denominator_bits=32)
    @test_throws ResourceLimitError exact_fraction(1, big(2)^40; policy=tight)
    @test_throws ResourceLimitError exact_fraction(big(2)^40, 1; policy=tight)
    # Check every candidate before summing: the second value has a 33-bit
    # denominator, although the reduced final sum would fit the 32-bit budget.
    @test_throws ResourceLimitError exact_rational_sum(
        [big(1)//big(2)^31, big(1)//big(2)^31 + big(1)//3]; policy=tight)
    # Independently, two reduced addends can each fit while their running total
    # cannot. The limit must therefore apply to the accumulated value as well.
    @test_throws ResourceLimitError exact_rational_sum(
        [big(1)//big(2)^31, big(1)//(big(2)^31 - 1)]; policy=tight)
end

@testset "Precision scopes and restoration" begin
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 128)
        value, achieved = with_precision(; bits=512) do
            @test precision(BigFloat) == 512
            BigFloat(1) / 3
        end
        @test achieved == 512
        @test precision(value) == 512
        @test precision(BigFloat) == 128
        nested, bits = with_precision(; bits=256) do
            inner, inner_bits = with_precision(() -> BigFloat(2); bits=384)
            @test precision(BigFloat) == 256
            @test precision(BigFloat(inner)) == 256
            inner, inner_bits
        end
        @test bits == 256
        @test precision(BigFloat) == 128
        @test_throws ArgumentError with_precision(() -> 1; bits=1)
        setprecision(BigFloat, 16)
        @test_throws UnsupportedRepresentationError require_count(2.0^53 + 2.0)
        @test require_count(2.0^53) == Int64(2)^53
    finally
        setprecision(BigFloat, old)
    end
end

@testset "Storage and display are distinct claims" begin
    huge_int = big(2)^100 + 1
    v = exact_value(huge_int)
    @test is_exact(v)
    @test to_storage(v) == string(huge_int)
    @test to_display(v) == string(huge_int)
    @test to_display(exact_value(big(1)//3); digits=10) ==
          "1/3 (exact; 10dp = 0.3333333333)"
    @test occursin("exact", to_display(exact_value(big(1)//3)))
    # Rounding is performed by integer quotient/remainder, not Float64/Printf.
    @test occursin("3dp = -0.013", to_display(exact_value(big(-25)//2000); digits=3))
    @test occursin("0dp = 1", to_display(exact_value(big(1)//2); digits=0))
    for T in (Int8, Int16, Int32, Int64, Int128)
        smallest = typemin(T)
        @test to_display(exact_value(smallest // one(T)); digits=0) ==
            "$smallest/1 (exact; 0dp = $smallest)"
        @test to_display(exact_value(smallest // T(3)); digits=6) ==
            to_display(exact_value(BigInt(smallest) // 3); digits=6)
    end
    @test endswith(to_display(exact_value(typemin(Int64) // 3); digits=6),
        "6dp = -3074457345618258602.666667)")
    @test to_storage(exact_value(3)) == 3
    @test to_storage(exact_value(big(2)^53)) == string(big(2)^53)
    @test to_storage(exact_value(big(1)//3)) == "1/3"
    @test_throws UnsupportedRepresentationError exact_value(0.5)
    @test_throws UnsupportedRepresentationError exact_value(1//0)

    approx = approximate_value(1.0/3)
    @test !is_exact(approx)
    @test approx.precision_bits == 53
    @test occursin("approximate", to_display(approx))
    rounded = rounded_value(1.0/3; digits=6)
    @test rounded.kind == ROUNDED
    @test to_display(rounded) == "0.333333 (rounded, 53-bit, rounded to 6dp)"
    @test_throws UnsupportedRepresentationError to_storage(approximate_value(NaN))
    @test_throws UnsupportedRepresentationError encode_json_number(Inf)

    huge, bits = with_precision(() -> BigFloat("1e400"); bits=256)
    @test bits == 256
    @test !occursin("Inf", to_display(approximate_value(huge); digits=3))
    @test to_storage(approximate_value(huge)) == string(huge)
    parsed = parse_approximate_decimal(to_storage(approximate_value(huge)); bits=256)
    @test parsed.value == huge
    @test parsed.precision_bits == 256
    # Every JCS record below contains only JSON-safe values.
    @test length(jcs_sha256(Dict("big" => to_storage(v),
                                "rational" => to_storage(exact_value(big(1)//3))))) == 64
end

@testset "Lossless numeric text parsing" begin
    @test parse_exact_count("42") == 42
    @test parse_exact_count("+42") == 42
    @test parse_exact_count("9007199254740993") == big"9007199254740993"
    @test parse_exact_count(string(big(2)^80)) == big(2)^80
    for bad in ("", "  ", "-1", "1.0", "1e3", "NaN", "0x10", "1/2")
        @test_throws UnsupportedRepresentationError parse_exact_count(bad)
    end
    @test parse_exact_rational("1/3") == big(1)//3
    @test parse_exact_rational("-2/6") == big(-1)//3
    @test parse_exact_rational("5") == big(5)//1
    @test parse_exact_rational(string(big(2)^80, "/", big(2)^80 + 1)) ==
          big(2)^80 // (big(2)^80 + 1)
    for bad in ("", "1/0", "0.5", "1e-3", "1/2/3", "NaN")
        @test_throws UnsupportedRepresentationError parse_exact_rational(bad)
    end
    @test_throws UnsupportedRepresentationError parse_approximate_decimal("NaN")
    @test_throws UnsupportedRepresentationError parse_approximate_decimal("1/3")
    direct = parse_approximate_decimal("0.1"; bits=512)
    via_float = with_precision(() -> BigFloat(0.1); bits=512)[1]
    @test direct.value != via_float
    @test parse_approximate_decimal(string(direct.value); bits=512).value == direct.value
end
