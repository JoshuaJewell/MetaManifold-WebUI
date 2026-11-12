# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Numeric integrity for the opt-in validated statistics layer (upstream issue #65).
# Exact integers/rationals, floating-point approximations and rounded displays
# are different claims and deliberately have different representations.
module NumericPolicy

using SHA

export NUMERIC_MODES, ORDINARY_PRECISION_BITS,
       DEFAULT_MAX_DENOMINATOR_BITS, MAX_DENOMINATOR_BITS,
       MAX_COUNT_BITS, MAX_CELLS, MAX_EXACT_FLOAT, MAX_REQUEST_BYTES,
       MAX_ROUND_DIGITS,
       NumericPolicySpec, numeric_policy, validate_policy, policy_fingerprint,
       policy_to_storage, assert_mode,
       CountOverflowError, ResourceLimitError, UnsupportedRepresentationError,
       OperationCancelledError, CancellationToken, cancel!, check_cancelled!,
       check_cancelled_periodically!,
       require_count, checked_count_sum, exact_fraction, exact_relative_abundance,
       exact_rational_sum,
       EXACT, APPROXIMATE, ROUNDED, NumericValue, exact_value, approximate_value,
       rounded_value, is_exact, to_storage, to_display, encode_json_number,
       parse_exact_count, parse_exact_rational

const NUMERIC_MODES = (:ordinary, :exact_counts)
const ORDINARY_PRECISION_BITS = 53
const DEFAULT_MAX_DENOMINATOR_BITS = 4096
# Configuration limits are deliberately finite. They can be changed in source
# after review, not raised per request until an unauthenticated request can OOM us.
const MAX_DENOMINATOR_BITS = 65536
const MAX_COUNT_BITS = 4096
const MAX_CELLS = 100_000
# The validated-statistics HTTP body limit. It lives here so the catalogue and
# the route cannot disagree about what an accepted request is.
const MAX_REQUEST_BYTES = 16 * 1024 * 1024
const MAX_ROUND_DIGITS = 100
# Every integer through 2^53 is representable in Float64. Only through 2^53-1
# is safe for the I-JSON/browser contract used for JSON numbers.
const MAX_EXACT_FLOAT = big(2)^53

struct CountOverflowError <: Exception
    left::BigInt
    right::BigInt
end
Base.showerror(io::IO, e::CountOverflowError) = print(io,
    "count overflow: $(e.left) + $(e.right) leaves the fixed-width integer range; ",
    "widen deliberately with on_overflow = :widen")

struct ResourceLimitError <: Exception
    what::String
    limit::Int
    observed::Int
end
Base.showerror(io::IO, e::ResourceLimitError) = print(io,
    "numeric resource limit exhausted: $(e.what) is $(e.observed), limit $(e.limit); ",
    "no floating-point fallback was performed")

struct UnsupportedRepresentationError <: Exception
    what::String
    value::String
    target::String
end
Base.showerror(io::IO, e::UnsupportedRepresentationError) = print(io,
    "$(e.target) cannot carry $(e.value) as $(e.what) without loss or a false claim")

struct OperationCancelledError <: Exception
    reason::String
end
Base.showerror(io::IO, e::OperationCancelledError) = print(io,
    "validated statistics operation cancelled: ", e.reason,
    "; no partial result was returned")

"""Cooperative cancellation for a bounded validated-statistics request.

Cancellation is checked between cells, samples and DuckDB rows, not inside a
single arithmetic operation. A deadline is therefore an upper bound on when the
next check can observe it, not a hard real-time interrupt. The token is
deliberately small and explicit: it is never ambient process state, so
concurrent analyses cannot cancel one another.
"""
mutable struct CancellationToken
    cancelled::Threads.Atomic{Bool}
    deadline_ns::UInt64
end

function CancellationToken(; deadline_seconds::Union{Nothing,Real}=nothing)
    deadline = UInt64(0)
    if !isnothing(deadline_seconds)
        seconds = Float64(deadline_seconds)
        isfinite(seconds) && seconds > 0 || throw(ArgumentError(
            "deadline_seconds must be a finite positive number"))
        seconds <= 3600 || throw(ArgumentError(
            "deadline_seconds must not exceed 3600"))
        deadline = time_ns() + UInt64(round(Int64, seconds * 1e9))
    end
    CancellationToken(Threads.Atomic{Bool}(false), deadline)
end

function cancel!(token::CancellationToken)
    Threads.atomic_xchg!(token.cancelled, true)
    token
end

check_cancelled!(::Nothing) = nothing

function check_cancelled!(token::CancellationToken)
    # Cancellation is sticky: every later check must also refuse work.
    token.cancelled[] &&
        throw(OperationCancelledError("cancelled by caller"))
    if !iszero(token.deadline_ns) && time_ns() >= token.deadline_ns
        Threads.atomic_xchg!(token.cancelled, true)
        throw(OperationCancelledError("time limit exceeded"))
    end
    nothing
end

# Checking every cell would dominate the arithmetic for ordinary tables. 256
# cells is small enough to bound a large request and large enough to be free.
const _CANCELLATION_CHECK_INTERVAL = 256

@inline function check_cancelled_periodically!(token, counter::Int)
    if (counter % _CANCELLATION_CHECK_INTERVAL) == 0
        check_cancelled!(token)
    end
    nothing
end

"""Immutable numeric configuration for one analysis; never ambient mutable state."""
struct NumericPolicySpec
    mode::Symbol
    precision_bits::Int
    max_denominator_bits::Int
    round_digits::Int
end

function _configuration_int(value::Integer, name::String)
    typemin(Int) <= value <= typemax(Int) || throw(ArgumentError(
        "$name does not fit this platform's configuration integer range"))
    Int(value)
end

function numeric_policy(; mode::Symbol=:ordinary,
                        precision_bits::Union{Integer,Nothing}=nothing,
                        max_denominator_bits::Integer=DEFAULT_MAX_DENOMINATOR_BITS,
                        round_digits::Integer=6)
    # Both modes record 53 bits (binary64); there is no precision setting. The
    # field stays in the policy and its fingerprint, so stored policies keep
    # their identity.
    bits = isnothing(precision_bits) ? ORDINARY_PRECISION_BITS :
           _configuration_int(precision_bits, "precision_bits")
    spec = NumericPolicySpec(mode, Int(bits),
        _configuration_int(max_denominator_bits, "max_denominator_bits"),
        _configuration_int(round_digits, "round_digits"))
    validate_policy(spec)
end
numeric_policy(mode::Symbol) = numeric_policy(; mode)

function validate_policy(spec::NumericPolicySpec)
    spec.mode in NUMERIC_MODES || throw(ArgumentError(
        "mode must be one of $(join(NUMERIC_MODES, ", ")); got '$(spec.mode)'"))
    spec.mode === :ordinary && spec.precision_bits != ORDINARY_PRECISION_BITS &&
        throw(ArgumentError("ordinary mode is fixed at 53 bits"))
    spec.mode === :exact_counts && spec.precision_bits != ORDINARY_PRECISION_BITS &&
        throw(ArgumentError("exact_counts uses integers/rationals, not a floating precision setting"))
    32 <= spec.max_denominator_bits <= MAX_DENOMINATOR_BITS || throw(ArgumentError(
        "max_denominator_bits must be in 32:$MAX_DENOMINATOR_BITS"))
    0 <= spec.round_digits <= MAX_ROUND_DIGITS || throw(ArgumentError(
        "round_digits must be in 0:$MAX_ROUND_DIGITS"))
    spec
end

function assert_mode(spec::NumericPolicySpec, mode::Symbol)
    validate_policy(spec)
    spec.mode === mode || throw(ArgumentError(
        "this analysis requires numeric mode '$mode' but the policy is '$(spec.mode)'; " *
        "select the required mode explicitly or use the existing ordinary workflow"))
    spec
end

# Named fields and validated integral values make this byte representation stable;
# unlike a Dict's iteration order it does not depend on hash-table history.
function policy_fingerprint(spec::NumericPolicySpec)
    canonical = join(("metamanifold.numeric-policy.v1",
                      "max_denominator_bits=$(spec.max_denominator_bits)",
                      "mode=$(spec.mode)",
                      "precision_bits=$(spec.precision_bits)",
                      "round_digits=$(spec.round_digits)"), '\n')
    bytes2hex(sha256(codeunits(canonical)))
end

policy_to_storage(spec::NumericPolicySpec) = Dict{String,Any}(
    "schema_version" => 1,
    "mode" => String(spec.mode),
    "precision_bits" => spec.precision_bits,
    "max_denominator_bits" => spec.max_denominator_bits,
    "round_digits" => spec.round_digits,
    "rounding" => "none for exact arithmetic",
    "fingerprint" => policy_fingerprint(spec),
)

"""Return an integer count, refusing Booleans, fractions, nonfinite values and
floating-point values beyond the range in which every integer is representable."""
function require_count(x)
    x isa Bool && throw(UnsupportedRepresentationError(
        "a count", string(x), "counts; Boolean is not an integer count"))
    x isa Integer && return x
    if x isa Rational
        denominator(x) == 1 && return numerator(x)
        throw(UnsupportedRepresentationError(
            "a count", string(x), "counts; the value is not an integer"))
    end
    if x isa Base.IEEEFloat
        isfinite(x) || throw(UnsupportedRepresentationError(
            "a count", string(x), "counts; NaN and infinity are not counts"))
        isinteger(x) || throw(UnsupportedRepresentationError(
            "a count", string(x), "counts; the value is not an integer"))
        abs(x) <= 2.0^53 || throw(UnsupportedRepresentationError(
            "an exactly-known count", string(x),
            "a floating-point input; above 2^53 its integer history cannot be proven"))
        return Int64(x)
    end
    throw(UnsupportedRepresentationError(
        "a count", string(x), "counts; unsupported type $(typeof(x))"))
end

function checked_count_sum(values; on_overflow::Symbol=:error)
    on_overflow in (:error, :widen) || throw(ArgumentError(
        "on_overflow must be :error or :widen, got :$on_overflow"))
    state = iterate(values)
    isnothing(state) && return 0
    first_value, rest = state
    acc = require_count(first_value)
    while true
        nxt = iterate(values, rest)
        isnothing(nxt) && break
        value, rest = nxt
        count = require_count(value)
        if acc isa BigInt || count isa BigInt
            acc += count
            continue
        end
        try
            acc = Base.Checked.checked_add(promote(acc, count)...)
        catch err
            err isa OverflowError || rethrow()
            on_overflow === :widen || throw(CountOverflowError(BigInt(acc), BigInt(count)))
            acc = BigInt(acc) + BigInt(count)
        end
    end
    acc
end

_bits(x::Integer) = iszero(x) ? 1 : ndigits(abs(x), base=2)

function _check_rational_bits(value::Rational, spec::NumericPolicySpec)
    numerator_bits = _bits(numerator(value))
    numerator_bits <= spec.max_denominator_bits || throw(ResourceLimitError(
        "numerator bits", spec.max_denominator_bits, numerator_bits))
    denominator_bits = _bits(denominator(value))
    denominator_bits <= spec.max_denominator_bits || throw(ResourceLimitError(
        "denominator bits", spec.max_denominator_bits, denominator_bits))
    value
end

function exact_fraction(num::Integer, den::Integer;
                        policy::NumericPolicySpec=numeric_policy(:exact_counts))
    iszero(den) && throw(UnsupportedRepresentationError(
        "a finite rational", "$num/$den", "exact_fraction; zero denominators are refused"))
    _check_rational_bits(Rational{BigInt}(num, den), policy)
end

function exact_rational_sum(values; policy::NumericPolicySpec=numeric_policy(:exact_counts))
    acc = nothing
    for value in values
        value isa Rational || throw(UnsupportedRepresentationError(
            "an exact rational", string(value), "exact_rational_sum"))
        iszero(denominator(value)) && throw(UnsupportedRepresentationError(
            "a finite rational", string(value), "exact_rational_sum"))
        candidate = _check_rational_bits(Rational{BigInt}(value), policy)
        acc = isnothing(acc) ? candidate :
              _check_rational_bits(acc + candidate, policy)
    end
    acc
end

function exact_relative_abundance(count::Integer, total::Integer;
                                  policy::NumericPolicySpec=numeric_policy(:exact_counts))
    total >= 0 || throw(ArgumentError("total must not be negative, got $total"))
    count >= 0 || throw(ArgumentError("count must not be negative, got $count"))
    count <= total || throw(ArgumentError("count ($count) exceeds total ($total)"))
    iszero(total) && return nothing
    exact_fraction(count, total; policy)
end

const EXACT = :exact
const APPROXIMATE = :approximate
const ROUNDED = :rounded

struct NumericValue
    value::Any
    kind::Symbol
    precision_bits::Union{Int,Nothing}
    digits::Union{Int,Nothing}
end

function exact_value(x)
    (x isa Integer && !(x isa Bool)) || x isa Rational ||
        throw(UnsupportedRepresentationError(
            "an exact integer or rational", string(x), "exact_value"))
    if x isa Rational && iszero(denominator(x))
        throw(UnsupportedRepresentationError(
            "a finite rational", string(x), "exact_value"))
    end
    NumericValue(x, EXACT, nothing, nothing)
end

# Approximate values are IEEE binary floats; the type fixes their precision.
# Any other float type is refused, not stored or displayed as if it were one.
function approximate_value(x)
    x isa Base.IEEEFloat || throw(UnsupportedRepresentationError(
        "a floating-point approximation", string(x), "approximate_value"))
    NumericValue(x, APPROXIMATE, precision(x), nothing)
end

function rounded_value(x; digits::Integer)
    x isa Base.IEEEFloat || throw(UnsupportedRepresentationError(
        "a floating-point approximation", string(x), "rounded_value"))
    0 <= digits <= MAX_ROUND_DIGITS || throw(ArgumentError("digits out of range"))
    NumericValue(x, ROUNDED, precision(x), Int(digits))
end

is_exact(v::NumericValue) = v.kind === EXACT

function to_storage(v::NumericValue)
    x = v.value
    if v.kind === EXACT
        if x isa Integer
            return abs(big(x)) <= MAX_EXACT_FLOAT - 1 ? x : string(x)
        end
        return string(numerator(x), "/", denominator(x))
    end
    isfinite(x) || throw(UnsupportedRepresentationError(
        "a stored finite number", string(x), "JSON/storage; NaN and infinity are refused"))
    x isa Float64 && return x
    Float64(x)
end

# Decimal rendering is integer arithmetic. Converting a finite float to
# Rational{BigInt} is exact (its binary significand/exponent), so this does
# not introduce another floating rounding step or depend on Printf.
function _rounded_scaled(num::Integer, den::Integer, scale::Integer)
    n = big(num) * scale
    q, r = divrem(n, den)
    if 2 * abs(r) >= abs(den)
        q += sign(n) * sign(den) # ties away from zero
    end
    q
end

function _rendered_decimal(x::Rational, digits::Integer)
    digits <= MAX_ROUND_DIGITS || throw(ResourceLimitError(
        "display decimal places", MAX_ROUND_DIGITS, Int(digits)))
    scale = big(10)^digits
    # abs(typemin(T)) overflows for fixed-width signed integers. Widen first.
    scaled = _rounded_scaled(abs(big(numerator(x))), denominator(x), scale)
    sign_string = numerator(x) < 0 && !iszero(scaled) ? "-" : ""
    digits == 0 && return string(sign_string, scaled)
    whole, frac = divrem(scaled, scale)
    string(sign_string, whole, ".", lpad(string(frac), digits, '0'))
end

function to_display(v::NumericValue; digits::Integer=6)
    0 <= digits <= MAX_ROUND_DIGITS || throw(ArgumentError("digits out of range"))
    if v.kind === EXACT
        return v.value isa Integer ? string(v.value) :
               string(numerator(v.value), "/", denominator(v.value),
                      " (exact; ", digits, "dp = ",
                      _rendered_decimal(v.value, digits), ")")
    end
    isfinite(v.value) || throw(UnsupportedRepresentationError(
        "a finite display number", string(v.value), "to_display"))
    kept = v.kind === ROUNDED && !isnothing(v.digits) ? v.digits : Int(digits)
    rendered = _rendered_decimal(Rational{BigInt}(v.value), kept)
    bits = isnothing(v.precision_bits) ? "unknown precision" : "$(v.precision_bits)-bit"
    # The display is always a rounding, whatever the value's kind: an
    # approximate value shown at 6dp has still been rounded, and a reader
    # must not mistake the short string for the stored value.
    "$rendered ($(v.kind), $bits, rounded to $(kept)dp)"
end

function encode_json_number(x)
    if x isa Integer && !(x isa Bool)
        return abs(big(x)) <= MAX_EXACT_FLOAT - 1 ? x : string(x)
    end
    x isa Rational && return to_storage(exact_value(x))
    if x isa Base.IEEEFloat
        isfinite(x) || throw(UnsupportedRepresentationError(
            "a finite JSON number", string(x), "JSON"))
        return x isa Float64 ? x : encode_json_number(Float64(x))
    end
    throw(UnsupportedRepresentationError(
        "a JSON number", string(x), "JSON; unsupported type $(typeof(x))"))
end

function parse_exact_count(s::AbstractString; max_bits::Integer=MAX_COUNT_BITS)
    1 <= max_bits <= MAX_COUNT_BITS || throw(ArgumentError("max_bits out of range"))
    text = strip(s)
    isempty(text) && throw(UnsupportedRepresentationError(
        "a count", s, "parse_exact_count; empty input"))
    occursin(r"^[+]?[0-9]+$", text) || throw(UnsupportedRepresentationError(
        "a non-negative whole-number literal", s,
        "parse_exact_count; decimals, exponents and negative counts are refused"))
    # Reject an absurd allocation before parse(BigInt), not after it.
    length(text) <= div(Int(max_bits), 3) + 2 || throw(ResourceLimitError(
        "count decimal digits", div(Int(max_bits), 3) + 2, length(text)))
    value = parse(BigInt, text)
    _bits(value) <= max_bits || throw(ResourceLimitError(
        "count bits", Int(max_bits), _bits(value)))
    typemin(Int64) <= value <= typemax(Int64) ? Int64(value) : value
end

function parse_exact_rational(s::AbstractString;
                              policy::NumericPolicySpec=numeric_policy(:exact_counts))
    text = strip(s)
    isempty(text) && throw(UnsupportedRepresentationError(
        "a rational", s, "parse_exact_rational; empty input"))
    occursin(r"^[+-]?[0-9]+(?:/[+-]?[0-9]+)?$", text) ||
        throw(UnsupportedRepresentationError("an integer or n/d literal", s,
            "parse_exact_rational; decimal points and exponent notation are refused"))
    parts = split(text, '/'; limit=2)
    # Apply the configured bit budget to the textual allocation as well as to
    # the resulting rational; do not parse an unbounded integer first.
    digit_limit = div(policy.max_denominator_bits, 3) + 3
    length(parts[1]) <= digit_limit || throw(ResourceLimitError(
        "rational numerator digits", digit_limit, length(parts[1])))
    length(parts[end]) <= digit_limit || throw(ResourceLimitError(
        "rational denominator digits", digit_limit, length(parts[end])))
    den = length(parts) == 2 ? parse(BigInt, parts[2]) : big(1)
    iszero(den) && throw(UnsupportedRepresentationError(
        "a finite rational", s, "parse_exact_rational; zero denominator"))
    exact_fraction(parse(BigInt, parts[1]), den; policy)
end

end # module NumericPolicy
