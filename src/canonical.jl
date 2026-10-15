# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# JCS (RFC 8785) canonical JSON: deterministic bytes for hashing, signing and
# content-addressing. The same value always yields the same byte string, so a
# chart spec or an analysis result can be hashed where it is produced and the
# hash compared across machines, Julia versions and refactorings that build
# the same object in a different order.
#
# What this implements, and what it deliberately does not:
#
#   * Object keys are sorted by UTF-16 code unit, the rule that makes an
#     astral-plane key sort BEFORE a key such as U+FF21, because its high
#     surrogate (U+D800..U+DBFF) is a smaller code unit than U+FF21.
#   * Float64 values serialize per ECMAScript's Number::toString, the rules
#     RFC 8785 builds on: shortest round-trip digits, positional notation for
#     1e-6 <= |x| < 1e21, exponent notation outside it, a signed exponent
#     ("1e+30"), and -0.0 collapsing to 0.
#   * Integers outside I-JSON's exact range (-(2^53)+1 .. 2^53-1) are
#     refused, as the estate's ijson-jcs crate refuses them, because a
#     consumer cannot read them back as the same number.
#   * NaN and Inf are refused for the same reason.
#
# The implementation is verified against the estate's ijson-jcs crate
# (Rust, RFC 7493 + RFC 8785): the fixtures in
# test/fixtures/jcs_oracle.jsonl were produced BY that crate, and
# test/unit/test_canonical.jl asserts this module reproduces them
# byte-for-byte. The crate is the oracle; this is the second
# implementation, and any disagreement between them is a bug in one of
# them.

module Canonical


using SHA

export jcs_bytes, jcs_sha256

# The I-JSON exact-integer range.
const SAFE_INT_MAX = 2^53 - 1

"""
    jcs_bytes(x) -> Vector{UInt8}

The RFC 8785 canonical form of `x`, as bytes. Accepts the JSON-shaped
values the server produces: `Dict` and `JSON3.Object` (key order is
irrelevant; it is sorted away), `AbstractVector` and `JSON3.Array` (order
is preserved, as JCS preserves array order), `AbstractString`, `Integer`,
`AbstractFloat`, `Bool` and `nothing`. Floats serialize per ECMAScript's
Number::toString; integers outside the exact I-JSON range, and NaN or Inf
floats, throw `ArgumentError`.
"""
function jcs_bytes(x)
    io = IOBuffer()
    _write!(io, x)
    take!(io)
end

"""
    jcs_sha256(x) -> String

The SHA-256 of `jcs_bytes(x)`, as lowercase hex -- the content address of
a JSON value: two producers that build the same value hash identically
whatever order they built it in.
"""
jcs_sha256(x) = bytes2hex(sha256(jcs_bytes(x)))

# Object keys sort by UTF-16 code unit, which for astral-plane characters
# means their surrogate pair, a smaller first unit than any BMP code point
# at or above U+E000. transcoding the key to UInt16 gives exactly those
# units to compare.
_utf16_units(s::AbstractString) = transcode(UInt16, String(s))

_write!(io, ::Nothing) = print(io, "null")
_write!(io, b::Bool) = print(io, b ? "true" : "false")

function _write!(io, n::Integer)
    # I-JSON: only |n| <= 2^53-1 survives a JSON round trip exactly.
    -SAFE_INT_MAX <= n <= SAFE_INT_MAX || throw(ArgumentError(
        "canonical JSON refuses the integer $n: it is outside the exact " *
        "IEEE-754 binary64 range I-JSON allows"))
    print(io, n)
end

function _write!(io, x::AbstractFloat)
    _write_float!(io, Float64(x))
end

function _write!(io, s::AbstractString)
    print(io, '"')
    for c in s
        _write_char!(io, c)
    end
    print(io, '"')
end

# ECMAScript's JSON.stringify escapes exactly these: the two mandatory
# characters, the five with short forms, and every other control
# character as \u00xx with lowercase hex.
function _write_char!(io, c::Char)
    if c == '"'
        print(io, "\\\"")
    elseif c == '\\'
        print(io, "\\\\")
    elseif c == '\b'
        print(io, "\\b")
    elseif c == '\t'
        print(io, "\\t")
    elseif c == '\n'
        print(io, "\\n")
    elseif c == '\f'
        print(io, "\\f")
    elseif c == '\r'
        print(io, "\\r")
    elseif c < ' '
        print(io, "\\u", string(UInt32(c); base=16, pad=4))
    else
        print(io, c)
    end
end

# AbstractDict (Dict and JSON3.Object alike): sorted keys, then each pair.
function _write!(io, d::AbstractDict)
    print(io, '{')
    first = true
    for k in sort!(collect(String.(keys(d))); by = _utf16_units)
        first || print(io, ',')
        first = false
        _write!(io, string(k))
        print(io, ':')
        _write!(io, d[k])
    end
    print(io, '}')
end

# AbstractVector (Vector and JSON3.Array alike): order preserved.
function _write!(io, a::AbstractVector)
    print(io, '[')
    for (i, v) in enumerate(a)
        i > 1 && print(io, ',')
        _write!(io, v)
    end
    print(io, ']')
end

# RFC 8785 §3.2.2.3: ECMAScript's Number::toString. Julia's own repr gives
# the shortest round-trip digits (Ryu, the same family ES specifies), so
# the work is reformatting those digits under ES's notation rules.
function _write_float!(io, x::Float64)
    isnan(x) && throw(ArgumentError(
        "canonical JSON refuses NaN: no consumer reads it back as the same value"))
    isinf(x) && throw(ArgumentError(
        "canonical JSON refuses Inf: no consumer reads it back as the same value"))
    if x == 0
        print(io, "0")   # +0.0 and -0.0 both collapse to 0, as ECMAScript does.
        return
    end
    if signbit(x)
        print(io, '-')
        x = -x
    end
    digits, n = _shortest_digits(x)
    k = length(digits)
    if n > 21 || n <= -6
        # Exponential: one digit, then the fraction, then e with an
        # explicit sign (ES prints "1e+30", not "1e30").
        print(io, digits[1])
        k > 1 && (print(io, '.'); print(io, digits[2:end]))
        e = n - 1
        print(io, 'e', e < 0 ? '-' : '+', string(abs(e)))
    elseif n >= k
        # An integer too big for the fraction: digits followed by zeros.
        print(io, digits, "0"^(n - k))
    elseif n > 0
        # Digits with the point inside them.
        print(io, digits[1:n], '.', digits[n+1:end])
    else
        # Below one, positional as far as ES goes: 0.000001, not 1e-6.
        print(io, "0.", "0"^(-n), digits)
    end
end

# The shortest round-trip decimal form of x as (digits, n) with
# x = 0.digits * 10^n: repr yields either "123", "4.5" or "1.5e-10".
function _shortest_digits(x::Float64)
    s = repr(x)
    mantissa, exponent = occursin('e', s) ? (split(s, 'e') .|> String) : (s, "")
    e = isempty(exponent) ? 0 : parse(Int, exponent)
    if occursin('.', mantissa)
        intpart, fracpart = split(mantissa, '.') .|> String
    else
        intpart, fracpart = mantissa, ""
    end
    digits = string(intpart, fracpart)
    n = length(intpart) + e
    # Trailing zeros are notation, not value: 1.0e30 is digits "1", n 31.
    while length(digits) > 1 && endswith(digits, '0')
        digits = digits[1:end-1]
    end
    digits, n
end

end # module
