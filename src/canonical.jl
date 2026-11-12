# SPDX-License-Identifier: MPL-2.0
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
# The stored fixture strings were checked against ECMAScript JSON.stringify
# with UTF-16 key ordering, and the Julia tests compare this writer with them.
# The inherited claim of Rust ijson-jcs fixture provenance has not been
# independently reproduced; no cross-language Rust verification is claimed.

module Canonical


using SHA, JSON3

export jcs_bytes, jcs_sha256, ijson_read

# The I-JSON exact-integer range.
const SAFE_INT_MAX = 2^53 - 1

"""
    jcs_bytes(x) -> Vector{UInt8}

The RFC 8785 canonical form of `x`, as bytes. Accepts the JSON-shaped
values the server produces: `Dict`, `NamedTuple` and `JSON3.Object` (key order is
irrelevant; it is sorted away), `AbstractVector` and `JSON3.Array` (order
is preserved, as JCS preserves array order), `AbstractString`, `Integer`,
`Float16`, `Float32`, `Float64`, `Bool` and `nothing`. Floats serialize per ECMAScript's
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
    print(io, string(n; base=10))
end

function _write!(io, x::Union{Float16,Float32,Float64})
    _write_float!(io, Float64(x))
end

_write!(io, ::AbstractFloat) = throw(ArgumentError(
    "canonical JSON supports only Float16, Float32 and Float64"))

# I-JSON disallows Unicode noncharacters (including the last two code points
# of every plane), whether they occur in values or property names.
function _valid_ijson_string(s::AbstractString)
    isvalid(s) || return false
    all(s) do c
        cp = UInt32(c)
        !(0xfdd0 <= cp <= 0xfdef || (cp & 0xfffe) == 0xfffe)
    end
end

function _write!(io, s::AbstractString)
    _valid_ijson_string(s) || throw(ArgumentError("canonical JSON refuses invalid I-JSON Unicode"))
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

# NamedTuple is convenient for internal records. Convert to string-keyed
# entries so it follows the same canonical sorting and duplicate checks as
# dictionary objects.
_write!(io, nt::NamedTuple) = _write!(io, Dict{String,Any}(string(k) => v for (k, v) in pairs(nt)))

# AbstractDict (Dict and JSON3.Object alike): sorted keys, then each pair.
function _write!(io, d::AbstractDict)
    print(io, '{')
    first = true
    entries = Pair{String,Any}[]
    names = Set{String}()
    for (key, value) in pairs(d)
        key isa Union{AbstractString,Symbol} || throw(ArgumentError(
            "canonical JSON object keys must be strings or symbols"))
        name = String(key)
        _valid_ijson_string(name) || throw(ArgumentError("canonical JSON refuses invalid I-JSON Unicode keys"))
        name in names && throw(ArgumentError("canonical JSON refuses duplicate property names"))
        push!(names, name)
        push!(entries, name => value)
    end
    sort!(entries; by = entry -> _utf16_units(entry.first))
    for (name, value) in entries
        first || print(io, ',')
        first = false
        _write!(io, name)
        print(io, ':')
        _write!(io, value)
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

"""
    ijson_read(raw) -> JSON3.Object / JSON3.Array / scalar

Parse a JSON document in the strict I-JSON mode of hyperpolymath/ijson-jcs.
Inspect the original tokens *before* relying on JSON3's parsed object: parsed
objects can lose duplicate names, and integral-looking floats/oversized integer
tokens can be normalised. The return value retains JSON3's usual types so
existing callers need not change their field access.
"""
function ijson_read(raw::Union{AbstractString,AbstractVector{UInt8}})
    text = raw isa AbstractString ? String(raw) : String(copy(raw))
    isvalid(text) || throw(ArgumentError("I-JSON requires UTF-8"))
    # Bound nesting before either JSON3 or the recursive raw-token validation
    # sees it. Brackets inside strings and escaped quotes do not affect depth.
    depth = 0
    in_string = false
    escaped = false
    for byte in codeunits(text)
        if in_string
            if escaped
                escaped = false
            elseif byte == 0x5c
                escaped = true
            elseif byte == 0x22
                in_string = false
            end
        elseif byte == 0x22
            in_string = true
        elseif byte == 0x7b || byte == 0x5b
            depth += 1
            depth <= 64 || throw(ArgumentError("JSON nesting exceeds 64 levels"))
        elseif byte == 0x7d || byte == 0x5d
            depth -= 1
            depth >= 0 || throw(ArgumentError("unbalanced JSON delimiters"))
        end
    end
    # JSON3 parses the value; the raw-token walk also checks grammar and
    # trailing content, which JSON3 may otherwise ignore.
    parsed = try
        JSON3.read(text)
    catch err
        throw(ArgumentError("invalid JSON: $(sprint(showerror, err))"))
    end
    bytes = codeunits(text)
    n = length(bytes)
    pos = Ref(1)
    peek() = pos[] <= n ? bytes[pos[]] : 0x00
    function whitespace!()
        while peek() in (0x20, 0x09, 0x0a, 0x0d) && pos[] <= n
            pos[] += 1
        end
    end
    function string!()
        start = pos[]
        pos[] += 1 # opening quote
        while pos[] <= n
            c = bytes[pos[]]
            pos[] += 1
            c == 0x22 && break
            if c == 0x5c
                # JSON3 checked the escape grammar; additionally reject lone
                # surrogate escapes even if a parser replaces them with U+FFFD.
                if bytes[pos[]] == 0x75
                    high = parse(UInt32, String(bytes[pos[]+1:pos[]+4]); base=16)
                    if 0xd800 <= high <= 0xdbff
                        pos[] + 10 <= n && bytes[pos[]+5] == 0x5c &&
                            bytes[pos[]+6] == 0x75 || throw(ArgumentError(
                                "I-JSON refuses an unpaired high surrogate"))
                        low = parse(UInt32, String(bytes[pos[]+7:pos[]+10]); base=16)
                        0xdc00 <= low <= 0xdfff || throw(ArgumentError(
                            "I-JSON refuses an unpaired high surrogate"))
                        pos[] += 11
                        continue
                    elseif 0xdc00 <= high <= 0xdfff
                        throw(ArgumentError("I-JSON refuses an unpaired low surrogate"))
                    end
                    pos[] += 5
                else
                    pos[] += 1
                end
            end
        end
        value = String(JSON3.read(String(bytes[start:pos[]-1])))
        _valid_ijson_string(value) || throw(ArgumentError("I-JSON refuses Unicode noncharacters"))
        value
    end
    function value!()
        whitespace!()
        if peek() == 0x7b # object
            pos[] += 1
            whitespace!()
            seen = Set{String}()
            while peek() != 0x7d
                key = string!()
                key in seen && throw(ArgumentError("I-JSON refuses duplicate object key $(repr(key))"))
                push!(seen, key)
                whitespace!()
                pos[] += 1 # colon (JSON3 already verified the grammar)
                value!()
                whitespace!()
                peek() == 0x2c || break
                pos[] += 1
                whitespace!()
            end
            pos[] += 1 # closing brace
        elseif peek() == 0x5b # array
            pos[] += 1
            whitespace!()
            while peek() != 0x5d
                value!()
                whitespace!()
                peek() == 0x2c || break
                pos[] += 1
                whitespace!()
            end
            pos[] += 1 # closing bracket
        elseif peek() == 0x22
            string!()
        else
            start = pos[]
            while pos[] <= n && !(peek() in (0x2c, 0x5d, 0x7d, 0x20, 0x09, 0x0a, 0x0d))
                pos[] += 1
            end
            token = String(bytes[start:pos[]-1])
            if occursin(r"^-?(0|[1-9][0-9]*)$", token)
                # Bound the token before parsing: an HTTP body could contain
                # millions of digits, but a safe integer has at most 16.
                length(token) <= 17 &&
                    -SAFE_INT_MAX <= parse(BigInt, token) <= SAFE_INT_MAX ||
                    throw(ArgumentError("I-JSON refuses an unsafe integer token $token"))
            elseif occursin(r"^-?(0|[1-9][0-9]*)\.[0-9]+([eE][+-]?[0-9]+)?$|^-?(0|[1-9][0-9]*)[eE][+-]?[0-9]+$", token)
                x = tryparse(Float64, token)
                (!isnothing(x) && isfinite(x)) || throw(ArgumentError(
                    "I-JSON refuses a non-finite number $token"))
                # Exponent/fraction input may stringify as an integer in
                # JCS. Never accept a value whose canonical token would be
                # rejected as an unsafe integer on the next read.
                canonical = String(jcs_bytes(x))
                if occursin(r"^-?(0|[1-9][0-9]*)$", canonical) &&
                   !(abs(x) <= SAFE_INT_MAX)
                    throw(ArgumentError("I-JSON refuses a number canonicalizing to an unsafe integer token $canonical"))
                end
            elseif !(token in ("true", "false", "null"))
                throw(ArgumentError("invalid JSON token $token"))
            end
        end
    end
    value!()
    whitespace!()
    pos[] > n || throw(ArgumentError("invalid JSON: trailing content"))
    parsed
end

end # module
