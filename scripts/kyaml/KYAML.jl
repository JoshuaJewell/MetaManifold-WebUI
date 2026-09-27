# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# KYAML.jl — the YAML <-> KYAML switch for this repository.
#
# Authority: hyperpolymath/standards `3-practice/YAML-POLICY.adoc`, rules Y-2 (writing)
# and Y-3 (KYAML as the target authoring dialect), owner ruling 2026-09-26 making this
# repository the pilot. docs/pilots/kyaml-pilot.md is the operating manual; this file is the
# tool it names.
#
#   use-kyaml  rewrite every estate-owned YAML file as KYAML: `---` then one flow-style
#              document (`{ ... }`), two-space nesting, trailing commas, every string value
#              double-quoted, keys unquoted where that is unambiguous (schema-ambiguous keys
#              such as `on` are quoted, which is what the dialect is for).
#   use-yaml   rewrite them back as block-style YAML. Multi-line strings become block
#              scalars again, and a file that was never converted is re-emitted from its
#              original block-scalar text verbatim.
#   check      parse each file and re-emit it as KYAML; a file whose bytes are not exactly
#              what the emitter writes fails. Canonical-form checking is idempotence by
#              construction: running it twice cannot give two answers.
#
# The contract, stated so that nothing is claimed that is not implemented:
#
#   * COMMENTS ARE PRESERVED, and so is their association: a comment on its own line stays
#     on its own line above the same entry; an end-of-line comment stays on the same line as
#     the same entry. A comment that cannot be placed losslessly is a REFUSAL, never a silent
#     drop. `--report` prints the comment count read and written per file, and
#     test/unit/test_kyaml.jl plants a dropped comment to prove the check goes red
#     (standards §2.2: a check proves nothing until a mutant dies).
#   * BYTES ARE RECOVERABLE: YAML -> KYAML -> YAML reproduces the canonical block form of the
#     same document with the same comments, and `git revert` of the pilot commit reproduces
#     the pre-pilot file byte for byte — the escape hatch the owner asked for.
#   * MEANING IS NOT SILENTLY CHANGED. Every place KYAML forces a decision YAML left implicit
#     is counted and printed by `--report`: `~` becomes `null`; a plain scalar that is not a
#     canonical number, `true`, `false` or `null` gets quoted (so a bare `no` stops being a
#     boolean by accident — the Norway problem); a schema-ambiguous key is quoted; an
#     end-of-line comment on a key whose value is a collection moves to its own line above the
#     key, because a flow collection ends several lines later.
#   * WHAT IS REFUSED IS NAMED, WITH A LINE NUMBER: anchors, aliases, tags, merge keys,
#     multiple documents, directives, tabs in indentation, duplicate keys in one mapping, and
#     multi-line plain scalars. A refusal leaves the file untouched and exits non-zero.
#
# Boundary: this is a repository tool, not the estate-wide formatter that standards#1022
# asks for. It handles the YAML this repository actually contains (see the construct census
# in docs/pilots/kyaml-pilot.md) and refuses the rest. That refusal list is the honest boundary to
# draw until #1022 lands.

module KYAML

export KyamlError, convert, check, git_yaml_paths, render_kyaml, render_yaml

# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

struct KyamlError <: Exception
    path::String
    line::Int
    message::String
end

Base.showerror(io::IO, e::KyamlError) = print(io, e.path, ":", e.line, ": ", e.message)

# ---------------------------------------------------------------------------
# The document model
# ---------------------------------------------------------------------------
#
# A Scalar keeps three things a bare value would lose: the style it was written in (so a
# quoted "15" stays a string while a bare 15 stays a number), the chomping of a block
# scalar, and the block scalar's original lines (so YAML -> YAML is exact rather than nearly).

abstract type Node end

struct Scalar <: Node
    text::String
    style::Symbol        # :plain | :double | :single | :literal | :folded | :empty
    chomp::Symbol        # :clip | :strip | :keep | :none
    raw::Vector{String}  # the block scalar's source lines, verbatim; empty otherwise
end

mutable struct Entry
    key::String
    key_quoted::Bool
    value::Node
    comments::Vector{String}
    inline::String
    blanks::Int
end

mutable struct Mapping <: Node
    entries::Vector{Entry}
    comments::Vector{String}   # comments at the end of the mapping
end

mutable struct Item
    value::Node
    comments::Vector{String}
    inline::String
    blanks::Int
end

mutable struct Sequence <: Node
    items::Vector{Item}
    comments::Vector{String}   # comments at the end of the sequence
end

Mapping() = Mapping(Entry[], String[])
Sequence() = Sequence(Item[], String[])

scalar(text::AbstractString, style::Symbol, chomp::Symbol = :none) =
    Scalar(String(text), style, chomp, String[])

# ---------------------------------------------------------------------------
# Decisions the report counts
# ---------------------------------------------------------------------------

mutable struct Stats
    comments_read::Int
    comments_written::Int
    keys_quoted::Int
    scalars_quoted::Int
    nulls_canonicalised::Int
    inline_comments_moved::Int
    long_lines_kept::Int
end

Stats() = Stats(0, 0, 0, 0, 0, 0, 0)

# ---------------------------------------------------------------------------
# Character helpers (Vector{Char}: this repository's YAML holds non-ASCII prose)
# ---------------------------------------------------------------------------

pad(n::Int) = repeat(" ", max(n, 0))

is_blank_line(l::Vector{Char}) = all(isspace, l)

function leading_spaces(l::Vector{Char})::Int
    n = 0
    while n < length(l) && l[n + 1] == ' '
        n += 1
    end
    return n
end

function rstrip_chars(l::Vector{Char})::Vector{Char}
    j = length(l)
    while j >= 1 && isspace(l[j])
        j -= 1
    end
    return l[1:j]
end

function strip_chars(l::Vector{Char})::Vector{Char}
    return rstrip_chars(l[leading_spaces(l) + 1:end])
end

# The first `#` that starts a comment: outside quotes and preceded by whitespace (or at the
# start of the line). Nothing when the line carries no comment.
function find_comment_start(chars::Vector{Char})::Union{Nothing,Int}
    i = 1
    n = length(chars)
    quote_char = '\0'
    while i <= n
        c = chars[i]
        if quote_char != '\0'
            if quote_char == '"' && c == '\\'
                i += 2
                continue
            elseif c == quote_char
                quote_char = '\0'
            end
        elseif c == '"' || c == '\''
            quote_char = c
        elseif c == '#' && (i == 1 || isspace(chars[i - 1]))
            return i
        end
        i += 1
    end
    return nothing
end

# `key: value` -> (key chars, key was quoted, value chars); nothing when not an entry.
function split_key_value(content::Vector{Char}, path::String, lineno::Int)
    isempty(content) && return nothing
    if content[1] == '"' || content[1] == '\''
        q = content[1]
        j = 2
        while j <= length(content)
            if q == '"' && content[j] == '\\'
                j += 2
                continue
            end
            content[j] == q && break
            j += 1
        end
        j > length(content) && throw(KyamlError(path, lineno, "unterminated quoted key"))
        key = content[2:j - 1]
        k = j + 1
        while k <= length(content) && content[k] == ' '
            k += 1
        end
        if k > length(content) || content[k] != ':'
            return nothing
        end
        k += 1
        while k <= length(content) && content[k] == ' '
            k += 1
        end
        return (key, true, k <= length(content) ? content[k:end] : Char[])
    end
    n = length(content)
    k = 1
    while k <= n
        if content[k] == ':' && (k == n || content[k + 1] == ' ')
            key = content[1:k - 1]
            isempty(key) && return nothing
            c0 = content[1]
            (c0 == '&' || c0 == '*' || c0 == '!') && throw(KyamlError(path, lineno,
                "anchors, aliases and tags are not supported by this tool; keep the file in YAML"))
            return (key, false, content[min(k + 1, n):end])
        end
        k += 1
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Double-quoted scalar escapes
# ---------------------------------------------------------------------------

function unescape_double(chars::Vector{Char}, path::String, lineno::Int)::String
    out = IOBuffer()
    i = 1
    while i <= length(chars)
        c = chars[i]
        if c == '\\'
            i + 1 > length(chars) &&
                throw(KyamlError(path, lineno, "trailing backslash in a double-quoted scalar"))
            e = chars[i + 1]
            if e == 'n'
                print(out, '\n')
            elseif e == 't'
                print(out, '\t')
            elseif e == 'r'
                print(out, '\r')
            elseif e == '"'
                print(out, '"')
            elseif e == '\\'
                print(out, '\\')
            elseif e == '/'
                print(out, '/')
            elseif e == 'u'
                i + 5 > length(chars) && throw(KyamlError(path, lineno, "short \\u escape"))
                print(out, Char(parse(UInt32, String(chars[i + 2:i + 5]); base = 16)))
                i += 4
            else
                throw(KyamlError(path, lineno, "unsupported escape \\$e in a double-quoted scalar"))
            end
            i += 2
        else
            print(out, c)
            i += 1
        end
    end
    return String(take!(out))
end

function escape_double(s::AbstractString)::String
    out = IOBuffer()
    for c in s
        if c == '"'
            print(out, "\\\"")
        elseif c == '\\'
            print(out, "\\\\")
        elseif c == '\n'
            print(out, "\\n")
        elseif c == '\t'
            print(out, "\\t")
        elseif c == '\r'
            print(out, "\\r")
        elseif Int(c) < 0x20
            print(out, "\\u", string(UInt16(Int(c)); base = 16, pad = 4))
        else
            print(out, c)
        end
    end
    return String(take!(out))
end

# ---------------------------------------------------------------------------
# Block-style parsing
# ---------------------------------------------------------------------------

mutable struct BlockParser
    path::String
    lines::Vector{Vector{Char}}   # a mutable copy: `- ` becomes two spaces in place
    i::Int
    comments::Vector{String}
    blanks::Int
    stats::Stats
end

function prepare!(p::BlockParser)
    while p.i <= length(p.lines)
        l = p.lines[p.i]
        if is_blank_line(l)
            p.blanks += 1
            p.i += 1
            continue
        end
        ind = leading_spaces(l)
        if ind < length(l) && l[ind + 1] == '\t'
            throw(KyamlError(p.path, p.i, "a tab is used for indentation; YAML forbids that"))
        end
        rest = rstrip_chars(l[ind + 1:end])
        if !isempty(rest) && rest[1] == '#'
            push!(p.comments, strip(String(rest[2:end])))
            p.stats.comments_read += 1
            p.i += 1
            continue
        end
        return (ind, rest, p.i)
    end
    return nothing
end

function take_pending!(p::BlockParser)
    c = p.comments
    b = p.blanks
    p.comments = String[]
    p.blanks = 0
    return (c, b)
end

function parse_node!(p::BlockParser, minindent::Int)::Union{Node,Nothing}
    sig = prepare!(p)
    sig === nothing && return nothing
    (ind, rest, _) = sig
    ind < minindent && return nothing
    if rest[1] == '-' && (length(rest) == 1 || rest[2] == ' ')
        return parse_sequence!(p, ind)
    end
    split_key_value(rest, p.path, p.i) === nothing && throw(KyamlError(p.path, p.i,
        "a document whose root is a bare scalar is not supported; make it a mapping or a sequence"))
    return parse_mapping!(p, ind)
end

function parse_mapping!(p::BlockParser, indent::Int)::Mapping
    m = Mapping()
    while true
        sig = prepare!(p)
        if sig === nothing
            (c, _) = take_pending!(p)
            append!(m.comments, c)
            break
        end
        (ind, rest, _) = sig
        if ind < indent
            (c, _) = take_pending!(p)
            append!(m.comments, c)
            break
        end
        ind > indent && throw(KyamlError(p.path, p.i,
            "unexpected indentation: expected $indent spaces, found $ind. A multi-line plain scalar is the usual cause and is not supported — quote the value or use a block scalar"))
        (comments, blanks) = take_pending!(p)
        chopped = rest
        inline = ""
        pos = find_comment_start(rest)
        if pos !== nothing
            inline = strip(String(rest[pos + 1:end]))
            chopped = rstrip_chars(rest[1:pos - 1])
        end
        kv = split_key_value(chopped, p.path, p.i)
        kv === nothing && throw(KyamlError(p.path, p.i, "expected `key: value`"))
        (keychars, key_quoted, valuechars) = kv
        key = String(keychars)
        for e in m.entries
            e.key == key && throw(KyamlError(p.path, p.i, "duplicate key `$key` in one mapping"))
        end
        local value::Node
        if isempty(valuechars)
            p.i += 1
            value = nothing
            sig2 = prepare!(p)
            if sig2 !== nothing
                (ind2, rest2, _) = sig2
                if ind2 == indent && rest2[1] == '-' && (length(rest2) == 1 || rest2[2] == ' ')
                    # `key:` followed by `- item` at the key's own indentation: YAML reads that
                    # sequence as this key's value, and so does the rest of the world.
                    value = parse_sequence!(p, ind2)
                elseif ind2 > indent
                    value = parse_node!(p, ind2)
                end
            end
            value === nothing && (value = scalar("", :empty))
        else
            value = parse_value!(p, valuechars, indent)
        end
        if value isa Scalar && (value.style == :literal || value.style == :folded) && !isempty(inline)
            # A block scalar swallows the rest of its line: an inline comment cannot follow it.
            throw(KyamlError(p.path, p.i, "an end-of-line comment on a block scalar cannot be represented"))
        end
        push!(m.entries, Entry(key, key_quoted, value, comments, inline, blanks))
    end
    return m
end

function parse_sequence!(p::BlockParser, indent::Int)::Sequence
    s = Sequence()
    while true
        sig = prepare!(p)
        if sig === nothing
            (c, _) = take_pending!(p)
            append!(s.comments, c)
            break
        end
        (ind, rest, _) = sig
        if ind < indent || !(rest[1] == '-' && (length(rest) == 1 || rest[2] == ' '))
            (c, _) = take_pending!(p)
            append!(s.comments, c)
            break
        end
        ind > indent && throw(KyamlError(p.path, p.i, "unexpected indentation inside a sequence"))
        (comments, blanks) = take_pending!(p)
        after = length(rest) == 1 ? Char[] : rest[3:end]
        # `- ` becomes two spaces in place, so the remainder parses at indent + 2.
        p.lines[p.i] = vcat(fill(' ', indent + 2), after)
        local value::Node
        if isempty(after)
            p.i += 1
            child = parse_node!(p, indent + 1)
            value = child === nothing ? scalar("", :empty) : child
            push!(s.items, Item(value, comments, "", blanks))
        elseif after[1] == '-' && (length(after) == 1 || after[2] == ' ')
            value = parse_sequence!(p, indent + 2)
            push!(s.items, Item(value, comments, "", blanks))
        elseif split_key_value(after, p.path, p.i) !== nothing
            value = parse_mapping!(p, indent + 2)
            push!(s.items, Item(value, comments, "", blanks))
        else
            inline = ""
            chopped = after
            pos = find_comment_start(after)
            if pos !== nothing
                inline = strip(String(after[pos + 1:end]))
                chopped = rstrip_chars(after[1:pos - 1])
            end
            value = parse_value!(p, chopped, indent)
            push!(s.items, Item(value, comments, inline, blanks))
        end
    end
    return s
end

const BLOCK_HEADER = r"^([|>])([0-9]?)([+-]?)$"

function parse_value!(p::BlockParser, chars::Vector{Char}, key_indent::Int)::Node
    lineno = p.i
    head = String(chars)
    m = match(BLOCK_HEADER, head)
    if m !== nothing
        p.i += 1
        return parse_block_scalar!(p, m[1][1], m[2][1], m[3][1], key_indent, lineno)
    end
    if chars[1] == '[' || chars[1] == '{'
        f = FlowParser(p.path, vcat(chars, ['\n']), 1, lineno, p.stats)
        node = parse_flow_node!(f)
        skip_flow_trivia!(f)
        if f.i <= length(f.chars)
            throw(KyamlError(p.path, lineno,
                "a flow collection that continues past the end of its line is not supported; write it on one line or use block style"))
        end
        p.i += 1
        return node
    end
    if chars[1] == '&' || chars[1] == '*' || chars[1] == '!'
        throw(KyamlError(p.path, lineno,
            "anchors, aliases and tags are not supported by this tool; keep the file in YAML"))
    end
    if chars[1] == '"'
        j = 2
        while j <= length(chars)
            if chars[j] == '\\'
                j += 2
                continue
            end
            chars[j] == '"' && break
            j += 1
        end
        j > length(chars) && throw(KyamlError(p.path, lineno, "unterminated double-quoted scalar"))
        isempty(rstrip_chars(chars[j + 1:end])) ||
            throw(KyamlError(p.path, lineno, "text after a double-quoted scalar"))
        p.i += 1
        return scalar(unescape_double(chars[2:j - 1], p.path, lineno), :double)
    end
    if chars[1] == '\''
        j = 2
        while j <= length(chars)
            if chars[j] == '\''
                (j < length(chars) && chars[j + 1] == '\'') ? (j += 2; continue) : break
            end
            j += 1
        end
        j > length(chars) && throw(KyamlError(p.path, lineno, "unterminated single-quoted scalar"))
        isempty(rstrip_chars(chars[j + 1:end])) ||
            throw(KyamlError(p.path, lineno, "text after a single-quoted scalar"))
        p.i += 1
        return scalar(replace(String(chars[2:j - 1]), "''" => "'"), :single)
    end
    p.i += 1
    return scalar(String(chars), :plain)
end

function parse_block_scalar!(p::BlockParser, indent_char::Char, digit_char::Char,
                             chomp_char::Char, key_indent::Int, lineno::Int)::Scalar
    style = indent_char == '|' ? :literal : :folded
    chomp = chomp_char == '-' ? :strip : chomp_char == '+' ? :keep : :clip
    explicit = digit_char == '0' ? 0 : (digit_char == '\0' ? 0 : digit_char - '0')
    block_indent = explicit > 0 ? key_indent + explicit : -1
    raw = String[]
    content = String[]
    while p.i <= length(p.lines)
        l = p.lines[p.i]
        if is_blank_line(l)
            push!(raw, String(l))
            push!(content, "")
            p.i += 1
            continue
        end
        ind = leading_spaces(l)
        ind <= key_indent && break
        block_indent < 0 && (block_indent = ind)
        push!(raw, String(l))
        push!(content, String(l[min(ind, block_indent) + 1:end]))
        p.i += 1
    end
    while !isempty(content) && isempty(content[end]) && chomp != :keep
        pop!(content)
    end
    text = style == :literal ? literal_text(content, chomp) : folded_text(content, chomp)
    _ = lineno
    return Scalar(text, style, chomp, raw)
end

function literal_text(content::Vector{String}, chomp::Symbol)::String
    body = join(content, "\n")
    chomp == :strip && return body
    return body * "\n"
end

function folded_text(content::Vector{String}, chomp::Symbol)::String
    out = IOBuffer()
    first = true
    previous_blank = false
    for line in content
        if first
            print(out, line)
            first = false
            previous_blank = false
            continue
        end
        if isempty(line)
            print(out, '\n')
            previous_blank = true
            continue
        end
        print(out, previous_blank ? "" : " ")
        print(out, line)
        previous_blank = false
    end
    body = String(take!(out))
    chomp == :strip && return body
    return body * "\n"
end

# ---------------------------------------------------------------------------
# Flow-style parsing (this tool's own KYAML output, and hand-written KYAML)
# ---------------------------------------------------------------------------

mutable struct FlowParser
    path::String
    chars::Vector{Char}
    i::Int
    lineno::Int
    stats::Stats
end

function skip_flow_trivia!(f::FlowParser)
    while f.i <= length(f.chars)
        c = f.chars[f.i]
        if c == '\n'
            f.lineno += 1
            f.i += 1
        elseif isspace(c)
            f.i += 1
        else
            return
        end
    end
end

# Consumes whitespace and own-line comments before the root node, returning the comments.
function consume_leading_comments!(f::FlowParser)::Vector{String}
    comments = String[]
    while true
        while f.i <= length(f.chars) && (f.chars[f.i] == ' ' || f.chars[f.i] == '\t' || f.chars[f.i] == '\n')
            f.chars[f.i] == '\n' && (f.lineno += 1)
            f.i += 1
        end
        if f.i <= length(f.chars) && f.chars[f.i] == '#'
            start = f.i + 1
            while f.i <= length(f.chars) && f.chars[f.i] != '\n'
                f.i += 1
            end
            push!(comments, strip(String(f.chars[start:f.i - 1])))
            f.stats.comments_read += 1
        else
            return comments
        end
    end
end

function flow_peek!(f::FlowParser)::Char
    skip_flow_trivia!(f)
    f.i > length(f.chars) &&
        throw(KyamlError(f.path, f.lineno, "unexpected end of file inside a flow collection"))
    return f.chars[f.i]
end

# Own-line comments and blank lines ahead of the next entry or item.
function collect_flow_comments!(f::FlowParser)
    comments = String[]
    blanks = 0
    while true
        sawblank = false
        while f.i <= length(f.chars) && (f.chars[f.i] == ' ' || f.chars[f.i] == '\n' || f.chars[f.i] == '\t')
            f.chars[f.i] == '\n' && (sawblank = true; f.lineno += 1)
            f.i += 1
        end
        sawblank && (blanks += 1)
        if f.i <= length(f.chars) && f.chars[f.i] == '#'
            start = f.i + 1
            while f.i <= length(f.chars) && f.chars[f.i] != '\n'
                f.i += 1
            end
            push!(comments, strip(String(f.chars[start:f.i - 1])))
            f.stats.comments_read += 1
        else
            return (comments, blanks)
        end
    end
end

# An inline comment sits on the same line as the entry's terminating comma.
function same_line_comment!(f::FlowParser)::String
    save = f.i
    while f.i <= length(f.chars) && (f.chars[f.i] == ' ' || f.chars[f.i] == '\t')
        f.i += 1
    end
    if f.i <= length(f.chars) && f.chars[f.i] == '#'
        start = f.i + 1
        while f.i <= length(f.chars) && f.chars[f.i] != '\n'
            f.i += 1
        end
        f.stats.comments_read += 1
        return strip(String(f.chars[start:f.i - 1]))
    end
    f.i = save
    return ""
end

function parse_flow_node!(f::FlowParser)::Node
    c = flow_peek!(f)
    c == '{' && return parse_flow_mapping!(f)
    c == '[' && return parse_flow_sequence!(f)
    return parse_flow_scalar!(f)
end

function parse_flow_mapping!(f::FlowParser)::Mapping
    m = Mapping()
    f.i += 1
    while true
        (comments, blanks) = collect_flow_comments!(f)
        if flow_peek!(f) == '}'
            f.i += 1
            append!(m.comments, comments)
            return m
        end
        key_quoted = flow_peek!(f) == '"'
        key = parse_flow_key!(f)
        flow_peek!(f) == ':' ||
            throw(KyamlError(f.path, f.lineno, "expected `:` after the key `$key`"))
        f.i += 1
        for e in m.entries
            e.key == key && throw(KyamlError(f.path, f.lineno, "duplicate key `$key` in one mapping"))
        end
        value = parse_flow_node!(f)
        nextc = flow_peek!(f)
        inline = ""
        if nextc == ','
            f.i += 1
            inline = same_line_comment!(f)
        elseif nextc != '}'
            throw(KyamlError(f.path, f.lineno, "expected `,` or `}` after the value of `$key`"))
        end
        if value isa Mapping || value isa Sequence
            isempty(inline) || (push!(comments, inline); f.stats.inline_comments_moved += 1)
            inline = ""
        end
        push!(m.entries, Entry(key, key_quoted, value, comments, inline, blanks))
    end
end

function parse_flow_key!(f::FlowParser)::String
    c = flow_peek!(f)
    if c == '"'
        f.i += 1
        buf = Char[]
        while f.i <= length(f.chars)
            if f.chars[f.i] == '\\'
                push!(buf, f.chars[f.i]); push!(buf, f.chars[f.i + 1]); f.i += 2; continue
            end
            f.chars[f.i] == '"' && break
            push!(buf, f.chars[f.i])
            f.i += 1
        end
        f.i > length(f.chars) && throw(KyamlError(f.path, f.lineno, "unterminated quoted key"))
        f.i += 1
        return unescape_double(buf, f.path, f.lineno)
    end
    start = f.i
    while f.i <= length(f.chars) && !(f.chars[f.i] == ':' && f.i < length(f.chars) &&
                                       (f.chars[f.i + 1] == ' ' || f.chars[f.i + 1] == '\n'))
        f.chars[f.i] == '\n' && throw(KyamlError(f.path, f.lineno, "a newline inside a plain key"))
        f.i += 1
    end
    return strip(String(f.chars[start:f.i - 1]))
end

function parse_flow_sequence!(f::FlowParser)::Sequence
    s = Sequence()
    f.i += 1
    while true
        (comments, blanks) = collect_flow_comments!(f)
        if flow_peek!(f) == ']'
            f.i += 1
            append!(s.comments, comments)
            return s
        end
        value = parse_flow_node!(f)
        nextc = flow_peek!(f)
        inline = ""
        if nextc == ','
            f.i += 1
            inline = same_line_comment!(f)
        elseif nextc != ']'
            throw(KyamlError(f.path, f.lineno, "expected `,` or `]` after an item"))
        end
        if value isa Mapping || value isa Sequence
            isempty(inline) || (push!(comments, inline); f.stats.inline_comments_moved += 1)
            inline = ""
        end
        push!(s.items, Item(value, comments, inline, blanks))
    end
end

function parse_flow_scalar!(f::FlowParser)::Scalar
    if flow_peek!(f) == '"'
        f.i += 1
        buf = Char[]
        while true
            f.i > length(f.chars) &&
                throw(KyamlError(f.path, f.lineno, "unterminated double-quoted scalar"))
            c = f.chars[f.i]
            if c == '\\'
                if f.i < length(f.chars) && f.chars[f.i + 1] == '\n'
                    # an escaped line break: the break and the next line's indentation vanish
                    f.lineno += 1
                    f.i += 2
                    while f.i <= length(f.chars) && (f.chars[f.i] == ' ' || f.chars[f.i] == '\t')
                        f.i += 1
                    end
                    continue
                end
                f.i + 1 > length(f.chars) &&
                    throw(KyamlError(f.path, f.lineno, "trailing backslash in a double-quoted scalar"))
                push!(buf, '\\'); push!(buf, f.chars[f.i + 1])
                f.i += 2
                continue
            elseif c == '"'
                f.i += 1
                break
            elseif c == '\n'
                f.lineno += 1
                push!(buf, ' ')
                f.i += 1
                continue
            end
            push!(buf, c)
            f.i += 1
        end
        return scalar(unescape_double(buf, f.path, f.lineno), :double)
    end
    if flow_peek!(f) == '\''
        f.i += 1
        buf = Char[]
        while true
            f.i > length(f.chars) &&
                throw(KyamlError(f.path, f.lineno, "unterminated single-quoted scalar"))
            c = f.chars[f.i]
            if c == '\''
                if f.i < length(f.chars) && f.chars[f.i + 1] == '\''
                    push!(buf, '\'')
                    f.i += 2
                    continue
                end
                f.i += 1
                break
            end
            push!(buf, c)
            f.i += 1
        end
        return scalar(String(buf), :single)
    end
    start = f.i
    while f.i <= length(f.chars) && !(f.chars[f.i] in (',', '}', ']', ':'))
        f.chars[f.i] == '\n' && throw(KyamlError(f.path, f.lineno, "a newline inside a plain scalar"))
        f.i += 1
    end
    f.i == start && throw(KyamlError(f.path, f.lineno, "empty scalar"))
    text = strip(String(f.chars[start:f.i - 1]))
    (text == "null" || text == "~") && return scalar(text, :plain)
    return scalar(text, :plain)
end

# ---------------------------------------------------------------------------
# KYAML emission
# ---------------------------------------------------------------------------

function is_canonical_int(t::AbstractString)::Bool
    isempty(t) && return false
    s = startswith(t, "-") ? t[2:end] : t
    isempty(s) && return false
    all(isdigit, s) || return false
    (length(s) > 1 && s[1] == '0') && return false   # 007 is octal in YAML 1.1: not canonical
    return true
end

is_canonical_float(t::AbstractString)::Bool =
    occursin(r"^-?[0-9]+\.[0-9]+([eE][-+]?[0-9]+)?$", t) || occursin(r"^-?\.[0-9]+$", t)

# A plain scalar stays bare only where every YAML 1.1 and YAML 1.2 reader agrees what it is.
function scalar_is_bare(s::Scalar)::Bool
    s.style == :plain || return false
    t = s.text
    (isempty(t) || t == "~" || t == "null") && return false
    lowercase(t) in ("true", "false") && return true
    is_canonical_int(t) && return true
    is_canonical_float(t) && return true
    return false
end

function key_is_bare(key::AbstractString, was_quoted::Bool)::Bool
    was_quoted && return false
    isempty(key) && return false
    lowercase(key) in ("on", "off", "yes", "no", "y", "n", "true", "false", "null", "~") && return false
    return occursin(r"^[A-Za-z_][A-Za-z0-9_.-]*$", key)
end

function quoted_lines(text::AbstractString, indent::Int)::Vector{String}
    segments = split(text, '\n')
    escaped = [escape_double(s) for s in segments]
    n = length(escaped)
    if n == 1
        return ["\"" * escaped[1] * "\""]
    end
    out = String[]
    push!(out, "\"" * escaped[1] * "\\n\\")
    for k in 2:n - 1
        push!(out, pad(indent) * escaped[k] * "\\n\\")
    end
    push!(out, pad(indent) * escaped[n] * "\"")
    return out
end

function scalar_lines(s::Scalar, indent::Int, stats::Stats)::Vector{String}
    if s.style == :empty || (s.style == :plain && (isempty(s.text) || s.text == "~" || s.text == "null"))
        stats.nulls_canonicalised += 1
        return ["null"]
    end
    scalar_is_bare(s) && return [s.text]
    stats.scalars_quoted += 1
    return quoted_lines(s.text, indent)
end

# Appends the value lines. The first line carries no indentation (the caller puts it after
# `key: ` or `- `); continuation lines are indented to `indent + 2`.
function flow_value!(out::Vector{String}, node::Node, indent::Int, stats::Stats)
    if node isa Scalar
        append!(out, scalar_lines(node, indent + 2, stats))
        return
    end
    if node isa Mapping
        push!(out, pad(indent) * "{")
        for e in node.entries
            for _ in 1:e.blanks
                push!(out, "")
            end
            for c in e.comments
                push!(out, pad(indent + 2) * "# " * c)
                stats.comments_written += 1
            end
            bare = key_is_bare(e.key, e.key_quoted)
            bare || (stats.keys_quoted += 1)
            keytext = bare ? e.key : "\"" * escape_double(e.key) * "\""
            inner = String[]
            flow_value!(inner, e.value, indent + 2, stats)
            push!(out, pad(indent + 2) * keytext * ": " * lstrip(inner[1]))
            append!(out, inner[2:end])
            suffix = ","
            isempty(e.inline) || (suffix *= "  # " * e.inline)
            out[end] *= suffix
        end
        for c in node.comments
            push!(out, pad(indent + 2) * "# " * c)
            stats.comments_written += 1
        end
        push!(out, pad(indent) * "}")
        return
    end
    if node isa Sequence
        push!(out, pad(indent) * "[")
        for it in node.items
            for _ in 1:it.blanks
                push!(out, "")
            end
            for c in it.comments
                push!(out, pad(indent + 2) * "# " * c)
                stats.comments_written += 1
            end
            inner = String[]
            flow_value!(inner, it.value, indent + 2, stats)
            push!(out, inner[1])
            append!(out, inner[2:end])
            suffix = ","
            isempty(it.inline) || (suffix *= "  # " * it.inline)
            out[end] *= suffix
        end
        for c in node.comments
            push!(out, pad(indent + 2) * "# " * c)
            stats.comments_written += 1
        end
        push!(out, pad(indent) * "]")
        return
    end
    throw(KyamlError("<emitter>", 0, "unknown node type $(typeof(node))"))
end

function render_kyaml(doc::Node; stats::Stats = Stats())::String
    out = String["---"]
    flow_value!(out, doc, 0, stats)
    return join(out, "\n") * "\n"
end

# ---------------------------------------------------------------------------
# Block-style (YAML) emission
# ---------------------------------------------------------------------------

function block_scalar_header(s::Scalar)::String
    header = s.style == :folded ? ">" : "|"
    if !isempty(s.raw)
        return header * (s.chomp == :strip ? "-" : s.chomp == :keep ? "+" : "")
    end
    trailing = 0
    while length(s.text) > trailing && s.text[end - trailing] == '\n'
        trailing += 1
    end
    return header * (trailing == 0 ? "-" : trailing == 1 ? "" : "+")
end

function block_scalar_body(s::Scalar, indent::Int)::Vector{String}
    if !isempty(s.raw)
        return copy(s.raw)
    end
    trailing = 0
    while length(s.text) > trailing && s.text[end - trailing] == '\n'
        trailing += 1
    end
    body = trailing >= 1 ? s.text[1:end - trailing] : s.text
    isempty(body) && return String[]
    return [isempty(l) ? "" : pad(indent) * l for l in split(body, '\n')]
end

function scalar_inline_text(s::Scalar, stats::Stats)::String
    if s.style == :empty
        stats.nulls_canonicalised += 1
        return "null"
    elseif s.style == :plain
        return s.text
    elseif s.style == :single
        return "'" * replace(s.text, "'" => "''") * "'"
    end
    stats.scalars_quoted += 1
    return "\"" * escape_double(s.text) * "\""
end

function block_mapping!(out::Vector{String}, m::Mapping, indent::Int, stats::Stats)
    for e in m.entries
        for _ in 1:e.blanks
            push!(out, "")
        end
        for c in e.comments
            push!(out, pad(indent) * "# " * c)
            stats.comments_written += 1
        end
        bare = key_is_bare(e.key, e.key_quoted)
        bare || (stats.keys_quoted += 1)
        keytext = bare ? e.key : "\"" * escape_double(e.key) * "\""
        push!(out, pad(indent) * keytext * ":")
        value = e.value
        if value isa Scalar
            if value.style == :literal || value.style == :folded
                out[end] *= " " * block_scalar_header(value)
                append!(out, block_scalar_body(value, indent + 2))
            else
                out[end] *= " " * scalar_inline_text(value, stats)
            end
            isempty(e.inline) || (out[end] *= "  # " * e.inline)
        elseif value isa Sequence
            block_sequence!(out, value, indent + 2, stats)
            isempty(e.inline) || (out[end] *= "  # " * e.inline)
        elseif value isa Mapping
            if isempty(e.inline)
                block_mapping!(out, value, indent + 2, stats)
            else
                push!(out, pad(indent + 2) * "# " * e.inline)
                stats.comments_written += 1
                block_mapping!(out, value, indent + 2, stats)
            end
        end
    end
    for c in m.comments
        push!(out, pad(indent) * "# " * c)
        stats.comments_written += 1
    end
    return
end

function block_sequence!(out::Vector{String}, s::Sequence, indent::Int, stats::Stats)
    for it in s.items
        for _ in 1:it.blanks
            push!(out, "")
        end
        for c in it.comments
            push!(out, pad(indent) * "# " * c)
            stats.comments_written += 1
        end
        value = it.value
        if value isa Scalar
            if value.style == :literal || value.style == :folded
                push!(out, pad(indent) * "- " * block_scalar_header(value))
                append!(out, block_scalar_body(value, indent + 2))
            else
                push!(out, pad(indent) * "- " * scalar_inline_text(value, stats))
            end
            isempty(it.inline) || (out[end] *= "  # " * it.inline)
        elseif value isa Sequence
            push!(out, pad(indent) * "-")
            block_sequence!(out, value, indent + 2, stats)
            isempty(it.inline) || (out[end] *= "  # " * it.inline)
        elseif value isa Mapping
            before = length(out)
            block_mapping!(out, value, indent + 2, stats)
            if length(out) > before
                out[before + 1] = pad(indent) * "- " * lstrip(out[before + 1])
            else
                push!(out, pad(indent) * "- {}")
            end
            isempty(it.inline) || (out[end] *= "  # " * it.inline)
        end
    end
    for c in s.comments
        push!(out, pad(indent) * "# " * c)
        stats.comments_written += 1
    end
    return
end

function render_yaml(doc::Node; stats::Stats = Stats())::String
    out = String[]
    if doc isa Mapping
        block_mapping!(out, doc, 0, stats)
    elseif doc isa Sequence
        block_sequence!(out, doc, 0, stats)
    else
        push!(out, scalar_inline_text(doc::Scalar, stats))
    end
    return join(out, "\n") * "\n"
end

# ---------------------------------------------------------------------------
# Front door
# ---------------------------------------------------------------------------

function parse_document(path::String, source::String; stats::Stats = Stats())::Node
    lines = [collect(l) for l in split(chomp(source), "\n")]
    first_significant = 0
    for (k, l) in enumerate(lines)
        is_blank_line(l) && continue
        rest = rstrip_chars(l[leading_spaces(l) + 1:end])
        isempty(rest) && continue
        rest[1] == '#' && (stats.comments_read += 1; continue)
        first_significant = k
        break
    end
    first_significant == 0 && return scalar("", :empty)
    ind = leading_spaces(lines[first_significant])
    rest = rstrip_chars(lines[first_significant][ind + 1:end])
    if String(rest) == "---"
        chars = Char[]
        for (k, l) in enumerate(lines)
            if k == first_significant
                continue
            end
            append!(chars, l)
            push!(chars, '\n')
        end
        f = FlowParser(path, chars, 1, first_significant + 1, stats)
        # Comments above the root node belong to the document, not to nothing. They are
        # attached to the first entry (or item) so that they survive the round trip; the
        # emitter places them there, and the placement is stable from the first conversion on.
        leading = consume_leading_comments!(f)
        doc = parse_flow_node!(f)
        if !isempty(leading)
            if doc isa Mapping && !isempty(doc.entries)
                prepend!(doc.entries[1].comments, leading)
            elseif doc isa Sequence && !isempty(doc.items)
                prepend!(doc.items[1].comments, leading)
            else
                append!(doc.comments, leading)
            end
        end
        skip_flow_trivia!(f)
        if f.i <= length(f.chars)
            c = f.chars[f.i]
            c == '.' && throw(KyamlError(path, f.lineno, "a second document is not supported"))
            c == '-' && throw(KyamlError(path, f.lineno, "a second document is not supported"))
            throw(KyamlError(path, f.lineno, "unexpected trailing content at the end of the document"))
        end
        return doc
    end
    if String(rest) == "..." || startswith(String(rest), "%")
        throw(KyamlError(path, first_significant,
            "document markers and directives are not supported; keep the file in YAML"))
    end
    p = BlockParser(path, lines, 1, String[], 0, stats)
    doc = parse_node!(p, 0)
    doc === nothing && return scalar("", :empty)
    while p.i <= length(p.lines)
        l = p.lines[p.i]
        if is_blank_line(l)
            p.i += 1
            continue
        end
        rest2 = rstrip_chars(l[leading_spaces(l) + 1:end])
        if !isempty(rest2) && rest2[1] == '#'
            stats.comments_read += 1
            p.i += 1
            continue
        end
        throw(KyamlError(path, p.i,
            "trailing content after the document; a second document is not supported"))
    end
    return doc
end

function convert(path::String; to::Symbol, stats::Stats = Stats())::String
    source = read(path, String)
    doc = parse_document(path, source; stats = stats)
    return to == :kyaml ? render_kyaml(doc; stats = stats) : render_yaml(doc; stats = stats)
end

function check(path::String; stats::Stats = Stats())::Bool
    source = read(path, String)
    doc = parse_document(path, source; stats = stats)
    return source == render_kyaml(doc; stats = stats)
end

function git_yaml_paths(dir::AbstractString = ".")::Vector{String}
    out = try
        read(`git -C $dir ls-files -- "*.yml" "*.yaml"`, String)
    catch err
        throw(KyamlError(String(dir), 0, "git ls-files failed: $(sprint(showerror, err))"))
    end
    return [String(l) for l in split(chomp(out), "\n") if !isempty(l)]
end

end # module KYAML

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

function _kyaml_cli(argv::Vector{String})::Int
    mode = ""
    expect = :kyaml
    report = false
    paths = String[]
    skip_file = joinpath(dirname(@__DIR__), "config", "kyaml", "drift.txt")
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--to-kyaml"
            mode = "to-kyaml"
        elseif a == "--to-yaml"
            mode = "to-yaml"
        elseif a == "--check"
            mode = "check"
        elseif a == "--report"
            report = true
        elseif a == "--skip-file"
            i += 1
            skip_file = argv[i]
        elseif a == "--expect"
            i += 1
            expect = Symbol(argv[i])
        elseif a == "--help" || a == "-h"
            print("""
            KYAML.jl — switch this repository's YAML between KYAML and YAML.

              julia --project=no scripts/kyaml/KYAML.jl --to-kyaml [paths...]
              julia --project=no scripts/kyaml/KYAML.jl --to-yaml  [paths...]
              julia --project=no scripts/kyaml/KYAML.jl --check    [paths...]

            Options:
              --report              print the decisions taken per file
              --skip-file PATH      --check-only exemptions (default config/kyaml/drift.txt)
              --expect kyaml|yaml   style --check compares against (default kyaml)

            Exit codes: 0 ok, 1 a check failed, 2 a file was refused (nothing was written).
            """)
            return 0
        else
            push!(paths, a)
        end
        i += 1
    end
    mode == "" && (println(stderr, "KYAML.jl: pass --to-kyaml, --to-yaml or --check"); return 2)

    skip = String[]
    if isfile(skip_file)
        for line in eachline(skip_file)
            t = strip(line)
            (isempty(t) || startswith(t, "#")) && continue
            push!(skip, t)
        end
    end
    isempty(paths) && (paths = KYAML.git_yaml_paths())
    # Drift exemptions are a gate policy, not a migration boundary. Bot-owned
    # workflows are still converted (and can be rolled back); only --check omits
    # them, so an external formatter cannot make the canonicality gate flap.
    skipped = String[]
    kept = paths
    if mode == "check"
        skipped = [p for p in paths if any(s -> startswith(p, s), skip)]
        kept = [p for p in paths if !any(s -> startswith(p, s), skip)]
    end

    failures = String[]
    reports = Dict{String,KYAML.Stats}()
    rendered = Dict{String,String}()
    for p in kept
        stats = KYAML.Stats()
        try
            src = read(p, String)
            doc = KYAML.parse_document(p, src; stats = stats)
            out = KYAML.render_kyaml(doc; stats = stats)
            reports[p] = stats
            if mode == "check"
                target = expect == :yaml ? KYAML.render_yaml(doc; stats = stats) : out
                src == target || push!(failures, p)
            elseif mode == "to-yaml"
                rendered[p] = KYAML.render_yaml(doc; stats = stats)
            else
                rendered[p] = out
            end
        catch err
            err isa KYAML.KyamlError || rethrow()
            println(stderr, sprint(showerror, err))
            return 2
        end
    end
    if mode == "check"
        if isempty(failures)
            println("kyaml check: $(length(kept)) file(s) in canonical $expect form" *
                    (isempty(skipped) ? "" : ", $(length(skipped)) exempt"))
            return 0
        end
        println(stderr, "kyaml check: $(length(failures)) file(s) are not canonical $expect:")
        for p in failures
            println(stderr, "  ✗ ", p)
        end
        println(stderr, "  remedy: `just use-kyaml` (or `just use-yaml` to switch back)")
        return 1
    end
    written = 0
    for (p, text) in rendered
        read(p, String) == text && continue
        write(p, text)
        written += 1
    end
    println("kyaml $mode: $(written) of $(length(kept)) file(s) rewritten" *
            (isempty(skipped) ? "" : ", $(length(skipped)) exempt"))
    if report
        for p in sort(collect(keys(reports)))
            s = reports[p]
            println("  ", p, "  comments=", s.comments_read, "/", s.comments_written,
                    " keys_quoted=", s.keys_quoted,
                    " scalars_quoted=", s.scalars_quoted,
                    " nulls=", s.nulls_canonicalised,
                    " inline_moved=", s.inline_comments_moved)
        end
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_kyaml_cli(collect(String, ARGS)))
end
