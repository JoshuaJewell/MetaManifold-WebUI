# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Primers library
# Loads, normalises, serialises and validates config/primers.yml: the Forward
# and Reverse primer vocabularies and the Pairs composed from them. The
# validation rules themselves live in Validation (one source of truth shared
# with the environment validator); this module reuses them.
module PrimersLibrary

using YAML
using ..Validation

# An empty document, used when the file is absent. A malformed file raises in
# load, since a Save of an empty document would overwrite the real primers.
_empty() = Dict{String,Any}("Forward" => Dict{String,Any}(),
                            "Reverse" => Dict{String,Any}(),
                            "Pairs"   => Any[])

# Coerce a value that is stringified into the document. `string(nothing)` is the
# literal "nothing", so a YAML null (absent or `~`) becomes the empty string, which
# the rules below read as missing.
_str_or_empty(v) = isnothing(v) ? "" : string(v)

# Coerce a name-to-sequence map's keys to String. A numeric primer name parses
# as an Int64 under YAML; a later name-keyed lookup would dangle otherwise.
_strmap(m) = m isa AbstractDict ?
    Dict{String,Any}(_str_or_empty(k) => v for (k, v) in m) : Dict{String,Any}()

_empty_pair() = Dict{String,Any}("name" => "", "forward" => "", "reverse" => "")

# Convert one native pair entry into the canonical flat shape, returning EVERY
# pair it carries. An entry is meant to be a single-key mapping name => [fwd,
# rev], and a well-formed one yields exactly one pair.
#
# A multi-key entry is malformed: one missing "- " in primers.yml fuses the next
# pair into the same mapping. Every key is returned, so the user sees both pairs
# and a Save rewrites them one mapping per pair. Validation.primer_document_errors
# also iterates every key.
#
# Missing members become empty strings, so validate flags a malformed entry.
function _flatten_pairs(entry)
    entry isa AbstractDict || return Any[_empty_pair()]
    out = Any[]
    for (name, members) in entry
        fwd = members isa AbstractVector && length(members) >= 1 ? _str_or_empty(members[1]) : ""
        rev = members isa AbstractVector && length(members) >= 2 ? _str_or_empty(members[2]) : ""
        push!(out, Dict{String,Any}("name" => _str_or_empty(name), "forward" => fwd, "reverse" => rev))
    end
    isempty(out) && return Any[_empty_pair()]
    # The keys of one fused entry are sorted for a stable order, since a mapping
    # has none. A well-formed entry has one key; the order of Pairs is untouched.
    sort!(out; by = p -> p["name"])
    out
end

# Normalise a raw/native document (as parsed from YAML) to the canonical shape.
function normalise(raw::AbstractDict)
    pairs_raw = get(raw, "Pairs", Any[])
    pairs = Any[]
    if pairs_raw isa AbstractVector
        for e in pairs_raw
            append!(pairs, _flatten_pairs(e))
        end
    end
    Dict{String,Any}(
        "Forward" => _strmap(get(raw, "Forward", Dict{String,Any}())),
        "Reverse" => _strmap(get(raw, "Reverse", Dict{String,Any}())),
        "Pairs"   => pairs,
    )
end

# The sections a primers document carries, and the shape each must have. These
# are exactly the sections validate demands of a document, so a file that fails
# this check is one no Save could ever have written.
const _SECTIONS = (("Forward", AbstractDict, "a mapping of name to sequence"),
                   ("Reverse", AbstractDict, "a mapping of name to sequence"),
                   ("Pairs",   AbstractVector, "a list"))

# Why `raw` is not a primers document, or nothing when it is one.
#
# Presence is checked as well as type. `get` defaults only an ABSENT key, so a
# section misspelt (`Pares:`) or nulled (`Pairs: ~`) would otherwise read as no
# pairs at all.
#
# An empty-but-present section is legitimate and must stay so: `Pairs: []` and
# `Forward: {}` say "no pairs" and "no primers" unambiguously, and deleting
# everything is a real edit. Only absence and the wrong type are refused.
function _document_error(raw)
    raw isa AbstractDict || return "it is not a YAML mapping"
    for (name, T, shape) in _SECTIONS
        haskey(raw, name) || return "it has no '$name:' section"
        raw[name] isa T || return "its '$name:' section is not $shape"
    end
    nothing
end

# Load the document. A missing file yields an empty document (no primers defined).
#
# A file that exists but is not a primers document raises, whether it is
# unparseable YAML or structurally wrong. Degrading would be silent data loss: the editor would
# render the file as "no primers", an empty document breaks no validation rule,
# and the user's next Save would overwrite their real primers with nothing. A
# caller that must not throw catches this at its own boundary and reports the
# file as unreadable; see the primers routes, which turn it into a 400 naming
# the file.
#
# The structural case needs the guard as much as the unparseable one, and only
# the write gate's REQUEST path was ever covered: a malformed section submitted
# by a client is refused by validate, but the same shape arriving from the FILE
# was laundered into a valid empty document before the gate could see it, so the
# gate never fired.
function load(path::String)
    isfile(path) || return _empty()
    raw = YAML.load_file(path)
    problem = _document_error(raw)
    isnothing(problem) || error("$path is not a primers document: $problem")
    normalise(raw)
end

# Convert the canonical shape back to the native file shape for YAML.write:
# Pairs becomes a list of single-key mappings name => [forward, reverse].
function to_yaml_doc(doc::AbstractDict)
    pairs = Any[]
    for p in Validation._seq(get(doc, "Pairs", nothing))
        p isa AbstractDict || continue
        name = _str_or_empty(get(p, "name", ""))
        fwd  = _str_or_empty(get(p, "forward", ""))
        rev  = _str_or_empty(get(p, "reverse", ""))
        push!(pairs, Dict{String,Any}(name => Any[fwd, rev]))
    end
    Dict{String,Any}(
        "Forward" => _strmap(get(doc, "Forward", Dict{String,Any}())),
        "Reverse" => _strmap(get(doc, "Reverse", Dict{String,Any}())),
        "Pairs"   => pairs,
    )
end

# Pair names in document order.
function pair_names(doc::AbstractDict)
    names = String[]
    for p in Validation._seq(get(doc, "Pairs", nothing))
        p isa AbstractDict && push!(names, _str_or_empty(get(p, "name", "")))
    end
    names
end

# Validate a canonical document. Delegates to the shared rule function so the
# environment validator and this module cannot disagree on what a valid primers
# document is. Never throws.
#
# The rules below are added on top, and only here. Each is a WRITE-time rule, not
# a file-validity rule: each is meaningless or destructive in a new document, but
# adding it to the shared rules would fail a primers.yml that validated
# yesterday, so it is applied where new documents are written and nowhere else.
#
# `Pairs` is checked first, and a malformed `Pairs` returns immediately instead of
# falling through to _seq, which only guards iteration over a null inner sequence.
# An empty document breaks no validation rule, so coercing a malformed section to
# empty would let the write gate accept a document that wipes every pair on save.
# Deleting every pair is spelled as the empty list. `Forward` and `Reverse` are
# checked the same way.
function validate(doc::AbstractDict; native::AbstractDict=to_yaml_doc(doc))
    errors = String[]
    for section in ("Forward", "Reverse")
        sect = get(doc, section, nothing)
        sect isa AbstractDict ||
            push!(errors, "the $section section must be a mapping")
    end

    pairs = get(doc, "Pairs", nothing)
    pairs isa AbstractVector || begin
        push!(errors, "the Pairs section must be a list")
        return errors
    end

    append!(errors, Validation.primer_document_errors(native))
    for section in ("Forward", "Reverse")
        primers = get(doc, section, Dict{String,Any}())
        primers isa AbstractDict || continue
        for (name, seq) in primers
            seq isa AbstractString && isempty(strip(seq)) &&
                push!(errors, "primer '$name' has an empty sequence")
        end
    end
    errors
end

end # module
