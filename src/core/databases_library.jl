# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Databases library
# Loads, normalises, serialises and validates config/databases.yml: the shared
# cache directory and, per reference database, its source URIs, taxonomy rank
# list, vsearch parser format and taxonomy corrections. The structural rules live
# in Validation (one source of truth shared with the environment validator); this
# module reuses them and adds the write-time rules. Whether a configured file is
# presently on disk is neither: it is an environment rule, and it lives with the
# environment validator alone.
module DatabasesLibrary

using YAML
using ..Validation

# The formats a database entry carries. dada2 additionally supports remote_path,
# a pre-existing path on the remote taxonomy host; vsearch has no remote stage.
const _FORMATS = ("dada2", "vsearch")

# The only two values vsearch_format selects between. make_db_meta defaults an
# unrecognised value to "generic", and only the literal "pr2" selects
# pipe-separated parsing, so this is a switch and not an open vocabulary.
const _VSEARCH_FORMATS = ("pr2", "generic")

# The schemes a uri may name. Downloads.download is libcurl-backed and would
# otherwise honour file:// and its kin, turning a saved uri into a read of any
# file the server can reach.
_has_fetchable_scheme(uri::AbstractString) =
    startswith(uri, "http://") || startswith(uri, "https://")

# The remote-path gate lives in Validation, shared with the remote stage runner
# and the pipeline config check so that one answer covers every path that
# reaches an ssh command string.
const _is_shell_safe = Validation.is_shell_safe

# The shared cache directory a document that names none falls back to, and the
# vsearch parser a database that names none is read with. Both are real defaults
# that the rest of the system already applies, so both an ABSENT and a null
# spelling of the field mean exactly this value.
const _DEFAULT_DIR = "./databases"
const _DEFAULT_VSEARCH_FORMAT = "generic"

# An empty document, used when the file is absent.
_empty() = Dict{String,Any}("dir" => _DEFAULT_DIR, "databases" => Any[])

# Coerce a value that is stringified into the document. `string(nothing)` is the
# literal "nothing", so a YAML null (absent or `~`) becomes the empty string, which
# the rules below read as missing.
_str_or_empty(v) = isnothing(v) ? "" : string(v)

# `dir` and `vsearch_format` have real defaults, so a null means the default.
# `string` would turn `dir: ~` into the directory "./nothing".
_dir_or_default(v) = isnothing(v) ? _DEFAULT_DIR : string(v)
_vsearch_format_or_default(v) = isnothing(v) ? _DEFAULT_VSEARCH_FORMAT : string(v)

# A database entry's key, trimmed, so validate and to_yaml_doc see the same key
# ("pr2 " would otherwise pass validation as "pr2" and be written as `"pr2 ":`).
# The duplicate rule compares trimmed keys, so it catches collisions trimming creates.
_entry_key(e) = strip(_str_or_empty(get(e, "key", "")))

_str_or_nothing(v) = isnothing(v) ? nothing : string(v)

# Convert one native format entry into the canonical shape, keeping absent and
# null alike as nothing so the editor renders an empty field either way.
function _flatten_format(raw, format::String)
    raw isa AbstractDict || return Dict{String,Any}("uri" => "", "local" => nothing)
    out = Dict{String,Any}(
        "uri"   => _str_or_empty(get(raw, "uri", "")),
        "local" => _str_or_nothing(get(raw, "local", nothing)),
    )
    format == "dada2" && (out["remote_path"] = _str_or_nothing(get(raw, "remote_path", nothing)))
    out
end

# Corrections values are a native mapping of source value => target value. They
# become ordered from/to rows for the same reason Pairs did: a map keyed by
# edited text collapses two rows the moment one key is typed into the other.
function _flatten_values(raw)
    raw isa AbstractDict || return Any[]
    Any[Dict{String,Any}("from" => _str_or_empty(k), "to" => _str_or_empty(v)) for (k, v) in raw]
end

function _flatten_correction(raw)
    raw isa AbstractDict || return Dict{String,Any}("source" => "", "target" => "", "values" => Any[])
    Dict{String,Any}(
        "source" => _str_or_empty(get(raw, "source", "")),
        "target" => _str_or_empty(get(raw, "target", "")),
        "values" => _flatten_values(get(raw, "values", nothing)),
    )
end

function _flatten_entry(key, raw)
    raw isa AbstractDict || return Dict{String,Any}(
        "key" => _str_or_empty(key), "label" => "",
        "dada2" => _flatten_format(nothing, "dada2"), "vsearch" => _flatten_format(nothing, "vsearch"),
        "levels" => Any[], "vsearch_format" => _DEFAULT_VSEARCH_FORMAT, "corrections" => Any[])
    levels_raw = get(raw, "levels", Any[])
    corr_raw   = get(raw, "corrections", Any[])
    Dict{String,Any}(
        "key"            => _str_or_empty(key),
        # A null label is the empty string, which to_yaml_doc omits.
        "label"          => _str_or_empty(get(raw, "label", "")),
        "dada2"          => _flatten_format(get(raw, "dada2",   nothing), "dada2"),
        "vsearch"        => _flatten_format(get(raw, "vsearch", nothing), "vsearch"),
        # A null level becomes a blank one, which validate's blank-level rule rejects.
        "levels"         => levels_raw isa AbstractVector ? Any[_str_or_empty(l) for l in levels_raw] : Any[],
        "vsearch_format" => _vsearch_format_or_default(get(raw, "vsearch_format", _DEFAULT_VSEARCH_FORMAT)),
        "corrections"    => corr_raw isa AbstractVector ? Any[_flatten_correction(c) for c in corr_raw] : Any[],
    )
end

# Normalise a raw/native document (as parsed from YAML) to the canonical shape.
# `dir` lifts out of the databases: namespace, where it sits alongside the
# entries on disk, so no database may be called dir.
#
# Entries are sorted by key, since a YAML mapping has no meaningful order, so the
# editor shows a stable order.
function normalise(raw::AbstractDict)
    dbs = get(raw, "databases", nothing)
    dbs isa AbstractDict || return _empty()
    entries = Any[]
    for (k, v) in dbs
        string(k) == "dir" && continue
        push!(entries, _flatten_entry(k, v))
    end
    sort!(entries; by = e -> e["key"])
    Dict{String,Any}(
        "dir"       => _dir_or_default(get(dbs, "dir", _DEFAULT_DIR)),
        "databases" => entries,
    )
end

# Why `raw` is not a databases document, or nothing when it is one.
#
# Presence is checked as well as type. `get` defaults only an ABSENT key, so a
# section misspelt (`datbases:`) or nulled (`databases: ~`) would otherwise read as
# no databases at all.
#
# An empty-but-present section is legitimate and must stay so: `databases: {}`,
# and a `databases:` carrying only `dir`, say "no databases" unambiguously, and
# deleting every entry is a real edit. Only absence and the wrong type are
# refused. `dir` is not required: it has a real default, so a document that omits
# it is complete.
function _document_error(raw)
    raw isa AbstractDict || return "it is not a YAML mapping"
    haskey(raw, "databases") || return "it has no 'databases:' section"
    raw["databases"] isa AbstractDict || return "its 'databases:' section is not a mapping"
    nothing
end

# Load the document. A missing file yields an empty document (no databases defined).
#
# A file that exists but is not a databases document raises, whether it is
# unparseable YAML or structurally wrong. Degrading would be silent data loss: the editor would
# render the file as "no databases", an empty document breaks no validation rule,
# and the user's next Save would overwrite their real config with nothing. A
# caller that must not throw catches this at its own boundary and reports the
# file as unreadable; see the databases routes, which turn it into a 400 naming
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
    isnothing(problem) || error("$path is not a databases document: $problem")
    normalise(raw)
end

_native_format(fmt::AbstractDict, format::String) = begin
    out = Dict{String,Any}("uri" => _str_or_empty(get(fmt, "uri", "")),
                           "local" => get(fmt, "local", nothing))
    format == "dada2" && (out["remote_path"] = get(fmt, "remote_path", nothing))
    out
end

_native_values(rows) = Dict{String,Any}(
    _str_or_empty(get(r, "from", "")) => _str_or_empty(get(r, "to", ""))
    for r in (rows isa AbstractVector ? rows : Any[]) if r isa AbstractDict)

# Convert the canonical shape back to the native file shape for YAML.write: dir
# returns to the databases: namespace and each entry is keyed by its name.
function to_yaml_doc(doc::AbstractDict)
    dbs = Dict{String,Any}("dir" => _dir_or_default(get(doc, "dir", _DEFAULT_DIR)))
    for e in Validation._seq(get(doc, "databases", nothing))
        e isa AbstractDict || continue
        entry = Dict{String,Any}()
        label = _str_or_empty(get(e, "label", ""))
        isempty(label) || (entry["label"] = label)
        for format in _FORMATS
            fmt = get(e, format, nothing)
            fmt isa AbstractDict && (entry[format] = _native_format(fmt, format))
        end
        levels = get(e, "levels", Any[])
        entry["levels"] = levels isa AbstractVector ? Any[_str_or_empty(l) for l in levels] : Any[]
        entry["vsearch_format"] = _vsearch_format_or_default(get(e, "vsearch_format", _DEFAULT_VSEARCH_FORMAT))
        corr = get(e, "corrections", Any[])
        entry["corrections"] = Any[
            Dict{String,Any}("source" => _str_or_empty(get(c, "source", "")),
                             "target" => _str_or_empty(get(c, "target", "")),
                             "values" => _native_values(get(c, "values", Any[])))
            for c in (corr isa AbstractVector ? corr : Any[]) if c isa AbstractDict]
        dbs[_entry_key(e)] = entry
    end
    Dict{String,Any}("databases" => dbs)
end

# Database keys in document order.
function database_keys(doc::AbstractDict)
    keys_ = String[]
    for e in Validation._seq(get(doc, "databases", nothing))
        e isa AbstractDict && push!(keys_, _entry_key(e))
    end
    keys_
end

# The taxonomy ranks of one database, or an empty vector when it is absent.
function levels_of(doc::AbstractDict, key::String)
    for e in Validation._seq(get(doc, "databases", nothing))
        e isa AbstractDict || continue
        _entry_key(e) == key || continue
        levels = get(e, "levels", Any[])
        return levels isa AbstractVector ? String[_str_or_empty(l) for l in levels] : String[]
    end
    String[]
end

# Validate a canonical document. Delegates to the shared structural rule function
# so the environment validator and this module cannot disagree on the shape of a
# databases document. Never throws.
#
# It does not ask whether a configured local: file exists. That depends on the
# machine: a user may legitimately save a
# path before the file arrives, and _validate_databases reports the absence when
# it matters. Asked here it did real harm, because a save carries the WHOLE
# document: one local: whose file had since moved rejected every save of every
# entry, including entries the user never touched, and it answered any client's
# question of whether an arbitrary path exists on the server.
#
# The rules below are added on top, and only here. Each is a WRITE-time rule, not
# a file-validity rule: each is meaningless or destructive in a new document, but
# adding it to the shared rules would fail a databases.yml that validated
# yesterday, so it is applied where new documents are written and nowhere else.
#
# `dir` and `databases` are checked first, and a malformed `databases` returns
# immediately. _seq only guards iteration over a null inner sequence; a malformed
# top-level section is an error.
# Coercing it here would be the same mistake load's docstring above refuses to
# make: an empty document breaks no validation rule, so the write gate would
# wave through a document that wipes every database on save. Deleting every
# database already has an unambiguous spelling, the empty list, so nothing
# legitimate is lost by rejecting anything else.
function validate(doc::AbstractDict; native::AbstractDict=to_yaml_doc(doc))
    errors = String[]
    # An ABSENT or null dir is not an error: normalise and to_yaml_doc both
    # default both spellings to ./databases, and `~` is how this file already
    # writes every unset field. A blank one is. normalise defaults only an
    # absent or null key, the editor always sends the field, and to_yaml_doc
    # writes the blank verbatim,
    # whereupon _db_info's abspath("") resolves to the working directory: the
    # cache silently relocates to the repository root, and the placeholder the
    # editor shows is a default the field never actually gets.
    dir = get(doc, "dir", nothing)
    if !(dir === nothing || dir isa AbstractString)
        push!(errors, "the dir must be a string")
    elseif dir isa AbstractString && isempty(strip(dir))
        push!(errors, "the dir must not be blank")
    end

    dbs = get(doc, "databases", nothing)
    dbs isa AbstractVector || begin
        push!(errors, "the databases section must be a list")
        return errors
    end

    append!(errors, Validation.database_document_errors(native))
    seen = Set{String}()
    for e in dbs
        e isa AbstractDict || continue
        # The TRIMMED key, exactly as to_yaml_doc now writes it, so the rules
        # below and the file on disk cannot describe different documents.
        key = _entry_key(e)
        if isempty(key)
            push!(errors, "a database has no name")
            continue
        end
        # `dir` is the shared cache directory and shares the entries' namespace on
        # disk, so a database of that name would silently become the cache path.
        key == "dir" && push!(errors, "a database may not be named 'dir'")
        key in seen && push!(errors, "duplicate database name '$key'")
        push!(seen, key)

        levels = get(e, "levels", Any[])
        levels = levels isa AbstractVector ? String[_str_or_empty(l) for l in levels] : String[]
        isempty(levels) && push!(errors, "database '$key' has no taxonomy levels")
        # Names are compared trimmed and each must be non-blank (a null arrives
        # as blank). Two names differing only in surrounding space are one rank
        # spelt twice. `isempty` above tests the list itself.
        trimmed = String[strip(l) for l in levels]
        any(isempty, trimmed) &&
            push!(errors, "database '$key' has a taxonomy level with no name")
        length(unique(trimmed)) == length(trimmed) ||
            push!(errors, "database '$key' has duplicate taxonomy levels")

        fmt = _str_or_empty(get(e, "vsearch_format", ""))
        fmt in _VSEARCH_FORMATS ||
            push!(errors, "database '$key' has vsearch_format '$fmt'; expected one of $(join(_VSEARCH_FORMATS, ", "))")

        for format in _FORMATS
            f = get(e, format, nothing)
            f isa AbstractDict || begin
                push!(errors, "database '$key' has no $format section")
                continue
            end
            uri   = strip(_str_or_empty(get(f, "uri", "")))
            local_ = get(f, "local", nothing)
            (isempty(uri) && (isnothing(local_) || isempty(strip(string(local_))))) &&
                push!(errors, "database '$key' $format has neither a uri nor a local path")
            # A uri names something to FETCH, and the download route hands it to
            # Downloads.download, which is libcurl-backed and so honours file://
            # among others. That would make a save able to read any file the
            # server can read and cache it under the databases directory, whence
            # the files route can serve it back. A pre-downloaded file already
            # has its own spelling here, `local:`, so restricting the scheme
            # costs nothing legitimate.
            isempty(uri) || _has_fetchable_scheme(uri) ||
                push!(errors, "database '$key' $format uri must begin with http:// or https://")
            # remote_path is interpolated into the command string that
            # _assign_taxonomy_remote hands to ssh, which the remote sshd runs
            # through a login shell, so shell metacharacters are refused.
            if format == "dada2"
                rp = get(f, "remote_path", nothing)
                isnothing(rp) || _is_shell_safe(string(rp)) ||
                    push!(errors, "database '$key' remote_path may not contain shell metacharacters")
            end
        end

        for c in Validation._seq(get(e, "corrections", nothing))
            c isa AbstractDict || continue
            for side in ("source", "target")
                name = _str_or_empty(get(c, side, ""))
                name in levels ||
                    push!(errors, "database '$key' has a correction whose $side '$name' is not one of its levels")
            end

            # _native_values re-keys these rows by `from` on the way to disk, so the
            # same hazard that drove the ordered-rows shape reappears here: a blank
            # `from` has nowhere to go and two rows sharing one collapse into a
            # single entry, silently discarding a correction.
            seen_from = Set{String}()
            for v in Validation._seq(get(c, "values", nothing))
                v isa AbstractDict || continue
                from = strip(_str_or_empty(get(v, "from", "")))
                if isempty(from)
                    push!(errors, "database '$key' has a correction value with no source label")
                    continue
                end
                from in seen_from &&
                    push!(errors, "database '$key' has a correction with duplicate source label '$from'")
                push!(seen_from, from)
                # A correction RENAMES `from` to `to`, so a target carrying no
                # name renames the taxon to nothing: the same defect the blank
                # level rule above catches, one field along. A null is the same
                # thing differently spelt and arrives here as a blank, so one
                # rule covers both; left alone it reached disk as the literal
                # target label "nothing".
                to = strip(_str_or_empty(get(v, "to", "")))
                isempty(to) &&
                    push!(errors, "database '$key' has a correction value '$from' with no target label")
            end
        end
    end
    errors
end

end # module
