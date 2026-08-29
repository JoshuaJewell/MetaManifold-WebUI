module Databases

# Ensures all databases declared in the config are available locally,
# downloading from their configured URIs as needed.
#
# Usage:
#   dbs = ensure_databases("config/databases.yml")
#   dbs["pr2_dada2"]   # -> resolved local path for DADA2 assignTaxonomy
#   dbs["pr2_vsearch"] # -> resolved local path for VSEARCH --db
#
# Keys follow the pattern  "<database_name>_<format>",  e.g. "pr2_dada2",
# "pr2_vsearch".  Format names map to the sub-keys under each database entry
# in the databases: config section.
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

import Downloads
using YAML, Logging, CodecZlib, SHA, OrderedCollections
using ..PipelineTypes

export ensure_databases, resolve_db, make_db_meta, verify_db_file

    ## Reference file integrity
    #
    # A reference database that arrives incomplete does not fail loudly: R's
    # gzip reader hands back the readable prefix of a truncated stream, so
    # assignTaxonomy trains on whatever survived and reports confident-looking
    # bootstraps against an amputated reference set. Nothing downstream can
    # detect that. So a file is not usable here until it has been shown to
    # decompress to its end.

    # Stream the whole gzip member through the decompressor. Returns
    # (ok, uncompressed_bytes, message); a truncated stream throws inside
    # CodecZlib and is reported in the message.
    function _gzip_intact(path::AbstractString)
        total = 0
        try
            open(path, "r") do raw
                stream = GzipDecompressorStream(raw)
                try
                    buf = Vector{UInt8}(undef, 1 << 20)
                    while !eof(stream)
                        n = readbytes!(stream, buf)
                        total += n
                    end
                finally
                    close(stream)
                end
            end
        catch err
            return (false, total, sprint(showerror, err))
        end
        return (true, total, "")
    end

    _sha256_of(path::AbstractString) = open(io -> bytes2hex(sha256(io)), path, "r")

    # Cache the verdict beside the file, keyed by size and mtime, so a
    # multi-hundred-megabyte reference is not re-read on every pipeline run.
    # Any change to the file misses the cache and forces a re-check.
    _verify_sidecar(path::AbstractString) = path * ".verified.yml"

    # `checks` records which expectations the stored verdict actually covered, so
    # adding a sha256 to databases.yml after a file was verified without one
    # forces a re-check.
    function _cached_verdict(path::AbstractString, checks::AbstractString)
        sidecar = _verify_sidecar(path)
        isfile(sidecar) || return nothing
        st = stat(path)
        try
            rec = YAML.load_file(sidecar)
            get(rec, "size", nothing) == Int(st.size) || return nothing
            get(rec, "mtime_us", nothing) == round(Int, st.mtime * 1_000_000) || return nothing
            get(rec, "checks", nothing) == checks || return nothing
            return get(rec, "ok", nothing) === true
        catch err
            @warn "Databases: ignoring unreadable verification sidecar $sidecar" exception=err
            return nothing
        end
    end

    function _store_verdict(path::AbstractString, ok::Bool, detail::AbstractString,
                            checks::AbstractString)
        st = stat(path)
        try
            # Ordered so the sidecar is byte-identical for an identical verdict.
            YAML.write_file(_verify_sidecar(path), OrderedDict{String,Any}(
                "file"     => basename(path),
                "ok"       => ok,
                "detail"   => detail,
                "checks"   => checks,
                "size"     => Int(st.size),
                "mtime_us" => round(Int, st.mtime * 1_000_000),
            ))
        catch err
            @warn "Databases: could not write verification sidecar" path exception=err
        end
    end

    """
        verify_db_file(path; expected_sha256=nothing, expected_size=nothing,
                       key="", recheck=false) -> Nothing

    Throw unless `path` is a complete, usable reference file.

    Checks, in order: the file exists and is non-empty; its size matches
    `expected_size` when one is configured; a `.gz` file decompresses cleanly to
    the end of the stream; and its SHA-256 matches `expected_sha256` when one is
    configured. The verdict is cached in a sidecar keyed by size and mtime, so
    the expensive checks run once per version of the file rather than per run.
    """
    function verify_db_file(path::AbstractString;
                            expected_sha256=nothing, expected_size=nothing,
                            key::AbstractString="", recheck::Bool=false)
        tag = isempty(key) ? basename(path) : key
        isfile(path) || error("Database '$tag': file not found: $path")
        sz = filesize(path)
        sz > 0 || error("Database '$tag': file is empty: $path")

        if !isnothing(expected_size) && sz != Int(expected_size)
            error("Database '$tag': size mismatch for $path - expected " *
                  "$(expected_size) bytes, found $sz. The file is incomplete or " *
                  "is not the configured release.")
        end

        is_gz = endswith(lowercase(path), ".gz")
        # The key records the expected *values* as well as which checks ran: a
        # digest changed in databases.yml must invalidate a verdict reached
        # against the old one. (`size` is compared above, ahead of the cache.)
        checks = (is_gz ? "gzip" : "") *
                 (isnothing(expected_sha256) ? "" : "|sha256=" * String(expected_sha256))

        # A cached pass short-circuits the expensive checks; a cached failure is
        # always re-run, so a repaired file is picked up without hand-editing.
        if !recheck && _cached_verdict(path, checks) === true
            return nothing
        end

        if is_gz
            ok, nbytes, msg = _gzip_intact(path)
            if !ok
                _store_verdict(path, false, msg, checks)
                error("Database '$tag': $path is a truncated or corrupt gzip stream " *
                      "($(nbytes) bytes decompressed before it failed: $msg). " *
                      "Delete the file and let it download again.")
            end
        end

        if !isnothing(expected_sha256)
            actual = _sha256_of(path)
            if actual != String(expected_sha256)
                _store_verdict(path, false, "sha256 $actual", checks)
                error("Database '$tag': SHA-256 mismatch for $path - expected " *
                      "$(expected_sha256), found $actual.")
            end
        end

        _store_verdict(path, true, "verified", checks)
        return nothing
    end

    """
        ensure_databases(config_path) -> Dict{String,String}

    Reads the `databases:` section of `config_path`, ensures every declared
    database file is present in `databases.dir` (downloading from `uri` if the
    file is absent), and returns a Dict mapping `"<name>_<format>"` keys to
    resolved absolute local paths.

    Set a `local:` path under any entry to use a pre-existing file directly.
    If the `local:` path does not exist, the function warns and falls back to
    downloading from `uri`.
    """
    function ensure_databases(config_path::String; only::Union{Nothing,AbstractSet}=nothing)
        if !isfile(config_path)
            template_path = joinpath(dirname(config_path), "defaults", "databases.yml")
            if isfile(template_path)
                cp(template_path, config_path)
                @warn "ensure_databases: $config_path not found - copied from $template_path. " *
                    "Set local: paths for any pre-downloaded databases."
            else
                @warn "ensure_databases: $config_path not found. " *
                    "Copy config/defaults/databases.yml to config/databases.yml and set local: paths."
                return Dict{String,String}()
            end
        end

        cfg    = YAML.load_file(config_path)
        db_cfg = get(cfg, "databases", nothing)

        if isnothing(db_cfg) || isempty(db_cfg)
            @warn "ensure_databases: no `databases:` section found in $config_path"
            return Dict{String,String}()
        end

        db_dir = abspath(get(db_cfg, "dir", "./databases"))
        mkpath(db_dir)

        resolved = Dict{String,String}()
        for (db_name, db_info) in db_cfg
            db_name == "dir" && continue
            isnothing(only) || db_name in only || continue
            !(db_info isa AbstractDict) && continue
            for (fmt, fmt_info) in db_info
                !(fmt_info isa AbstractDict) && continue
                # Skip entries where both uri and local are null/missing.
                uri_val   = get(fmt_info, "uri", nothing)
                local_val = get(fmt_info, "local", nothing)
                if isnothing(uri_val) && isnothing(local_val)
                    continue
                end
                key = "$(db_name)_$(fmt)"
                resolved[key] = _resolve_entry(key, fmt_info, db_dir)
            end
        end
        return resolved
    end

    """
        resolve_db(config_path, db_name, fmt; emit=nothing) -> String

    Resolve a single database entry from a databases.yml file.

    `config_path` is the path to a databases.yml-style config, `db_name` is the
    key under `databases:` (e.g., `"pr2"`), and `fmt` is the format sub-key
    (e.g., `"dada2"` or `"vsearch"`). Returns the resolved absolute local path.

    Pass an `emit` function (e.g., from `_emitter`) to route log messages through
    the pipeline's progress channel instead of the default `@info` logger.
    """
    function resolve_db(config_path::String, db_name::String, fmt::String; emit=nothing)
        cfg    = YAML.load_file(config_path)
        db_cfg = get(cfg, "databases", Dict())
        db_dir = abspath(get(db_cfg, "dir", "./databases"))
        mkpath(db_dir)

        haskey(db_cfg, db_name) ||
            error("Database '$db_name' not found in $config_path")
        fmt_cfg = get(db_cfg[db_name], fmt, nothing)
        isnothing(fmt_cfg) &&
            error("databases.$db_name.$fmt is not configured in $config_path")

        _resolve_entry("$(db_name)_$(fmt)", fmt_cfg, db_dir; emit)
    end

    # Fixed column names that are never sample counts, regardless of database.
    const _FIXED_NONCOUNTS = Set([
        "SeqName", "ASV", "ASVs", "OTU",
        "Pident", "Accession", "rRNA", "Organellum", "specimen",
        "Sequence", "sequence", "", "Column1",
    ])

    """
        make_db_meta(config_path, db_name) -> DatabaseMeta

    Construct a `DatabaseMeta` from the database entry in `config_path`.
    Pre-computes `noncounts` as the union of taxonomy levels (including
    `_dada2` suffixed variants) and fixed metadata column names.
    """
    function make_db_meta(config_path::String, db_name::String)
        cfg    = YAML.load_file(config_path)
        db_cfg = get(cfg, "databases", Dict())
        haskey(db_cfg, db_name) ||
            error("Database '$db_name' not found in $config_path")
        entry = db_cfg[db_name]

        levels     = String[string(l) for l in get(entry, "levels", String[])]
        vsformat   = string(get(entry, "vsearch_format", "generic"))
        raw_corr   = get(entry, "corrections", [])
        corrections = Dict{String,Any}[]
        if raw_corr isa Vector
            for c in raw_corr
                c isa AbstractDict && push!(corrections, Dict{String,Any}(string(k) => v for (k,v) in c))
            end
        end

        noncounts = copy(_FIXED_NONCOUNTS)
        for l in levels
            push!(noncounts, l)
            push!(noncounts, l * "_dada2")
            push!(noncounts, l * "_vsearch")
            push!(noncounts, l * "_boot")
        end

        return DatabaseMeta(db_name, levels, vsformat, corrections, noncounts)
    end

    function _resolve_entry(key, fmt_info, db_dir; emit=nothing)
        log = isnothing(emit) ? msg -> @info(msg) : emit
        want_sha  = get(fmt_info, "sha256", nothing)
        want_size = get(fmt_info, "size", nothing)
        want_sha  = isnothing(want_sha)  ? nothing : String(string(want_sha))
        want_size = isnothing(want_size) ? nothing : Int(want_size)
        verify(p) = verify_db_file(p; expected_sha256=want_sha,
                                   expected_size=want_size, key=string(key))

        local_p = get(fmt_info, "local", nothing)
        if !isnothing(local_p)
            local_p = string(local_p)
            if !isempty(local_p)
                if isfile(local_p)
                    log("[$key] Using local file: $local_p")
                    # A hand-placed file gets the same scrutiny as a downloaded
                    # one: it can be truncated too, and nothing downstream notices.
                    verify(local_p)
                    return local_p
                end
                @warn "[$key] Configured local path not found: $local_p - falling back to uri"
            end
        end

        uri = get(fmt_info, "uri", nothing)
        isnothing(uri) &&
            error("databases entry '$key': no valid local path and no uri is configured")

        cached = joinpath(db_dir, basename(uri))
        if isfile(cached)
            log("[$key] Using cached: $cached")
        else
            log("[$key] Downloading: $uri")
            # Download to a temporary name and move it into place only once it
            # verifies. Writing straight to `cached` means an interrupted
            # transfer leaves a partial file under the real name, and every later
            # run takes the `isfile` branch above and calls it cached.
            # A unique name, so two jobs fetching the same file cannot share a partial download.
            part = tempname(dirname(cached)) * ".part"
            try
                Downloads.download(string(uri), part)
                verify_db_file(part; expected_sha256=want_sha, expected_size=want_size,
                               key=string(key))
                mv(part, cached; force=true)
            finally
                # The verdict written while checking the temporary file would
                # otherwise be orphaned beside it; the real file gets its own
                # below.
                isfile(part) && rm(part; force=true)
                isfile(_verify_sidecar(part)) && rm(_verify_sidecar(part); force=true)
            end
            log("[$key] Saved to: $cached")
        end

        # Verify on every resolve: the file may
        # predate this check, or have been damaged since.
        verify(cached)
        return cached
    end

end