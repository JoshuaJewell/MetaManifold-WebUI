module Databases

# Ensures all databases declared in the config are available locally,
# downloading from their configured URIs as needed.
#
# Usage:
#   dbs = ensure_databases("config/dada2.yml")
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
using YAML, Logging, CodecZlib
using ..PipelineTypes

export ensure_databases, resolve_db, make_db_meta

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
    function ensure_databases(config_path::String)
        if !isfile(config_path)
            @warn "ensure_databases: $config_path not found. " *
                "Copy config/databases.example.yml to config/databases.yml and set local: paths."
            return Dict{String,String}()
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
        "SeqName", "Pident", "Accession", "rRNA", "Organellum", "specimen",
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
        local_p = get(fmt_info, "local", nothing)
        if !isnothing(local_p)
            local_p = string(local_p)
            if !isempty(local_p)
                if isfile(local_p)
                    log("[$key] Using local file: $local_p")
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
            Downloads.download(uri, cached)
            log("[$key] Saved to: $cached")
        end

        # Post-processing: reformat database if requested.
        reformat = get(fmt_info, "reformat", nothing)
        if !isnothing(reformat) && !isempty(string(reformat))
            cached = _apply_reformat(key, string(reformat), cached, db_dir; log)
        end

        return cached
    end

    """
    Apply a named reformat step to a downloaded database file.
    Returns the path to the (possibly new) reformatted file.
    """
    function _apply_reformat(key, reformat, src_path, db_dir; log=msg->@info(msg))
        if reformat == "silva_vsearch"
            return _reformat_silva_vsearch(key, src_path, db_dir; log)
        else
            @warn "[$key] Unknown reformat '$reformat' - skipping"
            return src_path
        end
    end

    """
    Reformat a SILVA vsearch FASTA so taxonomy is embedded in the sequence ID.

    SILVA headers are `>Accession Kingdom;Phylum;...;Genus` (taxonomy after a space).
    vsearch only captures the sequence ID (before the first space), so taxonomy is lost.

    This function rewrites headers to `>Accession;Kingdom;Phylum;...;Genus` with spaces
    within taxon names replaced by underscores, so vsearch returns the full taxonomy
    string in the `target` field.

    The reformatted file is cached alongside the original; the original is kept intact.
    """
    function _reformat_silva_vsearch(key, src_path, db_dir; log=msg->@info(msg))
        # Derive output filename from source
        src_base = basename(src_path)
        # Strip .gz if present to insert _reformatted before the extension
        if endswith(src_base, ".fasta.gz")
            out_base = src_base[1:end-9] * "_reformatted.fasta.gz"
        elseif endswith(src_base, ".fa.gz")
            out_base = src_base[1:end-6] * "_reformatted.fa.gz"
        elseif endswith(src_base, ".fasta")
            out_base = src_base[1:end-6] * "_reformatted.fasta"
        else
            out_base = src_base * "_reformatted"
        end
        out_path = joinpath(db_dir, out_base)

        if isfile(out_path)
            log("[$key] Using cached reformatted SILVA: $out_path")
            return out_path
        end

        log("[$key] Reformatting SILVA FASTA for vsearch compatibility: $src_path -> $out_path")

        # Determine if gzipped
        is_gz = endswith(src_path, ".gz")
        is_out_gz = endswith(out_path, ".gz")

        open_in  = is_gz     ? () -> GzipDecompressorStream(open(src_path, "r"))  : () -> open(src_path, "r")
        open_out = is_out_gz ? () -> GzipCompressorStream(open(out_path, "w"))    : () -> open(out_path, "w")

        in_io  = open_in()
        out_io = open_out()
        n_seq = 0
        try
            for line in eachline(in_io)
                if startswith(line, '>')
                    # Header: ">Accession.start.end Kingdom;Phylum;...;Genus"
                    rest = line[2:end]
                    sp = findfirst(' ', rest)
                    if isnothing(sp)
                        # No space - already no taxonomy; pass through
                        write(out_io, line, '\n')
                    else
                        acc = rest[1:sp-1]
                        tax = rest[sp+1:end]
                        # Replace spaces within taxonomy with underscores,
                        # then join accession + taxonomy with semicolon.
                        tax_clean = replace(tax, ' ' => '_')
                        write(out_io, '>', acc, ';', tax_clean, '\n')
                    end
                    n_seq += 1
                else
                    write(out_io, line, '\n')
                end
            end
        finally
            close(out_io)
            close(in_io)
        end

        log("[$key] Reformatted $n_seq sequences -> $out_path")
        return out_path
    end

end