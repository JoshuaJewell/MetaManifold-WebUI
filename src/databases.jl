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
using YAML, Logging

export ensure_databases, resolve_db

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
        return cached
    end

end