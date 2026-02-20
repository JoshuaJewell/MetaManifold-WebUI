module Databases

# Ensures all databases declared in the config are available locally,
# downloading from their configured URIs as needed.
#
# Usage:
#   dbs = ensure_databases("config/dada2.yml")
#   dbs["pr2_dada2"]   # → resolved local path for DADA2 assignTaxonomy
#   dbs["pr2_vsearch"] # → resolved local path for VSEARCH --db
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

export ensure_databases

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

function _resolve_entry(key, fmt_info, db_dir)
    local_p = get(fmt_info, "local", nothing)
    if !isnothing(local_p)
        local_p = string(local_p)
        if !isempty(local_p)
            if isfile(local_p)
                @info "[$key] Using local file: $local_p"
                return local_p
            end
            @warn "[$key] Configured local path not found: $local_p — falling back to uri"
        end
    end

    uri = get(fmt_info, "uri", nothing)
    isnothing(uri) &&
        error("databases entry '$key': no valid local path and no uri is configured")

    cached = joinpath(db_dir, basename(uri))
    if isfile(cached)
        @info "[$key] Using cached: $cached"
    else
        @info "[$key] Downloading: $uri"
        Downloads.download(uri, cached)
        @info "[$key] Saved to: $cached"
    end
    return cached
end

end