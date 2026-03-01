# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/studies/{study}/config
#         /api/v1/studies/{study}/runs/{run}/config

using JSON3, YAML

## Config resolution

# Returns the merged config cascade for a run as a flat dict of dotted keys,
# each annotated with its effective value and source level.
#
# Cascade order (lowest wins): default -> study -> group -> run
#
# "source" is the deepest level that explicitly sets the key.
# Keys not set at any user level have source "default".

function _load_yml(path::String)
    isfile(path) ? something(YAML.load_file(path), Dict()) : Dict()
end

function _flatten(d::Dict, prefix::String="")
    out = Dict{String,Any}()
    for (k, v) in d
        key = prefix == "" ? string(k) : "$(prefix).$(k)"
        if v isa Dict
            merge!(out, _flatten(v, key))
        else
            out[key] = v
        end
    end
    out
end

function _resolve_config(study::String, run::Union{String,Nothing}=nothing,
                         group::Union{String,Nothing}=nothing)
    root        = dirname(ServerState.data_dir())
    default_cfg = _flatten(_load_yml(joinpath(root, "config", "defaults", "pipeline.yml")))
    study_cfg   = _flatten(_load_yml(joinpath(ServerState.data_dir(), study, "pipeline.yml")))
    group_cfg   = isnothing(group) ? Dict{String,Any}() :
                  _flatten(_load_yml(joinpath(ServerState.data_dir(), study, group, "pipeline.yml")))
    run_cfg     = isnothing(run) ? Dict{String,Any}() :
                  _flatten(_load_yml(joinpath(ServerState.data_dir(), study,
                                              isnothing(group) ? run : joinpath(group, run),
                                              "pipeline.yml")))

    all_keys = union(keys(default_cfg), keys(study_cfg), keys(group_cfg), keys(run_cfg))

    Dict(k => begin
        if haskey(run_cfg, k)
            (; value=run_cfg[k],   source="run")
        elseif haskey(group_cfg, k)
            (; value=group_cfg[k], source="group")
        elseif haskey(study_cfg, k)
            (; value=study_cfg[k], source="study")
        elseif haskey(default_cfg, k)
            (; value=default_cfg[k], source="default")
        end
    end for k in all_keys)
end

function _write_override(path::String, dotted_key::String, value)
    cfg = _load_yml(path)
    parts = split(dotted_key, ".")
    d = cfg
    for p in parts[1:end-1]
        d = get!(d, p, Dict())
    end
    d[parts[end]] = value
    open(path, "w") do io
        YAML.write(io, cfg)
    end
end

function _delete_override(path::String, dotted_key::String)
    isfile(path) || return
    cfg   = _load_yml(path)
    parts = split(dotted_key, ".")
    d     = cfg
    for p in parts[1:end-1]
        haskey(d, p) || return
        d = d[p]
    end
    delete!(d, parts[end])
    open(path, "w") do io
        YAML.write(io, cfg)
    end
end

## Routes

@get "/api/v1/studies/{study}/config" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    json(_resolve_config(study))
end

@patch "/api/v1/studies/{study}/config" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    body = JSON3.read(req.body)
    path = joinpath(ServerState.data_dir(), study, "pipeline.yml")
    for (k, v) in body
        _write_override(path, string(k), v)
    end
    json(_resolve_config(study))
end

@delete "/api/v1/studies/{study}/config/{key}" function(req, study::String, key::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    path = joinpath(ServerState.data_dir(), study, "pipeline.yml")
    _delete_override(path, key)
    json(_resolve_config(study))
end

@get "/api/v1/studies/{study}/runs/{run}/config" function(req, study::String, run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    json(_resolve_config(study, run))
end

@patch "/api/v1/studies/{study}/runs/{run}/config" function(req, study::String, run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    body = JSON3.read(req.body)
    path = joinpath(ServerState.data_dir(), study, run, "pipeline.yml")
    for (k, v) in body
        _write_override(path, string(k), v)
    end
    json(_resolve_config(study, run))
end

@delete "/api/v1/studies/{study}/runs/{run}/config/{key}" function(req, study::String,
                                                                        run::String,
                                                                        key::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    path = joinpath(ServerState.data_dir(), study, run, "pipeline.yml")
    _delete_override(path, key)
    json(_resolve_config(study, run))
end
