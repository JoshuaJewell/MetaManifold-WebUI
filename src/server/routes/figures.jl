# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/studies/{study}/figures
#
# Multi-panel figure layouts. Each is a JSON document in
# projects/{study}/figures/{id}.json holding the page, the lettered groups and
# the chart definition of every pane. The frontend renders the charts from the
# data each time, so a figure follows the results it is built on, and saves a
# PNG of the page beside the layout whenever it changes.
using JSON3, Dates, UUIDs

const _FIGURE_MAX_BYTES = 2 * 1024 * 1024
const _FIGURE_RASTER_MAX_BYTES = 200 * 1024 * 1024

_figures_dir(study::String) = joinpath(ServerState.projects_dir(), study, "figures")
_figure_path(study::String, id::String) = joinpath(_figures_dir(study), id * ".json")
_valid_figure_id(id::String) = occursin(r"^[0-9a-f]{8}$", id)

function _figure_request_error(study::String, id::Union{String,Nothing}=nothing)
    _study_exists(study) || return json_error(404, "study_not_found", "Study '$study' not found")
    isnothing(id) && return nothing
    _valid_figure_id(id) || return json_error(400, "invalid_id", "Figure ids are 8 hex characters")
    isfile(_figure_path(study, id)) || return json_error(404, "not_found", "No figure '$id'")
    nothing
end

function _read_figure(study::String, id::String)
    try JSON3.read(read(_figure_path(study, id), String)) catch; nothing end
end

@get "/api/v1/studies/{study}/figures" function(req, study::String)
    err = _figure_request_error(study)
    isnothing(err) || return err
    dir = _figures_dir(study)
    isdir(dir) || return json([])
    ids = [splitext(f)[1] for f in readdir(dir) if endswith(f, ".json") && _valid_figure_id(splitext(f)[1])]
    list = map(ids) do id
        doc = _read_figure(study, id)
        (; id, title = doc isa JSON3.Object ? string(get(doc, :title, id)) : id,
           modified = string(unix2datetime(mtime(_figure_path(study, id)))))
    end
    json(sort(list; by = x -> x.modified, rev=true))
end

@get "/api/v1/studies/{study}/figures/{id}" function(req, study::String, id::String)
    err = _figure_request_error(study, id)
    isnothing(err) || return err
    HTTP.Response(200, ["Content-Type" => "application/json"]; body=read(_figure_path(study, id)))
end

# The body is the whole document; the server assigns the id.
@post "/api/v1/studies/{study}/figures" function(req, study::String)
    err = _figure_request_error(study)
    isnothing(err) || return err
    _save_figure(req, study, string(uuid4())[1:8])
end

@put "/api/v1/studies/{study}/figures/{id}" function(req, study::String, id::String)
    err = _figure_request_error(study, id)
    isnothing(err) || return err
    _save_figure(req, study, id)
end

function _save_figure(req, study::String, id::String)
    length(req.body) <= _FIGURE_MAX_BYTES || return json_error(413, "too_large", "Figure documents are limited to 2 MB")
    doc = try JSON3.read(String(req.body)) catch; nothing end
    doc isa JSON3.Object || return json_error(400, "invalid_body", "Body must be a JSON object")
    out = Dict{String,Any}(String(k) => v for (k, v) in doc)
    out["id"] = id
    _write_atomic(_figure_path(study, id), JSON3.write(out))
    json(out)
end

# Rendered copies are named after the figure's title and end in _{id}.png.
_figure_rasters(study::String, id::String) =
    [joinpath(_figures_dir(study), f) for f in readdir(_figures_dir(study)) if endswith(f, "_$(id).png")]

# The page rendered by the browser, saved beside the layout. The body is the PNG.
@put "/api/v1/studies/{study}/figures/{id}/raster" function(req, study::String, id::String)
    err = _figure_request_error(study, id)
    isnothing(err) || return err
    bytes = req.body
    length(bytes) <= _FIGURE_RASTER_MAX_BYTES || return json_error(413, "too_large", "Rendered figures are limited to 200 MB")
    (length(bytes) > 8 && bytes[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]) ||
        return json_error(400, "invalid_png", "Body must be a PNG image")
    doc = _read_figure(study, id)
    title = doc isa JSON3.Object ? string(get(doc, :title, "figure")) : "figure"
    slug = strip(first(replace(strip(title), r"[^A-Za-z0-9._-]+" => "_"), 60), '_')
    path = joinpath(_figures_dir(study), "$(isempty(slug) ? "figure" : slug)_$(id).png")
    tmp = path * ".tmp"
    write(tmp, bytes)
    foreach(f -> f == path || rm(f; force=true), _figure_rasters(study, id))
    mv(tmp, path; force=true)
    json((; file=basename(path)))
end

@delete "/api/v1/studies/{study}/figures/{id}" function(req, study::String, id::String)
    err = _figure_request_error(study, id)
    isnothing(err) || return err
    foreach(f -> rm(f; force=true), _figure_rasters(study, id))
    rm(_figure_path(study, id); force=true)
    json((; deleted=id))
end
