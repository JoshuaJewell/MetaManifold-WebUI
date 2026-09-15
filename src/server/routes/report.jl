# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/studies/{study}/report
#
# A per-study basket of figures, tables and tree exports gathered for a
# manuscript. Items live as files in projects/{study}/report/ beside an
# ordered index, and export together as one ZIP with a captions file.
using JSON3, Dates, UUIDs, ZipArchives

const _REPORT_MAX_BYTES = 50 * 1024 * 1024
const _REPORT_KINDS = ("figure", "table", "tree")
const _REPORT_EXTENSIONS = (".svg", ".png", ".tif", ".json", ".csv", ".xlsx", ".nwk", ".txt", ".pdf")
const _report_lock = ReentrantLock()

_report_dir(study::String) = joinpath(ServerState.projects_dir(), study, "report")
_report_index(study::String) = joinpath(_report_dir(study), "items.json")

function _report_items(study::String)::Vector{Dict{String,Any}}
    path = _report_index(study)
    isfile(path) || return Dict{String,Any}[]
    raw = try JSON3.read(read(path, String)) catch; return Dict{String,Any}[] end
    [Dict{String,Any}(String(k) => v for (k, v) in item) for item in raw]
end

function _save_report_items(study::String, items)
    mkpath(_report_dir(study))
    tmp = _report_index(study) * ".tmp"
    write(tmp, JSON3.write(items))
    mv(tmp, _report_index(study); force=true)
end

_report_edit(f::Function, study::String) = lock(_report_lock) do
    items = _report_items(study)
    result = f(items)
    _save_report_items(study, items)
    result
end

function _report_study_error(study::String)
    (_valid_name(study) && study in _study_names()) || return json_error(404, "study_not_found", "Study '$study' not found")
    nothing
end

# A file name safe to store and to put in a ZIP, from a caption.
_report_slug(s::AbstractString) = let t = replace(strip(s), r"[^A-Za-z0-9._-]+" => "_")
    isempty(t) ? "item" : first(t, 60)
end

@get "/api/v1/studies/{study}/report" function(req, study::String)
    err = _report_study_error(study)
    isnothing(err) || return err
    json(_report_items(study))
end

# The file is the raw request body; kind, title and extension come as query parameters.
@post "/api/v1/studies/{study}/report" function(req, study::String)
    err = _report_study_error(study)
    isnothing(err) || return err
    q = queryparams(req)
    kind = get(q, "kind", "")
    kind in _REPORT_KINDS || return json_error(400, "invalid_kind", "kind must be one of $(join(_REPORT_KINDS, ", "))")
    ext = lowercase(get(q, "ext", ""))
    ext in _REPORT_EXTENSIONS || return json_error(400, "invalid_ext", "ext must be one of $(join(_REPORT_EXTENSIONS, ", "))")
    title = strip(get(q, "title", ""))
    isempty(title) && (title = "Untitled $kind")
    bytes = req.body
    isempty(bytes) && return json_error(400, "empty", "The item has no content")
    length(bytes) <= _REPORT_MAX_BYTES || return json_error(413, "too_large", "Report items are limited to 50 MB")

    id = string(uuid4())[1:8]
    file = "$(id)$(ext)"
    mkpath(_report_dir(study))
    write(joinpath(_report_dir(study), file), bytes)
    item = Dict{String,Any}("id" => id, "kind" => kind, "title" => title, "file" => file,
                            "added" => string(now(UTC)))
    _report_edit(items -> push!(items, item), study)
    json(item)
end

@patch "/api/v1/studies/{study}/report/{id}" function(req, study::String, id::String)
    err = _report_study_error(study)
    isnothing(err) || return err
    body = try JSON3.read(String(req.body)) catch; nothing end
    title = body isa JSON3.Object ? get(body, :title, nothing) : nothing
    title isa AbstractString && !isempty(strip(title)) ||
        return json_error(400, "invalid_title", "Body must include a non-empty 'title'")
    found = _report_edit(study) do items
        i = findfirst(it -> it["id"] == id, items)
        isnothing(i) && return nothing
        items[i]["title"] = strip(title)
        items[i]
    end
    isnothing(found) ? json_error(404, "not_found", "No report item '$id'") : json(found)
end

@put "/api/v1/studies/{study}/report/order" function(req, study::String)
    err = _report_study_error(study)
    isnothing(err) || return err
    body = try JSON3.read(String(req.body)) catch; nothing end
    ids = body isa JSON3.Object ? get(body, :ids, nothing) : nothing
    ids isa AbstractVector || return json_error(400, "invalid_body", "Body must include 'ids'")
    ordered = _report_edit(study) do items
        rank = Dict(string(id) => k for (k, id) in enumerate(ids))
        sort!(items; by = it -> get(rank, it["id"], typemax(Int)))
        copy(items)
    end
    json(ordered)
end

@delete "/api/v1/studies/{study}/report/{id}" function(req, study::String, id::String)
    err = _report_study_error(study)
    isnothing(err) || return err
    removed = _report_edit(study) do items
        i = findfirst(it -> it["id"] == id, items)
        isnothing(i) && return nothing
        popat!(items, i)
    end
    isnothing(removed) && return json_error(404, "not_found", "No report item '$id'")
    rm(joinpath(_report_dir(study), removed["file"]); force=true)
    json((; deleted=id))
end

const _REPORT_MIME = Dict(".svg" => "image/svg+xml", ".png" => "image/png", ".json" => "application/json",
                          ".csv" => "text/csv", ".nwk" => "text/plain", ".txt" => "text/plain",
                          ".pdf" => "application/pdf", ".tif" => "image/tiff",
                          ".xlsx" => "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

@get "/api/v1/studies/{study}/report/{id}/file" function(req, study::String, id::String)
    err = _report_study_error(study)
    isnothing(err) || return err
    item = findfirst(it -> it["id"] == id, _report_items(study))
    isnothing(item) && return json_error(404, "not_found", "No report item '$id'")
    it = _report_items(study)[item]
    path = joinpath(_report_dir(study), it["file"])
    isfile(path) || return json_error(404, "not_found", "The file for '$id' is missing")
    HTTP.Response(200, ["Content-Type" => get(_REPORT_MIME, lowercase(splitext(path)[2]), "application/octet-stream")];
                  body=read(path))
end

# Every item numbered in report order, plus captions.md listing them.
@get "/api/v1/studies/{study}/report/export" function(req, study::String)
    err = _report_study_error(study)
    isnothing(err) || return err
    items = _report_items(study)
    isempty(items) && return json_error(400, "empty", "The report has no items")
    io = IOBuffer()
    captions = IOBuffer()
    println(captions, "# $study report\n")
    ZipWriter(io) do w
        for (k, it) in enumerate(items)
            path = joinpath(_report_dir(study), it["file"])
            isfile(path) || continue
            name = lpad(k, 2, '0') * "_" * _report_slug(it["title"]) * splitext(it["file"])[2]
            zip_writefile(w, name, read(path))
            println(captions, "$k. `$name` - $(it["title"])")
        end
        zip_writefile(w, "captions.md", take!(captions))
    end
    HTTP.Response(200, ["Content-Type" => "application/zip",
                        "Content-Disposition" => "attachment; filename=\"$(study)_report.zip\""];
                  body=take!(io))
end
