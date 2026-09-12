# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/studies/{study}/trees
#
# Tree files live in projects/{study}/trees/ and are never modified by the
# viewer. Each file's view state (reroot, collapsed clades, renames, display
# settings) is a separate JSON document in trees/.views/{file}.json, written
# by the frontend and replayed over the tree on load.
using JSON3, SHA, Dates

const _TREE_EXTENSIONS = (".jplace", ".nwk", ".newick", ".tre", ".tree", ".treefile")
const _TREE_MAX_BYTES  = 50 * 1024 * 1024
const _VIEW_MAX_BYTES  = 5 * 1024 * 1024
const _tree_view_lock  = ReentrantLock()

_trees_dir(study::String) = joinpath(ServerState.projects_dir(), study, "trees")
_tree_view_path(study::String, file::String) = joinpath(_trees_dir(study), ".views", file * ".json")

_tree_format(file::String) = endswith(lowercase(file), ".jplace") ? "jplace" : "newick"

function _valid_tree_file(file::String)::Bool
    _valid_name(file) && any(ext -> endswith(lowercase(file), ext), _TREE_EXTENSIONS)
end

function _study_exists(study::String)::Bool
    _valid_name(study) &&
        (isdir(joinpath(ServerState.projects_dir(), study)) || study in _study_names())
end

function _tree_request_error(study::String, file::Union{String,Nothing}=nothing)
    _study_exists(study) || return json_error(404, "study_not_found", "Study '$study' not found")
    isnothing(file) && return nothing
    _valid_tree_file(file) || return json_error(400, "invalid_name",
        "Tree files need a plain name ending in one of: $(join(_TREE_EXTENSIONS, ", "))")
    nothing
end

function _tree_listing(study::String)
    dir = _trees_dir(study)
    isdir(dir) || return NamedTuple[]
    files = sort(filter(f -> _valid_tree_file(f) && isfile(joinpath(dir, f)), readdir(dir)))
    map(files) do f
        st = stat(joinpath(dir, f))
        (; file = f,
           format   = _tree_format(f),
           size     = st.size,
           modified = string(unix2datetime(st.mtime)),
           has_view = isfile(_tree_view_path(study, f)))
    end
end

function _write_atomic(path::String, bytes::AbstractString)
    mkpath(dirname(path))
    lock(_tree_view_lock) do
        tmp = path * ".tmp"
        write(tmp, bytes)
        mv(tmp, path; force=true)
    end
    nothing
end

# Reject content that cannot be a tree before it lands in the study.
function _check_tree_content(file::String, content::String)
    if _tree_format(file) == "jplace"
        doc = try JSON3.read(content) catch
            return "Not valid JSON"
        end
        doc isa JSON3.Object || return "A jplace file must be a JSON object"
        haskey(doc, :tree) && doc.tree isa AbstractString || return "Missing \"tree\" string"
        haskey(doc, :placements) || return "Missing \"placements\""
        haskey(doc, :fields) || return "Missing \"fields\""
        return nothing
    end
    # Leading comments such as FigTree's [&R] come before the tree itself.
    s = strip(replace(content, r"\A(\s*\[[^\]]*\])*" => ""))
    (startswith(s, "(") && endswith(s, ";")) || return "A Newick tree starts with '(' and ends with ';'"
    nothing
end

@get "/api/v1/studies/{study}/trees" function(req, study::String)
    err = _tree_request_error(study)
    isnothing(err) || return err
    json(_tree_listing(study))
end

@get "/api/v1/studies/{study}/trees/{file}" function(req, study::String, file::String)
    err = _tree_request_error(study, file)
    isnothing(err) || return err
    path = joinpath(_trees_dir(study), file)
    isfile(path) || return json_error(404, "not_found", "Tree '$file' not found")
    content = read(path, String)
    vpath = _tree_view_path(study, file)
    view = isfile(vpath) ? try JSON3.read(read(vpath, String)) catch; nothing end : nothing
    json((; file,
            format   = _tree_format(file),
            sha256   = bytes2hex(sha256(content)),
            modified = string(unix2datetime(stat(path).mtime)),
            content,
            view))
end

@post "/api/v1/studies/{study}/trees" function(req, study::String)
    err = _tree_request_error(study)
    isnothing(err) || return err
    body = try JSON3.read(String(req.body)) catch
        return json_error(400, "invalid_body", "Expected a JSON body")
    end
    file    = string(get(body, :file, ""))
    content = get(body, :content, nothing)
    err = _tree_request_error(study, file)
    isnothing(err) || return err
    content isa AbstractString || return json_error(400, "invalid_body", "\"content\" must be the file text")
    sizeof(content) <= _TREE_MAX_BYTES || return json_error(413, "too_large", "Tree files are limited to 50 MB")
    problem = _check_tree_content(file, String(content))
    isnothing(problem) || return json_error(400, "invalid_tree", "$file: $problem")
    path = joinpath(_trees_dir(study), file)
    if isfile(path) && get(body, :overwrite, false) != true
        return json_error(409, "exists", "Tree '$file' already exists")
    end
    _write_atomic(path, String(content))
    json((; file, format = _tree_format(file)))
end

@put "/api/v1/studies/{study}/trees/{file}/view" function(req, study::String, file::String)
    err = _tree_request_error(study, file)
    isnothing(err) || return err
    isfile(joinpath(_trees_dir(study), file)) ||
        return json_error(404, "not_found", "Tree '$file' not found")
    text = String(req.body)
    sizeof(text) <= _VIEW_MAX_BYTES || return json_error(413, "too_large", "View state is limited to 5 MB")
    doc = try JSON3.read(text) catch
        return json_error(400, "invalid_body", "Expected a JSON object")
    end
    doc isa JSON3.Object || return json_error(400, "invalid_body", "Expected a JSON object")
    _write_atomic(_tree_view_path(study, file), text)
    json((; saved = true))
end

@delete "/api/v1/studies/{study}/trees/{file}" function(req, study::String, file::String)
    err = _tree_request_error(study, file)
    isnothing(err) || return err
    path = joinpath(_trees_dir(study), file)
    isfile(path) || return json_error(404, "not_found", "Tree '$file' not found")
    rm(path)
    vpath = _tree_view_path(study, file)
    isfile(vpath) && rm(vpath)
    json((; deleted = file))
end
