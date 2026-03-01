module Server

    # © 2026 Joshua Benjamin Jewell. All rights reserved.
    #
    # This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

    # MetaManifold WebUI backend
    #
    # Start with:
    #   julia --project=. src/server/server.jl
    # or from a Julia session:
    #   include("src/server/server.jl"); Server.start()

    using Oxygen, HTTP, JSON3, YAML

    # Pipeline modules
    include(joinpath(@__DIR__, "..", "core",     "types.jl"))
    include(joinpath(@__DIR__, "..", "core",     "graph.jl"))
    include(joinpath(@__DIR__, "..", "core",     "log.jl"))
    include(joinpath(@__DIR__, "..", "core",     "config.jl"))
    include(joinpath(@__DIR__, "..", "core",     "databases.jl"))
    include(joinpath(@__DIR__, "..", "core",     "validate.jl"))
    include(joinpath(@__DIR__, "..", "core",     "project.jl"))
    include(joinpath(@__DIR__, "..", "pipeline", "tools.jl"))
    include(joinpath(@__DIR__, "..", "pipeline", "merge_taxa.jl"))
    include(joinpath(@__DIR__, "..", "pipeline", "dada2.jl"))
    include(joinpath(@__DIR__, "..", "pipeline", "swarm.jl"))
    include(joinpath(@__DIR__, "..", "analysis", "diversity.jl"))
    include(joinpath(@__DIR__, "..", "analysis", "plots.jl"))
    include(joinpath(@__DIR__, "..", "analysis", "plotly.jl"))
    include(joinpath(@__DIR__, "..", "analysis", "analysis.jl"))

    using .PipelineTypes, .PipelineGraph, .PipelineLog, .Config, .Databases
    using .Tools, .TaxonomyTableTools, .ProjectSetup
    using .DADA2, .OTUPipeline, .Analysis

    ## Server state

    module ServerState
        const _root = Ref{String}("")
        data_dir()     = joinpath(_root[], "data")
        projects_dir() = joinpath(_root[], "projects")
        set_root!(path::String) = (_root[] = abspath(path))
    end

    ## Job queue

    include(joinpath(@__DIR__, "jobs.jl"))
    using .JobQueue

    ## Shared helpers (available to all included route files)

    function json_error(status::Int, code::String, message::String; detail=nothing)
        body = isnothing(detail) ? (; error=code, message) : (; error=code, message, detail)
        HTTP.Response(status,
            ["Content-Type" => "application/json"],
            body = JSON3.write(body))
    end

    ## Routes

    include(joinpath(@__DIR__, "routes", "studies.jl"))
    include(joinpath(@__DIR__, "routes", "runs.jl"))
    include(joinpath(@__DIR__, "routes", "config.jl"))
    include(joinpath(@__DIR__, "routes", "pipeline.jl"))
    include(joinpath(@__DIR__, "routes", "jobs.jl"))
    include(joinpath(@__DIR__, "routes", "results.jl"))
    include(joinpath(@__DIR__, "routes", "databases.jl"))
    include(joinpath(@__DIR__, "routes", "events.jl"))

    ## Static file serving (frontend build + on-demand PDFs)

    const _frontend_dir = joinpath(@__DIR__, "..", "..", "web", "dist")
    const _mime_map = Dict(
        ".html" => "text/html", ".js" => "application/javascript",
        ".css"  => "text/css",  ".png" => "image/png",
        ".svg"  => "image/svg+xml", ".pdf" => "application/pdf",
        ".csv"  => "text/csv",  ".json" => "application/json",
    )

    # Static file serving middleware (Oxygen middleware signature: handler -> req -> response)
    function _file_middleware(next)
        function(req::HTTP.Request)
            uri = req.target

            # /files/{study}/runs/{run}/... -> serve from projects/
            m = match(r"^/files/([^/]+)/runs/([^/]+)/(.+)$", uri)
            if !isnothing(m)
                full     = abspath(joinpath(ServerState.projects_dir(), m[1], m[2], m[3]))
                projects = abspath(ServerState.projects_dir())
                startswith(full, projects * Base.Filesystem.path_separator) ||
                    return HTTP.Response(403, "Forbidden")
                isfile(full) || return HTTP.Response(404, "File not found")
                ext  = last(splitext(full))
                mime = get(_mime_map, ext, "application/octet-stream")
                return HTTP.Response(200, ["Content-Type" => mime]; body=read(full))
            end

            # Pass API requests through to Oxygen
            startswith(uri, "/api/") && return next(req)

            # SPA catch-all: serve frontend build or index.html
            rel    = lstrip(uri, '/')
            target = joinpath(_frontend_dir, rel)
            if isfile(target)
                ext  = last(splitext(target))
                mime = get(_mime_map, ext, "application/octet-stream")
                return HTTP.Response(200, ["Content-Type" => mime]; body=read(target))
            end
            index = joinpath(_frontend_dir, "index.html")
            isfile(index) && return HTTP.Response(200, ["Content-Type"=>"text/html"]; body=read(index))
            # No frontend build yet - fall through to Oxygen (404)
            next(req)
        end
    end

    ## Entry point

    function start(; root=pwd(), host="127.0.0.1", port=8080)
        ServerState.set_root!(root)
        @info "MetaManifold server starting" root host port

        # Initialise all projects on startup
        initialised = _ensure_all_projects()
        isempty(initialised) || @info "Initialised projects" initialised

        serve(; host, port, access_log=nothing, middleware=[_file_middleware], show_errors=true)
    end

end

# Allow running directly: julia src/server/server.jl
# Respects env vars JULIA_METAMANIFOLD_ROOT and JULIA_METAMANIFOLD_PORT for
# headless/subprocess startup (e.g. integration tests).
if abspath(PROGRAM_FILE) == @__FILE__
    root = get(ENV, "JULIA_METAMANIFOLD_ROOT", pwd())
    port = parse(Int, get(ENV, "JULIA_METAMANIFOLD_PORT", "8080"))
    Server.start(; root, port)
end
