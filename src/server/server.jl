module Server

    # © 2026 Joshua Benjamin Jewell. All rights reserved.
    #
    # This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

    # MetaManifold WebUI backend, a submodule of the package so that its
    # compiled code is cached with it.
    #
    # Start with start.sh, `julia --project=. scripts/serve.jl`, or from a
    # Julia session: `using MetaManifold; MetaManifold.Server.start()`.

    import ..MetaManifold
    using Oxygen, HTTP, JSON3, YAML, Logging
    using PrecompileTools: @setup_workload, @compile_workload

    using ..PipelineTypes, ..PipelineLog, ..Config
    using ..Databases, ..DuckDBStore, ..Validation
    using ..Tools, ..TaxonomyTableTools, ..ProjectSetup
    using ..DADA2, ..OTUPipeline
    using ..DiversityMetrics, ..Analysis

    ## EPIPE log filter
    # HTTP.jl logs every broken-pipe error from SSE streams as @error
    # "handle_connection handler error". These are harmless (client closed the
    # connection) but very noisy. Filter them before they reach the console.
    struct _SuppressEpipe{L<:AbstractLogger} <: AbstractLogger
        inner::L
    end
    Logging.min_enabled_level(l::_SuppressEpipe) = Logging.min_enabled_level(l.inner)
    Logging.shouldlog(l::_SuppressEpipe, args...) = Logging.shouldlog(l.inner, args...)
    function Logging.handle_message(l::_SuppressEpipe, level, msg, _module, group, id, file, line; kwargs...)
        # HTTP.jl embeds the exception as text inside the message string via
        # current_exceptions_to_string() - it is NOT passed as kwargs[:exception].
        # Match on the message text itself.
        msg_str = string(msg)
        if occursin("handle_connection handler error", msg_str) &&
           (occursin("EPIPE", msg_str) || occursin("broken pipe", msg_str))
            return
        end
        try
            Logging.handle_message(l.inner, level, msg, _module, group, id, file, line; kwargs...)
        catch e
            # When a job is stopped the TTY handle becomes invalid. uv_write()
            # then returns Nothing instead of Int32, triggering a TypeError
            # before Julia can turn it into an IOError/EPIPE. Suppress both so
            # "Exception while generating log record" spam is silenced on cancel.
            (e isa Base.IOError || e isa TypeError) && return
            rethrow()
        end
    end

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

    ## Route registration
    # Oxygen's routing macros register with its global router as the file
    # loads. A precompiled module does not load its files again, so these
    # macros shadow Oxygen's and collect the routes, and `__init__` registers
    # them each time the package is loaded.
    const ROUTES = Tuple{String,Any,Function}[]
    _route(method, path, func) = let (p, f) = Oxygen.adjustparams(path, func)
        :(push!(ROUTES, ($method, $(esc(p)), $(esc(f)))))
    end
    macro get(path, func = nothing)    _route("GET", path, func) end
    macro post(path, func = nothing)   _route("POST", path, func) end
    macro put(path, func = nothing)    _route("PUT", path, func) end
    macro patch(path, func = nothing)  _route("PATCH", path, func) end
    macro delete(path, func = nothing) _route("DELETE", path, func) end
    macro stream(path, func = nothing) _route("STREAM", path, func) end

    register_routes!() = foreach(((m, p, f),) -> Oxygen.route([m], p, f), ROUTES)
    __init__() = register_routes!()

    ## Routes
    include(joinpath(@__DIR__, "routes", "duckdb_helpers.jl"))
    include(joinpath(@__DIR__, "routes", "heatmap.jl"))
    include(joinpath(@__DIR__, "routes", "studies.jl"))
    include(joinpath(@__DIR__, "routes", "runs.jl"))
    include(joinpath(@__DIR__, "routes", "config.jl"))
    include(joinpath(@__DIR__, "routes", "pipeline.jl"))
    include(joinpath(@__DIR__, "routes", "jobs.jl"))
    include(joinpath(@__DIR__, "routes", "results.jl"))
    include(joinpath(@__DIR__, "routes", "sources.jl"))
    include(joinpath(@__DIR__, "routes", "databases.jl"))
    include(joinpath(@__DIR__, "routes", "events.jl"))
    include(joinpath(@__DIR__, "routes", "analysis.jl"))
    include(joinpath(@__DIR__, "routes", "read_funnel.jl"))
    include(joinpath(@__DIR__, "routes", "composition.jl"))
    include(joinpath(@__DIR__, "routes", "publication_tables.jl"))
    include(joinpath(@__DIR__, "routes", "trees.jl"))
    include(joinpath(@__DIR__, "routes", "report.jl"))
    include(joinpath(@__DIR__, "routes", "figures.jl"))
    include(joinpath(@__DIR__, "routes", "phylogeny.jl"))

    ## Bad request middleware
    # Helpers deep in a route (e.g. `_req_group`) throw this to reject a request.
    struct BadRequest <: Exception
        message::String
    end

    # Oxygen URL-decodes path parameters, so an encoded '/' would reach a handler
    # as a real separator. No study, run, table, column or file name needs one.
    function _unsafe_path_segment(target::AbstractString)
        path = first(split(target, '?'; limit=2))
        startswith(path, "/api/") || return nothing
        for seg in split(path, '/')
            d = HTTP.URIs.unescapeuri(seg)
            (d in (".", "..") || any(c -> c in ('/', '\\', '"', '\0'), d)) && return d
        end
        nothing
    end

    function _bad_request_middleware(next)
        function(req::HTTP.Request)
            bad = _unsafe_path_segment(req.target)
            isnothing(bad) || return json_error(400, "bad_request", "Invalid path segment: $(repr(bad))")
            if startswith(req.target, "/api/")
                g = get(HTTP.queryparams(HTTP.URI(req.target)), "group", "")
                isempty(g) || _valid_name(g) ||
                    return json_error(400, "bad_request", "Invalid group name: $(repr(g))")
            end
            try
                next(req)
            catch e
                e isa BadRequest || rethrow()
                json_error(400, "bad_request", e.message)
            end
        end
    end

    ## Analysis source
    # The classifier an analysis request names (body or query `source`), held for
    # the request so column resolution can follow it; nothing keeps the configured one.
    _analysis_source() = get(task_local_storage(), :mm_analysis_source, nothing)

    function _request_source(req::HTTP.Request)
        occursin(r"/(analysis|composition)(/|\?|$)", req.target) || return nothing
        q = get(HTTP.queryparams(HTTP.URI(req.target)), "source", nothing)
        if isnothing(q) && req.method == "POST" && !isempty(req.body)
            body = try JSON3.read(String(copy(req.body))) catch; nothing end
            if body isa JSON3.Object
                q = get(body, :source, nothing)
                # Cross-run requests carry it on their run specs.
                runs = get(body, :runs, nothing)
                if isnothing(q) && runs isa AbstractVector && !isempty(runs) && first(runs) isa JSON3.Object
                    q = get(first(runs), :source, nothing)
                end
            end
        end
        q isa AbstractString || return nothing
        s = uppercase(q)
        s in ("VSEARCH", "DADA2") ? s : nothing
    end

    function _analysis_source_middleware(next)
        function(req::HTTP.Request)
            src = _request_source(req)
            isnothing(src) && return next(req)
            task_local_storage(() -> next(req), :mm_analysis_source, src)
        end
    end

    # Rank and category columns for the request's source.
    _rank_column(columns, rank::AbstractString) = let src = _analysis_source()
        isnothing(src) ? MetaManifold.Analysis.taxon_column(columns, String(rank)) :
                         MetaManifold.Categories.rank_col(rank, src, Set(columns))
    end
    _category_column(set_name::AbstractString) = MetaManifold.Categories.column_name(set_name, _analysis_source())
    _suffixed() = !isnothing(_analysis_source())

    ## R-runtime busy middleware
    # The embedded R interpreter is shared between the pipeline and the analysis
    # endpoints (see `RRuntime`). A pipeline run can hold it for hours, so an
    # analysis waits only briefly and then throws `RBusyError` rather than hanging
    # the response. Report that as a 503, matching the existing `r_unavailable`
    # shape, so the frontend can tell the user to retry once the run completes.
    function _r_busy_middleware(next)
        function(req::HTTP.Request)
            try
                next(req)
            catch e
                e isa MetaManifold.RRuntime.RBusyError || rethrow()
                json_error(503, "r_busy",
                    "The R runtime is busy with a pipeline run - retry once it completes")
            end
        end
    end

    ## CORS middleware (needed when the frontend is served from a different origin)
    # Loopback binds also check Host, so a DNS-rebinding page cannot read the API.
    const _bound_host = Ref("127.0.0.1")
    _is_loopback(h::AbstractString) = lowercase(h) in ("localhost", "127.0.0.1", "::1", "[::1]")

    function _cors_middleware(next)
        function(req::HTTP.Request)
            if _is_loopback(_bound_host[])
                host = replace(HTTP.header(req, "Host", ""), r":\d+\z" => "")
                isempty(host) || _is_loopback(host) ||
                    return HTTP.Response(403, "Forbidden: unexpected Host")
            end
            origin = HTTP.header(req, "Origin", "")
            if isempty(origin)
                # Same-origin request - no CORS headers needed
                return next(req)
            end
            # Only allow localhost origins
            origin_url = try HTTP.URIs.URI(origin) catch; nothing end
            if isnothing(origin_url) || !(lowercase(origin_url.host) in ("localhost", "127.0.0.1", "::1"))
                return HTTP.Response(403, "Forbidden: non-localhost origin")
            end
            cors_headers = [
                "Access-Control-Allow-Origin"  => origin,
                "Access-Control-Allow-Methods" => "GET, POST, PUT, PATCH, DELETE, OPTIONS",
                "Access-Control-Allow-Headers" => "Content-Type",
            ]
            # Handle preflight
            if req.method == "OPTIONS"
                return HTTP.Response(204, cors_headers)
            end
            resp = next(req)
            for (k, v) in cors_headers
                HTTP.setheader(resp, k => v)
            end
            resp
        end
    end

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
            uri = first(split(req.target, '?'; limit=2))

            # /files/{study}/runs/{rest...} -> serve from projects/{study}/{rest}
            # rest may be {run}/... or {group}/{run}/... (group paths have an extra segment)
            m = match(r"^/files/([^/]+)/runs/(.+)$", uri)
            if !isnothing(m)
                candidate = abspath(joinpath(ServerState.projects_dir(), HTTP.URIs.unescapeuri(m[1]), HTTP.URIs.unescapeuri(m[2])))
                isfile(candidate) || return HTTP.Response(404, "File not found")
                full     = realpath(candidate)
                projects = realpath(ServerState.projects_dir())
                startswith(full, projects * Base.Filesystem.path_separator) ||
                    return HTTP.Response(403, "Forbidden")
                ext  = last(splitext(full))
                mime = get(_mime_map, ext, "application/octet-stream")
                return HTTP.Response(200, ["Content-Type" => mime]; body=read(full))
            end

            # Pass API requests through to Oxygen
            startswith(uri, "/api/") && return next(req)

            # SPA catch-all: serve frontend build or index.html
            rel    = lstrip(HTTP.URIs.unescapeuri(uri), '/')
            root   = abspath(_frontend_dir)
            target = abspath(joinpath(root, rel))
            inside = startswith(target, root * Base.Filesystem.path_separator)
            if inside && isfile(target)
                ext  = last(splitext(target))
                mime = get(_mime_map, ext, "application/octet-stream")
                # Content-hashed assets (js/css in assets/) are immutable.
                # Everything else (index.html, config.json) must not be cached so
                # the browser always loads the latest bundle after a rebuild.
                cache = (ext in (".js", ".css") && occursin("/assets/", uri)) ?
                    "public, max-age=31536000, immutable" : "no-store"
                return HTTP.Response(200,
                    ["Content-Type" => mime, "Cache-Control" => cache];
                    body=read(target))
            end
            index = joinpath(_frontend_dir, "index.html")
            isfile(index) && return HTTP.Response(200,
                ["Content-Type" => "text/html", "Cache-Control" => "no-store"];
                body=read(index))
            # No frontend build yet - fall through to Oxygen (404)
            next(req)
        end
    end

    ## Startup warm-up
    # The precompile workload (precompile.jl) compiles the routes a first visit
    # uses into the package cache. A developer who turns it off in
    # LocalPreferences.toml gets the same requests made once in the background
    # as the server starts, so the first page visit does not wait on compilation.
    _workload_enabled() =
        get(Base.get_preferences(Base.PkgId(MetaManifold).uuid), "precompile_workload", true) !== false
    function _warm_routes(host::String, port::Int)
        base = "http://$(host == "0.0.0.0" ? "127.0.0.1" : host):$port/api/v1"
        fetch(path) = try
            # Closed after each request: an idle keep-alive connection back to this
            # server would hold its shutdown open.
            r = HTTP.get(base * path, ["Connection" => "close"];
                         status_exception=false, retry=false, readtimeout=600)
            r.status == 200 ? JSON3.read(r.body) : nothing
        catch
            nothing
        end
        for _ in 1:240
            isnothing(fetch("/studies")) || break
            sleep(0.5)
        end
        t0 = time()
        fetch("/jobs")
        for s in something(fetch("/studies"), [])
            study = HTTP.escapeuri(String(s.name))
            detail = fetch("/studies/$study")
            isnothing(detail) && continue
            foreach(fetch, ("/studies/$study/config", "/studies/$study/config/overrides",
                            "/studies/$study/runs"))
            runs = get(detail, :runs, [])
            isempty(runs) && continue
            run = "/studies/$study/runs/$(HTTP.escapeuri(String(first(runs))))"
            foreach(fetch, (run, "$run/config", "$run/results/qc", "$run/results/dada2",
                            "$run/results/tables", "$run/analysis/read-funnel"))
            break
        end
        @info "Routes warmed" seconds = round(time() - t0; digits=1)
    end

    ## Entry point
    function start(; root=pwd(), host="127.0.0.1", port=8080, listen::Bool=true)
        global_logger(_SuppressEpipe(global_logger()))
        ServerState.set_root!(root)
        _bound_host[] = host
        @info "MetaManifold server starting" root host port

        # Initialise all projects on startup
        initialised = _ensure_all_projects()
        isempty(initialised) || @info "Initialised projects" initialised

        listen || return nothing
        _workload_enabled() || errormonitor(@async _warm_routes(host, port))
        serve(; host, port, access_log=nothing,
              middleware=[_bad_request_middleware, _r_busy_middleware, _cors_middleware,
                          _analysis_source_middleware, _file_middleware],
              show_errors=false)
    end

    """
        main()

    Start the server with the root and port from `JULIA_METAMANIFOLD_ROOT` and
    `JULIA_METAMANIFOLD_PORT`, as start.sh and the integration tests do.
    """
    main() = start(; root=get(ENV, "JULIA_METAMANIFOLD_ROOT", pwd()),
                     port=parse(Int, get(ENV, "JULIA_METAMANIFOLD_PORT", "8080")))

    include(joinpath(@__DIR__, "precompile.jl"))

end
