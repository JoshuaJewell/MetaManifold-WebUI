# Unit tests for the server's request guards: the Origin allow-list, Oxygen's
# docs and metrics being off, and the headers report uploads are served with.
# Requests go through the app middleware in-process; no socket is opened.
using MetaManifold
SV = MetaManifold.Server
using JSON3

# A request through the middleware `start` installs.
_ss_request(method, path, headers=Pair{String,String}[], body="") =
    SV.Oxygen.internalrequest(SV.HTTP.Request(method, path, headers, body); middleware=SV._MIDDLEWARE)

# A request through the handler Oxygen's `serve` builds from `_SERVE_OPTIONS`,
# docs and metrics included, without starting the listener.
function _ss_served(path)
    ctx = SV.Oxygen.CONTEXT[]
    opts = SV._SERVE_OPTIONS
    saved = ctx.docs.router[]
    try
        ctx.docs.router[] = SV.HTTP.Router()
        opts.docs && SV.Oxygen.Core.setupdocs(ctx)
        opts.metrics && SV.Oxygen.Core.setupmetrics(ctx)
        handler = SV.Oxygen.Core.setupmiddleware(ctx; middleware=SV._MIDDLEWARE,
            docs=opts.docs, metrics=opts.metrics, show_errors=opts.show_errors)
        handler(SV.HTTP.Request("GET", path, ["Host" => "127.0.0.1:8080"]))
    finally
        ctx.docs.router[] = saved
    end
end

@testset "Server request guards" begin
    tmp = mktempdir()
    old_root = SV.ServerState._root[]
    old_host, old_port = SV._bound_host[], SV._bound_port[]
    SV.ServerState.set_root!(tmp)
    SV._bound_host[] = "127.0.0.1"
    SV._bound_port[] = 8080
    mkpath(joinpath(tmp, "data", "StudyS", "run1"))
    touch(joinpath(tmp, "data", "StudyS", "run1", "s1_R1.fastq.gz"))
    mkpath(joinpath(tmp, "projects", "StudyS", "run1"))
    try
        @testset "Origin allow-list" begin
            studies = "/api/v1/studies"
            host = "Host" => "127.0.0.1:8080"
            # No Origin (same-origin GET, curl): unchanged.
            @test _ss_request("GET", studies, [host]).status == 200
            # The served origin, in both spellings, and the Vite dev server.
            for o in ("http://127.0.0.1:8080", "http://localhost:8080", "http://localhost:5173",
                      "http://127.0.0.1:5173", "http://LOCALHOST:8080/")
                r = _ss_request("GET", studies, [host, "Origin" => o])
                @test r.status == 200
                @test SV.HTTP.header(r, "Access-Control-Allow-Origin") == o
                @test _ss_request("OPTIONS", studies, [host, "Origin" => o]).status == 204
            end
            # Another loopback page: other port, other scheme, or no origin at all.
            for o in ("http://localhost:3000", "http://127.0.0.1:9999", "https://localhost:8080",
                      "http://[::1]:8080", "null", "http://evil.example")
                @test _ss_request("GET", studies, [host, "Origin" => o]).status == 403
                @test _ss_request("OPTIONS", studies, [host, "Origin" => o]).status == 403
                # A mutating request is refused before it reaches the route.
                r = _ss_request("POST", "/api/v1/studies/StudyS/report?kind=figure&ext=.txt",
                                [host, "Origin" => o, "Content-Type" => "text/plain"], "x")
                @test r.status == 403
            end
            @test !isdir(joinpath(tmp, "projects", "StudyS", "report"))
            # The allowed set follows the bound port.
            SV._bound_port[] = 9090
            @test _ss_request("GET", studies, [host, "Origin" => "http://localhost:9090"]).status == 200
            @test _ss_request("GET", studies, [host, "Origin" => "http://localhost:8080"]).status == 403
            SV._bound_port[] = 8080
            # Extra origins come from JULIA_METAMANIFOLD_ALLOWED_ORIGINS.
            withenv("JULIA_METAMANIFOLD_ALLOWED_ORIGINS" => " http://lab-box:8080 , http://localhost:5174/") do
                @test _ss_request("GET", studies, [host, "Origin" => "http://lab-box:8080"]).status == 200
                @test _ss_request("GET", studies, [host, "Origin" => "http://localhost:5174"]).status == 200
                @test _ss_request("GET", studies, [host, "Origin" => "http://localhost:3000"]).status == 403
            end
            # The Host check still applies first.
            @test _ss_request("GET", studies, ["Host" => "evil.example:8080"]).status == 403
        end

        @testset "Oxygen docs and metrics are not served" begin
            @test SV._SERVE_OPTIONS.docs === false
            @test SV._SERVE_OPTIONS.metrics === false
            for path in ("/docs", "/docs/schema", "/docs/swagger", "/docs/metrics")
                r = _ss_served(path)
                body = lowercase(String(r.body))
                @test !occursin("openapi", body)
                @test !occursin("swagger", body)
                @test !occursin("redoc", body)
                @test r.status != 200
            end
            # The served handler still answers the API.
            @test _ss_served("/api/v1/studies").status == 200
        end

        @testset "Report uploads are served inert" begin
            host = "Host" => "127.0.0.1:8080"
            base = "/api/v1/studies/StudyS/report"
            svg = "<svg xmlns='http://www.w3.org/2000/svg'><script>alert(1)</script></svg>"
            r = _ss_request("POST", "$base?kind=figure&title=X&ext=.svg", [host], svg)
            @test r.status == 200
            item = JSON3.read(String(r.body))
            for path in ("$base/$(item.id)/file", "/files/StudyS/runs/report/$(item.file)")
                f = _ss_request("GET", path, [host])
                @test f.status == 200
                @test String(f.body) == svg
                @test SV.HTTP.header(f, "Content-Type") == "image/svg+xml"
                @test SV.HTTP.header(f, "X-Content-Type-Options") == "nosniff"
                @test SV.HTTP.header(f, "Content-Security-Policy") ==
                      "default-src 'none'; style-src 'unsafe-inline'; sandbox"
            end
            # Pipeline output under /files/ is sniff-proof but keeps working inline.
            write(joinpath(tmp, "projects", "StudyS", "run1", "plot.pdf"), "%PDF-1.4")
            f = _ss_request("GET", "/files/StudyS/runs/run1/plot.pdf", [host])
            @test f.status == 200
            @test SV.HTTP.header(f, "X-Content-Type-Options") == "nosniff"
            @test SV.HTTP.header(f, "Content-Security-Policy") == ""
        end
    finally
        SV.ServerState._root[] = old_root
        SV._bound_host[], SV._bound_port[] = old_host, old_port
        rm(tmp; recursive=true, force=true)
    end
end
