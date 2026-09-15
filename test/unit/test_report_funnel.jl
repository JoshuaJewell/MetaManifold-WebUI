# Unit tests for the report basket, read funnel, config value checks and
# stage-hash snapshots.
using MetaManifold
SV = MetaManifold.Server
using JSON3, YAML

_rf_request(method, path, body="", headers=["Content-Type" => "application/json"]) =
    SV.Oxygen.internalrequest(SV.HTTP.Request(method, path, headers, body))

@testset "Report basket" begin
    tmp = mktempdir()
    old_root = SV.ServerState._root[]
    SV.ServerState.set_root!(tmp)
    mkpath(joinpath(tmp, "data", "StudyR", "run1"))
    touch(joinpath(tmp, "data", "StudyR", "run1", "s1_R1.fastq.gz"))
    mkpath(joinpath(tmp, "projects", "StudyR"))
    try
        base = "/api/v1/studies/StudyR/report"
        @test isempty(JSON3.read(String(_rf_request("GET", base).body)))
        svg = "<svg xmlns='http://www.w3.org/2000/svg'/>"
        r = _rf_request("POST", "$base?kind=figure&title=Alpha&ext=.svg", svg,
                        ["Content-Type" => "application/octet-stream"])
        @test r.status == 200
        a = JSON3.read(String(r.body))
        r = _rf_request("POST", "$base?kind=table&title=Tables&ext=.csv", "a,b\n1,2\n",
                        ["Content-Type" => "application/octet-stream"])
        b = JSON3.read(String(r.body))
        @test _rf_request("POST", "$base?kind=bad&ext=.svg", svg).status == 400
        @test _rf_request("POST", "$base?kind=figure&ext=.exe", svg).status == 400

        @test String(_rf_request("GET", "$base/$(a.id)/file").body) == svg
        @test _rf_request("PATCH", "$base/$(a.id)", JSON3.write(Dict("title" => "Figure 1"))).status == 200
        order = JSON3.read(String(_rf_request("PUT", "$base/order", JSON3.write(Dict("ids" => [b.id, a.id]))).body))
        @test [o.id for o in order] == [b.id, a.id]

        z = _rf_request("GET", "$base/export")
        @test z.status == 200
        @test z.body[1:2] == UInt8['P', 'K']

        @test _rf_request("DELETE", "$base/$(a.id)").status == 200
        @test length(JSON3.read(String(_rf_request("GET", base).body))) == 1
    finally
        SV.ServerState.set_root!(old_root)
        rm(tmp; recursive=true, force=true)
    end
end

@testset "Read funnel" begin
    tmp = mktempdir()
    try
        log = joinpath(tmp, "stats.txt")
        write(log, """
        This is cutadapt 4.4
        Command line parameters: -g X -o /p/cutadapt/S1_R1_trimmed.fastq.gz -p /p/cutadapt/S1_R2_trimmed.fastq.gz a b
        Total read pairs processed:             42,196
        Pairs written (passing filters):        25,719 (61.0%)
        This is cutadapt 4.4
        Command line parameters: -g X -o /p/cutadapt/S2_R1_trimmed.fastq.gz a
        Total reads processed:             1,000
        Reads written (passing filters):   900 (90.0%)
        """)
        counts = SV._cutadapt_counts(log)
        @test counts["S1"] == (42196, 25719)
        @test counts["S2"] == (1000, 900)
        @test isempty(SV._cutadapt_counts(joinpath(tmp, "missing.txt")))
    finally
        rm(tmp; recursive=true, force=true)
    end
end

@testset "Config value checks" begin
    @test isnothing(SV._config_value_error("vsearch.strand", "both"))
    @test !isnothing(SV._config_value_error("vsearch.strand", "both; rm -rf x"))
    @test !isnothing(SV._config_value_error("cutadapt.optional_args", "-m 1\ntouch x"))
    @test isnothing(SV._config_value_error("cutadapt.optional_args", "--nextseq-trim 20"))
    @test !isnothing(SV._config_value_error("vsearch.identity", "0.97"))
    @test isnothing(SV._config_value_error("vsearch.identity", 0.97))
    @test !isnothing(SV._config_value_error("remote.host", "-oProxyCommand=x"))
    @test !isnothing(SV._config_value_error("analysis.beta.normalisation", "log"))
    @test isnothing(SV._config_value_error("analysis.alpha.normalisation", "srs"))
end

@testset "Stage hash snapshots" begin
    tmp = mktempdir()
    try
        cfg = joinpath(tmp, "run_config.yml")
        hash_file = joinpath(tmp, "stage.hash")
        write(cfg, "vsearch:\n  identity: 0.8\n")
        MetaManifold.Config._write_section_hash(cfg, "vsearch", hash_file)
        before = MetaManifold.Config.stage_run_count(cfg)
        snap = MetaManifold.Config._begin_section(cfg, "vsearch", hash_file)
        @test !isfile(hash_file)
        @test MetaManifold.Config.stage_run_count(cfg) == before + 1
        # A config edit made while the stage runs is not recorded as what it ran with.
        write(cfg, "vsearch:\n  identity: 0.9\n")
        MetaManifold.Config._write_section_hash(cfg, "vsearch", hash_file; snapshot=snap)
        @test MetaManifold.Config._section_stale(cfg, "vsearch", hash_file)
        @test occursin("identity=0.8", read(hash_file * ".values", String))
    finally
        rm(tmp; recursive=true, force=true)
    end
end

@testset "Repeated rarefaction" begin
    D = MetaManifold.DiversityMetrics
    mat = Float64[5 0 0; 50 30 20; 40 40 20]
    n = D.normalise_counts(mat; method="rarefy", depth=10, seed=1)
    @test n.kept == [2, 3]
    @test all(sum(n.mat; dims=2) .== 10)
    @test D.normalise_counts(mat; method="none", depth=0, seed=1).kept == [1, 2, 3]
    a = D.alpha_diversity(mat; method="rarefy", depth=10, seed=1, iterations=20)
    @test a.kept == [2, 3]
    @test all(1 .<= a.richness .<= 3)
    @test eltype(a.richness) == Float64
end

@testset "Source-suffixed category columns" begin
    C = MetaManifold.Categories
    @test C.column_name("default") == "Category__default"
    @test C.column_name("default", "DADA2") == "Category__default__DADA2"
end

@testset "Figure documents" begin
    tmp = mktempdir()
    old_root = SV.ServerState._root[]
    SV.ServerState.set_root!(tmp)
    mkpath(joinpath(tmp, "data", "StudyF", "run1"))
    touch(joinpath(tmp, "data", "StudyF", "run1", "s1_R1.fastq.gz"))
    mkpath(joinpath(tmp, "projects", "StudyF"))
    try
        base = "/api/v1/studies/StudyF/figures"
        @test isempty(JSON3.read(String(_rf_request("GET", base).body)))
        r = _rf_request("POST", base, JSON3.write(Dict("title" => "Fig 1", "groups" => [])))
        @test r.status == 200
        id = JSON3.read(String(r.body)).id
        @test occursin(r"^[0-9a-f]{8}$", id)
        r = _rf_request("PUT", "$base/$id", JSON3.write(Dict("id" => "ffffffff", "title" => "Renamed")))
        @test JSON3.read(String(r.body)).id == id
        @test JSON3.read(String(_rf_request("GET", "$base/$id").body)).title == "Renamed"
        @test only(JSON3.read(String(_rf_request("GET", base).body))).title == "Renamed"
        @test _rf_request("PUT", "$base/$id", "[1,2]").status == 400
        @test _rf_request("GET", "$base/nothex12").status == 400
        @test _rf_request("GET", "$base/00000000").status == 404
        png = vcat(UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a], zeros(UInt8, 16))
        dir = joinpath(tmp, "projects", "StudyF", "figures")
        @test _rf_request("PUT", "$base/$id/raster", png, ["Content-Type" => "image/png"]).status == 200
        @test isfile(joinpath(dir, "Renamed_$(id).png"))
        @test _rf_request("PUT", "$base/$id/raster", "not a png", ["Content-Type" => "image/png"]).status == 400
        _rf_request("PUT", "$base/$id", JSON3.write(Dict("title" => "Final figure")))
        _rf_request("PUT", "$base/$id/raster", png, ["Content-Type" => "image/png"])
        @test filter(f -> endswith(f, ".png"), readdir(dir)) == ["Final_figure_$(id).png"]
        @test _rf_request("DELETE", "$base/$id").status == 200
        @test _rf_request("GET", "$base/$id").status == 404
        @test isempty(readdir(dir))
    finally
        SV.ServerState.set_root!(old_root)
        rm(tmp; recursive=true, force=true)
    end
end
