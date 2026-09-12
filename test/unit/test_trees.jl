# Unit tests for the tree file and view-state routes.
using MetaManifold
SV = MetaManifold.Server
using JSON3

function _tr_request(method::String, path::String, body=nothing)
    payload = isnothing(body) ? "" : body isa AbstractString ? body : JSON3.write(body)
    req = SV.HTTP.Request(method, path, ["Content-Type" => "application/json"], payload)
    SV.Oxygen.internalrequest(req)
end
_tr_json(res) = JSON3.read(String(res.body))

@testset "Tree routes" begin
    @test SV._valid_tree_file("a-b_1.jplace")
    @test SV._valid_tree_file("x.treefile")
    @test !SV._valid_tree_file("x.txt")
    @test !SV._valid_tree_file("../x.nwk")
    @test !SV._valid_tree_file(".hidden.nwk")
    @test SV._tree_format("A.JPLACE") == "jplace"
    @test SV._tree_format("a.nwk") == "newick"
    @test isnothing(SV._check_tree_content("a.nwk", "((A,B),C);"))
    @test !isnothing(SV._check_tree_content("a.nwk", "hello"))
    @test !isnothing(SV._check_tree_content("a.jplace", "{\"tree\": \"(A,B);\"}"))

    tmp = mktempdir()
    old_root = SV.ServerState._root[]
    SV.ServerState.set_root!(tmp)
    mkpath(joinpath(tmp, "projects", "StudyT"))
    try
        base = "/api/v1/studies/StudyT/trees"
        @test isempty(_tr_json(_tr_request("GET", base)))
        @test _tr_request("GET", "/api/v1/studies/Nope/trees").status == 404

        jp = JSON3.write(Dict("version" => 3, "tree" => "((A:1{0},B:1{1}):1{2},C:1{3}){4};",
                              "placements" => [], "fields" => ["edge_num"]))
        @test _tr_request("POST", base, Dict("file" => "t.jplace", "content" => jp)).status == 200
        @test _tr_request("POST", base, Dict("file" => "t.jplace", "content" => jp)).status == 409
        @test _tr_request("POST", base, Dict("file" => "t.jplace", "content" => jp, "overwrite" => true)).status == 200
        @test _tr_request("POST", base, Dict("file" => "bad.nwk", "content" => "nope")).status == 400
        @test _tr_request("POST", base, Dict("file" => "evil.sh", "content" => "(A);")).status == 400

        listing = _tr_json(_tr_request("GET", base))
        @test length(listing) == 1
        @test listing[1].file == "t.jplace" && listing[1].format == "jplace" && !listing[1].has_view

        got = _tr_json(_tr_request("GET", "$base/t.jplace"))
        @test got.content == jp
        @test isnothing(got.view)
        @test length(got.sha256) == 64

        view = Dict("version" => 1, "ops" => [Dict("op" => "collapse", "clade" => "abc")])
        @test _tr_request("PUT", "$base/t.jplace/view", view).status == 200
        @test _tr_request("PUT", "$base/t.jplace/view", "[1,2]").status == 400
        @test _tr_request("PUT", "$base/missing.nwk/view", view).status == 404
        got = _tr_json(_tr_request("GET", "$base/t.jplace"))
        @test got.view.ops[1].clade == "abc"
        @test _tr_json(_tr_request("GET", base))[1].has_view
        # The tree file itself is untouched by view saves.
        @test read(joinpath(tmp, "projects", "StudyT", "trees", "t.jplace"), String) == jp

        @test _tr_request("DELETE", "$base/t.jplace").status == 200
        @test !isfile(joinpath(tmp, "projects", "StudyT", "trees", ".views", "t.jplace.json"))
        @test isempty(_tr_json(_tr_request("GET", base)))
    finally
        SV.ServerState.set_root!(old_root)
        rm(tmp; recursive=true, force=true)
    end
end

@testset "Request guards" begin
    @test !SV._valid_name("evil\n")
    @test !SV._valid_name("..")
    @test SV._valid_name("run_1")

    # Oxygen's internalrequest skips the app middleware unless it is passed in.
    guarded(req) = SV.Oxygen.internalrequest(req; middleware=[SV._bad_request_middleware])
    req = SV.HTTP.Request("POST", "/api/v1/studies/S/runs/R/results/tables/..%2F..%2Fx/query",
                          ["Content-Type" => "application/json"], "{}")
    @test guarded(req).status == 400

    req = SV.HTTP.Request("GET", "/api/v1/studies/S/runs/R/results/tables?group=..%2Fx")
    @test_throws SV.BadRequest SV._req_group(req)
    @test isnothing(SV._req_group(SV.HTTP.Request("GET", "/api/v1/studies/S/runs/R")))
end
