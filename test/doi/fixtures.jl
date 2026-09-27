# SPDX-License-Identifier: MPL-2.0
# Explicit synthetic contract fixtures, never scientific estimates or live deposits.
const FAKE_TOKEN = "not-a-real-token-DO-NOT-LOG"
const CONFIG_ID = "12345678-1234-4234-8234-123456789012"
const RESULT_ID = "12345678-1234-4234-8234-123456789013"

function metadata_fixture()
    Dict{String,Any}("title" => "Reproducible analysis fixture", "description" => "Synthetic publication protocol fixture, not a scientific result.",
        "creators" => [Dict("name" => "Example, Ada", "orcid" => "0000-0002-1825-0097")],
        "license" => "CC-BY-4.0", "publication_date" => "2026-01-01", "version" => "1.0.0",
        "github_release_url" => "https://github.com/example/research/releases/tag/v1.0.0",
        "github_project_url" => "https://github.com/users/example/projects/1")
end

function bundle_fixture(root; result=false, mock=false, dangerous=false)
    dir = joinpath(root, "bundle")
    mkpath(dir)
    cfg = Dict("id" => CONFIG_ID, "hash" => repeat("a", 64), "dangerous" => dangerous,
        "method" => "nb_glm", "created_at" => "2026-01-01T00:00:00", "created_by" => "Example, Ada")
    write(joinpath(dir, "analysis_config.json"), JSON3.write(cfg))
    write(joinpath(dir, "analysis_config.ncl"), "{ method = \"nb_glm\" }\n")
    write(joinpath(dir, "analysis_config_chora.deed"), "(repo-deed :schema-version \"1.0.0\" :dangerous #f)\n")
    write(joinpath(dir, "provenance.json"), "{\"fixture\":true}")
    write(joinpath(dir, "datacite.json"), "{}")
    write(joinpath(dir, "README.md"), "# Synthetic protocol fixture\n")
    write(joinpath(dir, "content_hash.txt"), cfg["hash"])
    dangerous && write(joinpath(dir, "DANGER_BANNER.txt"), "DANGER: synthetic unsafe configuration fixture")
    if result
        write(joinpath(dir, "analysis_result.json"), JSON3.write(Dict("id" => RESULT_ID,
            "config_id" => CONFIG_ID, "config_hash" => cfg["hash"], "hash" => repeat("b", 64),
            "method" => "nb_glm", "provenance" => Dict("mock" => mock),
            "results" => Dict("synthetic_fixture" => Dict("status" => "fixture")))))
    end
    Target.DOIBundles.write_checksums!(dir)
    return dir
end

mutable struct FakeZenodo
    deposit::Dict{String,Any}
    calls::Vector{Tuple{String,String}}
    uploaded::Vector{UInt8}
    faults::Vector{Function}
    publish_done::Bool
    environment::String
end
FakeZenodo(; environment="sandbox") = FakeZenodo(Dict{String,Any}(), Tuple{String,String}[], UInt8[], Function[], true, environment)

function fake_response(value; status=200, headers=Pair{String,String}[])
    HTTP.Response(status, headers; body=JSON3.write(value))
end

function transport(fake::FakeZenodo)
    return function(method, url, headers, body)
        @test get(Dict(headers), "Authorization", "") == "Bearer " * FAKE_TOKEN
        @test !occursin(FAKE_TOKEN, url)
        @test !occursin("access_token", url)
        @test startswith(url, Z.ORIGINS[fake.environment] * "/api/")
        push!(fake.calls, (method, url))
        if !isempty(fake.faults)
            fault = popfirst!(fake.faults)
            response = fault(method, url, body)
            !isnothing(response) && return response
        end
        if method == "POST" && endswith(url, "/deposit/depositions")
            @test isempty(fake.deposit) # a duplicate create is a test failure
            metadata = JSON3.read(String(body), Dict{String,Any})["metadata"]
            prefix = fake.environment == "sandbox" ? "10.5072/zenodo." : "10.5281/zenodo."
            metadata["prereserve_doi"] = Dict("doi" => prefix * "101", "recid" => 101)
            fake.deposit = Dict{String,Any}("id" => 101, "state" => "unsubmitted", "submitted" => false,
                "metadata" => metadata, "files" => Any[], "links" => Dict("bucket" => Z.ORIGINS[fake.environment] * "/api/files/11111111-1111-4111-8111-111111111111"))
            return fake_response(fake.deposit; status=201)
        elseif method == "GET" && endswith(url, "/101")
            return fake_response(fake.deposit)
        elseif method == "PUT" && occursin("/api/files/", url)
            @test body isa IO # streaming upload, not an in-memory ZIP body
            fake.uploaded = read(body)
            file = Dict{String,Any}("name" => last(split(url, '/')), "checksum" => bytes2hex(md5(fake.uploaded)), "filesize" => string(length(fake.uploaded)))
            fake.deposit["files"] = [file]
            return fake_response(Dict("key" => file["name"], "checksum" => "md5:" * file["checksum"], "size" => length(fake.uploaded)); status=201)
        elseif method == "POST" && endswith(url, "/101/actions/publish")
            if fake.publish_done
                mark_published!(fake)
            end
            return fake_response(fake.deposit; status=202)
        end
        error("Unexpected synthetic Zenodo request: $method $url")
    end
end

function mark_published!(fake)
    fake.deposit["state"] = "done"
    fake.deposit["submitted"] = true
    fake.deposit["doi"] = fake.deposit["metadata"]["prereserve_doi"]["doi"]
    fake.deposit["record_id"] = 101
end
client(fake) = Z.Client(FAKE_TOKEN; environment=fake.environment, transport=transport(fake), sleeper=_ -> nothing)
count_calls(fake, method, suffix) = count(c -> c[1] == method && endswith(c[2], suffix), fake.calls)

function prepared_fixture(f; result=false, dangerous=false, environment="sandbox")
    mktempdir() do tmp
        source = bundle_fixture(tmp; result, dangerous)
        root = joinpath(tmp, "publications")
        fake = FakeZenodo(; environment)
        c = client(fake)
        prepared = Target.DOIPublications.prepare!(root, source, metadata_fixture(), c)
        f(tmp, root, source, fake, c, prepared)
    end
end

function publish_fixture(root, prepared, c)
    Target.DOIPublications.publish!(root, prepared["id"], c; confirmation=prepared["confirmation_phrase"],
        bundle_sha256=prepared["bundle_sha256"], acknowledge_public=true)
end

function captured_error(f)
    try f(); nothing catch e; e end
end
