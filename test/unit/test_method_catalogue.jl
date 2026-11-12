# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The catalogue is a published claim. These tests make a claim that the code
# cannot back fail here rather than in a reader's inference.
using Test
using MetaManifold
using MetaManifold.MethodCatalogue
using MetaManifold.NumericPolicy

const _CATALOGUE_SV = MetaManifold.Server

function _deep_copy(value)
    value isa AbstractDict && return Dict{String,Any}(k => _deep_copy(v) for (k, v) in pairs(value))
    value isa AbstractVector && return Any[_deep_copy(v) for v in value]
    value
end

@testset "Method catalogue with explicit pending review" begin
    loaded = load_method_catalogue()
    document = loaded.document
    @test isfile(loaded.path)
    @test basename(loaded.path) == "validated_statistics.yml"
    @test document["schema"] == CATALOGUE_SCHEMA
    @test document["catalogue"]["status"] == "pending_review"
    @test occursin("pending", document["catalogue"]["review_status"])
    @test endswith(document["catalogue"]["review_reference"], "/issues/65")
    @test length(loaded.fingerprint) == 64
    @test occursin(r"^[0-9a-f]{64}$", loaded.fingerprint)

    supported = supported_method_ids(document)
    @test method_by_id(document, "exact-descriptive-summaries")["status"] == "supported"
    @test method_by_id(document, "not-a-method") === nothing
    @test supported == ["exact-descriptive-summaries"]
    # Only implemented methods are catalogued.
    @test [m["id"] for m in document["methods"]] == supported
    @test CATALOGUE_STATUSES == ("supported",)
    unimplemented = deepcopy(document)
    unimplemented["methods"][1]["status"] = "unsupported"
    @test_throws ArgumentError validate_catalogue(unimplemented)

    for method in document["methods"]
        @test !isempty(method["warnings"])
        @test !isempty(method["references"])
        @test !isempty(method["diagnostics"])
        @test !isempty(method["evidence"])
        @test !isempty(method["implementation"])
        @test isfile(joinpath(dirname(loaded.path), "..", "..", method["implementation"]))
    end

    storage = catalogue_to_storage(document)
    @test storage["catalogue_sha256"] == loaded.fingerprint
    @test storage["limits"]["max_cells"] == MAX_CELLS
    @test storage["limits"]["max_denominator_bits"] == MAX_DENOMINATOR_BITS
    @test storage["limits"]["max_count_bits"] == MAX_COUNT_BITS
    @test storage["limits"]["max_round_digits"] == MAX_ROUND_DIGITS
    @test storage["limits"]["max_request_bytes"] == MAX_REQUEST_BYTES
    @test MAX_REQUEST_BYTES == 16_777_216
    @test MAX_REQUEST_BYTES == _CATALOGUE_SV._STATISTICS_REQUEST_LIMIT

    text = read(loaded.path, String)
    @test occursin("---\n", text)
    @test occursin("methods: [", text)
    @test occursin("limits: {", text)
    @test !occursin('\t', text)
    @test endswith(text, "\n")
end

@testset "Catalogue mutations are rejected" begin
    document = load_method_catalogue().document

    function rejected(mutator)
        copy = _deep_copy(document)
        mutator(copy)
        try
            validate_catalogue(copy)
            false
        catch err
            err isa ArgumentError
        end
    end

    @test rejected(doc -> doc["schema"] = "something-else")
    @test rejected(doc -> doc["catalogue"]["status"] = "draft")
    @test rejected(doc -> doc["catalogue"]["review_reference"] = "")
    @test rejected(doc -> doc["methods"] = Any[])
    @test rejected(doc -> doc["methods"][1]["status"] = "maybe")
    @test rejected(doc -> doc["methods"][1]["evidence"] = Any[])
    # Every entry is implemented, so one without an implementation or without
    # evidence is refused, as is any status other than supported.
    @test rejected(doc -> delete!(doc["methods"][1], "implementation"))
    @test rejected(doc -> doc["methods"][1]["status"] = "deferred")
    @test rejected(doc -> push!(doc["methods"], _deep_copy(doc["methods"][1])))
    @test rejected(doc -> doc["limits"]["max_cells"] = MAX_CELLS + 1)
    @test rejected(doc -> doc["limits"]["max_request_bytes"] = 1)
    @test rejected(doc -> doc["methods"][1]["effect_sizes"] = "")
    @test rejected(doc -> doc["methods"][1]["references"] = Any[""])
end

@testset "Catalogue endpoint returns canonical JCS" begin
    root = mktempdir()
    old = _CATALOGUE_SV.ServerState._root[]
    try
        _CATALOGUE_SV.ServerState.set_root!(root)
        req = _CATALOGUE_SV.HTTP.Request("GET", "/api/v1/statistics/catalogue", String[], "")
        response = _CATALOGUE_SV.Oxygen.internalrequest(req)
        @test response.status == 200
        bytes = response.body isa AbstractString ? Vector{UInt8}(codeunits(response.body)) :
                                                    Vector{UInt8}(response.body)
        parsed = MetaManifold.Canonical.ijson_read(bytes)
        @test parsed.schema == CATALOGUE_SCHEMA
        @test parsed.catalogue_sha256 == load_method_catalogue().fingerprint
        @test Int(parsed.method_count) == length(load_method_catalogue().document["methods"])
        @test String(MetaManifold.Canonical.jcs_bytes(parsed)) == String(copy(bytes))
    finally
        _CATALOGUE_SV.ServerState._root[] = old
        rm(root; recursive=true, force=true)
    end
end
