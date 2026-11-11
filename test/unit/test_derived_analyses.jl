# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

using Test
using JSON3
using Sockets
using MetaManifold.Canonical
using MetaManifold.DerivedAnalyses
using MetaManifold.ExactSummaries
using MetaManifold.NumericPolicy

@testset "Derived records fail closed on invalid JSON and failed publication" begin
    DA = MetaManifold.DerivedAnalyses
    @test sprint(showerror, DerivedAnalysisError("invalid_record", "invalid JSON")) == "invalid JSON"
    @test_throws DerivedAnalysisError DA._validate_bytes(Vector{UInt8}(codeunits("{")))
    mktempdir() do root
        source, config, result = Dict("run" => "r"), Dict("mode" => "exact_counts"), Dict("value" => 1)
        id = derived_analysis_id(source, config, result)
        # A non-file destination blocks publication. The error must escape and
        # the temporary file must be removed, without changing that destination.
        destination = joinpath(root, "$id.json")
        mkdir(destination)
        @test_throws Base.IOError persist_derived_analysis(root; source, config, result)
        @test isdir(destination)
        @test readdir(root) == ["$id.json"]
    end
end

@testset "Derived-record I/O and storage limits fail explicitly" begin
    DA = MetaManifold.DerivedAnalyses
    mktempdir() do root
        # `read` on a directory raises nothing in Julia: it returns no bytes. So
        # a directory cannot stand in for an unreadable file. Opening a
        # Unix-domain socket path fails for every user, root included, and so
        # does not depend on filesystem permissions.
        unreadable = joinpath(root, "not-a-file")
        server = Sockets.listen(unreadable)
        error = try
            DA._existing(unreadable, repeat("0", 64))
            nothing
        catch err
            err
        finally
            close(server)
        end
        @test error isa DerivedAnalysisError
        @test error.code == "conflict"
        @test occursin("cannot read existing derived analysis", error.reason)

        # The limit applies to the complete canonical record, not just payload.
        # A payload at the limit necessarily makes the enclosing record exceed it.
        directory = joinpath(root, "oversized")
        error = try
            persist_derived_analysis(directory;
                source=Dict("run" => "r"), config=Dict("mode" => "exact_counts"),
                result=Dict("payload" => repeat("x", MAX_DERIVED_BYTES)))
            nothing
        catch err
            err
        end
        @test error isa DerivedAnalysisError
        @test error.code == "resource_limit"
        @test !isdir(directory) # refuse before allocating a publication directory
    end
end

@testset "Immutable content-addressed derived analyses" begin
    mktempdir() do root
        directory = joinpath(root, "derived")
        summary = exact_summary([1 2; 3 4]; sample_labels=["A", "B"],
                                feature_labels=["x", "y"])
        result = summary_to_storage(summary)
        source = Dict{String,Any}("study_id" => "s1", "run_id" => "r1",
                                  "taxonomic_level" => "genus")
        config = Dict{String,Any}("numeric_policy" => policy_to_storage(numeric_policy(:exact_counts)),
                                  "groups" => nothing, "source_format" => "direct-json")
        saved = persist_derived_analysis(directory; source, config, result)
        @test saved.created
        @test saved.id == derived_analysis_id(source, config, result)
        @test occursin(r"^[0-9a-f]{64}$", saved.id)
        @test isfile(saved.path)
        saved_bytes = read(saved.path)
        @test jcs_bytes(JSON3.read(saved_bytes)) == saved_bytes

        repeated = persist_derived_analysis(directory; source, config, result)
        @test repeated.id == saved.id
        @test !repeated.created
        @test read(repeated.path) == read(saved.path)
        loaded = load_derived_analysis(directory, saved.id)
        @test loaded["analysis_id"] == saved.id
        @test loaded["result"]["input_sha256"] == summary.input_sha256

        other_config = merge(config, Dict{String,Any}("groups" => ["g1", "g2"]))
        @test derived_analysis_id(source, other_config, result) != saved.id
        @test_throws DerivedAnalysisError load_derived_analysis(directory, "../escape")
        @test_throws DerivedAnalysisError load_derived_analysis(directory, repeat("0", 64))

        # Failure before publication leaves neither a final record nor a temp file.
        @test_throws ArgumentError persist_derived_analysis(directory;
            source, config, result=Dict("bad" => big(2)^80))
        @test length(readdir(directory)) == 1
    end
end

@testset "Concurrent publication and tamper recovery" begin
    mktempdir() do root
        directory = joinpath(root, "derived")
        source = Dict("run" => "r")
        config = Dict("mode" => "exact_counts")
        result = Dict("status" => "success", "values" => ["1/3", nothing])
        outcomes = Vector{Any}(undef, 16)
        @sync for i in eachindex(outcomes)
            Threads.@spawn outcomes[i] = persist_derived_analysis(
                directory; source, config, result)
        end
        @test count(x -> x.created, outcomes) == 1
        @test length(unique([x.id for x in outcomes])) == 1
        @test length(readdir(directory)) == 1
        id = outcomes[1].id

        tampered = joinpath(root, "tampered")
        mkpath(tampered)
        record = JSON3.read(read(joinpath(directory, "$id.json"), String))
        object = JSON3.read(jcs_bytes(record))
        # Write a same-ID record with altered content: loading must fail closed.
        bad = Dict{String,Any}("schema" => record["schema"],
            "analysis_id" => id, "source" => source, "config" => config,
            "result" => Dict("status" => "success", "values" => ["2/3", nothing]))
        write(joinpath(tampered, "$id.json"), jcs_bytes(bad))
        @test_throws DerivedAnalysisError load_derived_analysis(tampered, id)
        @test_throws DerivedAnalysisError persist_derived_analysis(
            tampered; source, config, result)
        # The verified original is still available after the failed conflicting use.
        @test load_derived_analysis(directory, id)["analysis_id"] == id
        @test object["analysis_id"] == id

        # A record with the right semantic ID but noncanonical bytes is not a
        # valid content-addressed I-JSON/JCS document.
        write(joinpath(tampered, "$id.json"), vcat(UInt8[0x20], read(joinpath(directory, "$id.json"))))
        @test_throws DerivedAnalysisError load_derived_analysis(tampered, id)
    end
end
