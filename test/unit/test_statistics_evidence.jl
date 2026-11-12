# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

using Test, JSON3

module StatisticsMutationScript
include(joinpath(@__DIR__, "../../scripts/statistics/mutation_check.jl"))
end

module StatisticsPipelineEvidenceScript
include(joinpath(@__DIR__, "../../scripts/statistics/collect_pipeline_evidence.jl"))
end

@testset "VSEARCH platform suffixes retain strict release matching" begin
    pin = Dict("version" => "2.30.5", "archives" => Dict("linux-x86_64" => Dict()))
    check = StatisticsPipelineEvidenceScript._tool_pin_matches
    @test check("vsearch", "2.30.5", pin)
    @test check("vsearch", "2.30.5_linux_x86_64", pin)
    @test !check("vsearch", "2.30.6_linux_x86_64", pin)
    @test !check("vsearch", "2.30.5_unreviewed", pin)
    @test !check("cutadapt", "2.30.5_linux_x86_64", pin)
end

@testset "Pipeline evidence refuses stale and absent matrices" begin
    mktempdir() do directory
        path = joinpath(directory, "counts.csv")
        check = StatisticsPipelineEvidenceScript._require_fresh_matrix
        @test_throws ErrorException check(path, 1)
        write(path, "feature,sample\na,1\n")
        recorded = mtime(path)
        @test check(path, recorded) === nothing
        @test_throws ErrorException check(path, recorded + 1)
        @test_throws ErrorException check(path, 0)
        @test_throws ErrorException check(path, NaN)
        write(path, "")
        @test_throws ErrorException check(path, recorded)
    end
end

# Known answers for counts 1 and 2 (total 3), worked by hand:
#   total 3; richness 2; proportions 1/3 and 2/3, exactly.
@testset "R row-name CSVs give the hand-worked descriptive results" begin
    mktempdir() do directory
        csv = joinpath(directory, "counts.csv")
        provenance = joinpath(directory, "provenance.json")
        result = joinpath(directory, "result.json")
        # The pipeline's R writer leaves the feature-column heading empty.
        write(csv, ",sample,empty\nASV1,1,0\nASV2,2,0\n")
        write(provenance, "{\"source\":\"known unit fixture\"}")
        StatisticsPipelineEvidenceScript.biological_results(csv, provenance, result)
        report = JSON3.read(read(result, String))
        sample, empty = report.exact.samples
        @test sample.label == "sample" && sample.status == "success"
        @test sample.total == 3 && sample.counts == [1, 2] && sample.richness == 2
        @test sample.proportions == ["1/3", "2/3"]
        @test empty.label == "empty" && empty.status == "undefined_zero_total"
        @test empty.total == 0 && all(isnothing, empty.proportions)
        # A blank sample heading still must fail, even though the feature
        # heading may legitimately be blank.
        write(csv, ",\nASV1,1\n")
        @test_throws ErrorException StatisticsPipelineEvidenceScript.biological_results(csv, provenance, result)
    end
end

@testset "Mutation evidence rejects infrastructure failures" begin
    for mutant in StatisticsMutationScript.MUTANTS
        source = read(joinpath(StatisticsMutationScript.ROOT, mutant.path), String)
        @test length(split(source, mutant.old)) - 1 == 1
    end
    files = ["first test", "second"]
    @test collect(StatisticsMutationScript._test_command("julia", files))[end-1:end] == files
    @test collect(StatisticsMutationScript._git_status_command(files))[end-1:end] == files
    detected = "case: Test Failed at test.jl:1\nTest Summary: | Pass Fail Total\n"
    check = StatisticsMutationScript._expected_test_failure
    @test check(1, detected)
    @test !check(0, detected)
    @test !check(137, detected)
    @test !check(1, "ERROR: syntax error")
    @test !check(1, "case: Test Failed at test.jl:1")
    @test !check(1, detected * "case: Error During Test at test.jl:2")
    rejected = "case: Error During Test at test.jl:1\nArgumentError: reviewed limit changed\nTest Summary: | Error Total\n"
    @test !check(1, rejected)
    @test check(1, rejected; expected_rejection="reviewed limit changed")
    @test !check(1, rejected; expected_rejection="a different error")
    @test !check(137, rejected; expected_rejection="reviewed limit changed")
end
