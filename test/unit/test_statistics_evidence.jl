# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

using Test, JSON3

module StatisticsSimulationScript
include(joinpath(@__DIR__, "../../scripts/statistics/simulation_study.jl"))
end

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

@testset "R row-name CSVs reach the independent reference checker" begin
    mktempdir() do directory
        csv = joinpath(directory, "counts.csv")
        provenance = joinpath(directory, "provenance.json")
        result = joinpath(directory, "result.json")
        receipt = joinpath(directory, "check.json")
        # The pipeline's R writer leaves the feature-column heading empty.
        write(csv, ",sample,empty\nASV1,1,0\nASV2,2,0\n")
        write(provenance, "{\"source\":\"known unit fixture\"}")
        StatisticsPipelineEvidenceScript.biological_results(csv, provenance, result)
        python = something(Sys.which("python3"), Sys.which("python"))
        checker = joinpath(@__DIR__, "../references/check_biological_reference.py")
        @test success(`$python $checker $csv $result $receipt $provenance`)
        checked = JSON3.read(read(receipt, String))
        @test checked.status == "passed"
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

@testset "Simulation errors retain reference precision and restore MPFR state" begin
    NP = MetaManifold.NumericPolicy
    reference, _ = NP.with_precision(; bits=4096) do
        BigFloat(1) + BigFloat(2)^(-600) + BigFloat(2)^(-900)
    end
    observed, _ = NP.with_precision(; bits=512) do
        BigFloat(1)
    end
    expected, _ = NP.with_precision(; bits=4096) do
        difference = BigFloat(2)^(-600) + BigFloat(2)^(-900)
        (absolute=difference, relative=difference / reference)
    end
    # The operands agree at 256 bits, while the difference has two terms more
    # than 256 bits apart. Comparing at ambient precision loses the lower term.
    for ambient in (128, 256, 512)
        NP.with_precision(; bits=ambient) do
            errors = StatisticsSimulationScript.metric_errors(observed, reference)
            @test precision(errors.absolute) == 4096
            @test precision(errors.relative) == 4096
            @test errors == expected
            @test errors.absolute > 0
            @test precision(BigFloat) == ambient
            @test rounding(BigFloat) == RoundNearest
        end
    end
end

@testset "Simulation evidence is retained as readable JSON" begin
    mktempdir() do directory
        output = joinpath(directory, "simulation.json")
        @test StatisticsSimulationScript.main(["--replicates", "1", "--output", output]) == 0
        report = JSON3.read(read(output, String))
        @test report.complete === true
        @test report.comparison_bits == report.reference_bits == 4096
        @test isempty(report.exact_failures)
        @test isempty(report.metric_failures)
    end
end
