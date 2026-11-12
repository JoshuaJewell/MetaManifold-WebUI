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

# Known answers for counts 1 and 2 (total 3), worked by hand:
#   proportions 1/3 and 2/3; richness 2;
#   Gini-Simpson 1 - (1^2 + 2^2) / 3^2 = 4/9, exactly;
#   Shannon -(1/3 ln 1/3 + 2/3 ln 2/3) = ln 3 - (2/3) ln 2.
# The Shannon digits come from `bc -l` with `scale=100; l(3) - (2/3)*l(2)`,
# not from Julia or MPFR. The 256-bit result must agree to 1e-65.
const HAND_SHANNON_1_2 = "0.6365141682948128184504238226170746592638238015825792823209076606418985452388124898030402034827862971"

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
        high, high_empty = report.high_precision.samples
        @test high.shannon.kind == "approximate" && high.shannon.precision_bits == 256
        setprecision(BigFloat, 512) do
            @test abs(parse(BigFloat, high.shannon.value) - parse(BigFloat, HAND_SHANNON_1_2)) <= big"1e-65"
            @test abs(parse(BigFloat, high.gini_simpson.value) - big(4) / 9) <= big"1e-65"
        end
        @test high_empty.status == "undefined_zero_total"
        @test isnothing(high_empty.shannon) && isnothing(high_empty.gini_simpson)
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

@testset "Reference error distinguishes a threshold-straddling excess" begin
    NP = MetaManifold.NumericPolicy
    observed = first(NP.with_precision(; bits=512) do
        BigFloat(1) + BigFloat(2)^(-496)
    end)
    reference = first(NP.with_precision(; bits=4096) do
        BigFloat(1) - BigFloat(2)^(-800)
    end)
    before = precision(BigFloat)
    comparison = StatisticsSimulationScript.metric_errors(observed, reference)
    exact_observed = first(NP.with_precision(; bits=4096) do
        Rational{BigInt}(comparison.absolute)
    end)
    @test exact_observed == big(1)//(big(2)^496) + big(1)//(big(2)^800)
    @test comparison.absolute > StatisticsSimulationScript.tolerance(512)
    @test precision(comparison.absolute) == 4096
    @test precision(comparison.relative) == 4096
    @test precision(BigFloat) == before
    low_error = first(NP.with_precision(; bits=256) do
        abs(observed-reference)
    end)
    @test low_error == StatisticsSimulationScript.tolerance(512)
end
