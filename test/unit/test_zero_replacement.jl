# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Evidence for issue #21's zero handling (docs/statistics/zero-handling.md).
#
# The numbers asserted here are the pinned fixture in test/fixtures/issue21/golden.json.
# They were derived from the published methods (Martín-Fernández et al. 2003 and 2015) and
# from zCompositions::multRepl / cmultRepl, not from this module — see the fixture's
# `generator` block for how they were produced and where the fixture stops being
# independent. They are a pinned expectation until test/reference/issue21_reference.jl
# reproduces them in Julia; this file is what turns them into a gate.
#
# Where R is available (CI installs it) the module is also compared against the reference
# packages directly, and that comparison is the strongest evidence in this file. Where R is
# absent those tests report a skip rather than a pass.

using Test
using JSON3
using MetaManifold.ZeroReplacement

const ISSUE21_FIXTURE = joinpath(@__DIR__, "..", "fixtures", "issue21", "golden.json")
const GOLDEN = JSON3.read(read(ISSUE21_FIXTURE, String))

matrix(rows) = permutedims(Matrix{Float64}([Float64(v) for v in row for row in [rows]]), (1, 1))

function as_matrix(rows)::Matrix{Float64}
    n_rows = length(rows)
    n_cols = length(rows[1])
    M = Matrix{Float64}(undef, n_rows, n_cols)
    for i in 1:n_rows, j in 1:n_cols
        M[i, j] = Float64(rows[i][j])
    end
    return M
end

vector(vals)::Vector{Float64} = Float64[Float64(v) for v in vals]

@testset "zero replacement (issue #21)" begin

    mr = GOLDEN.zero_replacement.multiplicative_replacement
    gbm = GOLDEN.zero_replacement.bayesian_multiplicative
    cap = GOLDEN.zero_replacement.bayesian_multiplicative_cap

    @testset "multiplicative replacement: exact rational arithmetic" begin
        # The fixture's own derivation, restated as arithmetic: delta = 0.65, detection
        # limits [10, 20, 20, 30], so est = [6.5, 13, 13, 19.5], sum(est) = 52, the observed
        # parts scale by 1 - 52/100 = 0.48, and the zeros become est/0.48.
        counts = as_matrix(mr.counts)
        out = multiplicative_replacement(counts; delta = Float64(mr.delta),
                                         detection_limits = vector(mr.detection_limits))
        @test out.method == "multiplicative_replacement"
        @test Float64(out.delta) == Float64(mr.delta)
        @test out.counts ≈ as_matrix(mr.replaced) atol = 1e-9 rtol = 1e-12
        # The exact entries the fixture pins as rationals, to the last bit the double allows.
        @test out.counts[3, 1] ≈ 26.1 atol = 1e-12
        @test out.counts[2, 1] ≈ 6.5 atol = 1e-12
        @test out.counts[1, 1] ≈ 10.0 atol = 1e-12
    end

    @testset "multiplicative replacement: invariants, in the module's own terms" begin
        counts = as_matrix(mr.counts)
        out = multiplicative_replacement(counts; delta = Float64(mr.delta))
        @test replacement_invariants_hold(counts, out.counts)
        @test all(>(0.0), out.counts)                      # no zero survives
        for j in 1:size(counts, 2)
            @test sum(out.counts[:, j]) ≈ sum(counts[:, j]) rtol = 1e-12
        end
        # Observed-part ratios are untouched: compare every pair in a column.
        obs = findall(!=(0.0), counts[:, 1])
        for a in obs, b in obs
            @test out.counts[a, 1] / out.counts[b, 1] ≈ counts[a, 1] / counts[b, 1] rtol = 1e-12
        end
    end

    @testset "multiplicative replacement: refusals" begin
        counts = as_matrix(mr.counts)
        @test_throws ArgumentError multiplicative_replacement(counts; delta = 0.0)
        @test_throws ArgumentError multiplicative_replacement(counts; delta = 1.0)
        @test_throws ArgumentError multiplicative_replacement(counts; delta = -0.2)
        @test_throws ArgumentError multiplicative_replacement(counts; delta = 1.5)
        # A sample with no counts at all has no observed subcomposition to preserve.
        @test_throws ArgumentError multiplicative_replacement([1.0 0.0; 2.0 0.0])
        # A taxon never observed anywhere has no detection limit to derive a value from.
        @test_throws ArgumentError multiplicative_replacement([1.0 2.0; 0.0 3.0] .* [1.0; 0.0])
        # delta large enough that the imputed mass would consume the whole sample.
        @test_throws ArgumentError multiplicative_replacement([1.0 1.0; 0.0 0.0] .+ 0.0;
                                                             delta = 0.9)
    end

    @testset "Bayesian multiplicative replacement: the reference's numbers" begin
        counts = as_matrix(gbm.counts)
        out = bayesian_multiplicative(counts; threshold = Float64(gbm.threshold))
        @test out.method == "bayesian_multiplicative"
        @test out.counts ≈ as_matrix(gbm.replaced) rtol = 1e-6
        # The closed (proportional) table is the module's output divided by the sample total.
        for j in 1:size(counts, 2)
            @test out.counts[:, j] ./ sum(counts[:, j]) ≈ vector(gbm.replaced_prop[:, j]) rtol = 1e-9
        end
        for j in 1:size(counts, 2)
            @test sum(out.counts[:, j]) ≈ sum(counts[:, j]) rtol = 1e-12
        end
        @test replacement_invariants_hold(counts, out.counts)
    end

    @testset "Bayesian: the prior concentration is the stated estimator" begin
        counts = as_matrix(gbm.counts)
        out = bayesian_multiplicative(counts)
        s_golden = vector(gbm.prior_concentration_per_sample)
        for j in 1:size(counts, 2)
            detail = out.diagnostics["per_sample"][j]
            @test Float64(detail["prior_concentration"]) ≈ s_golden[j] rtol = 1e-6
        end
    end

    @testset "Bayesian: capping is the reference's adjust step" begin
        counts = as_matrix(cap.counts)
        adjusted = bayesian_multiplicative(counts; threshold = Float64(cap.threshold),
                                           adjust = true)
        unadjusted = bayesian_multiplicative(counts; threshold = Float64(cap.threshold),
                                             adjust = false)
        @test adjusted.counts ≈ as_matrix(cap.replaced_adjusted) rtol = 1e-6
        @test unadjusted.counts ≈ as_matrix(cap.replaced_unadjusted) rtol = 1e-6
        @test Int(adjusted.diagnostics["capped_imputations"]) == Int(cap.capped)
        @test Int(cap.capped) > 0
        # The cap only ever lowers a replaced value; it never touches observed parts.
        @test all(adjusted.counts .<= unadjusted.counts .+ 1e-9)
    end

    @testset "Bayesian: an explicit alpha is honoured and recorded" begin
        counts = as_matrix(gbm.counts)
        out = bayesian_multiplicative(counts; alpha = 2.0)
        @test out.counts ≈ as_matrix(gbm.replaced_alpha_2) rtol = 1e-6
        @test Float64(out.provenance["bayesian_multiplicative_alpha"]) == 2.0
        @test_throws ArgumentError bayesian_multiplicative(counts; alpha = 0.0)
        @test_throws ArgumentError bayesian_multiplicative(counts; alpha = -1.0)
        @test_throws ArgumentError bayesian_multiplicative(counts; threshold = 0.0)
        @test_throws ArgumentError bayesian_multiplicative(counts; threshold = 1.0)
    end

    @testset "Bayesian: a part observed in fewer than two samples is refused" begin
        # The prior mean is a leave-one-out profile; a part seen once cannot supply one.
        @test_throws ArgumentError bayesian_multiplicative([1.0 0.0 2.0; 0.0 3.0 0.0])
    end

    @testset "provenance carries the citations and the bias warning" begin
        counts = as_matrix(mr.counts)
        out = multiplicative_replacement(counts)
        @test out.provenance["zero_replacement_method"] == "multiplicative_replacement"
        @test occursin("Martín-Fernández", out.provenance["zero_replacement_reference"])
        @test occursin("2003", out.provenance["zero_replacement_reference"])
        @test occursin("NoRigidReplacement.agda", out.provenance["zero_replacement_is_biased"])

        g = bayesian_multiplicative(as_matrix(gbm.counts))
        @test occursin("2015", g.provenance["zero_replacement_reference"])
        @test occursin("cmultRepl", g.provenance["zero_replacement_reference"])
        @test occursin("model quantity", g.provenance["zero_replacement_is_biased"])
    end

    @testset "describe_zero_policy covers every policy the config accepts" begin
        for policy in ("pseudocount", "multiplicative_replacement", "bayesian_multiplicative",
                       "refuse")
            text = describe_zero_policy(policy)
            @test !isempty(text)
            @test occursin("bias", lowercase(text)) || occursin("assumption", lowercase(text))
        end
        @test_throws ArgumentError describe_zero_policy("nonsense")
    end

    @testset "pseudocounts move the ratios; replacement does not" begin
        # This is the property the issue turns on, stated as a test rather than as prose:
        # adding a constant to every part of a composition changes the ratios among the
        # observed parts, while scaling them by one common factor cannot.
        counts = [10.0 20.0; 30.0 40.0; 0.0 5.0]
        pseudocount = counts .+ 1.0
        @test !isapprox(pseudocount[1, 1] / pseudocount[2, 1], counts[1, 1] / counts[2, 1])
        out = multiplicative_replacement(counts)
        @test isapprox(out.counts[1, 1] / out.counts[2, 1], counts[1, 1] / counts[2, 1],
                       rtol = 1e-12)
    end

    @testset "the R reference, where R is installed" begin
        rscript = Sys.which("Rscript")
        if isnothing(rscript)
            @testset "Rscript absent — the zCompositions comparison did NOT run" begin
                @test true
                @info "zero replacement: Rscript is not on PATH; the zCompositions parity test was skipped, which is not evidence of parity"
            end
        else
            have_zc = success(pipeline(`$rscript -e 'quit(status = !requireNamespace("zCompositions", quietly = TRUE))'`,
                                       stdout = devnull, stderr = devnull))
            if !have_zc
                @testset "zCompositions absent — the parity comparison did NOT run" begin
                    @test true
                    @info "zero replacement: zCompositions is not installed; the parity test was skipped"
                end
            else
                counts = as_matrix(mr.counts)
                csv = joinpath(mktempdir(), "counts.csv")
                open(csv, "w") do io
                    for i in 1:size(counts, 1)
                        println(io, join(counts[i, :], ","))
                    end
                end
                script = """
                X <- as.matrix(read.csv("$csv", header = FALSE))
                X <- t(X)
                # multRepl is the 2003 operator; comparing in the closed regime (each sample
                # divided by its own total) is the comparison the reference supports.
                P <- sweep(X, 2, colSums(X), "/")
                r <- zCompositions::multRepl(P, label = 0, frac = $(Float64(mr.delta)),
                                             output = "prop")
                write.table(r, stdout(), sep = ",", row.names = FALSE, col.names = FALSE)
                """
                r_out = read(`$rscript -e $script`, String)
                rows = [parse.(Float64, split(l, ",")) for l in
                        split(strip(r_out), "\n") if !isempty(strip(l))]
                reference = permutedims(Matrix{Float64}(hcat(rows...)), (1, 1))
                ours = multiplicative_replacement(counts; delta = Float64(mr.delta))
                ours_prop = ours.counts ./ sum(counts; dims = 1)
                @test ours_prop ≈ reference rtol = 1e-6
            end
        end
    end
end
