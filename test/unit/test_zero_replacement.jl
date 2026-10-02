# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Multiplicative zero replacement: a hand-computed known answer, the
# invariants (totals kept, observed ratios kept, no zero left), every refusal,
# and a control showing the invariant check can fail. Pure Julia; no R needed.

using MetaManifold
using MetaManifold.ZeroReplacement
using Test

@testset "Zero replacement" begin

    @testset "Known answer, computed by hand" begin
        # Detection limits: taxon 1 -> 4, taxon 2 -> 2, taxon 3 -> 6. With
        # delta = 0.5, sample 1 replaces one zero worth 0.5*4 = 2 of a total of
        # 10 (Δ = 0.2) and sample 2 one worth 0.5*2 = 1 (Δ = 0.1).
        x = [0 2 8; 4 0 6]
        out = multiplicative_replacement(x; delta=0.5)
        @test out.counts ≈ [2.0 1.6 6.4; 3.6 1.0 5.4]
        @test out.detection_limits == [4.0, 2.0, 6.0]
        @test out.zeros_replaced == 2
        @test out.imputed_mass ≈ [0.2, 0.1]
        @test out.delta == 0.5
    end

    @testset "Invariants on a larger table" begin
        x = [0 3 0 12 5; 7 0 1 0 9; 2 2 2 2 2; 0 0 4 30 1]
        out = multiplicative_replacement(x)
        @test replacement_invariants_hold(x, out.counts)
        @test vec(sum(out.counts; dims=2)) ≈ vec(sum(x; dims=2))
        @test all(>(0), out.counts)
        # Every replaced value is below its taxon's detection limit.
        for i in axes(x, 1), j in axes(x, 2)
            x[i, j] == 0 && @test out.counts[i, j] < out.detection_limits[j]
        end
        # A sample without zeros is returned unchanged.
        @test out.counts[3, :] == Float64.(x[3, :])
        @test out.imputed_mass[3] == 0.0
        @test out.zeros_replaced == 6
    end

    @testset "The invariant check can fail (control)" begin
        x = [0 2 8; 4 0 6]
        good = multiplicative_replacement(x; delta=0.5).counts
        bad_total = copy(good); bad_total[1, 1] += 1
        @test !replacement_invariants_hold(x, bad_total)
        bad_ratio = copy(good); bad_ratio[1, 2] += 0.1; bad_ratio[1, 3] -= 0.1
        @test !replacement_invariants_hold(x, bad_ratio)
        @test !replacement_invariants_hold(x, Float64.(x))
        @test !replacement_invariants_hold(x, good[:, 1:2])
    end

    @testset "Refusals" begin
        x = [0 2 8; 4 0 6]
        for d in (0, 1, -0.5, 1.5, NaN)
            @test_throws ArgumentError multiplicative_replacement(x; delta=d)
        end
        @test_throws ArgumentError multiplicative_replacement([-1 2; 3 4])
        @test_throws ArgumentError multiplicative_replacement([Inf 2; 3 4])
        @test_throws ArgumentError multiplicative_replacement([0 2; 0 4])
        @test_throws ArgumentError multiplicative_replacement([0 0; 3 4])
        @test_throws ArgumentError multiplicative_replacement(x; samples=["a"], taxa=["p", "q", "r"])

        msg = try
            multiplicative_replacement([0 0 1; 50 50 1]; samples=["s1", "s2"], taxa=["p", "q", "r"])
            ""
        catch e
            e.msg
        end
        # Sample s1: limits 50 + 50 against a total of 1, so delta must stay below 0.01.
        @test occursin("'s1'", msg)
        @test occursin("below 0.01", msg)

        msg = try multiplicative_replacement([0 2; 0 4]; taxa=["ghost", "real"]); "" catch e; e.msg end
        @test occursin("'ghost'", msg)
    end
end
