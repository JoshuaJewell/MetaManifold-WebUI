@testset "DiversityMetrics" begin

    @testset "richness" begin
        @test richness([1, 2, 3])    == 3
        @test richness([0, 1, 0, 2]) == 2
        @test richness([0, 0, 0])    == 0
        @test richness(Int[])        == 0
        @test richness([5])          == 1
    end

    @testset "shannon" begin
        @test shannon([10]) ≈ 0.0
        @test shannon([0, 0, 0]) ≈ 0.0
        @test shannon([50, 50]) ≈ log(2) atol=1e-10
        @test shannon([25, 25, 25, 25]) ≈ log(4) atol=1e-10
        @test shannon([50, 0, 50]) ≈ log(2) atol=1e-10

        expected = -(0.9*log(0.9) + 0.1*log(0.1))
        @test shannon([90, 10]) ≈ expected atol=1e-10
    end

    @testset "simpson" begin
        @test simpson([100]) ≈ 0.0
        @test simpson([0, 0]) ≈ 0.0
        @test simpson([50, 50]) ≈ 0.5 atol=1e-10
        @test simpson([25, 25, 25, 25]) ≈ 0.75 atol=1e-10
        @test simpson([50, 0, 50]) ≈ 0.5 atol=1e-10

        for counts in ([1,1,1,1,1], [100,1], [1,2,3,4,5])
            d = simpson(counts)
            @test 0 <= d < 1
        end
    end

    @testset "Normalisation" begin

        @testset "rarefy" begin
            mat = [100.0 200.0 300.0;
                   150.0 150.0 200.0]
            depth = 50

            result = rarefy(mat; depth, seed=42)
            @test all(sum(result; dims=2) .≈ Float64(depth))
            @test all(result .== floor.(result))
            @test rarefy(mat; depth, seed=42) == rarefy(mat; depth, seed=42)
            @test rarefy(mat; depth, seed=42) != rarefy(mat; depth, seed=99)
        end

        @testset "normalise_counts" begin
            # 3 samples by 2 features; lib sizes: 300, 400, 100
            mat = [100.0 200.0;
                   300.0 100.0;
                    50.0  50.0]

            # "none": input unchanged, all row indices returned
            r = normalise_counts(mat; method="none", depth=0, seed=1)
            @test r.mat == Matrix{Float64}(mat)
            @test r.kept == [1, 2, 3]

            # "rarefy" depth=0 (auto): resolved depth = min positive lib size = 100
            r_rar = normalise_counts(mat; method="rarefy", depth=0, seed=42)
            @test r_rar.kept == [1, 2, 3]
            @test all(sum(r_rar.mat; dims=2) .≈ 100.0)

            # Fixed depth=150: row 3 (lib_size=100) is below threshold, dropped
            r_drop = normalise_counts(mat; method="rarefy", depth=150, seed=42)
            @test r_drop.kept == [1, 2]
            @test size(r_drop.mat, 1) == 2
            @test all(sum(r_drop.mat; dims=2) .≈ 150.0)

            # Edge: single sample - no crash
            single = reshape([10.0, 20.0, 30.0], 1, 3)
            r_single = normalise_counts(single; method="rarefy", depth=0, seed=42)
            @test r_single.kept == [1]
            @test size(r_single.mat, 1) == 1

            # Edge: all-zero sample is dropped (lib_size=0 < resolved_depth of positive min)
            with_zero = [0.0 0.0; 50.0 50.0]
            r_zero = normalise_counts(with_zero; method="rarefy", depth=0, seed=42)
            @test r_zero.kept == [2]
            @test size(r_zero.mat, 1) == 1
        end

        @testset "srs" begin
            mat = [100.0 200.0 300.0;
                   150.0 150.0 200.0]
            depth = 60

            result = srs(mat; depth, seed=42)
            @test all(sum(result; dims=2) .≈ Float64(depth))
            @test all(result .== floor.(result))
            @test all(result .>= 0)

            # Scaling, not resampling: an exactly divisible library lands on its
            # scaled counts with nothing left to rank, whatever the seed.
            @test result[1, :] == [10.0, 20.0, 30.0]
            @test srs(mat; depth, seed=42) == srs(mat; depth, seed=99)

            # SRS never invents reads for a taxon the sample did not observe.
            sparse = [0.0 7.0 3.0]
            @test srs(sparse; depth=5, seed=1)[1, 1] == 0.0

            # The rounding shortfall goes to the largest discarded fractions.
            # Scaling [5,3,2] (total 10) to 7 gives [3.5, 2.1, 1.4] -> floors
            # [3,2,1] with one read left; fractions .5 > .4 > .1 award it to
            # feature 1.
            ranked = srs([5.0 3.0 2.0]; depth=7, seed=1)
            @test ranked[1, :] == [4.0, 2.0, 1.0]

            # An all-zero library has nothing to scale and stays zero.
            @test srs([0.0 0.0]; depth=5, seed=1) == [0.0 0.0]

            # Reproducible when a tie does force the RNG.
            tied = [4.0 4.0 4.0]
            @test srs(tied; depth=5, seed=7) == srs(tied; depth=5, seed=7)
            @test sum(srs(tied; depth=5, seed=7)) == 5.0
        end

        @testset "normalise_counts with srs" begin
            mat = [100.0 200.0;
                   300.0 100.0;
                    50.0  50.0]

            r_auto = normalise_counts(mat; method="srs", depth=0, seed=42)
            @test r_auto.kept == [1, 2, 3]
            @test all(sum(r_auto.mat; dims=2) .≈ 100.0)

            # Fixed depth drops the short library, as rarefaction does: neither
            # method can scale a library up to a depth it never reached.
            r_drop = normalise_counts(mat; method="srs", depth=150, seed=42)
            @test r_drop.kept == [1, 2]
            @test all(sum(r_drop.mat; dims=2) .≈ 150.0)

            @test_throws ErrorException normalise_counts(mat; method="nope", depth=0, seed=1)
        end

    end  # Normalisation

    @testset "Transforms" begin

        @testset "hellinger" begin
            mat = [1.0 0.0 3.0;
                   4.0 4.0 0.0]
            h = hellinger(mat)

            # Each row is a unit vector: the squared entries are the relative
            # abundances, which sum to one.
            @test all(sum(h .^ 2; dims=2) .≈ 1.0)
            @test h[1, :] ≈ [0.5, 0.0, sqrt(0.75)]
            # A zero entry stays zero, so absences survive the transform.
            @test h[1, 2] == 0.0
            # An empty library has no relative abundances to take; it stays zero
            # rather than dividing by zero.
            @test hellinger([0.0 0.0]) == [0.0 0.0]

            # Scale invariance: doubling a library leaves its profile unchanged.
            @test hellinger(mat) ≈ hellinger(2 .* mat)
        end

        @testset "transform_counts" begin
            mat = [1.0 3.0; 4.0 4.0]
            @test transform_counts(mat; method="none") == mat
            @test transform_counts(mat; method="hellinger") == hellinger(mat)
            @test_throws ErrorException transform_counts(mat; method="nope")
        end

    end  # Transforms

end
