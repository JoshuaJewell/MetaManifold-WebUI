@testset "DiversityMetrics" begin

    @testset "richness" begin
        @test richness([1, 2, 3])    == 3
        @test richness([0, 1, 0, 2]) == 2
        @test richness([0, 0, 0])    == 0
        @test richness(Int[])        == 0
        @test richness([5])          == 1
    end

    @testset "shannon" begin
        # Single species: no diversity, H = 0
        @test shannon([10]) ≈ 0.0

        # All zeros: defined as 0
        @test shannon([0, 0, 0]) ≈ 0.0

        # Two equal species: H = ln(2)
        @test shannon([50, 50]) ≈ log(2) atol=1e-10

        # Four equal species: H = ln(4)
        @test shannon([25, 25, 25, 25]) ≈ log(4) atol=1e-10

        # Zeros ignored
        @test shannon([50, 0, 50]) ≈ log(2) atol=1e-10

        # Known asymmetric value: p=[0.9, 0.1], H = -(0.9*ln(0.9) + 0.1*ln(0.1))
        expected = -(0.9*log(0.9) + 0.1*log(0.1))
        @test shannon([90, 10]) ≈ expected atol=1e-10
    end

    @testset "simpson" begin
        # Single species: all dominance, D = 0
        @test simpson([100]) ≈ 0.0

        # All zeros
        @test simpson([0, 0]) ≈ 0.0

        # Two equal species: D = 1 - 2*(0.5^2) = 0.5
        @test simpson([50, 50]) ≈ 0.5 atol=1e-10

        # Four equal species: D = 1 - 4*(0.25^2) = 0.75
        @test simpson([25, 25, 25, 25]) ≈ 0.75 atol=1e-10

        # Zeros ignored
        @test simpson([50, 0, 50]) ≈ 0.5 atol=1e-10

        # Range: always 0 <= D < 1
        for counts in ([1,1,1,1,1], [100,1], [1,2,3,4,5])
            d = simpson(counts)
            @test 0 <= d < 1
        end
    end

end
