# Unit tests for the shared heatmap helpers and the count-column maxima.
using MetaManifold
SV = MetaManifold.Server
using DataFrames, DuckDB, DBInterface

@testset "Heatmap helpers" begin
    @testset "colour" begin
        @test SV._heat_colour(10, 10) == "#6fa8dc"
        @test SV._heat_colour(5, 10) == SV._heat_blend(SV._HEAT_HIGH, 0.5)
        @test SV._heat_colour(20, 10) == "#6fa8dc"               # clamped
        @test isnothing(SV._heat_colour(0, 10))
        @test isnothing(SV._heat_colour(3, 0))
        @test isnothing(SV._heat_colour("x", 10))
        @test isnothing(SV._heat_colour(missing, 10))
        @test SV._heat_colour(-5, 10; diverging=true) == SV._heat_blend(SV._HEAT_LOW, 0.5)
        @test SV._heat_blend(SV._HEAT_HIGH, 0.0) == "#ffffff"
    end

    @testset "count fills" begin
        df = DataFrame(Genus=["A", "B", "C"], s1=[10, 5, 0], s2=[1, 2, missing])
        @test isempty(SV._count_fills(df, ["s1", "s2"], "none"))

        col = SV._count_fills(df, ["s1", "s2"], "column")
        @test (1, 2, "#6fa8dc") in col                           # s1 max
        @test (2, 3, "#6fa8dc") in col                           # s2 max
        @test !any(f -> f[1] == 3, col)                          # zero and missing unfilled
        @test !any(f -> f[2] == 1, col)                          # label column untouched

        whole = SV._count_fills(df, ["s1", "s2"], "table")
        @test (2, 3, SV._heat_blend(SV._HEAT_HIGH, 0.2)) in whole
        @test (1, 2, "#6fa8dc") in whole
    end

    @testset "count maxima" begin
        con = DBInterface.connect(DuckDB.DB)
        DBInterface.execute(con, "CREATE TABLE t (Genus VARCHAR, s1 INTEGER, s2 INTEGER)")
        DBInterface.execute(con, "INSERT INTO t VALUES ('A', 10, 1), ('B', 5, 7), ('C', NULL, NULL)")
        @test SV._count_maxima(con, "t", ["s1", "s2"]) == Dict("s1" => 10.0, "s2" => 7.0)
        @test SV._count_maxima(con, "t", ["s1", "s2"], "WHERE Genus = ?", Any["B"]) ==
              Dict("s1" => 5.0, "s2" => 7.0)
        @test SV._count_maxima(con, "t", ["s1"], "WHERE Genus = ?", Any["C"]) == Dict("s1" => 0.0)
        @test isempty(SV._count_maxima(con, "t", String[]))
        DBInterface.close!(con)
    end
end
