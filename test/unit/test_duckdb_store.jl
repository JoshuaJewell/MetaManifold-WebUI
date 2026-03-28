@testset "DuckDBStore" begin

    @testset "load_results_db loads CSVs into DuckDB" begin
        dir = mktempdir()
        write(joinpath(dir, "merged.csv"), "SeqName,Domain,s1,s2\nseq1,Eukaryota,100,50\nseq2,Bacteria,30,20\n")
        write(joinpath(dir, "protist.csv"), "SeqName,Domain,s1,s2\nseq1,Eukaryota,100,50\n")

        db_path = load_results_db(dir)
        @test isfile(db_path)
        @test endswith(db_path, "results.duckdb")

        # Verify tables were created
        with_results_db(dir) do con
            result = DataFrame(DBInterface.execute(con,
                "SELECT table_name FROM information_schema.tables ORDER BY table_name"))
            tables = string.(result.table_name)
            @test "merged" in tables
            @test "protist" in tables

            # Verify data integrity
            df = DataFrame(DBInterface.execute(con, "SELECT * FROM merged"))
            @test nrow(df) == 2
            @test "SeqName" in names(df)
        end
        rm(dir; recursive=true)
    end

    @testset "load_results_db with swarm_dir" begin
        merge_dir = mktempdir()
        swarm_dir = mktempdir()
        write(joinpath(merge_dir, "merged.csv"), "SeqName,s1\nseq1,100\n")
        write(joinpath(swarm_dir, "cluster_membership.csv"), "OTU,ASV\notu1,seq1\n")

        load_results_db(merge_dir; swarm_dir)

        with_results_db(merge_dir) do con
            result = DataFrame(DBInterface.execute(con,
                "SELECT table_name FROM information_schema.tables ORDER BY table_name"))
            tables = string.(result.table_name)
            @test "merged" in tables
            @test "cluster_membership" in tables
        end
        rm(merge_dir; recursive=true)
        rm(swarm_dir; recursive=true)
    end

    @testset "load_results_db overwrites stale DB" begin
        dir = mktempdir()
        write(joinpath(dir, "merged.csv"), "SeqName,s1\nseq1,100\n")
        load_results_db(dir)
        # Modify CSV and reload
        write(joinpath(dir, "merged.csv"), "SeqName,s1\nseq1,100\nseq2,200\n")
        load_results_db(dir)

        with_results_db(dir) do con
            df = DataFrame(DBInterface.execute(con, "SELECT * FROM merged"))
            @test nrow(df) == 2
        end
        rm(dir; recursive=true)
    end

    @testset "with_results_db errors on missing file" begin
        dir = mktempdir()
        @test_throws ErrorException with_results_db(identity, dir)
        rm(dir; recursive=true)
    end

end
