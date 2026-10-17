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

    @testset "load_results_db leaves no WAL behind, so tagging survives the swap" begin
        # A database renamed into place with its data still in the WAL arrived
        # with neither its tables nor its category columns: the WAL kept its
        # .building name beside it.
        dir = mktempdir()
        ranks = ["Domain", "Supergroup", "Division", "Subdivision", "Class",
                 "Order", "Family", "Genus", "Species"]
        header = join(vcat(["SeqName"], ranks, ranks .* "_boot", ["s1", "s2"]), ",")
        row(name, taxa, boot, a, b) = join(vcat([name], taxa, fill(string(boot), length(ranks)),
                                                [string(a), string(b)]), ",")
        euk = ["Eukaryota", "Obazoa", "Opisthokonta", "Fungi", "Ascomycota",
               "Saccharomycetales", "Saccharomycetaceae", "Saccharomyces", "Saccharomyces_cerevisiae"]
        bac = ["Bacteria", "Bacteria_X", "Pseudomonadota", "Pseudomonadota_X", "Gammaproteobacteria",
               "Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia_coli"]
        write(joinpath(dir, "merged.csv"),
              join([header, row("seq1", euk, 100, 10, 5), row("seq2", bac, 80, 3, 7)], "\n") * "\n")
        tagging = Dict("source" => "DADA2", "max_x" => -1, "category_sets" => ["default"],
                       "library_path" => joinpath(@__DIR__, "..", "..", "config", "defaults", "composition.yml"))

        load_results_db(dir; tagging)
        @test !isfile(joinpath(dir, "results.duckdb.building.wal"))
        @test !isfile(joinpath(dir, "results.duckdb.building"))
        with_results_db(dir) do con
            cols = string.(DataFrame(DBInterface.execute(con,
                "SELECT column_name FROM information_schema.columns WHERE table_name = 'merged'")).column_name)
            @test "Category__default" in cols
            @test nrow(DataFrame(DBInterface.execute(con, "SELECT * FROM merged"))) == 2
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
