@testset "Databases" begin

    @testset "make_db_meta parses database config" begin
        dir = mktempdir()
        cfg_path = joinpath(dir, "databases.yml")
        write(cfg_path, """
databases:
  dir: ./databases
  pr2:
    levels: [Domain, Supergroup, Division, Class, Order, Family, Genus, Species]
    vsearch_format: pr2
    corrections:
      - source: Division
        target: Supergroup
        values:
          Rhizaria: Rhizaria
    dada2:
      uri: "https://example.com/pr2_dada2.fasta"
    vsearch:
      uri: "https://example.com/pr2_vsearch.fasta"
""")
        meta = Databases.make_db_meta(cfg_path, "pr2")
        @test meta.name == "pr2"
        @test meta.levels == ["Domain", "Supergroup", "Division", "Class",
                              "Order", "Family", "Genus", "Species"]
        @test meta.vsearch_format == "pr2"
        @test length(meta.corrections) == 1
        @test meta.corrections[1]["source"] == "Division"

        # noncounts includes taxonomy levels and their suffixed variants
        @test "Domain" in meta.noncounts
        @test "Domain_dada2" in meta.noncounts
        @test "Domain_vsearch" in meta.noncounts
        @test "Domain_boot" in meta.noncounts
        @test "SeqName" in meta.noncounts  # fixed noncount

        rm(dir; recursive=true)
    end

    @testset "make_db_meta errors on missing database" begin
        dir = mktempdir()
        cfg_path = joinpath(dir, "databases.yml")
        write(cfg_path, "databases:\n  dir: ./databases\n")
        @test_throws ErrorException Databases.make_db_meta(cfg_path, "nonexistent")
        rm(dir; recursive=true)
    end

    @testset "make_db_meta with no corrections" begin
        dir = mktempdir()
        cfg_path = joinpath(dir, "databases.yml")
        write(cfg_path, """
databases:
  dir: ./databases
  silva:
    levels: [Domain, Phylum, Class]
    vsearch_format: silva
    dada2:
      uri: "https://example.com/silva.fasta"
""")
        meta = Databases.make_db_meta(cfg_path, "silva")
        @test meta.name == "silva"
        @test meta.vsearch_format == "silva"
        @test isempty(meta.corrections)
        rm(dir; recursive=true)
    end

    @testset "ensure_databases with local paths" begin
        dir = mktempdir()
        db_dir = joinpath(dir, "databases")
        mkpath(db_dir)

        # Create fake local database files
        local_dada2 = joinpath(db_dir, "pr2_dada2.fasta")
        local_vsearch = joinpath(db_dir, "pr2_vsearch.fasta")
        write(local_dada2, ">seq1\nACGT\n")
        write(local_vsearch, ">seq1\nACGT\n")

        cfg_path = joinpath(dir, "databases.yml")
        write(cfg_path, """
databases:
  dir: $db_dir
  pr2:
    dada2:
      local: "$local_dada2"
    vsearch:
      local: "$local_vsearch"
""")

        resolved = Databases.ensure_databases(cfg_path)
        @test haskey(resolved, "pr2_dada2")
        @test haskey(resolved, "pr2_vsearch")
        @test resolved["pr2_dada2"] == local_dada2
        @test resolved["pr2_vsearch"] == local_vsearch

        rm(dir; recursive=true)
    end

    @testset "ensure_databases with missing config" begin
        dir = mktempdir()
        result = Databases.ensure_databases(joinpath(dir, "nonexistent.yml"))
        @test isempty(result)
        rm(dir; recursive=true)
    end

    @testset "ensure_databases skips null uri and local" begin
        dir = mktempdir()
        cfg_path = joinpath(dir, "databases.yml")
        write(cfg_path, """
databases:
  dir: ./databases
  pr2:
    dada2:
      uri: ~
      local: ~
""")
        resolved = Databases.ensure_databases(cfg_path)
        @test !haskey(resolved, "pr2_dada2")
        rm(dir; recursive=true)
    end

end
