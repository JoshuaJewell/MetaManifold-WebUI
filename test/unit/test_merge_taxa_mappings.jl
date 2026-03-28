@testset "merge_taxa mappings and corrections" begin

    @testset "_apply_single_mapping! remaps values" begin
        df = DataFrame(
            Division   = ["Rhizaria", "Alveolata", "Fungi"],
            Supergroup = ["SAR",      "SAR",       "Opisthokonta"],
        )
        mapping = Dict("Rhizaria" => "Rhizaria", "Alveolata" => "Alveolata")
        TaxonomyTableTools._apply_single_mapping!(df, "Division", "Supergroup", mapping)

        # Mapped values override target column
        @test df.Supergroup[1] == "Rhizaria"    # was SAR, Division=Rhizaria -> mapped
        @test df.Supergroup[2] == "Alveolata"   # was SAR, Division=Alveolata -> mapped
        @test df.Supergroup[3] == "Opisthokonta" # Fungi not in mapping -> unchanged
    end

    @testset "_apply_single_mapping! handles missing values" begin
        df = DataFrame(
            Division   = [missing, "Rhizaria"],
            Supergroup = ["SAR",   "SAR"],
        )
        mapping = Dict("Rhizaria" => "Rhizaria")
        TaxonomyTableTools._apply_single_mapping!(df, "Division", "Supergroup", mapping)
        @test df.Supergroup[1] == "SAR"       # missing Division -> keeps target
        @test df.Supergroup[2] == "Rhizaria"  # mapped
    end

    @testset "_apply_single_mapping! skips when columns missing" begin
        df = DataFrame(Domain = ["Eukaryota"], sample1 = [10])
        # Should warn and not error when columns don't exist
        mapping = Dict("x" => "y")
        @test_logs (:warn, r"Column remapping skipped") TaxonomyTableTools._apply_single_mapping!(
            df, "Nonexistent", "AlsoMissing", mapping)
        @test nrow(df) == 1  # unchanged
    end

    @testset "_apply_single_mapping! empty mapping is no-op" begin
        df = DataFrame(Division = ["Rhizaria"], Supergroup = ["SAR"])
        TaxonomyTableTools._apply_single_mapping!(df, "Division", "Supergroup", Dict())
        @test df.Supergroup[1] == "SAR"  # unchanged
    end

    @testset "_apply_mappings! with new format (vector of entries)" begin
        df = DataFrame(
            Division   = ["Rhizaria", "Fungi"],
            Supergroup = ["SAR",      "Opisthokonta"],
            Class      = ["Cercozoa",  "Ascomycota"],
            Order      = ["Euglyphida", "Saccharomycetales"],
        )
        mapping_config = [
            Dict("source_column" => "Division", "target_column" => "Supergroup",
                 "values" => Dict("Rhizaria" => "Rhizaria")),
            Dict("source_column" => "Class", "target_column" => "Order",
                 "values" => Dict("Cercozoa" => "Cercozoa_order")),
        ]
        TaxonomyTableTools._apply_mappings!(df, mapping_config)
        @test df.Supergroup[1] == "Rhizaria"
        @test df.Order[1] == "Cercozoa_order"
        @test df.Order[2] == "Saccharomycetales"  # not in mapping
    end

    @testset "_apply_mappings! with legacy format (flat dict)" begin
        df = DataFrame(
            Division   = ["Rhizaria"],
            Supergroup = ["SAR"],
        )
        legacy_mapping = Dict("Rhizaria" => "Rhizaria")
        TaxonomyTableTools._apply_mappings!(df, legacy_mapping)
        @test df.Supergroup[1] == "Rhizaria"
    end

    @testset "_apply_mappings! skips incomplete entries" begin
        df = DataFrame(Division = ["Rhizaria"], Supergroup = ["SAR"])
        # Missing target_column -> should skip
        mapping_config = [
            Dict("source_column" => "Division", "values" => Dict("Rhizaria" => "X")),
        ]
        TaxonomyTableTools._apply_mappings!(df, mapping_config)
        @test df.Supergroup[1] == "SAR"  # unchanged
    end

    @testset "filter_table applies mappings before filters" begin
        df = DataFrame(
            SeqName    = ["seq1", "seq2"],
            Division   = ["Rhizaria", "Fungi"],
            Supergroup = ["SAR",      "Opisthokonta"],
            sample1    = [10, 20],
        )
        tmp = tempname() * ".yml"
        write(tmp, """
mappings:
  - source_column: Division
    target_column: Supergroup
    values:
      Rhizaria: Rhizaria

filters:
  - column: Supergroup
    pattern: "Rhizaria"
    action: keep
""")
        result = filter_table(df, tmp)
        rm(tmp)
        @test nrow(result) == 1
        @test result.SeqName[1] == "seq1"
        @test result.Supergroup[1] == "Rhizaria"
    end

    @testset "merge_taxonomy_counts applies database corrections" begin
        # Database with a correction rule
        db = DatabaseMeta("test_db", ["Domain", "Phylum"], "generic",
                          [Dict{String,Any}(
                              "source" => "Domain",
                              "target" => "Phylum",
                              "values" => Dict("Bacteria" => "Bacteria_corrected")
                          )],
                          Set{String}())

        vsearch_tmp = tempname() * ".tsv"
        write(vsearch_tmp, "seq1\tBacteria;Firmicutes\t95.0\n")

        counts_tmp = tempname() * ".csv"
        write(counts_tmp, "SeqName,s1\nseq1,100\n")

        df = TaxonomyTableTools.merge_taxonomy_counts(vsearch_tmp, counts_tmp, db)
        rm(vsearch_tmp); rm(counts_tmp)

        # Correction should have been applied: Domain=Bacteria maps Phylum to Bacteria_corrected
        @test df[1, :Phylum] == "Bacteria_corrected"
    end

end
