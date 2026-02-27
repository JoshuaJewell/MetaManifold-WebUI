@testset "merge_taxa helpers" begin

    @testset "seqnum" begin
        @test TaxonomyTableTools.seqnum("seq1")   == 1
        @test TaxonomyTableTools.seqnum("seq10")  == 10
        @test TaxonomyTableTools.seqnum("seq100") == 100
        @test TaxonomyTableTools.seqnum("otu1")   == 1
        @test TaxonomyTableTools.seqnum("otu42")  == 42
        # Non-matching strings sort to end
        @test TaxonomyTableTools.seqnum("other")  == typemax(Int)
        @test TaxonomyTableTools.seqnum(missing)  == typemax(Int)
        # Ordering is consistent with sort
        ids = ["seq3", "seq1", "seq10", "seq2"]
        @test sort(ids, by=TaxonomyTableTools.seqnum) == ["seq1", "seq2", "seq3", "seq10"]
    end

    @testset "vsearch_to_df" begin
        rows = [
            ["seq1", "Bacteria;Firmicutes;Clostridia", "95.5"],
            ["seq2", "Bacteria;Proteobacteria;Gammaproteobacteria", "88.0"],
        ]
        df = TaxonomyTableTools.vsearch_to_df(rows)
        @test nrow(df) == 2
        @test "SeqName" in names(df)
        @test "Pident"  in names(df)
        @test "sseqtax" in names(df)
        @test df.SeqName == ["seq1", "seq2"]
        @test df.Pident  ≈ [95.5, 88.0]

        # Short rows are skipped
        df2 = TaxonomyTableTools.vsearch_to_df([["seq1", "tax"]])  # only 2 fields
        @test nrow(df2) == 0
    end

    @testset "import_vsearch" begin
        tmp = tempname() * ".tsv"
        write(tmp, "seq1\tBacteria;Firmicutes\t95.0\nseq2\tBacteria;Proteobacteria\t88.0\n")
        rows = TaxonomyTableTools.import_vsearch(tmp)
        rm(tmp)
        @test length(rows) == 2
        @test rows[1] == ["seq1", "Bacteria;Firmicutes", "95.0"]
        @test rows[2] == ["seq2", "Bacteria;Proteobacteria", "88.0"]
    end

    @testset "filter_table - exclude" begin
        df = DataFrame(
            SeqName = ["seq1", "seq2", "seq3"],
            Domain  = ["Bacteria", "Eukaryota", "Bacteria"],
            sample1 = [10, 20, 30],
        )
        tmp = tempname() * ".yml"
        write(tmp, """
filters:
  - column: Domain
    pattern: "Bacteria"
    action: exclude
""")
        result = filter_table(df, tmp)
        rm(tmp)
        @test nrow(result) == 1
        @test result.SeqName == ["seq2"]
    end

    @testset "filter_table - keep" begin
        df = DataFrame(
            SeqName = ["seq1", "seq2", "seq3"],
            Domain  = ["Bacteria", "Eukaryota", "Bacteria"],
            sample1 = [10, 20, 30],
        )
        tmp = tempname() * ".yml"
        write(tmp, """
filters:
  - column: Domain
    pattern: "Bacteria"
    action: keep
""")
        result = filter_table(df, tmp)
        rm(tmp)
        @test nrow(result) == 2
        @test all(result.Domain .== "Bacteria")
    end

    @testset "filter_table - remove_empty" begin
        df = DataFrame(
            SeqName = ["seq1", "seq2", "seq3"],
            Domain  = ["Bacteria", "", missing],
            sample1 = [10, 20, 30],
        )
        tmp = tempname() * ".yml"
        write(tmp, "remove_empty:\n  - Domain\n")
        result = filter_table(df, tmp)
        rm(tmp)
        @test nrow(result) == 1
        @test result.SeqName == ["seq1"]
    end

    @testset "filter_table - regex" begin
        df = DataFrame(
            SeqName = ["seq1", "seq2", "seq3"],
            Species = ["Homo_sapiens", "Mus_musculus", "Canis_lupus"],
            sample1 = [1, 2, 3],
        )
        tmp = tempname() * ".yml"
        write(tmp, """
filters:
  - column: Species
    pattern: "^Homo"
    action: keep
    regex: true
""")
        result = filter_table(df, tmp)
        rm(tmp)
        @test nrow(result) == 1
        @test result.Species[1] == "Homo_sapiens"
    end

end
