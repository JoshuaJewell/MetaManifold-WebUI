@testset "PipelinePlotsPlotly" begin

    @testset "_palette_hex" begin
        c3 = PipelinePlotsPlotly._palette_hex(3)
        @test length(c3) == 3
        @test all(s -> startswith(s, "#") && length(s) == 7, c3)

        c10 = PipelinePlotsPlotly._palette_hex(10)
        @test length(c10) == 10
        @test allunique(c10)

        # Base palette returned unchanged for small n
        c7 = PipelinePlotsPlotly._palette_hex(7)
        @test c7[1] == "#E69F00"
    end

    @testset "taxa_bar_chart" begin
        df = DataFrame(
            SeqName = ["seq1", "seq2", "seq3"],
            Domain  = ["Eukaryota", "Eukaryota", "Bacteria"],
            s1      = [100, 200, 50],
            s2      = [30,   80, 20],
        )
        out = tempname() * ".json"
        PipelinePlotsPlotly.taxa_bar_chart(df, ["s1", "s2"], out;
            rank="Domain", rank_order=["Domain"])
        @test isfile(out)
        fig = JSON3.read(read(out, String))
        @test haskey(fig, :data) && haskey(fig, :layout)
        @test fig[:layout][:barmode] == "stack"
        # 2 unique domain labels -> 2 traces
        @test length(fig[:data]) == 2
        rm(out)

        # Absolute mode
        out_abs = tempname() * ".json"
        PipelinePlotsPlotly.taxa_bar_chart(df, ["s1", "s2"], out_abs;
            rank="Domain", relative=false, rank_order=["Domain"])
        fig_abs = JSON3.read(read(out_abs, String))
        @test !haskey(fig_abs[:layout][:yaxis], :range)
        rm(out_abs)

        # Subtitle appended
        out_sub = tempname() * ".json"
        PipelinePlotsPlotly.taxa_bar_chart(df, ["s1", "s2"], out_sub;
            rank="Domain", subtitle="run A", rank_order=["Domain"])
        fig_sub = JSON3.read(read(out_sub, String))
        @test occursin("<br><sup>run A</sup>", fig_sub[:layout][:title][:text])
        rm(out_sub)

        # Empty df -> no file written, returns nothing
        out_empty = tempname() * ".json"
        ret = PipelinePlotsPlotly.taxa_bar_chart(DataFrame(), String[], out_empty;
            rank_order=["Domain"])
        @test ret === nothing
        @test !isfile(out_empty)
    end

    @testset "filter_composition_plot" begin
        totals = Dict("Eukaryota" => [100.0, 200.0], "Bacteria" => [50.0, 30.0])
        out = tempname() * ".json"
        PipelinePlotsPlotly.filter_composition_plot(totals, ["s1", "s2"], out)
        @test isfile(out)
        fig = JSON3.read(read(out, String))
        @test haskey(fig, :data) && haskey(fig, :layout)
        @test length(fig[:data]) == 2
        rm(out)

        # Unclassified sorted to end
        totals2 = Dict("Unclassified" => [10.0], "Eukaryota" => [90.0])
        out2 = tempname() * ".json"
        PipelinePlotsPlotly.filter_composition_plot(totals2, ["s1"], out2)
        fig2 = JSON3.read(read(out2, String))
        names = [t[:name] for t in fig2[:data]]
        @test last(names) == "Unclassified"
        rm(out2)

        # Empty -> no file
        out3 = tempname() * ".json"
        PipelinePlotsPlotly.filter_composition_plot(Dict{String,Vector{Float64}}(), String[], out3)
        @test !isfile(out3)
    end

    @testset "alpha_diversity_plot" begin
        adf = DataFrame(sample=["s1","s2"], richness=[10,20],
                        shannon=[1.2,2.1], simpson=[0.8,0.9])
        out = tempname() * ".json"
        PipelinePlotsPlotly.alpha_diversity_plot(adf, out; subtitle="test run")
        @test isfile(out)
        fig = JSON3.read(read(out, String))
        @test haskey(fig, :data) && haskey(fig, :layout)
        @test length(fig[:data]) == 3  # one trace per metric (richness, shannon, simpson)
        @test occursin("<br><sup>test run</sup>", fig[:layout][:title][:text])
        # 3-panel layout keys
        @test haskey(fig[:layout], :yaxis2) && haskey(fig[:layout], :yaxis3)
        rm(out)

        # Empty -> no file
        out2 = tempname() * ".json"
        PipelinePlotsPlotly.alpha_diversity_plot(DataFrame(), out2)
        @test !isfile(out2)
    end

    @testset "pipeline_stats_plot" begin
        sdf = DataFrame(sample=["s1","s2"],
                        input=[1000,800], filtered=[900,700], nochim=[850,650])
        out = tempname() * ".json"
        PipelinePlotsPlotly.pipeline_stats_plot(sdf, out)
        @test isfile(out)
        fig = JSON3.read(read(out, String))
        @test haskey(fig, :data) && haskey(fig, :layout)
        @test fig[:layout][:barmode] == "group"
        @test length(fig[:data]) == 2  # one bar series per sample
        rm(out)

        # Empty -> no file
        out2 = tempname() * ".json"
        PipelinePlotsPlotly.pipeline_stats_plot(DataFrame(), out2)
        @test !isfile(out2)
    end

    @testset "nmds_plot" begin
        coords = [0.1 0.2; -0.3 0.4; 0.5 -0.1]
        labels = ["s1", "s2", "s3"]

        out = tempname() * ".json"
        PipelinePlotsPlotly.nmds_plot(coords, labels, out;
            colour_by=["A","B","A"], stress=0.12)
        @test isfile(out)
        fig = JSON3.read(read(out, String))
        @test haskey(fig, :data) && haskey(fig, :layout)
        @test length(fig[:data]) == 2  # 2 colour groups: A, B
        # Stress annotation present
        @test !isempty(fig[:layout][:annotations])
        @test occursin("stress", fig[:layout][:annotations][1][:text])
        rm(out)

        # No colour_by -> single trace
        out2 = tempname() * ".json"
        PipelinePlotsPlotly.nmds_plot(coords, labels, out2)
        fig2 = JSON3.read(read(out2, String))
        @test length(fig2[:data]) == 1
        rm(out2)

        # Empty -> no file
        out3 = tempname() * ".json"
        PipelinePlotsPlotly.nmds_plot(zeros(0,2), String[], out3)
        @test !isfile(out3)
    end

    @testset "alpha_boxplot" begin
        adf1 = DataFrame(richness=[10,20], shannon=[1.0,2.0], simpson=[0.7,0.8])
        adf2 = DataFrame(richness=[5,15],  shannon=[0.8,1.8], simpson=[0.6,0.75])

        out = tempname() * ".json"
        PipelinePlotsPlotly.alpha_boxplot([adf1, adf2], ["GroupA","GroupB"], out)
        @test isfile(out)
        fig = JSON3.read(read(out, String))
        @test haskey(fig, :data) && haskey(fig, :layout)
        @test length(fig[:data]) == 6  # 3 panels * 2 groups
        # Only first panel traces show legend
        legend_traces = [t for t in fig[:data] if get(t, :showlegend, false)]
        @test length(legend_traces) == 2
        rm(out)

        # Empty -> no file
        out2 = tempname() * ".json"
        PipelinePlotsPlotly.alpha_boxplot(DataFrame[], String[], out2)
        @test !isfile(out2)
    end

    @testset "group_comparison_chart" begin
        df1 = DataFrame(Domain=["Eukaryota","Bacteria"], s1=[100,50])
        df2 = DataFrame(Domain=["Eukaryota","Archaea"],  s1=[200,30])

        out = tempname() * ".json"
        PipelinePlotsPlotly.group_comparison_chart(
            [df1, df2], [["s1"],["s1"]], ["G1","G2"], out;
            rank="Domain", rank_order=["Domain"])
        @test isfile(out)
        fig = JSON3.read(read(out, String))
        @test haskey(fig, :data) && haskey(fig, :layout)
        @test fig[:layout][:barmode] == "stack"
        # 3 unique labels (Eukaryota, Bacteria, Archaea)
        @test length(fig[:data]) == 3
        rm(out)

        # Absolute mode
        out_abs = tempname() * ".json"
        PipelinePlotsPlotly.group_comparison_chart(
            [df1, df2], [["s1"],["s1"]], ["G1","G2"], out_abs;
            rank="Domain", relative=false, rank_order=["Domain"])
        fig_abs = JSON3.read(read(out_abs, String))
        @test occursin("absolute", fig_abs[:layout][:title][:text])
        rm(out_abs)

        # Empty -> no file
        out2 = tempname() * ".json"
        PipelinePlotsPlotly.group_comparison_chart(
            DataFrame[], Vector{String}[], String[], out2; rank_order=["Domain"])
        @test !isfile(out2)
    end

end
