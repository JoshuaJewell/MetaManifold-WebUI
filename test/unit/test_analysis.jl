@testset "Analysis chart builders" begin

    @testset "_palette_hex" begin
        c3 = Analysis._palette_hex(3)
        @test length(c3) == 3
        @test all(s -> startswith(s, "#") && length(s) == 7, c3)

        c10 = Analysis._palette_hex(10)
        @test length(c10) == 10
        @test allunique(c10)

        c7 = Analysis._palette_hex(7)
        @test c7[1] == "#E69F00"
    end

    @testset "alpha_chart" begin
        samples = ["s1", "s2", "s3"]
        r = [10, 20, 30]
        h = [1.0, 2.0, 2.5]
        s = [0.7, 0.8, 0.9]

        fig = Analysis.alpha_chart(samples, r, h, s)
        @test haskey(fig, "data") && haskey(fig, "layout")
        @test length(fig["data"]) == 3  # one trace per metric
        # 3-panel layout keys
        @test haskey(fig["layout"], "yaxis2") && haskey(fig["layout"], "yaxis3")

        fig_empty = Analysis.alpha_chart(String[], Int[], Float64[], Float64[])
        @test haskey(fig_empty, "data")
    end

    @testset "taxa_bar_chart" begin
        labels = ["Eukaryota", "Bacteria"]
        samples = ["s1", "s2"]
        counts = [100.0 30.0; 50.0 20.0]  # 2 taxa x 2 samples

        fig = Analysis.taxa_bar_chart(labels, samples, counts)
        @test haskey(fig, "data") && haskey(fig, "layout")
        @test fig["layout"]["barmode"] == "stack"
        @test length(fig["data"]) == 2

        # Absolute mode
        fig_abs = Analysis.taxa_bar_chart(labels, samples, counts; relative=false)
        @test !haskey(fig_abs["layout"]["yaxis"], "range")

        # Relative mode has range [0, 1]
        fig_rel = Analysis.taxa_bar_chart(labels, samples, counts; relative=true)
        @test fig_rel["layout"]["yaxis"]["range"] == [0, 1]
    end

    @testset "taxa_bar_chart top_n collapsing" begin
        labels = ["T$i" for i in 1:20]
        samples = ["s1"]
        counts = Float64[i for i in 1:20] |> c -> reshape(c, 20, 1)

        fig = Analysis.taxa_bar_chart(labels, samples, counts; top_n=5)
        trace_names = [t["name"] for t in fig["data"]]
        @test "Other" in trace_names
        @test length(trace_names) == 6  # top 5 + Other
    end

    @testset "nmds_chart" begin
        coords = [0.1 0.2; -0.3 0.4; 0.5 -0.1]
        labels = ["s1", "s2", "s3"]

        fig = Analysis.nmds_chart(coords, labels;
            colour_by=["A", "B", "A"], stress=0.12)
        @test haskey(fig, "data") && haskey(fig, "layout")
        @test length(fig["data"]) == 2  # 2 colour groups: A, B
        # Stress annotation present
        @test !isempty(fig["layout"]["annotations"])
        @test occursin("stress", fig["layout"]["annotations"][1]["text"])

        fig2 = Analysis.nmds_chart(coords, labels)
        @test length(fig2["data"]) == 1

        fig3 = Analysis.nmds_chart(zeros(0, 2), String[])
        @test haskey(fig3, "data")
    end

    @testset "alpha_boxplot" begin
        groups = [
            ("GroupA", [10, 20], [1.0, 2.0], [0.7, 0.8]),
            ("GroupB", [5, 15],  [0.8, 1.8], [0.6, 0.75]),
        ]

        fig = Analysis.alpha_boxplot(groups)
        @test haskey(fig, "data") && haskey(fig, "layout")
        @test length(fig["data"]) == 6  # 3 panels * 2 groups
        # Only first panel traces show legend
        legend_traces = [t for t in fig["data"] if get(t, "showlegend", false)]
        @test length(legend_traces) == 2

        # The overall-significance label must anchor to its panel's axis domain,
        # not to paper: a paper-referenced annotation cannot survive the
        # frontend's per-metric axis renumbering and spawns phantom axes.
        annotated = Analysis.alpha_boxplot(groups; annotate_significance=true)
        sig_anns = get(annotated["layout"], "annotations", Any[])
        @test !isempty(sig_anns)
        @test all(a -> get(a, "xref", "") != "paper" && get(a, "yref", "") != "paper",
                  sig_anns)
        @test all(a -> endswith(String(get(a, "yref", "")), " domain"), sig_anns)
    end

    @testset "alpha_boxplot paired test needs shared samples" begin
        label(fig) = first(get(fig["layout"], "annotations", Any[]))["text"]
        unshared = [
            ("GroupA", ["s1", "s2"], [10, 20], [1.0, 2.0], [0.7, 0.8]),
            ("GroupB", ["s3"],       [15],     [1.5],      [0.75]),
        ]
        fig = Analysis.alpha_boxplot(unshared; annotate_significance=true, paired_samples=true)
        @test startswith(label(fig), "KW")
        shared = [
            ("GroupA", ["s1", "s2"], [10, 20], [1.0, 2.0], [0.7, 0.8]),
            ("GroupB", ["s1", "s2"], [15, 25], [1.5, 2.5], [0.75, 0.85]),
        ]
        fig = Analysis.alpha_boxplot(shared; annotate_significance=true, paired_samples=true)
        @test startswith(label(fig), "Paired Wilcoxon")
    end

    @testset "alpha_boxplot paired lines" begin
        # "solo" appears in one group only, so nothing can be paired with it.
        groups = [
            ("GroupA", ["p1", "p2", "solo"], [10, 12, 9], [1.0, 1.2, 0.9], [0.5, 0.55, 0.45]),
            ("GroupB", ["p1", "p2"],         [20, 18],    [2.0, 1.8],      [0.7, 0.65]),
        ]

        plain = Analysis.alpha_boxplot(groups; significance_test="none")
        @test all(t -> t["type"] == "box", plain["data"])

        fig = Analysis.alpha_boxplot(groups; paired_lines=true, significance_test="none")
        lines = [t for t in fig["data"] if t["type"] == "scatter"]
        boxes = [t for t in fig["data"] if t["type"] == "box"]
        @test length(boxes) == 6                 # 3 panels * 2 groups
        @test length(lines) == 6                 # 3 panels * 2 paired samples
        @test all(t -> t["text"][1] != "solo", lines)
        @test all(t -> !t["showlegend"], lines)

        # Each line spans the two groups on one panel's own axes.
        richness_lines = [t for t in lines if t["xaxis"] == "x"]
        @test length(richness_lines) == 2
        @test all(t -> t["x"] == ["GroupA", "GroupB"], richness_lines)
        p1 = first(filter(t -> t["text"][1] == "p1", richness_lines))
        @test p1["y"] == [10.0, 20.0]
        @test p1["yaxis"] == "y"
        @test Set(t["xaxis"] for t in lines) == Set(["x", "x2", "x3"])

        # Jitter would leave the lines hanging off their points, so it is off.
        @test all(t -> t["jitter"] == 0.0, boxes)
        @test all(t -> t["jitter"] == 0.35, [t for t in plain["data"] if t["type"] == "box"])

        # A single group has no pairs to draw.
        one = Analysis.alpha_boxplot(groups[1:1]; paired_lines=true, significance_test="none")
        @test all(t -> t["type"] == "box", one["data"])
    end

    @testset "bar_chart keep_empty" begin
        labels = ["A", "B"]
        samples = ["s1", "s2", "s3"]
        counts = Float64[10 0 5; 4 0 6]

        dropped = Analysis.bar_chart(labels, samples, copy(counts))
        @test dropped["data"][1]["x"] == ["s1", "s3"]

        kept = Analysis.bar_chart(labels, samples, copy(counts); keep_empty=true)
        @test kept["data"][1]["x"] == ["s1", "s2", "s3"]
        # The retained sample is a zero-height bar, not a hole in the data.
        @test all(t -> t["y"][2] == 0.0, kept["data"])
        # Plotly reserves a slot for an all-zero bar only on a categorical axis.
        @test kept["layout"]["xaxis"]["type"] == "category"
        @test dropped["layout"]["xaxis"]["type"] == "category"
    end

    @testset "faceted_bar_chart" begin
        panels = [
            (; row="RunA", col="Ctrl",  segment_labels=["A", "B"],
               sample_names=["a1", "a2"], counts=Float64[6 4; 4 6]),
            (; row="RunA", col="Treat", segment_labels=["A", "C"],
               sample_names=["a3"], counts=reshape(Float64[7, 3], 2, 1)),
            (; row="RunB", col="Ctrl",  segment_labels=["B"],
               sample_names=["b1"], counts=reshape(Float64[9], 1, 1)),
        ]
        fig = Analysis.faceted_bar_chart(panels, ["RunA", "RunB"], ["Ctrl", "Treat"];
                                         row_title="run", col_title="group")

        @test fig["layout"]["grid"] == Dict("rows" => 2, "columns" => 2,
                                            "pattern" => "independent")
        # Every cell gets a trace per retained segment, so the legend is complete
        # however sparse an individual panel is: 4 cells * 3 segments.
        @test length(fig["data"]) == 12
        @test length([t for t in fig["data"] if t["showlegend"]]) == 3

        # The grid fills row-major, so RunB/Treat is panel 4 -- present but empty
        # rather than absent, keeping the hole visible.
        empty_cell = [t for t in fig["data"] if t["xaxis"] == "x4"]
        @test length(empty_cell) == 3
        @test all(t -> isempty(t["x"]), empty_cell)

        # A segment absent from a panel is zero there, not missing.
        ctrl_a = first(filter(t -> t["xaxis"] == "x" && t["name"] == "C", fig["data"]))
        @test ctrl_a["y"] == [0.0, 0.0]

        # One colour per label across every panel, so a taxon reads the same
        # wherever it appears in the grid.
        by_label = Dict{String, Set{String}}()
        for t in fig["data"]
            push!(get!(by_label, t["name"], Set{String}()), t["marker"]["color"])
        end
        @test length(by_label) == 3
        @test all(cs -> length(cs) == 1, values(by_label))

        # Row and column headers anchor to panel axis domains, never to paper,
        # so they travel with their panel.
        anns = fig["layout"]["annotations"]
        @test any(a -> occursin("group: Ctrl", a["text"]), anns)
        @test any(a -> occursin("run: RunB", a["text"]), anns)
        @test all(a -> endswith(a["xref"], " domain"), anns)

        # The top-N cut is global, so "Other" means the same thing in every panel.
        wide = [(; row="R", col="C", segment_labels=["t$i" for i in 1:5],
                   sample_names=["s1"], counts=reshape(Float64[50, 40, 30, 20, 10], 5, 1))]
        capped = Analysis.faceted_bar_chart(wide, ["R"], ["C"]; top_n=2, relative=false)
        @test [t["name"] for t in capped["data"]] == ["t1", "t2", "Other"]
        @test last(capped["data"])["y"] == [60.0]

        # An empty grid is a figure with nothing in it, not an error.
        empty_fig = Analysis.faceted_bar_chart(typeof(panels[1])[], String[], String[])
        @test isempty(empty_fig["data"])
    end

    @testset "natural sample order" begin
        names = ["Caecum_10c_m", "Caecum_2c_m", "Caecum_1c_m", "Caecum_12c_m",
                 "Large_Intestine_3l_m"]
        @test Analysis.natural_sort(names) ==
            ["Caecum_1c_m", "Caecum_2c_m", "Caecum_10c_m", "Caecum_12c_m",
             "Large_Intestine_3l_m"]
        # Leading zeros tie on value, then fall back to a total string order.
        @test Analysis.natural_sort(["s010", "s9", "s10"]) == ["s9", "s010", "s10"]

        fig = Analysis.bar_chart(["A"], ["x_10", "x_2", "x_1"],
                                 Float64[1 2 3]; relative=false)
        @test fig["data"][1]["x"] == ["x_1", "x_2", "x_10"]
        @test fig["data"][1]["y"] == [3.0, 2.0, 1.0]
    end

    @testset "short_sample_labels" begin
        # Within one compartment of one run only the individual number varies.
        @test Analysis.short_sample_labels(
            ["Large_Intestine_20l_m", "Large_Intestine_2l_m", "Large_Intestine_33l_m"]) ==
            ["20", "2", "33"]
        # Across a whole run only the run suffix is shared.
        @test Analysis.short_sample_labels(["Caecum_1c_m", "Large_Intestine_2l_m"]) ==
            ["Caecum_1c", "Large_Intestine_2l"]
        # Trimming works on whole tokens: "10" and "20" share no token.
        @test Analysis.short_sample_labels(["10c", "20c"]) == ["10", "20"]
        # Nothing to compare against, or nothing left: names are kept.
        @test Analysis.short_sample_labels(["Caecum_1c_m"]) == ["Caecum_1c_m"]
        @test Analysis.short_sample_labels(["a_1", "a_1_x"]) == ["a_1", "a_1_x"]

        fig = Analysis.bar_chart(["A"], ["Caecum_2c_m", "Caecum_10c_m"],
                                 Float64[1 1]; short_labels=true)
        @test fig["data"][1]["x"] == ["2", "10"]
        @test fig["data"][1]["hovertext"] == ["Caecum_2c_m", "Caecum_10c_m"]

        panels = [(; row="Multiplex", col="Caecum", segment_labels=["A"],
                     sample_names=["Caecum_10c_m", "Caecum_1c_m"],
                     counts=Float64[1 2])]
        grid = Analysis.faceted_bar_chart(panels, ["Multiplex"], ["Caecum"];
                                          relative=false, short_labels=true)
        @test grid["data"][1]["x"] == ["1", "10"]
        @test grid["data"][1]["y"] == [2.0, 1.0]
    end

    @testset "pinned segment order" begin
        labels = ["Unassigned", "Host", "Protozoa", "Fungi"]
        counts = Float64[40 40; 30 30; 25 25; 1 1]
        names_of(fig) = [t["name"] for t in fig["data"]]

        @test names_of(Analysis.bar_chart(labels, ["s1", "s2"], copy(counts))) ==
            ["Protozoa", "Host", "Fungi", "Unassigned"]
        # Unassigned stays last even behind a collapsed Other.
        @test names_of(Analysis.bar_chart(labels, ["s1", "s2"], copy(counts); top_n=3)) ==
            ["Protozoa", "Host", "Other", "Unassigned"]

        panels = [(; row="R", col="C", segment_labels=labels,
                     sample_names=["s1", "s2"], counts=copy(counts))]
        grid = Analysis.faceted_bar_chart(panels, ["R"], ["C"]; top_n=3, relative=false)
        @test names_of(grid) == ["Protozoa", "Host", "Other", "Unassigned"]
        # The collapsed Fungi reads land in Other, not in Unassigned.
        other = first(filter(t -> t["name"] == "Other", grid["data"]))
        @test other["y"] == [1.0, 1.0]
    end

    @testset "bar_chart modes and colour_for" begin
        fig = Analysis.bar_chart(["A", "B"], ["s1", "s2"],
            Float64[1 2; 3 4]; mode="group", relative=false,
            colour_for = l -> l == "A" ? "#111111" : "#222222")
        @test fig["layout"]["barmode"] == "group"
        # Traces are ordered by total descending: B(7) first, A(3) second.
        # Find the trace named "A" and check its colour.
        trace_a = first(filter(t -> t["name"] == "A", fig["data"]))
        @test trace_a["marker"]["color"] == "#111111"
        # Grouped absolute mode must not fix the y-axis range to [0, 1].
        @test !haskey(fig["layout"]["yaxis"], "range")

        # Stacked relative mode keeps the [0, 1] range.
        fig2 = Analysis.bar_chart(["A", "B"], ["s1", "s2"],
            Float64[1 2; 3 4]; mode="stacked", relative=true)
        @test fig2["layout"]["barmode"] == "stack"
        @test fig2["layout"]["yaxis"]["range"] == [0, 1]
    end

    @testset "pool_columns" begin
        counts = [10.0 20.0 30.0 40.0;
                   5.0 10.0 15.0 20.0]
        names = ["A_s1", "A_s2", "B_s1", "B_s2"]

        # Pool by prefix
        pn, pc = Analysis.pool_columns(names, counts, ["A", "B"])
        @test pn == ["A", "B"]
        @test pc[:, 1] == [30.0, 15.0]   # A_s1 + A_s2
        @test pc[:, 2] == [70.0, 35.0]   # B_s1 + B_s2

        # Pool with unmatched -> Other
        pn2, pc2 = Analysis.pool_columns(names, counts, ["A"])
        @test pn2 == ["A", "Other"]
        @test pc2[:, 1] == [30.0, 15.0]
        @test pc2[:, 2] == [70.0, 35.0]

        # Empty groups -> pool all
        pn3, pc3 = Analysis.pool_columns(names, counts, String[])
        @test pn3 == ["Total"]
        @test pc3[:, 1] == [100.0, 50.0]

        # Custom fallback label
        pn4, _ = Analysis.pool_columns(names, counts, String[]; fallback_label="MyRun")
        @test pn4 == ["MyRun"]
    end

end
