# Unit tests for the publication-table builders and their CSV/XLSX renderings.
using MetaManifold
SV = MetaManifold.Server
using XLSX, CSV, DataFrames

@testset "Publication tables" begin
    run_a = (; label="RunA", n=3,
               tally=Dict("Protozoa" => (60, 4), "Host" => (40, 1), "Unassigned" => (0, 0)))
    run_b = (; label="RunB", n=2,
               tally=Dict("Protozoa" => (1, 1), "Bacteria" => (999, 7)))
    all_values = [:asvs, :reads, :pct]
    # One column per run, headed by the run.
    by_run(units...; values=all_values, order=["Protozoa", "Host", "Bacteria"], title="T", notes=String[]) =
        SV._pub_table(; title, first_col="Category", groups=[(; label=u.label, columns=[u]) for u in units],
                      nested=false, order, values, diff=false, notes)

    @testset "category order" begin
        order = SV._pub_category_order(["Protozoa", "Host", "Fungi"], [run_a.tally, run_b.tally])
        # Set order first (including categories with no reads), extras after,
        # Unassigned only when it has reads.
        @test order == ["Protozoa", "Host", "Fungi", "Bacteria"]
        order = SV._pub_category_order(["Protozoa"], [Dict("Unassigned" => (5, 1))])
        @test order == ["Protozoa", "Unassigned"]
    end

    @testset "taxon order" begin
        t1 = Dict("A" => (10, 1), "B" => (90, 2), "Unclassified" => (50, 1), "Z" => (0, 0))
        t2 = Dict("A" => (100, 1))
        # Mean share: A = (10/150 + 1)/2, B = (90/150)/2; zero-read labels dropped.
        @test SV._pub_taxon_order([t1, t2]) == ["A", "B", "Unclassified"]
    end

    @testset "table by run" begin
        tbl = by_run(run_a, run_b)
        @test [h["label"] for h in only(tbl["header_rows"])] == ["", "RunA (n = 3)", "RunB (n = 2)"]
        @test [h["span"] for h in only(tbl["header_rows"])] == [1, 3, 3]
        @test length(tbl["columns"]) == 7
        @test tbl["rows"][1] == ["Protozoa", 4, 60, 60.0, 1, 1, 0.1]
        @test tbl["rows"][3] == ["Bacteria", 0, 0, 0.0, 7, 999, 99.9]
        @test only(tbl["footer"]) == ["Total", 5, 100, 100.0, 8, 1000, 100.0]
        @test all(v -> v isa Integer, only(tbl["footer"])[[2, 3, 5, 6]])
    end

    @testset "chosen values" begin
        tbl = by_run(run_a, run_b; values=[:pct])
        @test [c["label"] for c in tbl["columns"]] == ["Category", "%", "%"]
        @test tbl["rows"][1] == ["Protozoa", 60.0, 0.1]
        @test only(tbl["footer"]) == ["Total", 100.0, 100.0]
    end

    @testset "nested columns and difference" begin
        run_x = (; label="Run_X", columns=[
            (; label="Caecum", n=2, tally=Dict("Host" => (50, 1), "Protozoa" => (50, 2))),
            (; label="Large_Intestine", n=2, tally=Dict("Host" => (75, 1), "Protozoa" => (25, 1))),
        ])
        # A single sample is headed by its name alone.
        run_y = (; label="Run_Y", columns=[(; label="Y_s1", n=1, tally=Dict("Host" => (0, 0)))])
        tbl = SV._pub_table(; title="T", first_col="Category", groups=[run_x, run_y], nested=true,
                            order=["Host", "Protozoa"], values=[:reads, :pct], diff=true, notes=String[])
        top, sub = tbl["header_rows"]
        @test [(h["label"], h["span"]) for h in top] == [("", 1), ("Run X", 5), ("Run Y", 2)]
        @test [h["label"] for h in sub] == ["", "Caecum (n = 2)", "Large Intestine (n = 2)", "", "Y s1"]
        @test sum(h["span"] for h in top) == length(tbl["columns"])
        @test sum(h["span"] for h in sub) == length(tbl["columns"])
        # Difference column only for the two-column group; empty scopes give no share.
        @test tbl["rows"][1] == ["Host", 50, 50.0, 75, 75.0, 25.0, 0, nothing]
        @test only(tbl["footer"]) == ["Total", 100, 100.0, 100, 100.0, nothing, 0, nothing]
        # No percentages, no difference.
        plain = SV._pub_table(; title="T", first_col="Category", groups=[run_x], nested=true,
                              order=["Host"], values=[:reads], diff=true, notes=String[])
        @test plain["rows"][1] == ["Host", 50, 75]
    end

    @testset "flat headers and CSV" begin
        tbl = SV._pub_table(; title="T", first_col="Genus", groups=[(; label="RunA", columns=[run_a])],
                            nested=false, order=["Protozoa"], values=all_values, diff=false, notes=String[])
        @test SV._pub_flat_headers(tbl) == ["Genus", "RunA (n = 3) ASVs", "RunA (n = 3) Reads", "RunA (n = 3) %"]
        df = CSV.read(IOBuffer(SV._pub_csv(tbl)), DataFrame)
        @test df[1, "RunA (n = 3) Reads"] === 60
        @test df[2, "Genus"] == "Total"
    end

    @testset "heatmap and hidden zeros" begin
        tbl = by_run(run_a, run_b)
        plain = SV._pub_present(tbl)
        @test all(isnothing, reduce(vcat, plain["fills"]))
        @test plain["rows"] == tbl["rows"]

        # Per column: RunA reads peak at Protozoa (60), so it gets the full colour.
        col = SV._pub_present(tbl; heatmap="column")
        @test col["fills"][1][3] == "#6fa8dc"
        @test col["fills"][2][3] == SV._heat_blend(SV._HEAT_HIGH, 40 / 60)
        @test isnothing(col["fills"][1][1])                 # label column
        @test isnothing(col["fills"][3][3])                 # zero reads
        # Whole table: RunA reads share a scale with RunB reads (max 999).
        whole = SV._pub_present(tbl; heatmap="table")
        @test whole["fills"][1][3] == SV._heat_blend(SV._HEAT_HIGH, 60 / 999)
        @test whole["fills"][3][6] == "#6fa8dc"

        hidden = SV._pub_present(tbl; hide_zeros=true)
        @test hidden["rows"][3][2:4] == ["", "", ""]         # Bacteria absent from RunA
        @test hidden["rows"][3][1] == "Bacteria"
        @test hidden["footer"] == tbl["footer"]
        @test tbl["rows"][3][2] == 0                        # original untouched

        diff = Dict{String,Any}("columns" => [SV._pub_col("x", "label"), SV._pub_col("d", "pp")],
                                "rows" => [Any["a", 2.0], Any["b", -1.0]], "footer" => [])
        d = SV._pub_present(diff; heatmap="column")
        @test d["fills"][1][2] == "#6fa8dc"
        @test d["fills"][2][2] == SV._heat_blend(SV._HEAT_LOW, 0.5)

        # Fills reach the workbook.
        path = tempname() * ".xlsx"
        write(path, SV._pub_xlsx([SV._pub_present(tbl; heatmap="column", hide_zeros=true)]))
        xf = XLSX.readxlsx(path)
        @test ismissing(xf["Table 1"]["B6"])               # hidden zero written blank
        @test xf["Table 1"]["C4"] == 60
        rm(path)
    end

    @testset "display values" begin
        @test SV._pub_display_value("", "int") == ""
        @test SV._pub_display_value(nothing, "pct") == "–"
        @test SV._pub_display_value(0.004, "pct") == "<0.01"
        @test SV._pub_display_value(0.0, "pct") == 0.0
        @test SV._pub_display_value(12.5, "pct") == 12.5
        @test SV._pub_display_value(-0.004, "pp") === 0.0
        @test SV._pub_display_value(2.468, "pp") == 2.47
    end

    @testset "xlsx workbook" begin
        tbl = by_run(run_a, run_b; order=["Protozoa", "Host"], title="Composition.", notes=["A note."])
        path = tempname() * ".xlsx"
        write(path, SV._pub_xlsx([tbl]))
        xf = XLSX.readxlsx(path)
        @test XLSX.sheetnames(xf) == ["Table 1"]
        sh = xf["Table 1"]
        @test sh["A1"] == "Composition."
        @test sh["B2"] == "RunA (n = 3)"
        @test sh["A3"] == "Category"
        @test sh["A4"] == "Protozoa" && sh["C4"] == 60
        @test sh["A6"] == "Total" && sh["F6"] == 1000
        @test sh["A8"] == "A note."
        rm(path)
    end
end
