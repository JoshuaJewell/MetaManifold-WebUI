# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Differential abundance: Benjamini-Hochberg, TSS/RLE size factors, input
# validation, and the negative-binomial fit checked against a direct
# MASS::glm.nb call. Tests needing R or MASS skip loudly when either is absent.

using MetaManifold
using MetaManifold.Differential
using RCall
using Test
using DuckDB, DBInterface, JSON3

"""
    _r_has(pkg) -> Bool

Whether the embedded R session starts and can load `pkg`. Probed directly,
not through `Server.r_available`, which also requires vegan.
"""
function _r_has(pkg::AbstractString)
    try
        rcopy(RCall.reval("requireNamespace('$pkg', quietly = TRUE)"))
    catch
        false
    end
end

const HAVE_R = _r_has("stats")
const HAVE_MASS = HAVE_R && _r_has("MASS")
HAVE_MASS || @info "R or MASS unavailable - the negative-binomial fit tests are SKIPPED"

# Rows per run: Up is split over two ASVs so the route must pool them by genus.
const _DA_ROWS = Dict(
    "runA" => [("s1", "Up", [3, 5, 4, 4]), ("s2", "Up", [2, 3, 2, 3]),
               ("s3", "Flat", [30, 25, 35, 28]), ("s4", "Down", [50, 60, 45, 55])],
    "runB" => [("s1", "Up", [20, 30, 20, 30]), ("s2", "Up", [20, 25, 15, 30]),
               ("s3", "Flat", [27, 33, 29, 31]), ("s4", "Down", [6, 4, 9, 5])])

"""
    _da_make_run(root, study, run)

Write `<root>/projects/<study>/<run>/merged/results.duckdb` with four samples
named after the run, plus the data-side run directory the route guards need.
"""
function _da_make_run(root::String, study::String, run::String)
    data_run = joinpath(root, "data", study, run)
    mkpath(data_run)
    samples = ["$(run)_$i" for i in 1:4]
    for s in samples, d in ("R1", "R2")
        touch(joinpath(data_run, "$(s)_$d.fastq.gz"))
    end
    write(joinpath(data_run, "pipeline.yml"), "analysis:\n  exclude_categories: []\n")
    merge_dir = joinpath(root, "projects", study, run, "merged")
    mkpath(merge_dir)
    db = DuckDB.DB(joinpath(merge_dir, "results.duckdb"))
    con = DBInterface.connect(db)
    try
        cols = join(["\"$s\" BIGINT" for s in samples], ", ")
        DBInterface.execute(con, """CREATE TABLE merged ("SeqName" VARCHAR, "sequence" VARCHAR,
                                    "Domain" VARCHAR, "Genus" VARCHAR, $cols)""")
        for (seq, genus, n) in _DA_ROWS[run]
            DBInterface.execute(con, "INSERT INTO merged VALUES ('$seq', 'ACGT', 'Eukaryota', '$genus', $(join(n, ", ")))")
        end
    finally
        DBInterface.close!(con); close(db)
    end
end

"""
    _da_fixture(f)

Run `f()` with the server root pointed at a temporary study `studyD` holding
runs `runA` and `runB`, restoring the previous root afterwards.
"""
function _da_fixture(f::Function)
    SV = MetaManifold.Server
    root = mktempdir()
    old_root = SV.ServerState._root[]
    try
        SV.ServerState.set_root!(root)
        mkpath(joinpath(root, "config"))
        _da_make_run(root, "studyD", "runA")
        _da_make_run(root, "studyD", "runB")
        f()
    finally
        SV.ServerState._root[] = old_root
        rm(root; recursive=true, force=true)
    end
end

"""
    _da_post(study, body) -> HTTP.Response

POST `body` to the study's differential route through Oxygen's in-process router.
"""
function _da_post(study::String, body)
    SV = MetaManifold.Server
    req = SV.HTTP.Request("POST", "/api/v1/studies/$study/analysis/differential",
                          ["Content-Type" => "application/json"], JSON3.write(body))
    SV.Oxygen.internalrequest(req)
end

@testset "Differential abundance" begin

    @testset "BH adjustment: known answers" begin
        @test bh_adjust([0.01, 0.02, 0.03, 0.04]) ≈ fill(0.04, 4)
        @test bh_adjust([0.001, 0.5, 0.5]) ≈ [0.003, 0.5, 0.5]
        # Order of the input is preserved, monotonicity enforced by the step-up.
        @test bh_adjust([0.04, 0.01]) ≈ [0.04, 0.02]
        @test bh_adjust([1.0]) == [1.0]
        @test isempty(bh_adjust(Float64[]))
        @test_throws ArgumentError bh_adjust([0.1, NaN])
        @test_throws ArgumentError bh_adjust([1.4])
        @test_throws ArgumentError bh_adjust([-0.2])
    end

    @testset "BH adjustment: parity with R p.adjust" begin
        if !HAVE_R
            @info "R unavailable - p.adjust parity SKIPPED"
            @test_skip false
        else
            p = [0.001, 0.5, 0.5, 0.02, 0.7, 0.04, 0.0009, 0.31]
            ref = rcopy(R"p.adjust($p, method = 'BH')")
            @test maximum(abs.(bh_adjust(p) .- ref)) < 1e-12
        end
    end

    # samples x taxa
    fixture = [10 20 30; 20 20 60; 40 10 50]
    samples = ["s1", "s2", "s3"]

    @testset "TSS size factors" begin
        f = tss_factors(fixture, samples)
        @test f ≈ [0.711378660898012, 1.185631101496687, 1.185631101496687] rtol=1e-12
        @test exp(sum(log.(f)) / 3) ≈ 1.0
        @test tss_factors(7.5 .* fixture, samples) ≈ f
        e = try tss_factors([1 2; 0 0], ["a", "empty"]); nothing catch err; err end
        @test e isa ScalingRefusal
        @test occursin("empty", sprint(showerror, e))
    end

    @testset "RLE size factors" begin
        f = rle_factors(fixture, samples)
        @test f ≈ [0.683132584572384, 1.285704749170459, 1.13855430762064] rtol=1e-12
        raw = [0.6694329500821695, 1.259921049894873, 1.1157215834702825]
        @test f ./ f[1] ≈ raw ./ raw[1]
        @test exp(sum(log.(f)) / 3) ≈ 1.0
        @test rle_factors(7.5 .* fixture, samples) ≈ f
        doubled = copy(fixture); doubled[2, :] .*= 2
        g = rle_factors(doubled, samples)
        @test (g[2] / g[1]) ≈ 2 * (f[2] / f[1])
        e = try rle_factors([1 2; 0 0], ["a", "b"]); nothing catch err; err end
        @test e isa ScalingRefusal
        @test occursin("geometric mean", sprint(showerror, e))
        @test size_factors(fixture, samples, "rle") == f
        @test_throws ArgumentError size_factors(fixture, samples, "css")
    end

    @testset "Configuration" begin
        @test DifferentialConfig().offset == "tss"
        @test DifferentialConfig().min_prevalence == 0.0
        @test DifferentialConfig("rle", 0.5).offset == "rle"
        @test_throws ArgumentError DifferentialConfig("css")
        @test_throws ArgumentError DifferentialConfig("tss", 1.5)
        @test_throws ArgumentError DifferentialConfig("tss", NaN)
    end

    @testset "Input validation" begin
        taxa = ["a", "b", "c"]
        groups = ["A", "A", "B"]
        kw = (reference="A", contrast="B")
        @test_throws ArgumentError validate_counts([1.5 2 3; 1 2 3; 1 2 3], samples, taxa)
        @test_throws ArgumentError validate_counts([-1 2 3; 1 2 3; 1 2 3], samples, taxa)
        msg = try validate_counts([1 2 3; 1 2.5 3; 1 2 3], samples, taxa); "" catch e; e.msg end
        @test occursin("'b'", msg) && occursin("'s2'", msg)
        @test_throws ArgumentError differential_abundance(fixture, samples, taxa, groups;
                                                         reference="A", contrast="A")
        @test_throws ArgumentError differential_abundance(fixture, samples, taxa, ["A", "B", "C"]; kw...)
        @test_throws ArgumentError differential_abundance(fixture, samples, taxa, ["A", "A", "A"]; kw...)
        @test_throws ArgumentError differential_abundance(fixture[1:2, :], samples[1:2], taxa, ["A", "B"]; kw...)
        @test_throws ArgumentError differential_abundance(fixture, samples, String[], groups; kw...)
        @test_throws ArgumentError differential_abundance(Float64.(fixture) .+ 0.5, samples, taxa, groups; kw...)
    end

    # taxa x samples as written, transposed to samples x taxa.
    effect_groups = vcat(fill("A", 12), fill("B", 12))
    effect_samples = ["s$i" for i in 1:24]
    effect_taxa = ["t_up", "t_flat", "t_down"]
    effect = permutedims(Int.(
        [7 17 1 12 10 3 15 14 6 11 15 9 52 18 40 57 25 58 8 43 65 31 47 36;
         32 18 44 6 43 27 39 14 30 49 23 35 14 35 27 49 6 32 44 23 43 30 18 39;
         81 50 22 65 39 71 59 10 45 73 54 31 4 8 8 2 6 1 5 9 1 7 3 6]))

    @testset "Negative-binomial fit recovers an injected effect" begin
        if !HAVE_MASS
            @test_skip false
        else
            res = differential_abundance(effect, effect_samples, effect_taxa, effect_groups;
                                         reference="A", contrast="B")
            byt = Dict(r["taxon"] => r for r in res["rows"])
            @test res["status"] == "ok"
            @test res["diagnostics"]["n_tested"] == 3
            @test all(r -> r["status"] == "ok", res["rows"])
            @test isapprox(byt["t_up"]["estimate"], log(4); atol=0.5)
            @test isapprox(byt["t_up"]["log2_fold_change"], 2; atol=0.8)
            @test byt["t_up"]["pvalue"] < 0.01
            @test byt["t_down"]["estimate"] < 0
            @test abs(byt["t_flat"]["estimate"]) < 0.6
            @test [byt[t]["padj"] for t in effect_taxa] ≈
                  bh_adjust([byt[t]["pvalue"] for t in effect_taxa])
            # Swapping the reference flips the sign.
            swapped = differential_abundance(effect, effect_samples, effect_taxa, effect_groups;
                                             reference="B", contrast="A")
            sb = Dict(r["taxon"] => r for r in swapped["rows"])
            @test sb["t_up"]["estimate"] ≈ -byt["t_up"]["estimate"] rtol=1e-6

            @testset "matches a direct MASS::glm.nb call" begin
                off = log.(tss_factors(effect, effect_samples))
                for (j, t) in enumerate(effect_taxa)
                    y = effect[:, j]
                    g = effect_groups
                    direct = rcopy(R"""
                        local({
                          fit <- MASS::glm.nb($y ~ factor($g, levels = c("A", "B")) + offset($off),
                                              control = glm.control(maxit = 100))
                          co <- summary(fit)$coefficients[2, ]
                          c(co[[1]], co[[2]], co[[4]], fit$theta)
                        })
                    """)
                    r = byt[t]
                    @test r["estimate"] ≈ direct[1] rtol=1e-6
                    @test r["standard_error"] ≈ direct[2] rtol=1e-6
                    @test r["pvalue"] ≈ direct[3] rtol=1e-6
                    @test r["dispersion_theta"] ≈ direct[4] rtol=1e-6
                end
            end
        end
    end

    @testset "Degenerate and filtered taxa are reported, not dropped" begin
        if !HAVE_MASS
            @test_skip false
        else
            counts = hcat(effect, fill(5, 24), vcat(fill(0, 20), [3, 0, 1, 0]))
            taxa = vcat(effect_taxa, ["t_const", "t_rare"])
            res = differential_abundance(counts, effect_samples, taxa, effect_groups;
                                         reference="A", contrast="B",
                                         config=DifferentialConfig("tss", 0.2))
            byt = Dict(r["taxon"] => r for r in res["rows"])
            @test length(res["rows"]) == 5
            @test res["status"] == "partial"
            c = byt["t_const"]
            @test c["status"] == "failed"
            @test occursin("constant", c["note"])
            @test all(isnothing, (c["estimate"], c["standard_error"], c["pvalue"], c["padj"]))
            r = byt["t_rare"]
            @test r["status"] == "filtered"
            @test occursin("min_prevalence", r["note"])
            @test isnothing(r["padj"])
            @test r["prevalence"] ≈ 2 / 24
            d = res["diagnostics"]
            @test (d["n_taxa"], d["n_tested"], d["n_failed"], d["n_filtered"]) == (5, 3, 1, 1)
            # The BH family is the three fitted taxa only.
            @test [byt[t]["padj"] for t in effect_taxa] ≈
                  bh_adjust([byt[t]["pvalue"] for t in effect_taxa])
            # Unfitted rows sort last.
            @test [x["taxon"] for x in res["rows"][4:5]] ⊆ ["t_const", "t_rare"]
        end
    end

    @testset "No fittable taxon is an explicit error" begin
        if !HAVE_MASS
            @test_skip false
        else
            e = try
                differential_abundance(fill(5, 4, 2), ["a", "b", "c", "d"], ["x", "y"],
                                       ["A", "A", "B", "B"]; reference="A", contrast="B")
                nothing
            catch err
                err
            end
            @test e isa ErrorException
            @test occursin("no taxon could be fitted", e.msg)
        end
    end

    @testset "Volcano chart" begin
        res = Dict{String,Any}(
            "groups" => Dict("reference" => "A", "contrast" => "B"),
            "rows" => [
                Dict{String,Any}("taxon" => "x", "log2_fold_change" => 2.0, "pvalue" => 1e-4, "padj" => 1e-3),
                Dict{String,Any}("taxon" => "y", "log2_fold_change" => -0.1, "pvalue" => 0.6, "padj" => 0.6),
                Dict{String,Any}("taxon" => "z", "log2_fold_change" => nothing, "pvalue" => nothing, "padj" => nothing),
            ])
        fig = volcano_chart(res)
        @test length(fig["data"]) == 2
        xs = reduce(vcat, [collect(t["x"]) for t in fig["data"]])
        @test sort(xs) == [-0.1, 2.0]
        @test occursin("B vs A", fig["layout"]["title"]["text"])
    end

    @testset "Config value checks" begin
        SV = MetaManifold.Server
        @test isnothing(SV._config_value_error("analysis.differential.offset", "rle"))
        @test !isnothing(SV._config_value_error("analysis.differential.offset", "css"))
        @test isnothing(SV._config_value_error("analysis.differential.min_prevalence", 0.1))
    end

    @testset "Route: two runs through the in-process router" begin
        _da_fixture() do
            r = _da_post("studyD", Dict("runs" => [Dict("run" => "runA"), Dict("run" => "runB")],
                                        "table" => "merged", "rank" => "Genus"))
            body = JSON3.read(String(r.body))
            if !HAVE_MASS
                @test r.status == 503
                @test body.error == "r_unavailable"
            else
                @test r.status == 200
                @test body.rank == "Genus"
                @test body.groups.reference == "runA" && body.groups.contrast == "runB"
                @test body.n_samples.runA == 4 && body.n_samples.runB == 4
                @test length(body.size_factors) == 8
                byt = Dict(String(x.taxon) => x for x in body.rows)
                @test Set(keys(byt)) == Set(["Up", "Flat", "Down"])
                @test byt["Up"].estimate > 0 && byt["Down"].estimate < 0
                @test all(x -> x.padj !== nothing, body.rows)
                @test haskey(body, :figure) && length(body.figure.data) >= 1
                # Aggregation pooled the two Up rows: the route's estimate is the
                # library's on the summed counts.
                up_a = [5, 8, 6, 7]; up_b = [40, 55, 35, 60]
                mat = hcat(vcat(up_a, up_b), vcat([30, 25, 35, 28], [27, 33, 29, 31]),
                           vcat([50, 60, 45, 55], [6, 4, 9, 5]))
                direct = differential_abundance(mat, ["a$i" for i in 1:8], ["Up", "Flat", "Down"],
                                                vcat(fill("runA", 4), fill("runB", 4));
                                                reference="runA", contrast="runB")
                up = only(filter(x -> x["taxon"] == "Up", direct["rows"]))
                @test byt["Up"].estimate ≈ up["estimate"] rtol=1e-8
            end

            r = _da_post("studyD", Dict("runs" => [Dict("run" => "runA")], "table" => "merged"))
            @test r.status == 400
            @test JSON3.read(String(r.body)).error == "two_conditions_required"

            r = _da_post("noSuchStudy", Dict("runs" => [Dict("run" => "runA"), Dict("run" => "runB")]))
            @test r.status == 404
        end
    end
end
