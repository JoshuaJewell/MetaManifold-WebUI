# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: AGPL-3.0-or-later

# Native (CompositionalDA-backed) differential abundance. The configuration
# boundary fuzz and the refusals run everywhere; the adapter comparisons against
# the package run only where it is loadable.
#
# CompositionalDA is not in the General registry yet, so it is not a declared
# dependency (`Project.toml` records why) and `Pkg.instantiate()` does not
# install it: on such a checkout `compositional_da_available()` is false, the two
# native methods refuse by name, and the adapter tests below skip loudly.
# `Pkg.develop(path="/…/CompositionalDA.jl")` in a scratch environment activates
# them.

using MetaManifold, Test
using Random
using MetaManifold.Differential

# "Not installed" is an ArgumentError. Any other failure while loading the
# package (a broken install, say) is a real problem and must not read as a skip.
const HAVE_COMPOSITIONAL_DA = try
    @eval import CompositionalDA
    true
catch err
    err isa ArgumentError || err isa LoadError || rethrow()
    false
end
HAVE_COMPOSITIONAL_DA || @info "CompositionalDA unavailable - the native adapter tests are SKIPPED"

# The module's own probe must agree with the one above, or the refusal path and
# the tests disagree about what is installed.
@test HAVE_COMPOSITIONAL_DA == compositional_da_available()

@testset "Native configuration boundary fuzz" begin
    rng = MersenneTwister(23)
    SV = MetaManifold.Server
    for _ in 1:256
        mc = rand(rng,-128:4352)
        seed = rand(rng,(-1,0,1,17,Int64(typemax(Int32)),Int64(typemax(Int32))+1))
        @test isnothing(SV._config_value_error("analysis.differential.mc_samples",mc)) == (2 <= mc <= 4096)
        @test isnothing(SV._config_value_error("analysis.differential.seed",seed)) == (0 <= seed <= typemax(Int32))
        if 2 <= mc <= 4096 && 0 <= seed <= typemax(Int32)
            config = DifferentialConfig(;method="aldex2",mc_samples=mc,seed)
            @test config.mc_samples == mc && config.seed == seed
        else
            @test_throws ArgumentError DifferentialConfig(;method="aldex2",mc_samples=mc,seed)
        end
    end
    for invalid in (NaN,Inf,-Inf,nothing,missing,[],Dict(),true,false,"128",2.5)
        @test !isnothing(SV._config_value_error("analysis.differential.mc_samples",invalid))
        @test !isnothing(SV._config_value_error("analysis.differential.seed",invalid))
    end
end

# samples x taxa, two groups, four taxa: a valid input for both native methods.
const NATIVE_COUNTS = [10 20 30 8; 12 21 27 9; 9 18 31 7; 40 19 29 8; 42 20 30 9; 38 22 28 7]
const NATIVE_SAMPLES = ["s$i" for i in 1:6]
const NATIVE_TAXA = ["a", "b", "c", "d"]
# Reverse lexical order is deliberate: the user's reference is authoritative.
const NATIVE_GROUPS = ["z", "z", "z", "a", "a", "a"]

@testset "Native input boundaries do not depend on the package" begin
    config = DifferentialConfig(; method="aldex2", mc_samples=16, seed=19)
    # The oversize refusal is an input error raised before the dependency probe,
    # so it is reachable whether or not CompositionalDA is installed.
    for T in (BigInt, Float64)
        oversized = T.(NATIVE_COUNTS)
        oversized[1, 1] = T(big(typemax(Int)) + 1)
        @test_throws ArgumentError differential_abundance(oversized, NATIVE_SAMPLES, NATIVE_TAXA,
            NATIVE_GROUPS; reference="z", contrast="a", config)
    end
    @test_throws ArgumentError DifferentialConfig(; method="aldex2", mc_samples=1)
    @test_throws ArgumentError DifferentialConfig(; method="aldex2", seed=-1)
    @test_throws ArgumentError DifferentialConfig(; method="ancombc", dispersion="glmGamPoi")
end

@testset "Native methods refuse by name when the package is absent" begin
    if HAVE_COMPOSITIONAL_DA
        @test_skip false    # the adapter tests below cover the available case
    else
        for method in ("aldex2", "ancombc")
            config = DifferentialConfig(; method)
            err = try
                differential_abundance(NATIVE_COUNTS, NATIVE_SAMPLES, NATIVE_TAXA,
                    NATIVE_GROUPS; reference="z", contrast="a", config)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("CompositionalDA", sprint(showerror, err))
            @test occursin(method, sprint(showerror, err))
            # The message has to point at what does work, not just refuse.
            @test occursin("nb_glm", sprint(showerror, err))
        end
    end
end

@testset "Native compositional DA adapter" begin
    if !HAVE_COMPOSITIONAL_DA
        @test_skip false
    else
        counts = NATIVE_COUNTS
        samples = NATIVE_SAMPLES
        taxa = NATIVE_TAXA
        groups = NATIVE_GROUPS
        config = DifferentialConfig(; method="aldex2", mc_samples=16, seed=19)
        actual = differential_abundance(counts, samples, taxa, groups;
            reference="z", contrast="a", config)
        expected = CompositionalDA.aldex2(counts, taxa, ["0", "0", "0", "1", "1", "1"];
            mc_samples=16, seed=19)
        for row in actual["rows"]
            i = findfirst(==(row["taxon"]), expected.table.feature)
            @test row["estimate"] == expected.table.diff_btw[i]
            @test row["pvalue"] == expected.table.pvalue[i]
            @test row["padj"] == expected.table.qvalue[i]
            @test row["log2_fold_change"] === nothing
            @test row["aldex_effect"] == expected.table.effect[i]
            @test row["aldex_diff_within"] == expected.table.diff_win[i]
            @test row["aldex_abundance"] == expected.table.rab_all[i]
            @test row["wilcoxon_pvalue"] == expected.table.wi_pvalue[i]
            @test row["wilcoxon_padj"] == expected.table.wi_qvalue[i]
        end
        @test actual["config"]["seed"] == 19
        @test actual["config"]["mc_samples"] == 16
        @test isempty(actual["size_factors"])
        @test actual["groups"]["reference"] == "z"
        @test actual["diagnostics"]["native_method"]["package"] == "CompositionalDA"
    end
end

@testset "ANCOM-BC adapter preserves units and structural decisions" begin
    if !HAVE_COMPOSITIONAL_DA
        @test_skip false
    else
        rng = MersenneTwister(2026)
        counts = round.(Int, permutedims(collect(range(1000, 5000; length=40))) .*
            exp.(0.25 .* randn(rng, 20, 40)))
        counts[11:20, 1:3] .*= 8
        counts[11:20, 4:5] .= max.(1, counts[11:20, 4:5] .÷ 4)
        counts = hcat(counts, vcat(zeros(Int, 10), fill(200, 10)), zeros(Int, 20))
        taxa = vcat(["f$j" for j in 1:40], ["structural", "all_zero"])
        groups = vcat(fill("z", 10), fill("a", 10))
        actual = differential_abundance(counts, ["s$i" for i in 1:20], taxa, groups;
            reference="z", contrast="a", config=DifferentialConfig(; method="ancombc"))
        expected = CompositionalDA.ancombc(counts, taxa, vcat(fill("0", 10), fill("1", 10));
            prevalence=0.0, structural_zeros=true)
        by_taxon = Dict(row["taxon"] => row for row in actual["rows"])
        for i in 1:40
            row = by_taxon[taxa[i]]
            @test row["estimate"] == expected.table.effect[i]
            @test row["log2_fold_change"] == expected.table.effect[i] / log(2)
            @test row["standard_error"] == expected.table.se[i]
            @test row["padj"] == expected.table.qvalue[i]
        end
        @test by_taxon["structural"]["status"] == "structural_zero"
        @test by_taxon["structural"]["pvalue"] === nothing
        @test by_taxon["structural"]["padj"] === nothing
        @test by_taxon["all_zero"]["status"] == "filtered"
        @test actual["diagnostics"]["n_structural_zero"] == 1
        @test actual["diagnostics"]["n_filtered"] == 1
        @test actual["diagnostics"]["n_tested"] == 40
    end
end

## The adapter, exercised against a stand-in backend
#
# The two testsets above run only where CompositionalDA is installed, and CI
# cannot install it: the package is not in the General registry. Without the
# testsets below, nothing in CI runs the adapter that turns a DAResult into
# result rows, or the branch that keeps BH from being applied twice.
#
# `differential_abundance(...; native_backend=…)` makes the adapter reachable
# without changing what production does (the default is still the CompositionalDA
# this session loaded, or the refusal when it did not). The stand-in in
# `test/support/native_da_stub.jl` offers the real package's entry points,
# columns and verdict strings, so what runs below is the code that runs in
# production; only the estimator behind the boundary is a stand-in, and the
# result names it under `diagnostics.native_method.package`.

@testset "The native backend seam selects the module, and defaults to the package" begin
    @test MetaManifold.Differential._native_backend(NativeDAStub) === NativeDAStub
    @test isnothing(MetaManifold.Differential._native_backend(nothing)) == !HAVE_COMPOSITIONAL_DA
    # An available backend satisfies the guard; an absent one refuses by name.
    @test isnothing(MetaManifold.Differential._native_methods_available("aldex2", NativeDAStub))
    err = try
        MetaManifold.Differential._native_methods_available("ancombc", nothing)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("CompositionalDA", sprint(showerror, err))
end

@testset "Native ALDEx2 rows are wired from the backend's own columns" begin
    NativeDAStub.reset!()
    config = DifferentialConfig(; method="aldex2", mc_samples=32, seed=7)
    actual = differential_abundance(NATIVE_COUNTS, NATIVE_SAMPLES, NATIVE_TAXA, NATIVE_GROUPS;
        reference="z", contrast="a", config, native_backend=NativeDAStub)

    call = NativeDAStub.last_call()
    @test call.method === :aldex2
    @test call.settings[:mc_samples] == 32
    @test call.settings[:seed] == 7
    # The user's groups are relabelled "0"/"1" in the requested direction, so the
    # package's sorted contrast is the contrast the user asked for even though
    # the reference ("z") sorts after the contrast ("a").
    @test call.groups == ["0", "0", "0", "1", "1", "1"]
    @test call.features == NATIVE_TAXA
    @test size(call.counts) == size(NATIVE_COUNTS)
    @test eltype(call.counts) <: Integer

    rows = Dict(String(r["taxon"]) => r for r in actual["rows"])
    @test length(rows) == length(NATIVE_TAXA)
    for (j, taxon) in enumerate(NATIVE_TAXA)
        row = rows[taxon]
        @test row["status"] == "ok"
        @test row["note"] == ""
        @test row["estimate"] == 0.5 * j                 # diff_btw
        @test row["aldex_effect"] == 0.25 * j            # effect
        @test row["aldex_diff_within"] == 0.125 * j      # diff_win
        @test row["aldex_abundance"] == 1.0 * j          # rab_all
        @test row["wilcoxon_pvalue"] == 0.0075 * j       # wi_pvalue
        @test row["wilcoxon_padj"] == 0.00375 * j        # wi_qvalue
        @test row["pvalue"] == 0.01 * j
        # ALDEx2 averages BH within each posterior instance. Applying BH a second
        # time to its expected p-values would implement a different method, so
        # padj is the package's qvalue even where BH would disagree.
        @test row["padj"] == 0.005 * j
        @test row["log2_fold_change"] === nothing
        @test row["standard_error"] === nothing
        @test row["statistic"] === nothing
    end
    @test [rows[t]["padj"] for t in NATIVE_TAXA] != bh_adjust([0.01 * j for j in 1:4])
    @test actual["effect"]["key"] == "estimate"
    @test occursin("ALDEx2", actual["method"])
    @test occursin("32 Dirichlet instances", actual["method"])
    @test isempty(actual["size_factors"])
    @test actual["config"]["method"] == "aldex2"
    @test actual["config"]["mc_samples"] == 32
    @test actual["config"]["seed"] == 7
    @test actual["diagnostics"]["native_method"]["package"] == "NativeDAStub"
    @test actual["diagnostics"]["native_method"]["version"] == "0.0.0-stub"
    @test actual["diagnostics"]["native_method"]["parameters"]["mc_samples"] == 32
    @test actual["diagnostics"]["n_structural_zero"] == 0
    @test actual["diagnostics"]["n_tested"] == 4

    # A feature the caller has no row for is skipped, not invented.
    NativeDAStub.reset!()
    NativeDAStub.EXTRA_FEATURE[] = true
    extra = differential_abundance(NATIVE_COUNTS, NATIVE_SAMPLES, NATIVE_TAXA, NATIVE_GROUPS;
        reference="z", contrast="a", config, native_backend=NativeDAStub)
    @test length(extra["rows"]) == length(NATIVE_TAXA)
    @test !any(r -> String(r["taxon"]) == "not_requested", extra["rows"])
    NativeDAStub.reset!()
end

@testset "Native ANCOM-BC keeps units and reports structural decisions separately" begin
    NativeDAStub.reset!()
    samples = ["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8"]
    taxa = ["spike", "flat", "structural", "half", "structural_both", "absent"]
    groups = ["z", "z", "z", "z", "a", "a", "a", "a"]
    counts = [10 50 0  8   3 0;
              10 50 0  8   3 0;
              10 50 0  9   4 0;
              13 50 0  7   4 0;
              80 50 9  9   3 0;
              80 50 9  8   3 0;
              83 50 8  9   4 0;
              77 50 8  8   3 0]
    config = DifferentialConfig(; method="ancombc")
    actual = differential_abundance(counts, samples, taxa, groups;
        reference="z", contrast="a", config, native_backend=NativeDAStub)

    call = NativeDAStub.last_call()
    @test call.method === :ancombc
    @test call.settings[:prevalence] == 0.0
    @test call.settings[:structural_zeros] === true
    @test call.settings[:conserve] === true
    # "absent" has no reads at all: it is filtered before the backend is called.
    @test call.features == ["spike", "flat", "structural", "half", "structural_both"]
    @test call.groups == ["0", "0", "0", "0", "1", "1", "1", "1"]

    rows = Dict(String(r["taxon"]) => r for r in actual["rows"])
    @test rows["absent"]["status"] == "filtered"
    @test isnothing(rows["absent"]["pvalue"])

    for (j, taxon) in enumerate(call.features)
        row = rows[taxon]
        if row["status"] == "structural_zero"
            @test isnothing(row["estimate"])
            @test isnothing(row["pvalue"])
            @test isnothing(row["padj"])
            continue
        end
        @test row["status"] == "ok"
        @test row["estimate"] == log(2) * j / 2
        @test row["log2_fold_change"] ≈ j / 2
        @test row["standard_error"] == 0.1 * j
        @test row["statistic"] == (log(2) * j / 2) / (0.1 * j)
        @test row["pvalue"] == 0.02 * j
        @test row["padj"] == 0.005 * j
    end

    # Structural absences are decisions, not tests: named as such, and carrying
    # no stand-in p-value.
    @test rows["structural"]["status"] == "structural_zero"
    @test occursin("one condition", rows["structural"]["note"])
    @test rows["structural_both"]["status"] == "structural_zero"
    @test occursin("both conditions", rows["structural_both"]["note"])
    @test actual["diagnostics"]["n_structural_zero"] == 2
    @test actual["diagnostics"]["n_filtered"] == 1
    @test actual["diagnostics"]["n_tested"] == 3   # spike, flat, half
    @test actual["effect"]["key"] == "log2_fold_change"
    @test occursin("ANCOM-BC", actual["method"])
    @test actual["config"]["method"] == "ancombc"
    # ALDEx2's sampling settings have no meaning here and must not be echoed.
    @test !haskey(actual["config"], "mc_samples")
    @test actual["diagnostics"]["native_method"]["package"] == "NativeDAStub"
    @test actual["diagnostics"]["native_method"]["parameters"]["structural_zeros"] === true
    NativeDAStub.reset!()
end

@testset "Native refusals happen in the right order, whatever the backend" begin
    # Input problems are refused before the backend is consulted, so a malformed
    # matrix gets the same answer on a machine with the package and without it.
    one_taxon = [10 0; 12 0; 9 0; 40 0; 42 0; 38 0]
    for method in ("aldex2", "ancombc")
        config = DifferentialConfig(; method)
        err = try
            differential_abundance(one_taxon, NATIVE_SAMPLES, ["a", "b"], NATIVE_GROUPS;
                reference="z", contrast="a", config, native_backend=NativeDAStub)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("at least two retained taxa", sprint(showerror, err))
    end
    @test isempty(NativeDAStub.CALLS)

    config = DifferentialConfig(; method="aldex2")

    # The package's own refusal is translated into an ErrorException carrying
    # its message, which the route reports as 422 no_taxon_fitted.
    NativeDAStub.reset!()
    NativeDAStub.FAILURE[] = :composition
    err = try
        differential_abundance(NATIVE_COUNTS, NATIVE_SAMPLES, NATIVE_TAXA, NATIVE_GROUPS;
            reference="z", contrast="a", config, native_backend=NativeDAStub)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("CompositionError", err.msg)

    # Any other failure is the backend's own and is rethrown untouched: a broken
    # install must not be dressed up as a refusal of the data.
    NativeDAStub.reset!()
    NativeDAStub.FAILURE[] = :other
    err = try
        differential_abundance(NATIVE_COUNTS, NATIVE_SAMPLES, NATIVE_TAXA, NATIVE_GROUPS;
            reference="z", contrast="a", config, native_backend=NativeDAStub)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("a reason of its own", err.msg)
    @test !occursin("CompositionError", err.msg)

    # A statistic that is not finite is refused by name rather than published.
    NativeDAStub.reset!()
    NativeDAStub.FAILURE[] = :nonfinite
    err = try
        differential_abundance(NATIVE_COUNTS, NATIVE_SAMPLES, NATIVE_TAXA, NATIVE_GROUPS;
            reference="z", contrast="a", config, native_backend=NativeDAStub)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("non-finite statistics", err.msg)
    @test occursin("aldex2", err.msg)
    NativeDAStub.reset!()
end
