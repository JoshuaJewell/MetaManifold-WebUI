# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# PERMANOVA: each covariate is tested as its own sequential term, checked against
# a direct vegan::adonis2 call; the permutation count and p-value floor are
# reported; and a term within-block permutation cannot test is refused.
# Skips loudly when R/vegan is absent.

using MetaManifold
using MetaManifold.Analysis
using RCall
using DataFrames
using Random
using Test

const _PERM_R = MetaManifold.Server.r_available()
_PERM_R || @info "R/vegan unavailable - the PERMANOVA reference tests are SKIPPED"

"""
    _adonis2_reference(mat, meta, perm_expr) -> DataFrame

Run vegan::adonis2 directly on Bray-Curtis distances of `mat`, with
`group + run` terms tested sequentially and permutations given by the R
expression `perm_expr` (which may refer to `ref_meta`), seeded as
`run_permanova` seeds it. Returns the term rows: term, r2, f, p.
"""
function _adonis2_reference(mat, meta, perm_expr)
    MetaManifold.RRuntime.with_r_lock() do
        RCall.globalEnv[:ref_mat] = mat
        RCall.globalEnv[:ref_meta] = meta
        RCall.reval("""
            set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
            ref_d <- vegdist(ref_mat, method = "bray")
            ref_res <- adonis2(ref_d ~ group + run, data = ref_meta,
                               permutations = $perm_expr, by = "terms", parallel = 1)
            ref_rows <- which(!(rownames(ref_res) %in% c("Residual", "Total")))
            ref_out <- data.frame(term = rownames(ref_res)[ref_rows], r2 = ref_res\$R2[ref_rows],
                                  f = ref_res\$F[ref_rows], p = ref_res[["Pr(>F)"]][ref_rows])
        """)
        out = rcopy(RCall.reval("ref_out"))
        RCall.reval("rm(ref_mat, ref_meta, ref_d, ref_res, ref_rows, ref_out)")
        out
    end
end

@testset "PERMANOVA" begin
    @testset "_varies_within_any_block" begin
        blocks = ["i1", "i1", "i2", "i2"]
        @test !Analysis._varies_within_any_block(["a", "a", "b", "b"], blocks)
        @test Analysis._varies_within_any_block(["r1", "r2", "r1", "r1"], blocks)
    end

    if _PERM_R
        rng = MersenneTwister(7)
        # Six individuals, each sampled in two runs; three individuals per group.
        mat = Float64.(rand(rng, 0:60, 12, 40))
        individuals = ["i$(cld(k, 2))" for k in 1:12]
        meta = DataFrame(sample = ["s$k" for k in 1:12],
                         group = [cld(k, 2) <= 3 ? "A" : "B" for k in 1:12],
                         run = repeat(["r1", "r2"], 6))

        @testset "each term is tested on its own, as vegan by = \"terms\"" begin
            res = Analysis.run_permanova(mat, meta; seed=123)
            ref = _adonis2_reference(mat, meta, "999")
            @test [t.term for t in res.terms] == ["group", "run"] == ref.term
            @test [t.r2 for t in res.terms] ≈ ref.r2
            @test [t.f_statistic for t in res.terms] ≈ ref.f
            @test [t.p_value for t in res.terms] == ref.p
            # The headline result is the group term, not the whole model.
            @test res.term == "group"
            @test res.r2 ≈ ref.r2[1] && res.p_value == ref.p[1]
            @test res.permutations == 999
            @test res.min_p_value ≈ 1 / 1000
            @test all(t -> isnothing(t.untestable_reason), res.terms)
        end

        @testset "within-individual blocks: permutation count, floor and refusal" begin
            res = Analysis.run_permanova(mat, meta; seed=123, blocks=individuals)
            ref = _adonis2_reference(mat, meta,
                                     "permute::how(blocks = factor(rep(paste0('i', 1:6), each = 2)), nperm = 999)")
            @test res.blocked
            # Six blocks of two admit 2^6 orderings, the observed one excluded.
            @test res.permutations == 63
            @test res.min_p_value ≈ 1 / 64
            # group is constant within every individual, so within-block
            # permutation never exchanges its labels. vegan still prints an
            # ordinary-looking probability for it (the control: that is the
            # number withheld here), but it does not test a group difference.
            @test 0 < ref.p[1] <= 1
            group = res.terms[1]
            @test group.term == "group" && isnothing(group.p_value)
            @test occursin("constant within every permutation block", group.untestable_reason)
            @test isnothing(res.p_value) && res.untestable_reason == group.untestable_reason
            @test group.r2 ≈ ref.r2[1]
            @test occursin("withheld for group", res.text)
            # run varies within individuals and is tested as vegan tests it.
            run = res.terms[2]
            @test run.term == "run" && isnothing(run.untestable_reason)
            @test run.p_value == ref.p[2]
            @test run.f_statistic ≈ ref.f[2]
        end
    end
end
