# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## R runtime mutual exclusion
# RCall hosts one embedded R session per process, and R is neither thread-safe
# nor reentrant. The pipeline stages (on spawned job threads) and the NMDS,
# PERMANOVA, and alpha-significance analyses (on HTTP handler tasks) all evaluate
# against that one interpreter.
#
# These once guarded it with two separate `_r_lock` objects, one in `Analysis` and
# one in the pipeline route, which is a data race rather than mutual exclusion: a
# pipeline stage resets the workspace with `rm(list=ls())` on acquisition and
# would destroy globals an analysis was midway through evaluating. A single user
# was enough to provoke it, since a pipeline runs on a background thread while the
# browser stays live.

using MetaManifold
SV = MetaManifold.Server

using MetaManifold.RRuntime
const RR = MetaManifold.RRuntime

@testset "R runtime" begin

    ## The regression itself: there must be exactly one lock over the one
    # interpreter, and every caller must reach it by the same path.
    @testset "exactly one lock guards the interpreter" begin
        @test !isdefined(MetaManifold.Analysis, :_r_lock)
        @test !isdefined(SV, :_r_lock)
        @test MetaManifold.Analysis.with_r_lock === RR.with_r_lock
        @test SV.with_r_lock                    === RR.with_r_lock
    end

    @testset "the lock is exclusive across tasks" begin
        @test !RR.r_busy()

        held    = Channel{Bool}(1)
        release = Channel{Bool}(1)
        holder = Threads.@spawn RR.with_r_lock() do
            put!(held, true)
            take!(release)      # hold the runtime until the test says otherwise
            :done
        end

        take!(held)
        @test RR.r_busy()
        # A pipeline stage waits indefinitely, so it cannot be tested for
        # give-up behaviour; an interactive analysis passes a timeout and must
        # refuse to hang behind the run.
        @test_throws RR.RBusyError RR.with_r_lock(() -> :never; timeout=0.1)

        put!(release, true)
        @test fetch(holder) == :done
        @test !RR.r_busy()
    end

    @testset "an uncontended acquire runs and returns" begin
        @test RR.with_r_lock(() -> 6 * 7) == 42
        @test RR.with_r_lock(() -> 6 * 7; timeout=5) == 42
        @test !RR.r_busy()
    end

    ## The lock is reentrant, so a locked helper calling another locked helper on
    # the same task must not deadlock against itself.
    @testset "the lock is reentrant within a task" begin
        result = RR.with_r_lock() do
            RR.with_r_lock(() -> :inner)
        end
        @test result == :inner
        @test !RR.r_busy()
    end

    ## Alpha significance waits on the runtime like the other R analyses, and
    # refuses once the wait runs out.
    @testset "alpha significance refuses while the runtime is busy" begin
        held    = Channel{Bool}(1)
        release = Channel{Bool}(1)
        holder = Threads.@spawn RR.with_r_lock() do
            put!(held, true)
            take!(release)
            :done
        end
        take!(held)

        previous = MetaManifold.Analysis.R_WAIT_SECONDS[]
        MetaManifold.Analysis.R_WAIT_SECONDS[] = 0.1
        try
            @test_throws RR.RBusyError MetaManifold.Analysis._alpha_significance(
                [1.0, 2.0, 3.0, 4.0], ["a", "a", "b", "b"], ["s1", "s2", "s3", "s4"])
        finally
            MetaManifold.Analysis.R_WAIT_SECONDS[] = previous
            put!(release, true)
        end
        @test fetch(holder) == :done
    end

    ## A test that never ran must not compare equal to one that ran and found
    # nothing. Both used to be `(nothing, 0-row DataFrame)`, so this assertion
    # could not be written at all.
    @testset "a test that did not run is not a null result" begin
        AN = MetaManifold.Analysis
        not_run = AN._not_computed(:r_unavailable, "R/vegan is not available")
        @test !AN.was_computed(not_run)
        @test isnothing(not_run.omnibus)
        @test nrow(not_run.pairwise) == 0
        @test !isempty(not_run.reason)

        genuine_null = AN._computed(0.87, AN._empty_pairwise())
        @test AN.was_computed(genuine_null)
        @test not_run.status !== genuine_null.status
        @test not_run != genuine_null

        # R's own tryCatch yields NA for a degenerate group: a failed test keeps
        # the pairwise rows it got but never claims an omnibus result.
        failed = AN._not_computed(:test_failed, "the omnibus statistic could not be computed";
                                  pairwise=AN._empty_pairwise())
        @test !AN.was_computed(failed)
        @test isnothing(failed.omnibus)
    end

    ## The surface must say the test was not run.
    # Drawing nothing is how "no pair reached significance" looks. A run that was
    # never performed must therefore add something to the chart, not stay silent
    # and borrow that appearance.
    @testset "a chart states that pairwise tests were not run" begin
        AN = MetaManifold.Analysis
        groups = ["a", "b"]
        vals   = Dict("a" => [1.0, 2.0, 3.0], "b" => [4.0, 5.0, 6.0])

        not_run = AN._not_computed(:r_unavailable, "R/vegan is not available")
        layout_not_run = Dict{String,Any}()
        AN._add_pairwise_annotations!(layout_not_run, "x", "y", "yaxis",
                                      groups, vals, not_run)
        notices = [a for a in get(layout_not_run, "annotations", Any[])
                   if occursin("not run", String(a["text"]))]
        @test length(notices) == 1
        @test occursin("R/vegan", String(notices[1]["text"]))

        ## The positive control that makes the assertion above mean something.
        # A genuine result in which no pair reached significance must NOT gain the
        # notice, or the notice would be noise rather than a signal.
        genuine_null = AN._computed(0.87, AN._empty_pairwise())
        layout_null = Dict{String,Any}()
        AN._add_pairwise_annotations!(layout_null, "x", "y", "yaxis",
                                      groups, vals, genuine_null)
        @test isempty([a for a in get(layout_null, "annotations", Any[])
                       if occursin("not run", String(a["text"]))])
    end

    ## Each reason gets its own wording, because they call for different actions:
    # an absent runtime is a deployment fault, unpairable groups are a property
    # of the data, and a failed statistic is a property of these groups.
    @testset "the caption distinguishes the reasons a test did not run" begin
        AN = MetaManifold.Analysis
        failed  = AN._not_computed(:test_failed, "the omnibus statistic could not be computed")
        absent  = AN._not_computed(:r_unavailable, "R/vegan is not available")
        unpaired = AN._not_computed(:no_paired_samples, "the groups share no sample IDs")

        for r in (failed, absent, unpaired)
            @test occursin("not run", AN._significance_caption(r, "KW"))
            @test !occursin("n/a", AN._significance_caption(r, "KW"))
        end
        @test AN._significance_caption(failed, "KW") != AN._significance_caption(absent, "KW")
        @test AN._significance_caption(absent, "KW") != AN._significance_caption(unpaired, "KW")

        # A computed result still reads as a result.
        computed = AN._computed(0.02, AN._empty_pairwise())
        caption  = AN._significance_caption(computed, "KW")
        @test occursin("p = 0.02", caption)
        @test !occursin("not run", caption)
    end

end
