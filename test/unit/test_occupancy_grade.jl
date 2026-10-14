# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Grading a stage against the high-water-mark algebra (MMW-01, fork issue #105)
# The algebra is proved in occupancy-types' `OccupancyTypes.Reclaim`, over the HWM
# instance in resource-algebra; `scripts/fetch-siblings.sh` pins both and
# `scripts/check-occupancy.sh` type-checks what consumes them. Nothing here
# restates a sibling's definition: what is tested is the *reading* of a probe's
# samples as grades — that a task is sequenced from the previous one's drift,
# that the siblings' laws come out of that reading, and that the prediction
# reproduces the measurement the occupancy reading recorded on 2026-10-14.

using MetaManifold.MemoryProbe
using OrderedCollections, YAML

const MG = MetaManifold.MemoryProbe

const RECORDED_OCCUPANCY =
    normpath(joinpath(@__DIR__, "..", "..", "docs", "occupancy", "results-2026-10-14.yml"))
const RECORDED_SCHEDULES =
    normpath(joinpath(@__DIR__, "..", "..", "docs", "occupancy", "schedules.yml"))

# One task as the probe records it, with the sample window a schedule gives it.
# The window is what the grade reads: a task still live when the next one starts
# hands its peak on as drift, an exited one hands on nothing.
function _occ_task(peak; label = "t", first_sample = 1, last_sample = 2)
    OrderedDict{String,Any}("label" => label, "processes" => 1, "samples" => 3,
                            "peak_rss_kib" => peak, "hwm_kib" => peak,
                            "first_sample" => first_sample, "last_sample" => last_sample)
end

# One stage as `finish!` writes it: its tasks and the peak the samples saw.
_occ_stage(tasks; measured) = OrderedDict{String,Any}(
    "probe" => MG.PROBE, "interval_ms" => 10, "samples" => 10,
    "julia_allocated_bytes" => 1000, "julia_live_bytes_delta" => 500,
    "child_peak_rss_kib" => measured, "tasks" => Any[tasks...])

@testset "occupancy grades" begin
    @testset "a grade is (peak, drift), drift ≤ peak" begin
        g = MG.Grade(90, 40)
        @test g.peak == 90
        @test g.drift == 40
        @test MG.release(g) == 50                     # what the frontier reclaims
        # The invariant the HWM instance carries as its `bounded` field.
        @test_throws ArgumentError MG.Grade(40, 90)
    end

    @testset "sequencing starts the suffix from the drift (Reclaim.reclaim-law)" begin
        a = MG.Grade(100, 40)
        b = MG.Grade(60, 10)
        s = MG.seq(a, b)
        # peak (a ⊗ b) ≡ max (peak a) (drift a + peak b)
        @test s.peak == max(a.peak, a.drift + b.peak)
        @test s.drift == a.drift + b.drift
        @test MG.seq(a, MG.NO_GRADE) == a             # the unit holds nothing
        @test s.peak <= a.peak + b.peak               # Reclaim.occupancy≤cost
        # Alternatives keep the larger of each component.
        @test MG.alt(a, b) == MG.Grade(100, 40)
    end

    @testset "the two shapes a stage has: serial and concurrent" begin
        # One after the other: the first has exited (drift = 0), so the peak is
        # the larger of the two and the cut reclaims the difference.
        serial = _occ_stage([_occ_task(100; label = "asv", first_sample = 1, last_sample = 2),
                             _occ_task(60;  label = "otu", first_sample = 3, last_sample = 4)];
                            measured = 100)
        @test MG.predicted_peak_rss_kib(serial) == 100
        @test MG.cost_bound_rss_kib(serial) == 160
        @test MG.reclaimed_kib(serial) == 60          # Reclaim.reclaimed: strictly below cost
        @test MG.disagreement(serial) === nothing

        # Together: the first is still resident when the second starts
        # (drift = peak), so the peak is the sum and nothing is reclaimed.
        together = _occ_stage([_occ_task(100; label = "asv", first_sample = 1, last_sample = 4),
                               _occ_task(60;  label = "otu", first_sample = 2, last_sample = 5)];
                              measured = 160)
        @test MG.predicted_peak_rss_kib(together) == 160
        @test MG.cost_bound_rss_kib(together) == 160
        @test MG.reclaimed_kib(together) == 0         # Reclaim.no-reclaim: drift = peak
        @test MG.disagreement(together) === nothing
    end

    @testset "the graded record says what it is" begin
        rec = _occ_stage([_occ_task(100; label = "vsearch:asv", first_sample = 1, last_sample = 4),
                          _occ_task(60;  label = "vsearch:otu", first_sample = 2, last_sample = 5)];
                         measured = 160)
        g = MG.graded(rec)
        @test g["algebra"] == "ResourceAlgebra.Instances.HWM"
        @test g["predicted_peak_kib"] == 160
        @test g["cost_bound_kib"] == 160
        @test g["reclaimed_kib"] == 0
        @test g["measured_peak_kib"] == 160
        @test g["disagreement"] === nothing
        @test [t["label"] for t in g["tasks"]] == ["vsearch:asv", "vsearch:otu"]
        @test g["tasks"][1]["drift_kib"] == 100       # still resident at the handover
        @test g["tasks"][2]["release_kib"] == 60      # the last one hands on nothing
    end

    @testset "a measured peak above the prediction is a disagreement" begin
        # The prediction is an upper bound, so this is the one contradiction it
        # allows: a sample saw more than the timed timeline can hold.
        rec = _occ_stage([_occ_task(100; first_sample = 1, last_sample = 2)]; measured = 250)
        d = MG.disagreement(rec)
        @test d isa String
        @test occursin("250", d) && occursin("100", d)
        # At or below the bound there is nothing to report.
        @test MG.disagreement(_occ_stage([_occ_task(100; first_sample = 1, last_sample = 2)];
                                         measured = 100)) === nothing
    end

    @testset "a record with no sample windows is not graded" begin
        # Records written before the windows were sampled cannot be read as a
        # timeline, so nothing is predicted for them rather than an order assumed.
        rec = OrderedDict{String,Any}("child_peak_rss_kib" => 100,
            "tasks" => Any[OrderedDict{String,Any}("label" => "t", "peak_rss_kib" => 100)])
        @test MG.predicted_peak_rss_kib(rec) === nothing
        @test MG.graded(rec) === nothing
        @test MG.disagreement(rec) === nothing        # not measured is not a disagreement
    end

    @testset "finish! grades the record it writes" begin
        # A record the probe writes itself carries its tasks' sample windows, so
        # `finish!` attaches the grades read from them; the run with real vsearch
        # children is test/occupancy_disagreement_test.jl.
        if !Sys.islinux()
            @info "occupancy grades: not Linux, so there is no /proc to sample; the finish! grading check is SKIPPED"
        else
            probe = MG.start(Dict("occupancy_probe" => Dict("enabled" => true,
                                                            "interval_ms" => 10)))
            MG.within(probe) do
                MG.task("sleep") do
                    proc = run(`sleep 0.2`; wait = false)
                    pids = MG.watch!(proc)
                    try
                        wait(proc)
                    finally
                        MG.unwatch!(pids)
                    end
                end
            end
            rec = MG.finish!(probe)
            @test length(rec["tasks"]) == 1            # one task; no stage trace
            t = rec["tasks"][1]
            @test t["label"] == "sleep" && t["samples"] > 0
            @test t["peak_rss_kib"] > 0                # sampled from /proc
            @test t["last_sample"] >= t["first_sample"] # the window watch!/unwatch! set
            @test haskey(rec, "grades")
            @test rec["grades"]["predicted_peak_kib"] == MG.predicted_peak_rss_kib(rec)
            @test rec["grades"]["measured_peak_kib"] == rec["child_peak_rss_kib"]
            @test MG.disagreement(rec) === nothing
        end
    end

    @testset "the prediction reproduces the recorded measurement" begin
        # docs/occupancy/results-2026-10-14.yml is the probe's own record of the
        # run docs/occupancy-2026-10.adoc publishes. The schedules are the ones the
        # .sir files in docs/occupancy state — serial for two orders, together for
        # the third — and the peaks are the recorded ones: nothing here is invented.
        orders = YAML.load_file(RECORDED_OCCUPANCY)["orders"]
        # The interleavings the run used, from the same file the gate reads
        # (scripts/check-occupancy-grades.jl), so the two cannot disagree.
        schedules = YAML.load_file(RECORDED_SCHEDULES)
        for (name, rec) in orders
            windows = schedules[name]
            tasks = [_occ_task(t["peak_rss_kib"]; label = t["label"],
                               first_sample = windows[i][1], last_sample = windows[i][2])
                     for (i, t) in enumerate(rec["tasks"])]
            occ    = _occ_stage(tasks; measured = rec["child_peak_rss_kib"])
            graded = MG.graded(occ)
            @test MG.predicted_peak_rss_kib(occ) == rec["child_peak_rss_kib"]
            @test MG.disagreement(occ) === nothing
            peaks = [t["peak_rss_kib"] for t in rec["tasks"]]
            # Serial: the larger peak. Together: both resident, so nothing is reclaimed.
            @test graded["reclaimed_kib"] ==
                  (name == "concurrent" ? 0 : sum(peaks) - maximum(peaks))
            @test graded["cost_bound_kib"] == sum(peaks)
        end
    end
end
