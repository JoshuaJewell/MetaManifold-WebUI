# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# MMW-01: the occupancy probe, and the two-thread disagreement test.
#
# The vsearch stage runs two independent searches, ASV and OTU, concurrently in
# the full pipeline and one after the other in a single-stage re-run. The
# high-water-mark algebra (occupancy-types, resource-algebra `Instances/HWM`)
# predicts:
#
#   * cost commutes: the work done, measured as Julia's allocation total for the
#     stage, is the same in every order (choreographic-types A-2). This is
#     ASSERTED, and it can fail.
#   * each search's own peak is the same in every order;
#   * the stage's peak does not commute: serially it is the larger of the two
#     searches' peaks (each child exits holding nothing, so its drift is 0), while
#     concurrently it can approach their sum (A-3). These are RECORDED in the test
#     log, not asserted.
#
# If the stage peak came out the same in every order, the prediction would be
# falsified for this workload: the searches' memory would not overlap.
#
# The fixture is the smallest of a calibrated ladder (250, 500, 1000, 2000
# reference sequences of 400 bp; 300 ASV and 150 OTU queries) for which, in five
# repeats each on a 4-core host, both searches were resident together for at
# least 15 samples at 10 ms and the concurrent stage peak was at least 1.9 times
# the larger search's: 1000 sequences (17-20 samples, 2.0 times). 250 and 500
# also overlap, by 12-18 samples, too close to the edge for a slower runner. It is
# generated here, deterministically, so no data is checked in.
#
# vsearch runs with --threads 1 here. With several threads it writes its hits in
# whatever order its threads finish, and `_sort_hits_by_query!` then does a
# different amount of work (its sort key splits each line on every comparison),
# which moved the stage's allocation total by up to 12% between repeats of the
# SAME order. On one thread the output order is fixed, and the three orders
# agreed to within 0.08% over five repeats each. Run directly with
#   julia --project=. -t 2 test/occupancy_disagreement_test.jl
# and set OCCUPANCY_RESULTS=<path> to also write the numbers as YAML.

using Test, YAML, OrderedCollections
using MetaManifold

const MP = MetaManifold.MemoryProbe
const PV = MetaManifold.Provenance
const TL = MetaManifold.Tools
# Not const: other test files bind SV in Main the same way, and runtests.jl
# includes them all into one module.
SV = MetaManifold.Server

const PROBE_ON = Dict("occupancy_probe" => Dict("enabled" => true, "interval_ms" => 10))

# Relative spread of the stage's Julia allocation total across the three orders
# that still counts as "the same cost". Measured: under 0.08% (see above).
const COST_TOLERANCE = 0.01

## Fixture
# A 64-bit LCG, so the fixture is the same bytes on every Julia version (the
# streams of Random's generators are not guaranteed to be).
mutable struct _LCG
    s :: UInt64
end
_next!(g::_LCG) = (g.s = g.s * 0x5851f42d4c957f2d + 0x14057b7ef767814f; g.s >> 33)

const _NT = (UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T'))
_randseq(g, n) = String([_NT[_next!(g) % 4 + 1] for _ in 1:n])

# A copy of `s` with about `per_1000` substitutions per thousand bases: close
# enough to its reference (~96% identity) to pass vsearch's --id 0.75.
function _mutate(g, s, per_1000)
    b = Vector{UInt8}(s)
    for i in eachindex(b)
        _next!(g) % 1000 < per_1000 && (b[i] = _NT[_next!(g) % 4 + 1])
    end
    String(b)
end

function occupancy_fixture(dir; refs = 1000, asv = 300, otu = 150, len = 400)
    g = _LCG(0x4d4d57)
    mkpath(dir)
    seqs = [_randseq(g, len) for _ in 1:refs]
    open(joinpath(dir, "ref.fasta"), "w") do io
        for (i, s) in enumerate(seqs)
            println(io, ">ref$i;tax=k:K$(i % 7)")
            println(io, s)
        end
    end
    for (name, n) in (("asv", asv), ("otu", otu))
        open(joinpath(dir, "$name.fasta"), "w") do io
            for i in 1:n
                println(io, ">$name$i")
                println(io, _mutate(g, seqs[_next!(g) % refs + 1], 40))
            end
        end
    end
    dir
end

## Running the stage in one order
function _search(dir, out, name)
    TL.vsearch(joinpath(dir, "$name.fasta"), joinpath(dir, "ref.fasta"), joinpath(out, name);
               optional_args = "--id 0.75 --query_cov 0.8 --threads 1",
               vsearch_bin = something(Sys.which("vsearch"), "vsearch"))
end

const ORDERS = ("asv_then_otu", "otu_then_asv", "concurrent")

# The stage's work in `order`, under a fresh probe. Returns the probe's record.
# ASV-then-OTU and concurrent are the stage's own two shapes, run through the
# pipeline's _vsearch_searches; OTU-then-ASV has no caller and is spelled out.
function run_order(dir, out, order)
    mkpath(out)
    probe = MP.start(PROBE_ON)
    asv() = _search(dir, out, "asv")
    otu() = _search(dir, out, "otu")
    MP.within(probe) do
        if order == "otu_then_asv"
            MP.task(otu, "vsearch:otu")
            MP.task(asv, "vsearch:asv")
        else
            SV._vsearch_searches(asv, otu; concurrent = order == "concurrent")
        end
    end
    MP.finish!(probe)
end

_task(rec, label) = only(t for t in rec["tasks"] if t["label"] == label)

@testset "Occupancy probe (MMW-01)" begin

    @testset "off unless the run config enables it" begin
        @test !MP.enabled(Dict())
        @test !MP.enabled(Dict("occupancy_probe" => Dict()))
        @test !MP.enabled(Dict("occupancy_probe" => Dict("enabled" => "true")))
        @test !MP.enabled(Dict("occupancy_probe" => true))
        @test MP.enabled(PROBE_ON)
        @test MP.start(Dict()) === nothing
        @test MP.finish!(nothing) === nothing
        @test MP.within(() -> 42, nothing) == 42
        @test MP.task(() -> 7, "outside a probe") == 7
        @test MP.interval_seconds(Dict()) == 0.010
        @test MP.interval_seconds(Dict("occupancy_probe" => Dict("interval_ms" => 25))) == 0.025
        @test MP.interval_seconds(Dict("occupancy_probe" => Dict("interval_ms" => 0))) == 0.010
    end

    @testset "flag off: the stage record is unchanged" begin
        # The keys a stage record had before the probe existed. With the probe
        # off, `finish!` gives `nothing` and the record must be exactly these.
        before = ["name", "started", "finished", "status", "degraded",
                  "degraded_components", "runtimes", "tools", "databases",
                  "commands", "outputs"]
        att = PV.Attestation()
        PV.record_stage!(att, "vsearch", PV.CapturedEnvironment();
                         occupancy = MP.finish!(MP.start(Dict())))
        @test collect(keys(only(att.doc["stages"]))) == before

        # Negative control: with the probe on, the key is there.
        on = PV.Attestation()
        PV.record_stage!(on, "vsearch", PV.CapturedEnvironment();
                         occupancy = MP.finish!(MP.start(PROBE_ON)))
        @test collect(keys(only(on.doc["stages"]))) == [before; "occupancy"]
    end

    @testset "attestation round trip; older records still parse" begin
        mktempdir() do tmp
            att = PV.Attestation()
            rec = MP.finish!(MP.start(PROBE_ON))
            PV.record_stage!(att, "vsearch", PV.CapturedEnvironment(); occupancy = rec)
            PV.record_stage!(att, "swarm", PV.CapturedEnvironment())   # an unprobed stage
            path = joinpath(tmp, "attestation.yml")
            PV.write_attestation(att, path)
            back = PV.read_attestation(path)
            stages = Dict(s["name"] => s for s in back.doc["stages"])
            @test stages["vsearch"]["occupancy"]["probe"] == "MemoryProbe 1"
            @test stages["vsearch"]["occupancy"]["interval_ms"] == 10
            @test !haskey(stages["swarm"], "occupancy")

            # A record written before the probe existed merges with a new one.
            old = PV.Attestation()
            PV.record_stage!(old, "cutadapt", PV.CapturedEnvironment())
            PV.write_attestation(old, path; merge_existing = false)
            PV.write_attestation(att, path)
            merged = Dict(s["name"] => s for s in PV.read_attestation(path).doc["stages"])
            @test !haskey(merged["cutadapt"], "occupancy")
            @test haskey(merged["vsearch"], "occupancy")
        end
    end

    @testset "the vsearch stage's two shapes" begin
        # Serial and concurrent return the same pair, and record one task per
        # search; with no OTU search the pair ends in `nothing` and no OTU task
        # is recorded. A single-stage re-run attests both searches' outputs.
        for concurrent in (false, true)
            probe = MP.start(PROBE_ON)
            @test MP.within(() -> SV._vsearch_searches(() -> :asv, () -> :otu; concurrent), probe) ==
                  (:asv, :otu)
            @test MP.within(() -> SV._vsearch_searches(() -> :asv, nothing; concurrent), probe) ==
                  (:asv, nothing)
            @test sort([t["label"] for t in MP.finish!(probe)["tasks"]]) ==
                  ["vsearch:asv", "vsearch:asv", "vsearch:otu"]
        end

        # The pipeline's call: the same search on each result, against one database.
        stub(project, result, db) = (project, result, db)
        for concurrent in (false, true)
            @test SV._vsearch_stage(stub, :p, :asvs, :otus, :db; concurrent) ==
                  ((:p, :asvs, :db), (:p, :otus, :db))
            @test SV._vsearch_stage(stub, :p, :asvs, nothing, :db; concurrent) ==
                  ((:p, :asvs, :db), nothing)
        end

        # Concurrently, a failed ASV search is rethrown only once the OTU search
        # has finished, so it cannot be left running.
        otu_done = Threads.Atomic{Bool}(false)
        @test_throws TaskFailedException SV._vsearch_searches(() -> error("ASV search failed"),
                                                               () -> (sleep(0.3); otu_done[] = true; :otu);
                                                               concurrent = true)
        @test otu_done[]
    end

    @testset "a handle the probe cannot watch is ignored" begin
        # run(cmd; wait=false) returns a Process or a ProcessChain; anything else
        # must leave the stage alone rather than fail it.
        @test MP.watch!(:not_a_process) === nothing
        @test MP.unwatch!(MP.watch!(:not_a_process)) === nothing
    end

    if !Sys.islinux()
        @info "Occupancy probe: not Linux, so there is no /proc to sample; the child-memory tests are SKIPPED"
    else
        @testset "samples a child's resident memory" begin
            probe = MP.start(PROBE_ON)
            MP.within(probe) do
                # bash holds a 64 MiB string for half a second.
                MP.task("64MiB") do
                    TL._run_killable(`bash -c 'x=$(head -c 67108864 /dev/zero | tr "\0" a); sleep 0.5; : "${#x}"'`)
                end
                # Negative control: a child that allocates nothing.
                MP.task("sleep") do
                    TL._run_killable(`sleep 0.5`)
                end
            end
            rec = MP.finish!(probe)
            big, small = _task(rec, "64MiB"), _task(rec, "sleep")
            @test big["processes"] == 1 && big["samples"] > 0
            @test big["peak_rss_kib"] >= 64 * 1024
            @test big["hwm_kib"] >= big["peak_rss_kib"]
            @test small["samples"] > 0
            @test 0 < small["peak_rss_kib"] < 16 * 1024
            # Serial: the stage peak is the larger task's, not the sum.
            @test rec["child_peak_rss_kib"] < big["peak_rss_kib"] + small["peak_rss_kib"]
            @test !haskey(rec, "error")
        end

        @testset "through the pipeline's _attest, flag off and on" begin
            mktempdir() do tmp
                cfg_path = joinpath(tmp, "run_config.yml")
                for probed in (false, true)
                    cfg = Dict{String,Any}("vsearch" => Dict("enabled" => true))
                    probed && (cfg["occupancy_probe"] = Dict("enabled" => true))
                    YAML.write_file(cfg_path, cfg)
                    att = PV.Attestation(; config = cfg, config_sha256 = PV.file_sha256(cfg_path))
                    SV._attest(att, PV.CapturedEnvironment(), "vsearch"; config_path = cfg_path) do
                        MetaManifold.Config._begin_section(cfg_path, "vsearch", joinpath(tmp, "vsearch.hash"))
                        MP.task(() -> TL._run_killable(`sleep 0.2`), "vsearch:asv")
                        nothing
                    end
                    stage = only(att.doc["stages"])
                    @test stage["status"] == "ok"
                    @test haskey(stage, "occupancy") == probed
                    probed && @test only(stage["occupancy"]["tasks"])["label"] == "vsearch:asv"

                    # A stage that throws is recorded as failed, with its
                    # occupancy when probed, and the error still propagates.
                    failed = PV.Attestation(; config = cfg, config_sha256 = PV.file_sha256(cfg_path))
                    @test_throws ProcessFailedException SV._attest(failed, PV.CapturedEnvironment(), "swarm";
                                                                    config_path = cfg_path) do
                        TL._run_killable(`false`)
                    end
                    stage = only(failed.doc["stages"])
                    @test stage["status"] == "failed"
                    @test haskey(stage, "occupancy") == probed
                end
            end
        end

        @testset "a pipeline of commands counts each process" begin
            probe = MP.start(PROBE_ON)
            MP.within(() -> MP.task(() -> TL._run_killable(pipeline(`sleep 0.3`, `cat`)), "chain"), probe)
            t = _task(MP.finish!(probe), "chain")
            @test t["processes"] == 2
            @test t["samples"] > 0
        end

        @testset "a child that exits at once is counted, not sampled" begin
            probe = MP.start(PROBE_ON)
            MP.within(() -> MP.task(() -> TL._run_killable(`true`), "instant"), probe)
            t = _task(MP.finish!(probe), "instant")
            @test t["processes"] == 1
            @test t["samples"] == 0 || t["peak_rss_kib"] > 0
        end

        vsearch_bin = Sys.which("vsearch")
        if isnothing(vsearch_bin)
            @info "vsearch not on PATH - the two-thread disagreement test is SKIPPED"
        else
            @testset "two-thread disagreement: vsearch ASV and OTU in three orders" begin
                mktempdir() do tmp
                    dir = occupancy_fixture(joinpath(tmp, "fixture"))
                    # Warm-up: compile every path once, so no order pays for it.
                    for order in ORDERS
                        run_order(dir, joinpath(tmp, "warmup", order), order)
                    end
                    GC.gc()
                    records = OrderedDict{String,Any}()
                    for order in ORDERS
                        GC.gc()
                        records[order] = run_order(dir, joinpath(tmp, order), order)
                    end

                    # The same work: identical hits in every order.
                    for name in ("asv", "otu")
                        hits = [read(joinpath(tmp, o, name, "taxonomy.tsv"), String) for o in ORDERS]
                        @test !isempty(first(hits))
                        @test allequal(hits)
                    end

                    # Cost commutes (asserted).
                    cost = [records[o]["julia_allocated_bytes"] for o in ORDERS]
                    spread = (maximum(cost) - minimum(cost)) / maximum(cost)
                    @test spread <= COST_TOLERANCE

                    # Peaks (recorded, not asserted).
                    for o in ORDERS
                        r = records[o]
                        a, u = _task(r, "vsearch:asv"), _task(r, "vsearch:otu")
                        @info "occupancy" order = o julia_allocated_bytes = r["julia_allocated_bytes"] julia_live_bytes_delta = r["julia_live_bytes_delta"] child_peak_rss_kib = r["child_peak_rss_kib"] asv_peak_rss_kib = a["peak_rss_kib"] asv_hwm_kib = a["hwm_kib"] otu_peak_rss_kib = u["peak_rss_kib"] otu_hwm_kib = u["hwm_kib"] samples = r["samples"]
                    end
                    conc = records["concurrent"]
                    serial_max = maximum(records[o]["child_peak_rss_kib"] for o in ORDERS[1:2])
                    overlap = conc["child_peak_rss_kib"] > serial_max
                    @info "occupancy: concurrent stage peak vs the larger serial stage peak" concurrent = conc["child_peak_rss_kib"] serial = serial_max overlapped = overlap cost_spread = round(spread; sigdigits = 3)
                    overlap || @warn "occupancy: the concurrent searches did not overlap in memory on this runner, so this run cannot show the order dependence"

                    path = get(ENV, "OCCUPANCY_RESULTS", "")
                    if !isempty(path)
                        io = IOBuffer()
                        run(pipeline(`$vsearch_bin --version`; stdout = io, stderr = io))
                        version = first(split(String(take!(io)), '\n'))
                        YAML.write_file(path, OrderedDict{String,Any}(
                        "vsearch" => version,
                        "julia" => string(VERSION), "threads" => Threads.nthreads(),
                        "cpus" => Sys.CPU_THREADS, "cost_spread" => spread,
                        "orders" => records))
                    end
                end
            end
        end
    end
end
