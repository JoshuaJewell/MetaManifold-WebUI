# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Gate for the graded reading of the recorded run (fork issue #105): rebuild each
# order's timeline from its sample window (docs/occupancy/schedules.yml), grade it
# with the module the probe itself uses (src/core/memory_probe.jl — included
# rather than the whole package, so this needs OrderedCollections and YAML only),
# and compare the algebra's predicted peak with the measured one.
#
# The recorded run is fixed evidence, so the comparison is an equality: this is
# the reading docs/occupancy-2026-10.adoc publishes, and a change that moves it
# must be a deliberate one. Two further claims hold of every order and are
# checked: the prediction is never above the cost bound (occupancy≤cost) and the
# slack between them is never negative (what the cuts reclaim).
#
# Run from the repository root:  julia --project=. scripts/check-occupancy-grades.jl
include(joinpath(@__DIR__, "..", "src", "core", "memory_probe.jl"))
using .MemoryProbe
using YAML

const ROOT      = normpath(joinpath(@__DIR__, ".."))
const RESULTS   = joinpath(ROOT, "docs", "occupancy", "results-2026-10-14.yml")
const SCHEDULES = joinpath(ROOT, "docs", "occupancy", "schedules.yml")

orders    = YAML.load_file(RESULTS)["orders"]
schedules = YAML.load_file(SCHEDULES)

failures = 0
function fail(message)
    global failures += 1
    println("FAIL: ", message)
end

println(rpad("order", 16), rpad("predicted", 11), rpad("measured", 11),
        rpad("cost", 8), "reclaimed")
# An order's timeline, from the peaks the probe recorded and the window each task
# ran in. Nothing else is taken from the schedule: the grades are the
# measurement's.
function timeline(rec, windows)
    tasks = Any[Dict{String,Any}("label"        => t["label"],
                                 "peak_rss_kib" => t["peak_rss_kib"],
                                 "first_sample" => windows[i][1],
                                 "last_sample"  => windows[i][2])
                for (i, t) in enumerate(rec["tasks"])]
    Dict{String,Any}("child_peak_rss_kib" => rec["child_peak_rss_kib"], "tasks" => tasks)
end

for name in sort(collect(keys(orders)))
    rec       = orders[name]
    occ       = timeline(rec, schedules[name])
    predicted = MemoryProbe.predicted_peak_rss_kib(occ)
    cost      = MemoryProbe.cost_bound_rss_kib(occ)
    reclaimed = MemoryProbe.reclaimed_kib(occ)
    measured  = rec["child_peak_rss_kib"]
    println(rpad(name, 16), rpad(string(predicted), 11), rpad(string(measured), 11),
            rpad(string(cost), 8), string(reclaimed))

    predicted == measured ||
        fail("$name: predicted $predicted KiB != measured $measured KiB")
    predicted <= cost ||
        fail("$name: predicted $predicted KiB above the cost bound $cost KiB (occupancy≤cost)")
    reclaimed >= 0 ||
        fail("$name: reclaimed $reclaimed KiB is negative")
    d = MemoryProbe.disagreement(occ)
    d === nothing || fail("$name: $d")
end

if failures > 0
    println("check-occupancy-grades: $failures failure(s)")
    exit(1)
end
println("check-occupancy-grades: OK — the prediction reproduces the recorded measurement")
