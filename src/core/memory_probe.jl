# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

module MemoryProbe

## Occupancy probe
# Measures what an attested stage occupied, for the stage's Attestation record:
#
#   * Julia-side allocation, as the `Base.gc_num()` delta over the stage. This is
#     a cost: it only ever adds up, whatever order the stage's work ran in.
#   * Each child process's resident memory, sampled from /proc. The tools run as
#     child processes (vsearch, swarm, cd-hit, ...), so `Sys.maxrss()` of this
#     process never sees them.
#
# Two figures come from the samples. A task's `peak_rss_kib` is the largest VmRSS
# sampled for the processes that task started. The stage's `child_peak_rss_kib` is
# the largest SUM of VmRSS over every process of the stage at one sampling
# instant. For two tasks run one after the other it is the larger of their peaks;
# for two that run concurrently it can approach their sum. That order dependence,
# invisible in every output file, is what the measurement exists to show.
#
# `hwm_kib` is the kernel's own high-water mark (VmHWM) for a task's processes, as
# last read before they exited. It is exact up to the last sample, so it is at
# least `peak_rss_kib`; growth after the last sample is not seen. A process that
# lives for less than one interval may not be sampled at all, which `samples: 0`
# reports rather than a misleading zero peak.
#
# The probe is off unless the run config sets `occupancy_probe.enabled: true`.
# Off, `start` returns `nothing` and every other entry point is a pass-through,
# so the stage records and the outputs are exactly what they were without it.
#
# Linux only (/proc). Elsewhere, and for stages offloaded with `remote:` (where
# the local child is an ssh client), the child figures do not describe the work.
# The allocation figure is process-wide: while the study loop runs several runs'
# pipelines at once, a stage's figure includes the others' allocations too. The
# child figures are not affected, since each stage samples only its own children.

using Base.ScopedValues: ScopedValue, with
using OrderedCollections: OrderedDict

const PROBE = "MemoryProbe 1"
const DEFAULT_INTERVAL_MS = 10

# One labelled unit of a stage's work (one spawned search, say) and the child
# processes it started.
mutable struct TaskTrace
    label     :: String
    live      :: Vector{Int}   # pids started and not yet waited for
    processes :: Int           # every child started under this label
    samples   :: Int           # samples that found at least one live process
    peak_kib  :: Int
    hwm_kib   :: Int
    # The sample window the task's children were live in: the index of the first
    # sample taken after the first of them started, and the last taken before
    # the last was waited for. 0 means never. Two tasks overlap when the later
    # one's `first_sample` falls at or before the earlier one's `last_sample`,
    # which is the interleaving the grade algebra sequences over.
    first_sample :: Int
    last_sample  :: Int
end
TaskTrace(label::AbstractString) =
    TaskTrace(String(label), Int[], 0, 0, 0, 0, 0, 0)

mutable struct Probe
    interval    :: Float64
    lock        :: ReentrantLock
    tasks       :: Vector{TaskTrace}
    stage       :: TaskTrace           # children started outside any `task`
    samples     :: Int
    peak_kib    :: Int
    gc_before   :: Base.GC_Num
    live_before :: Int
    running     :: Bool
    sampler     :: Union{Task,Nothing}
    timer       :: Union{Timer,Nothing}
    record      :: Union{OrderedDict{String,Any},Nothing}
    error       :: Union{String,Nothing}   # why sampling stopped early, if it did
    buf         :: Vector{UInt8}       # /proc reads; preallocated so a sample
    path        :: Vector{UInt8}       # allocates nothing that would count as work
end

const STAGE = ScopedValue{Union{Probe,Nothing}}(nothing)
const TASK  = ScopedValue{Union{TaskTrace,Nothing}}(nothing)

## Configuration
_section(cfg) = cfg isa AbstractDict ? get(cfg, "occupancy_probe", nothing) : nothing

"""
    enabled(run_cfg) -> Bool

Whether the run config turns the probe on: `occupancy_probe.enabled: true`. Any
other value, or no `occupancy_probe` section at all, leaves it off.
"""
function enabled(cfg)
    s = _section(cfg)
    s isa AbstractDict && get(s, "enabled", false) === true
end

"""
    interval_seconds(run_cfg) -> Float64

The sampling interval, from `occupancy_probe.interval_ms` (default 10 ms, at
least 1 ms).
"""
function interval_seconds(cfg)
    s = _section(cfg)
    ms = s isa AbstractDict ? get(s, "interval_ms", DEFAULT_INTERVAL_MS) : DEFAULT_INTERVAL_MS
    (ms isa Real && !(ms isa Bool) && ms >= 1) || (ms = DEFAULT_INTERVAL_MS)
    Float64(ms) / 1000
end

## Reading /proc without allocating
# "/proc/<pid>/status" into `path`, NUL-terminated.
function _status_path!(path::Vector{UInt8}, pid::Int)
    k = 0
    for c in codeunits("/proc/"); k += 1; path[k] = c; end
    nd = ndigits(pid); d = pid
    for i in nd:-1:1
        path[k + i] = UInt8('0') + UInt8(d % 10); d ÷= 10
    end
    k += nd
    for c in codeunits("/status"); k += 1; path[k] = c; end
    path[k + 1] = 0x00
    nothing
end

# Bytes of /proc/<pid>/status read into `buf`, or -1 if the process is gone.
function _read_status!(buf::Vector{UInt8}, path::Vector{UInt8}, pid::Int)
    _status_path!(path, pid)
    fd = GC.@preserve path ccall(:open, Cint, (Ptr{UInt8}, Cint), path, 0)
    fd < 0 && return -1
    total = 0
    GC.@preserve buf while total < length(buf)
        r = ccall(:read, Cssize_t, (Cint, Ptr{UInt8}, Csize_t),
                  fd, pointer(buf, total + 1), length(buf) - total)
        r <= 0 && break
        total += r
    end
    ccall(:close, Cint, (Cint,), fd)
    total
end

# The kB value of the line starting `key` (6 bytes, e.g. "VmRSS:"), or -1.
function _field_kib(buf::Vector{UInt8}, n::Int, key::NTuple{6,UInt8})
    i = 1
    while i <= n - 6
        if (i == 1 || buf[i - 1] == UInt8('\n')) &&
           buf[i] == key[1] && buf[i + 1] == key[2] && buf[i + 2] == key[3] &&
           buf[i + 3] == key[4] && buf[i + 4] == key[5] && buf[i + 5] == key[6]
            j = i + 6
            while j <= n && (buf[j] == UInt8(' ') || buf[j] == UInt8('\t')); j += 1; end
            v = 0
            while j <= n && UInt8('0') <= buf[j] <= UInt8('9')
                v = 10v + Int(buf[j] - UInt8('0')); j += 1
            end
            return v
        end
        i += 1
    end
    -1
end

const VMRSS = Tuple(codeunits("VmRSS:"))
const VMHWM = Tuple(codeunits("VmHWM:"))

## Sampling
# One sample: each live child's VmRSS and VmHWM, summed per task and over the stage.
function _sample!(p::Probe)
    lock(p.lock) do
        total = 0
        for t in p.tasks
            total += _sample_task!(p, t)
        end
        total += _sample_task!(p, p.stage)
        p.samples += 1
        p.peak_kib = max(p.peak_kib, total)
    end
    nothing
end

function _sample_task!(p::Probe, t::TaskTrace)
    isempty(t.live) && return 0
    rss = 0; seen = false
    for pid in t.live
        n = _read_status!(p.buf, p.path, pid)
        n > 0 || continue
        r = _field_kib(p.buf, n, VMRSS)
        r < 0 && continue            # a zombie has no memory lines
        seen = true
        rss += r
        t.hwm_kib = max(t.hwm_kib, _field_kib(p.buf, n, VMHWM))
    end
    seen || return 0
    t.samples += 1
    t.peak_kib = max(t.peak_kib, rss)
    rss
end

"""
    start(run_cfg) -> Union{Probe,Nothing}

A running probe for one stage, or `nothing` when the run config does not enable
it. Pass the result to `within` and then `finish!`.
"""
function start(cfg)
    enabled(cfg) || return nothing
    p = Probe(interval_seconds(cfg), ReentrantLock(), TaskTrace[], TaskTrace("(stage)"),
              0, 0, Base.gc_num(), Base.gc_live_bytes(), true, nothing, nothing, nothing,
              nothing, Vector{UInt8}(undef, 8192), Vector{UInt8}(undef, 64))
    p.timer = Timer(0; interval = p.interval)
    timer = p.timer
    # A probe must never fail the stage it measures: a sampling error ends the
    # sampling and is reported in the record instead.
    p.sampler = Threads.@spawn begin
        try
            while p.running
                _sample!(p)
                isopen(timer) || break
                wait(timer)
            end
        catch err
            (err isa EOFError && !p.running) || (p.error = sprint(showerror, err))
        end
    end
    p
end

"""
    within(f, probe)

Run `f()` with `probe` as the current stage's probe, so the children it starts
(directly or from tasks it spawns) are sampled. With `probe === nothing` this is
just `f()`.
"""
within(f, ::Nothing) = f()
within(f, p::Probe) = with(f, STAGE => p)

"""
    task(f, label)

Run `f()` as one labelled unit of the current stage. Children it starts are
reported under `label`. Outside a probed stage this is just `f()`.
"""
function task(f, label::AbstractString)
    p = STAGE[]
    p === nothing && return f()
    t = TaskTrace(label)
    lock(() -> push!(p.tasks, t), p.lock)
    with(f, TASK => t)
end

# The trace a child started now belongs to, or nothing when nothing is probed.
function _current()
    p = STAGE[]
    p === nothing && return nothing, nothing
    t = TASK[]
    p, something(t, p.stage)
end

"""
    watch!(process) -> Union{Int,Nothing}

Register a just-started child process with the current stage's probe, sample it
at once, and return its pid for `unwatch!`. Called by `Tools._run_killable`.
Returns `nothing`, having done nothing, when the stage is not probed. A child that
has already exited is counted but cannot be sampled. Never throws: a probe must
not fail the stage it measures.
"""
function watch!(proc::Base.Process)
    p, t = _current()
    p === nothing && return nothing
    # `getpid` throws once the child has exited, so the pid is read now or never.
    pid = try
        Int(getpid(proc))
    catch
        nothing
    end
    lock(p.lock) do
        t.processes += 1
        t.processes == 1 && (t.first_sample = p.samples + 1)
        if pid !== nothing
            push!(t.live, pid)
        end
    end
    pid === nothing || _sample!(p)
    pid
end

# A pipeline of commands: each of its processes, as one task's children.
function watch!(chain::Base.ProcessChain)
    pids = Int[]
    for proc in chain.processes
        pid = watch!(proc)
        pid === nothing || push!(pids, pid)
    end
    isempty(pids) ? nothing : pids
end

# Anything else `run` can return is not sampled.
watch!(_) = nothing

"""
    unwatch!(pid)

Stop sampling a child that has been waited for (the pid `watch!` returned), so
the pid cannot be mistaken for a later, unrelated process. A no-op for `nothing`.
"""
unwatch!(::Nothing) = nothing
unwatch!(pids::Vector{Int}) = foreach(unwatch!, pids)
function unwatch!(pid::Int)
    p, t = _current()
    p === nothing && return nothing
    lock(p.lock) do
        filter!(!=(pid), t.live)
        isempty(t.live) && (t.last_sample = p.samples)
    end
    nothing
end

## Grading a stage against the high-water-mark algebra
# `docs/occupancy` reads a stage as a timeline of grades and checks the peak the
# algebra predicts against the one the probe measured. The grade is the pinned
# HWM instance's (`occupancy-types`' sibling `resource-algebra`,
# `ResourceAlgebra.Instances.HWM`; both are fetched by `scripts/fetch-siblings.sh`):
#
#     hwm peak drift        with   drift ≤ peak
#
#   * `peak`  — the high-water mark this interaction reached, in KiB;
#   * `drift` — the change in the live peak sum across the interaction's
#     boundary, so a boundary that also releases an earlier task is negative;
#   * `release` = peak − drift — everything freed between the interaction's start
#     and the next boundary, which can exceed the interaction's own peak.
#
# A task's drift subtracts earlier tasks that exited before the next start, so it
# can be negative — the signed-drift extension `Instances/HWM`'s header records
# as H-2. The folded drift is the sum of peaks still live at the stage's end.
#
# A task's drift subtracts earlier tasks that exited before the next start, so
# it can be negative. The folded drift is the sum of peaks still live there.
#
# Sequencing starts the next interaction from the previous one's drift:
#
#     hwm p₁ n₁ ⊗ hwm p₂ n₂ = hwm (max p₁ (n₁ + p₂)) (n₁ + n₂)
#
# and `OccupancyTypes.Reclaim` proves of that definition
#
#     peak (a ⊗ b) = max (peak a) (drift a + peak b)           (reclaim-law)
#     peak (a ⊗ b) ≤ peak a + peak b                           (occupancy≤cost)
#     drift a = peak a  ⇒  peak (a ⊗ c) = peak a + peak c       (no-reclaim)
#     drift a < peak a, 0 < peak c  ⇒  peak (a ⊗ c) < peak a + peak c (reclaimed)
#
# The last two are the two shapes a stage actually has. Run together, the first
# task is still resident when the second starts (drift = peak) and the stage peak
# is the sum; run one after the other it has exited (drift = 0) and the stage peak
# is the larger. The samples already tell the two apart, so the drift is read from
# the task windows and never assumed.
#
# The prediction is an upper bound on the measurement: at a handover the task
# still live holds at most its peak, so the drift the algebra adds is at least
# what a sample could have seen. `disagreement` therefore reports only what that
# bound forbids — a measured stage peak above the prediction — which is the "one
# of them is a bug" the occupancy reading asks about.

const GRADE_ALGEBRA = "ResourceAlgebra.Instances.HWM"

"A grade of the pinned high-water-mark algebra: a peak and a net memory change."
struct Grade
    peak  :: Int
    drift :: Int
    function Grade(peak::Integer, drift::Integer)
        drift <= peak || throw(ArgumentError(
            "a grade holds at most what it reached: drift $drift > peak $peak"))
        new(Int(peak), Int(drift))
    end
end

"`release(grade)` — `peak − drift`: everything freed by the interaction's boundary."
release(g::Grade) = g.peak - g.drift

"`seq(a, b)` — sequence `b` after `a`: the suffix starts from `a`'s drift."
seq(a::Grade, b::Grade) = Grade(max(a.peak, a.drift + b.peak), a.drift + b.drift)

"`alt(a, b)` — the larger grade, componentwise (what the timeline never needs, kept for the pair of laws)."
alt(a::Grade, b::Grade) = Grade(max(a.peak, b.peak), max(a.drift, b.drift))

"Nothing allocated and nothing held: the algebra's unit."
const NO_GRADE = Grade(0, 0)

# The trace's sample window, or `nothing` for a record that carries none (one
# written before the windows were sampled).
function _window(rec::AbstractDict)
    first_seen = get(rec, "first_sample", nothing)
    last_seen  = get(rec, "last_sample", nothing)
    (first_seen isa Integer && last_seen isa Integer) ? (Int(first_seen), Int(last_seen)) : nothing
end

_is_windowed(rec::AbstractDict) = _window(rec) !== nothing

# Whether `prev` was still live when `next` started, so the handover hands its
# peak on rather than reclaiming it. A started `prev` never unwatched (window
# end 0) was still live; a task that never started (window start 0) never overlapped.
_overlaps(prev::Tuple{Int,Int}, next::Tuple{Int,Int}) =
    prev[1] != 0 && next[1] != 0 && (prev[2] == 0 || next[1] <= prev[2])

# Grades and labels in sampled order, or `nothing` when the record carries no
# windows: a comparison is then skipped rather than made on an assumed order.
function _grades_and_labels(occ::AbstractDict)
    recs = get(occ, "tasks", nothing)
    recs isa AbstractVector || return Grade[], String[]
    all(_is_windowed, recs) || return nothing
    # Sampled order, ties broken by record order so the result never depends on
    # the sort's stability.
    order  = sort(collect(eachindex(recs)); by = i -> (_window(recs[i])[1], i))
    ordered = collect(recs)[order]
    windows = [_window(r) for r in ordered]
    peaks   = [Int(get(r, "peak_rss_kib", 0)) for r in ordered]
    grades  = Grade[]
    labels  = String[]
    live    = Int[]
    held    = 0
    for (i, rec) in enumerate(ordered)
        push!(live, i)
        if i == length(ordered)
            # The stage waits for its children: a task still resident at the end
            # keeps its peak in the folded drift, an exited one holds nothing.
            filter!(j -> windows[j][1] != 0 && windows[j][2] == 0, live)
        else
            filter!(j -> _overlaps(windows[j], windows[i + 1]), live)
        end
        remaining = sum(j -> peaks[j], live; init = 0)
        # seq adds drifts, so reclaim expired predecessors as well as this task.
        push!(grades, Grade(peaks[i], remaining - held))
        held = remaining
        push!(labels, String(get(rec, "label", "?")))
    end
    grades, labels
end

"""
    task_grades(occupancy) -> Union{Vector{Grade},Nothing}

One grade per task of a stage's `occupancy` record, in sampled order. Its drift
is the change in the live peak sum by the next task's start, including releases
of earlier tasks; the final boundary retains only unclosed windows. Folding
with `seq` carries exactly the peaks still live at each boundary. `nothing` when
the record carries no sample windows.
"""
function task_grades(occ::AbstractDict)
    graded = _grades_and_labels(occ)
    graded === nothing ? nothing : graded[1]
end

"""
    stage_grade(occupancy) -> Union{Grade,Nothing}

The stage's grade: its tasks' grades folded with `seq`, the algebra's sequencing.
`nothing` when the record carries no sample windows.
"""
function stage_grade(occ::AbstractDict)
    grades = task_grades(occ)
    grades === nothing && return nothing
    isempty(grades) ? NO_GRADE : foldl(seq, grades)
end

"`predicted_peak_rss_kib(occupancy)` — the peak the algebra predicts, or `nothing`."
function predicted_peak_rss_kib(occ::AbstractDict)
    g = stage_grade(occ)
    g === nothing ? nothing : g.peak
end

"`cost_bound_rss_kib(occupancy)` — the peaks added: the peak if no cut reclaimed anything (occupancy≤cost's right-hand side)."
function cost_bound_rss_kib(occ::AbstractDict)
    grades = task_grades(occ)
    grades === nothing ? nothing : sum(g -> g.peak, grades; init = 0)
end

"`reclaimed_kib(occupancy)` — the slack between the cost bound and the predicted peak: what the cuts reclaim."
function reclaimed_kib(occ::AbstractDict)
    cost = cost_bound_rss_kib(occ)
    predicted = predicted_peak_rss_kib(occ)
    (cost === nothing || predicted === nothing) ? nothing : cost - predicted
end

"""
    disagreement(occupancy) -> Union{String,Nothing}

Why the measured stage peak contradicts the grade algebra, or `nothing` when it
does not. The prediction is an upper bound, so the one contradiction is a
measured peak above it: a sample then saw more resident memory than the timed
timeline of tasks can hold, and one of the two is wrong. A record with no
windows, or none the probe graded, gives `nothing`: not measured is not a
disagreement.
"""
function disagreement(occ::AbstractDict)
    predicted = predicted_peak_rss_kib(occ)
    measured  = get(occ, "child_peak_rss_kib", nothing)
    (predicted === nothing || !(measured isa Integer)) && return nothing
    measured > predicted && return "measured $(measured) KiB above the predicted $(predicted) KiB"
    nothing
end

"""
    graded(occupancy) -> Union{OrderedDict,Nothing}

The grades of a stage's `occupancy` record, ready to sit beside it in the
attestation: per task `peak_kib`/`drift_kib`/`release_kib`, the stage's
`predicted_peak_kib`, the `cost_bound_kib` (its peaks added) and the
`reclaimed_kib` slack between the two, the `measured_peak_kib` it is checked
against, and the `disagreement` between them or `nothing`. `nothing` when the
record carries no sample windows.
"""
function graded(occ::AbstractDict)
    graded_pair = _grades_and_labels(occ)
    graded_pair === nothing && return nothing
    grades, labels = graded_pair
    stage = isempty(grades) ? NO_GRADE : foldl(seq, grades)
    cost  = sum(g -> g.peak, grades; init = 0)
    OrderedDict{String,Any}(
        "algebra"          => GRADE_ALGEBRA,
        "tasks"            => Any[OrderedDict{String,Any}(
                                      "label"   => labels[i],
                                      "peak_kib"    => grades[i].peak,
                                      "drift_kib"   => grades[i].drift,
                                      "release_kib" => release(grades[i]))
                                  for i in eachindex(grades)],
        "predicted_peak_kib" => stage.peak,
        "drift_kib"          => stage.drift,
        "release_kib"        => release(stage),
        "cost_bound_kib"     => cost,
        "reclaimed_kib"      => cost - stage.peak,
        "measured_peak_kib"  => get(occ, "child_peak_rss_kib", nothing),
        "disagreement"       => disagreement(occ),
    )
end

_task_record(t::TaskTrace) = OrderedDict{String,Any}(
    "label"        => t.label,
    "processes"    => t.processes,
    "samples"      => t.samples,
    "peak_rss_kib" => t.peak_kib,
    "hwm_kib"      => t.hwm_kib,
    "first_sample" => t.first_sample,
    "last_sample"  => t.last_sample,
)

"""
    finish!(probe) -> Union{OrderedDict,Nothing}

Stop sampling and return the stage's `occupancy` record (the same record on every
later call). `nothing` for `nothing`.
"""
finish!(::Nothing) = nothing
function finish!(p::Probe)
    p.record === nothing || return p.record
    p.running = false
    close(p.timer)
    wait(p.sampler)
    allocated = Base.GC_Diff(Base.gc_num(), p.gc_before).allocd
    traces = [p.tasks; p.stage.processes > 0 ? [p.stage] : TaskTrace[]]
    p.record = OrderedDict{String,Any}(
        "probe"                  => PROBE,
        "interval_ms"            => round(Int, p.interval * 1000),
        "samples"                => p.samples,
        "julia_allocated_bytes"  => allocated,
        "julia_live_bytes_delta" => Base.gc_live_bytes() - p.live_before,
        "child_peak_rss_kib"     => p.peak_kib,
        "tasks"                  => Any[_task_record(t) for t in traces],
    )
    p.error === nothing || (p.record["error"] = p.error)
    # The grades the sample windows imply, beside the samples they came from. A
    # record whose tasks carry no windows cannot be graded, and is left as it is
    # rather than graded on an assumed order.
    g = graded(p.record)
    g === nothing || (p.record["grades"] = g)
    p.record
end

end # module MemoryProbe
