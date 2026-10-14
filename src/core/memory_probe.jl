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
end
TaskTrace(label::AbstractString) = TaskTrace(String(label), Int[], 0, 0, 0, 0)

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
        pid === nothing || push!(t.live, pid)
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
    end
    nothing
end

_task_record(t::TaskTrace) = OrderedDict{String,Any}(
    "label"        => t.label,
    "processes"    => t.processes,
    "samples"      => t.samples,
    "peak_rss_kib" => t.peak_kib,
    "hwm_kib"      => t.hwm_kib,
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
    p.record
end

end # module MemoryProbe
