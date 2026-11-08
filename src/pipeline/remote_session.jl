# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

# Included inside RemoteExec. The journal belongs to a workflow's stage log,
# not to its content digest: two workflows with identical inputs stay separate.
using SHA, JSON3, UUIDs

_sq(s::AbstractString) = "'" * replace(s, "'" => "'\\''") * "'"

const _CONTROL = ".metamanifold"
const _WORKER = joinpath(@__DIR__, "remote_worker.sh")

_digest(path) = open(sha256, path) |> bytes2hex

function _relative(path)
    s = string(path)
    parts = split(s, '/')
    (isempty(s) || occursin('\\', s) || any(p -> p in ("", ".", ".."), parts) ||
     first(parts) == _CONTROL || any(iscntrl, s)) &&
        error("Remote path must stay inside its staging directory: $s")
    s
end

function _identity(target, stage, command_for, inputs, outputs, optional_outputs, context)
    uploads = sort!([(_relative(rel), _digest(path)) for (path, rel) in inputs])
    length(unique(first.(uploads))) == length(uploads) || error("Duplicate remote input path")
    destinations = sort!([_relative(rel) for (rel, _) in [outputs; optional_outputs]])
    length(unique(destinations)) == length(destinations) || error("Duplicate remote output path")
    # Arrays have an unambiguous, ordered encoding; neither path concatenation
    # nor Dict iteration order is used as an identity. Callers include external
    # (not uploaded) data checksums in context.
    JSON3.write([1, target.host, target.base, stage, target.threads,
                 command_for("/METAMANIFOLD_STAGE"), _digest(_WORKER), uploads,
                 sort!([_relative(rel) for (rel, _) in outputs]),
                 sort!([_relative(rel) for (rel, _) in optional_outputs]), context])
end

function _save_journal(path, journal)
    tmp, io = mktemp(dirname(path))
    try
        write(io, JSON3.write(journal))
        flush(io)
        Sys.isunix() && ccall(:fsync, Cint, (Cint,), fd(io)) != 0 &&
            error("Could not sync remote execution journal")
        close(io)
        # mv(...; force=true) can unlink its destination first. A direct rename
        # within this directory gives readers either complete journal version.
        Base.Filesystem.rename(tmp, path)
        if Sys.isunix()
            parent = ccall(:open, Cint, (Cstring, Cint), dirname(path), 0)
            parent >= 0 || error("Could not open remote journal directory for sync")
            try
                ccall(:fsync, Cint, (Cint,), parent) == 0 || error("Could not sync remote journal directory")
            finally
                ccall(:close, Cint, (Cint,), parent)
            end
        end
    finally
        isopen(io) && close(io)
        isfile(tmp) && rm(tmp; force=true)
    end
end

# flock is released by the kernel even after SIGKILL or a local reboot. Keep the
# inode: unlinking a lock file lets another process lock a different inode at the
# same path. MetaManifold runs on Unix (including WSL); remote workers are Linux.
function _with_owner(f, path; on_interrupt=nothing)
    Sys.isunix() || error("Remote execution requires Unix or WSL")
    open(path, "a") do io
        ccall(:flock, Cint, (Cint, Cint), fd(io), 2 | 4) == 0 ||
            error("This workflow's remote stage already has an active caller")
        try
            f()
        catch e
            e isa InterruptException && !isnothing(on_interrupt) && on_interrupt()
            rethrow()
        finally
            ccall(:flock, Cint, (Cint, Cint), fd(io), 8)
        end
    end
end

# Every foreground transport is interruptible, including uploads and downloads.
# Killing one of these SSH processes never signals the detached workload.
function _capture(cmd)
    out = Pipe()
    p = run(pipeline(cmd; stdout=out); wait=false)
    close(out.in)
    try
        bytes = read(out)
        wait(p)
        success(p) || throw(ProcessFailedException(p))
        bytes
    catch
        process_running(p) && kill(p, Base.SIGKILL)
        wait(p)
        rethrow()
    finally
        close(out)
    end
end

function _transport(target, ctl)
    id_opt = isnothing(target.identity_file) ? `` : `-i $(target.identity_file)`
    opts = `$id_opt -o ControlMaster=auto -o ControlPath=$ctl -o ControlPersist=60 -o ConnectTimeout=15 -o NumberOfPasswordPrompts=1 -o ServerAliveInterval=10 -o ServerAliveCountMax=3`
    ((args...) -> `ssh $opts $args`, (args...) -> `scp $opts $args`)
end

function _remote_state(ssh, target, dir)
    worker = _sq("$dir/$_CONTROL/worker.sh")
    strip(String(_capture(ssh(target.host,
        "if test -f $worker; then bash $worker status; else echo missing; fi"))))
end

function _control(ssh, target, dir, action)
    _capture(ssh(target.host, "bash $(_sq("$dir/$_CONTROL/worker.sh")) $action"))
end

# `reclaim` distinguishes an abandoned attempt from a requested stop. A stop
# goes through the worker, which fences a launch, kills the process group and
# removes the directory only once the execution lock is free. An attempt whose
# worker was never installed never launched, and one that already published a
# terminal marker has no live command, so both are reclaimed directly: the
# directory is journal-owned, and no wildcard or user-supplied path is removed.
function _cancel_remote(ssh, target, dir; reclaim::Bool=false)
    state = _remote_state(ssh, target, dir)
    if state == "missing" || (reclaim && !(state in ("running", "stopping")))
        _capture(ssh(target.host, "rm -rf -- $(_sq(dir))"))
    else
        _control(ssh, target, dir, "cancel")
        _control(ssh, target, dir, "cleanup")
    end
end

function _new_journal(identity, target)
    key = bytes2hex(sha256(identity))
    Dict{String,Any}("identity" => identity,
                    "directory" => joinpath(target.base, "mm_$(key)_$(uuid4())"),
                    "identity_file" => get(target, :identity_file, nothing),
                    "failed" => false, "cancel_requested" => false)
end

function _check_directory(journal, target)
    dir = string(journal["directory"])
    parent, name = dirname(dir), basename(dir)
    normpath(parent) == normpath(target.base) &&
        occursin(r"^mm_[0-9a-f]{64}_[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", name) ||
        error("Invalid directory in remote execution journal")
    dir
end

# Cancellation covers the whole owner-locked operation, including hashing and
# command construction before a reconnect. Read the durable journal so the stop
# always addresses the attempt's original host/key, even after settings change.
# `connection` is the transport of the operation being interrupted, if there is
# one: an interactive password prompt cannot be answered from an interrupt
# handler, so a live master is reused and only a transport failure (typically a
# ControlPath whose master has already gone) falls back to a fresh connection.
function _request_cancellation(journal_path, connection)
    isfile(journal_path) || return
    journal = JSON3.read(read(journal_path, String), Dict{String,Any})
    identity = JSON3.read(journal["identity"])
    target = (; host=String(identity[2]), base=String(identity[3]), threads=identity[5],
        identity_file=get(journal, "identity_file", nothing))
    dir = _check_directory(journal, target)
    journal["cancel_requested"] = true
    _save_journal(journal_path, journal)
    ctl = joinpath(tempdir(), "mm_$(uuid4())")
    own = isnothing(connection)
    ssh, _ = own ? _transport(target, ctl) : connection
    try
        try
            _cancel_remote(ssh, target, dir)
        catch
            # A master that is already gone leaves the attempt running; try
            # once more over a transport this handler opens and owns.
            own && rethrow()
            own = true
            ssh, _ = _transport(target, ctl)
            _cancel_remote(ssh, target, dir)
        end
        rm(journal_path; force=true)
    catch cleanup_error
        @warn "Remote cancellation not acknowledged; retry the stage to finish cancelling" host=target.host directory=dir exception=cleanup_error
    finally
        if own
            try _capture(`ssh -o ControlPath=$ctl -O exit $(target.host)`) catch end
        end
    end
end

function _prepare(ssh, scp, target, dir, command, identity, inputs, outputs, optional_outputs, before_run, emit)
    subdirs = unique([dirname(_relative(rel)) for (_, rel) in inputs])
    append!(subdirs, [dirname(_relative(rel)) for (rel, _) in [outputs; optional_outputs]])
    paths = join(_sq.(unique(["$dir/$_CONTROL"; ["$dir/$d" for d in subdirs]])), " ")
    _capture(ssh(target.host, "umask 077; mkdir -p -- $paths"))
    # worker.sh is installed last: its presence makes the status probe usable.
    emit("  Transferring $(length(inputs)) input file(s) to $(target.host)")
    for (path, rel) in inputs
        _capture(scp(path, "$(target.host):$dir/$rel"))
    end
    # Check bytes, not just sizes. A local input changed during transfer is a
    # failed upload, never a command launched under the old content identity.
    uploads = JSON3.read(identity)[8]
    sums = join(["$(entry[2])  $(entry[1])\n" for entry in uploads])
    isempty(uploads) || _capture(ssh(target.host,
        "cd $(_sq(dir)) && printf %s $(_sq(sums)) | sha256sum -c --strict -"))
    isnothing(before_run) || before_run(ssh, dir)
    # Unique temporary names keep an interrupted transfer from exposing or
    # later modifying the installed control scripts. Publish only checked bytes.
    command_part = "$_CONTROL/command.sh.part-$(uuid4())"
    worker_part = "$_CONTROL/worker.sh.part-$(uuid4())"
    mktemp() do path, io
        write(io, command * "\n")
        close(io)
        _capture(scp(path, "$(target.host):$dir/$command_part"))
    end
    _capture(scp(_WORKER, "$(target.host):$dir/$worker_part"))
    control_sums = "$(JSON3.read(identity)[7])  $worker_part\n" *
                   "$(bytes2hex(sha256(command * "\n")))  $command_part\n"
    _capture(ssh(target.host,
        "command -v bash flock setsid nohup sha256sum dd >/dev/null && " *
        "test -r /proc/sys/kernel/random/boot_id && cd $(_sq(dir)) && " *
        "printf %s $(_sq(control_sums)) | sha256sum -c --strict - && " *
        "mv -- $(_sq(command_part)) $_CONTROL/command.sh && " *
        "mv -- $(_sq(worker_part)) $_CONTROL/worker.sh && " *
        "touch $(_sq("$dir/$_CONTROL/ready"))"))
end

# A zero-exit attempt can still have produced nothing: a stage that left a
# required file unwritten must be diagnosable and retrievable, not a fetch that
# fails again on every run. The check asks the same question as the fetch, so a
# required output is reported before any transfer starts.
function _missing_outputs(ssh, target, dir, outputs)
    missing = String[]
    for (rel, _) in outputs
        remote_path = _sq("$dir/$rel")
        present = strip(String(_capture(ssh(target.host,
            "if test -f $remote_path; then echo present; else echo missing; fi"))))
        present == "missing" && push!(missing, _relative(rel))
    end
    missing
end

function _fetch_outputs(ssh, scp, target, dir, outputs, optional_outputs)
    fetched = Pair{String,String}[]
    absent = String[]
    try
        for (rel, path) in [outputs; optional_outputs]
            remote_path = _sq("$dir/$rel")
            # A negative test and an SSH failure are distinct outcomes.
            hash = strip(String(_capture(ssh(target.host,
                "if test -f $remote_path; then sha256sum -- $remote_path; else echo missing; fi"))))
            if hash == "missing"
                if (rel => path) in optional_outputs
                    push!(absent, path)
                    continue
                end
                error("Remote stage finished without required output: $rel")
            end
            expected = first(split(hash))
            occursin(r"^[0-9a-f]{64}$", expected) || error("Invalid remote output checksum: $rel")
            mkpath(dirname(path))
            part = path * ".part-$(uuid4())"
            push!(fetched, part => path)
            _capture(scp("$(target.host):$dir/$rel", part))
            _digest(part) == expected || error("Remote output transfer checksum mismatch: $rel")
        end
        # No partial download can appear under a checkpoint filename. Publish
        # only after *all* required files have arrived and passed their checks.
        for (part, path) in fetched
            mv(part, path; force=true)
        end
        for path in absent
            rm(path; force=true)
        end
    finally
        for (part, _) in fetched
            isfile(part) && rm(part; force=true)
        end
    end
end

function _verified_receipt(journal, outputs, optional_outputs)
    receipt = get(journal, "retrieved", nothing)
    !isnothing(receipt) &&
        all(get(receipt, abspath(p), nothing) isa AbstractString for (_, p) in outputs) &&
        all(haskey(receipt, abspath(p)) for (_, p) in optional_outputs) &&
        all(isnothing(h) ? !ispath(p) : isfile(p) && _digest(p) == h for (p, h) in receipt)
end

function _sweep(ssh, target)
    # Only terminal attempts older than seven days are eligible. In particular,
    # the age of a directory says nothing about whether its command is alive.
    script = """
    for d in $(_sq(target.base))/mm_*; do
        test ! -L "\$d" && test -f "\$d/$_CONTROL/worker.sh" || continue
        if test -n "\$(find "\$d/$_CONTROL" -maxdepth 1 -name exit -mtime +7 -print)"; then
            (
                exec 9>"\$d/$_CONTROL/run.lock"
                flock -n 9 || exit 0
                exec 8>"\$d/$_CONTROL/control.lock"
                flock -n 8 || exit 0
                # Recheck after acquiring locks; never execute code found in an
                # abandoned directory merely because its name matches a glob.
                test -n "\$(find "\$d/$_CONTROL" -maxdepth 1 -name exit -mtime +7 -print)" &&
                    rm -rf -- "\$d"
            ) || true
        fi
    done
    """
    _capture(ssh(target.host, script))
end

# Recovery paths reach the stage without relaunching it: a verified receipt, or
# a completed attempt retrieved for the first time. Stage entry points open a
# fresh log before calling the runner, so without this the retry would keep only
# the header and lose the command and diagnostics that the remote attempt has
# already produced. Anything the local side already wrote is left untouched, and
# the remote text is merged only when it is not already there: a staged,
# interrupted stage can be reached more than once with the same header.
function _replay_remote_log(ssh, target, dir, log_path)
    _remote_state(ssh, target, dir) == "missing" && return
    emitted = _control(ssh, target, dir, "log 0")
    isempty(emitted) && return
    # A tool's output is not necessarily valid UTF-8, so this compares bytes.
    existing = isfile(log_path) ? read(log_path) : UInt8[]
    isnothing(findfirst(emitted, existing)) || return
    open(io -> write(io, emitted), log_path, "a")
end

# `on_receipt` and `before_run` run inside the owner lock, and a stage callback
# may probe a remote database or tool. A callback that catches every exception
# turns a cancellation into an ordinary error: no stop is requested, and the
# detached workload keeps running while the local stage reports itself
# cancelled. An interrupt is therefore surfaced once the callback finishes
# whatever it is doing, so `_with_owner` sees exactly one stop request.
function _interruptible(f, args...)
    interrupted = InterruptException[]
    result = try
        f(args...)
    catch e
        e isa InterruptException || e isa CompositeException || rethrow()
        e isa InterruptException && push!(interrupted, e)
        e isa CompositeException && append!(interrupted,
            [x for x in e.exceptions if x isa InterruptException])
        nothing
    finally
    end
    isempty(interrupted) || throw(first(interrupted))
    result
end

"""
    run_remote(emit, target, command_for; inputs, outputs, optional_outputs,
               before_run, on_receipt, prepare_local, log_path, stage, context)

Run a detached command and retrieve verified outputs. A durable journal beside
`log_path` reattaches the same workflow after a restart. Content identity includes
the command, input checksums, output names, target and caller-supplied `context`.
Different workflows have different attempt directories even for equal content.
Transport failures retain the attempt; InterruptException explicitly cancels it.
`prepare_local(verified_receipt)` runs under the owner lock before network access;
callers may retain published outputs only for a matching, verified receipt.
See docs/guide/remote-execution.md for the protocol, assumptions and recovery limits.
"""
function run_remote(emit, target, command_for::Function;
                    inputs::AbstractVector{<:Pair}, outputs::AbstractVector{<:Pair},
                    optional_outputs::AbstractVector{<:Pair}=Pair{String,String}[],
                    before_run::Union{Function,Nothing}=nothing,
                    on_receipt::Union{Function,Nothing}=nothing,
                    prepare_local::Union{Function,Nothing}=nothing,
                    log_path::AbstractString, stage::AbstractString, context=String[],
                    _connection=nothing, _poll_interval=1.0)
    mkpath(dirname(log_path))
    journal_path = abspath(log_path) * ".remote.json"
    _with_owner(journal_path * ".lock";
                on_interrupt=() -> _request_cancellation(journal_path, _connection)) do
        identity = _identity(target, stage, command_for, inputs, outputs, optional_outputs, context)
        ctl = joinpath(tempdir(), "mm_$(uuid4())")
        ssh, scp = isnothing(_connection) ? _transport(target, ctl) : _connection
        journal = isfile(journal_path) ? JSON3.read(read(journal_path, String), Dict{String,Any}) : nothing
        verified_receipt = !isnothing(journal) && journal["identity"] == identity &&
            !journal["cancel_requested"] && _verified_receipt(journal, outputs, optional_outputs)
        isnothing(prepare_local) || prepare_local(verified_receipt)
        # The host/root are part of the identity. Retain an obsolete journal for
        # its original destination rather than deleting a path on the new host.
        if !isnothing(journal) && journal["identity"] != identity
            old_identity = JSON3.read(journal["identity"])
            old_target = (; host=String(old_identity[2]), base=String(old_identity[3]),
                threads=old_identity[5], identity_file=get(journal, "identity_file", target.identity_file))
            old_dir = _check_directory(journal, old_target)
            old_ctl = joinpath(tempdir(), "mm_$(uuid4())")
            old_ssh, _ = isnothing(_connection) ? _transport(old_target, old_ctl) : _connection
            if journal["cancel_requested"]
                # A settings change cannot erase a pending stop request. Use
                # the journal's original host/root/key, not the new target.
                try
                    _cancel_remote(old_ssh, old_target, old_dir)
                finally
                    if isnothing(_connection)
                        try _capture(`ssh -o ControlPath=$old_ctl -O exit $(old_target.host)`) catch end
                    end
                end
            else
                # The superseded attempt is abandoned, not awaited: the new
                # identity must be able to launch immediately. It is reclaimed
                # on a best-effort basis -- a denial (an unreachable host, a
                # newer worker that refuses) must not stop the fresh attempt.
                try
                    _cancel_remote(old_ssh, old_target, old_dir; reclaim=true)
                catch superseded_error
                    superseded_error isa InterruptException && rethrow()
                    @warn "Superseded remote attempt not reclaimed; it will be swept once its marker is old enough" host=old_target.host directory=old_dir exception=superseded_error
                finally
                    if isnothing(_connection)
                        try _capture(`ssh -o ControlPath=$old_ctl -O exit $(old_target.host)`) catch end
                    end
                end
            end
            mv(journal_path, journal_path * ".superseded-$(uuid4())")
            journal = nothing
        end
        if isnothing(journal)
            journal = _new_journal(identity, target)
            _save_journal(journal_path, journal)
        end
        dir = _check_directory(journal, target)
        command = command_for(dir)
        try
            if journal["cancel_requested"]
                _cancel_remote(ssh, target, dir)
                journal = _new_journal(identity, target)
                _save_journal(journal_path, journal)
                dir = _check_directory(journal, target)
                command = command_for(dir)
            end
            # The receipt closes the crash window between remote cleanup and
            # the caller committing its own stage checkpoint. Keep it until
            # content changes; reuse requires the actual local output bytes.
            receipt = get(journal, "retrieved", nothing)
            if _verified_receipt(journal, outputs, optional_outputs)
                # Preparation checks may read uploaded files that successful
                # cleanup has already removed (notably taxonomy databases).
                # Receipt recovery can restore caller metadata, but must not
                # repeat checks that require the old staging inputs.
                isnothing(on_receipt) || _interruptible(on_receipt, ssh, dir)
                # The stage's own log is replayed before cleanup, because the
                # attempt's log lives in the directory cleanup removes.
                _replay_remote_log(ssh, target, dir, log_path)
                if _remote_state(ssh, target, dir) != "missing"
                    _control(ssh, target, dir, "cleanup")
                end
                emit("  Recovered verified outputs for $stage")
                return nothing
            end
            state = _remote_state(ssh, target, dir)
            if state == "missing" && !isnothing(receipt)
                journal = _new_journal(identity, target)
                _save_journal(journal_path, journal)
                dir = _check_directory(journal, target)
                command = command_for(dir)
            end
            if state in ("missing", "preparing")
                _prepare(ssh, scp, target, dir, command, identity, inputs, outputs, optional_outputs, before_run, emit)
                # Rehash after upload; changing data or code must never label the
                # remote copy with a fingerprint of earlier local bytes.
                _identity(target, stage, command_for, inputs, outputs, optional_outputs, context) == identity ||
                    error("Remote stage inputs changed during upload; run the stage again")
                command_for(dir) # restore the caller's resolved command record
            else
                emit("  Reattaching $stage on $(target.host):$dir (no upload)")
                isnothing(before_run) || _interruptible(before_run, ssh, dir)
            end
            if journal["failed"] || state == "interrupted" ||
                    (startswith(state, "exit ") && state != "exit 0")
                _control(ssh, target, dir, "retry")
                journal["failed"] = false
                _save_journal(journal_path, journal)
            elseif state == "exit 0"
                # The workload finished and its exit marker is retained, but no
                # receipt was stored: either an earlier fetch did not find a
                # required output, or the local process died between the marker
                # and the receipt. `retry` refuses a successful attempt and
                # `start` refuses a terminal one, so the stage would be wedged
                # until its state was deleted by hand. `replay` prepares the
                # retained attempt for retrieval alone, keeping the marker for
                # the next run if this fetch fails too. Every byte is still
                # checksum-verified against the remote copy.
                emit("  Retrieving completed $stage outputs from $(target.host) without relaunching")
                _control(ssh, target, dir, "replay")
                journal["failed"] = false
                _save_journal(journal_path, journal)
            end
            state == "cancelled" && error("Remote attempt was cancelled; cancel it locally before retrying")
            log_command(command, log_path)
            _control(ssh, target, dir, "start")
            offset = 0
            while true
                state = _remote_state(ssh, target, dir)
                bytes = _control(ssh, target, dir, "log $offset")
                open(io -> write(io, bytes), log_path, "a")
                offset += length(bytes)
                if startswith(state, "exit ")
                    # Drain the bounded log before returning a terminal result.
                    length(bytes) == 65536 && continue
                    code = tryparse(Int, split(state)[2])
                    isnothing(code) && error("Invalid remote exit status")
                    if code != 0
                        journal["failed"] = true
                        _save_journal(journal_path, journal)
                        error("Remote $stage exited with status $code; its checkpoint is retained")
                    end
                    break
                elseif state in ("interrupted", "cancelled", "missing")
                    error("Remote $stage is $state; run the stage again to recover")
                end
                sleep(_poll_interval)
            end
            emit("  Retrieving $(length(outputs)) output file(s) from $(target.host)")
            # A failed fetch (a missing required output, a truncated transfer)
            # throws with the attempt's exit marker still in place. The next run
            # takes the `exit 0` barrier above, which replays the retained
            # attempt instead of relaunching it, so a short read never wedges
            # the stage and its outputs stay available to retrieve again.
            missing = _missing_outputs(ssh, target, dir, outputs)
            if !isempty(missing)
                journal["failed"] = false
                _save_journal(journal_path, journal)
                error("Remote $stage finished with status 0 but did not write required output(s): " *
                      join(missing, ", ") * ". The attempt is retained; run the stage again to " *
                      "retrieve it, or change its settings to run it afresh.")
            end
            _fetch_outputs(ssh, scp, target, dir, outputs, optional_outputs)
            journal["retrieved"] = Dict(abspath(path) => (isfile(path) ? _digest(path) : nothing)
                for (_, path) in [outputs; optional_outputs])
            _save_journal(journal_path, journal)
            # A receipt survives cleanup, including a lost cleanup acknowledgement.
            _control(ssh, target, dir, "cleanup")
            try
                _sweep(ssh, target)
            catch e
                e isa InterruptException && rethrow()
                @debug "Remote terminal-attempt sweep failed" exception=e
            end
        finally
            if isnothing(_connection)
                try
                    _capture(`ssh -o ControlPath=$ctl -O exit $(target.host)`)
                catch
                end
            end
        end
    end
    nothing
end
