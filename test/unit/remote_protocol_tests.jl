# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

# Exercise the actual shell worker, process groups, files, transfers and locks.
# Only the transport is local: no account, SSH daemon or bioinformatics tool is
# needed. The real-SSH subprocess tests use the same driver (see docs).
using JSON3
const _RX = MetaManifold.RemoteExec

function _remote_fixture(f)
    @info "Remote execution regression" case=Test.get_testset().description
    mktempdir() do d
        target = (; host="fixture", base=joinpath(d, "remote"), identity_file=nothing, threads=2)
        mkpath(target.base)
        input = joinpath(d, "input")
        write(input, "known input\n")
        uploads = Ref(0)
        downloads = Ref(0)
        fail_download = Ref(false)
        fail_control = Ref(:none)
        probe_launch = Ref(false)
        drop_status = Ref(false)
        drop_start_reply = Ref(false)
        fail_cancel = Ref(false)
        retries = Ref(0)
        ssh = function (_, command)
            if probe_launch[] && occursin("worker.sh' start", command)
                probe_launch[] = false
                dir = JSON3.read(read(joinpath(d, "stage.log.remote.json"), String))["directory"]
                lockfile = joinpath(dir, ".metamanifold", "run.lock")
                held = joinpath(dir, "probe-held")
                return `bash -c '(exec 9>"$1"; flock -x 9; touch "$2"; sleep 0.5) & while ! test -f "$2"; do sleep 0.01; done; bash -c "$3"' _ $lockfile $held $command`
            end
            occursin("worker.sh' retry", command) && (retries[] += 1)
            if drop_status[] && occursin("worker.sh' status", command)
                return `bash -c 'exit 255'`
            end
            if drop_start_reply[] && occursin("worker.sh' start", command)
                drop_start_reply[] = false
                return `bash -c $(command * "; exit 255")`
            end
            if fail_cancel[] && occursin("worker.sh' cancel", command)
                return `bash -c 'exit 255'`
            end
            `bash -c $command`
        end
        scp = function (src, dst)
            if startswith(src, "fixture:")
                downloads[] += 1
                source = src[length("fixture:")+1:end]
                if fail_download[]
                    return `bash -c 'head -c 2 -- "$1" > "$2"; exit 255' _ $source $dst`
                end
                `cp -- $source $dst`
            else
                uploads[] += 1
                destination = dst[length("fixture:")+1:end]
                if (fail_control[] == :worker && src == _RX._WORKER) ||
                   (fail_control[] == :command && occursin("/.metamanifold/command.sh", destination))
                    return `bash -c 'head -c 2 -- "$1" > "$2"; exit 255' _ $src $destination`
                end
                `cp -- $src $destination`
            end
        end
        invoke = function (command; log="stage", value=input, extra_outputs=Pair{String,String}[],
                           optional_outputs=Pair{String,String}[], context=String[],
                           before_run=nothing, on_receipt=nothing)
            _RX.run_remote(_ -> nothing, target, command;
                inputs=[value => "input"], outputs=["output" => joinpath(d, log * ".out"); extra_outputs],
                optional_outputs, before_run, on_receipt,
                log_path=joinpath(d, log * ".log"), stage="denoise", context,
                _connection=(ssh, scp), _poll_interval=0.02)
        end
        try
            f((; d, target, input, uploads, downloads, fail_download, fail_control, probe_launch, drop_status,
                drop_start_reply, fail_cancel, retries, invoke, ssh, scp))
        finally
            # Never let a failed assertion strand a detached test process.
            for dir in readdir(target.base; join=true)
                isfile(joinpath(dir, ".metamanifold", "worker.sh")) || continue
                try _RX._cancel_remote(ssh, target, dir) catch end
            end
        end
    end
end

function _wait_remote(predicate; timeout=120)
    timedwait(predicate, timeout; pollint=0.02) == :ok || error("Timed out waiting for remote fixture")
end

_journal(f, log="stage") = JSON3.read(read(joinpath(f.d, "$log.log.remote.json"), String), Dict{String,Any})
_dir(f, log="stage") = _journal(f, log)["directory"]
_task_error(t) = try fetch(t); nothing catch e; e end

@testset "Durable remote execution protocol" begin
    @testset "a transient status probe cannot lose a detached launch" begin
        _remote_fixture() do f
            f.probe_launch[] = true
            task = @async f.invoke(d -> "cd '$d' && cat input > output")
            try
                finished = timedwait(() -> istaskdone(task), 120; pollint=0.02)
                @test finished == :ok
                if finished == :ok
                    fetch(task)
                    @test read(joinpath(f.d, "stage.out"), String) == "known input\n"
                    @test isempty(readdir(f.target.base))
                end
            finally
                if !istaskdone(task)
                    schedule(task, InterruptException(); error=true)
                    _task_error(task)
                end
            end
        end
    end

    @testset "interrupted control uploads remain retryable and cancellable" begin
        for part in (:command, :worker), cancel in (false, true)
            _remote_fixture() do f
                f.fail_control[] = part
                command = d -> "cd '$d' && cat input > output"
                @test_throws ProcessFailedException f.invoke(command)
                dir = _dir(f)
                @test !isfile(joinpath(dir, ".metamanifold", "worker.sh"))
                @test _RX._remote_state(f.ssh, f.target, dir) == "missing"
                if cancel
                    _RX._request_cancellation(joinpath(f.d, "stage.log.remote.json"), (f.ssh, f.scp))
                    @test !isfile(joinpath(f.d, "stage.log.remote.json"))
                else
                    f.fail_control[] = :none
                    f.invoke(command)
                    @test read(joinpath(f.d, "stage.out"), String) == "known input\n"
                end
                @test !isdir(dir)
            end
        end
    end

    @testset "detached command, output verification, receipt and log" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && echo started && cat input > output"
            f.invoke(command)
            @test read(joinpath(f.d, "stage.out"), String) == "known input\n"
            @test occursin("started", read(joinpath(f.d, "stage.log"), String))
            @test isempty(readdir(f.target.base))
            uploads = f.uploads[]
            downloads = f.downloads[]
            f.invoke(command)
            @test f.uploads[] == uploads
            @test f.downloads[] == downloads
            @test haskey(_journal(f), "retrieved")
        end
    end

    @testset "lost connection reattaches without upload or a second command" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && echo start >> starts && while ! test -f release; do sleep 0.05; done; cat input > output"
            t = @async f.invoke(command)
            _wait_remote(() -> isfile(joinpath(f.d, "stage.log.remote.json")) &&
                isfile(joinpath(_dir(f), "starts")))
            dir = _dir(f)
            f.drop_status[] = true
            @test !isnothing(_task_error(t))
            @test isdir(dir)
            @test !isfile(joinpath(dir, ".metamanifold", "exit"))
            uploads = f.uploads[]
            f.drop_status[] = false
            resumed = @async f.invoke(command)
            sleep(0.1)
            @test f.uploads[] == uploads
            @test readlines(joinpath(dir, "starts")) == ["start"]
            touch(joinpath(dir, "release"))
            fetch(resumed)
            @test f.uploads[] == uploads
            @test !isdir(dir)
        end
    end

    @testset "receipt recovery does not validate deleted staging inputs" begin
        _remote_fixture() do f
            prepared, recovered = Ref(0), Ref(0)
            before_run = function (_, dir)
                isfile(joinpath(dir, "input")) || error("Staged input is missing")
                prepared[] += 1
            end
            on_receipt = (_, _) -> (recovered[] += 1)
            command = d -> "cd '$d' && cat input > output"
            f.invoke(command; before_run, on_receipt)
            @test isempty(readdir(f.target.base))
            uploads = f.uploads[]
            f.invoke(command; before_run, on_receipt)
            @test prepared[] == 1
            @test recovered[] == 1
            @test f.uploads[] == uploads
        end
    end

    @testset "a new Julia process recovers after SIGKILL" begin
        _remote_fixture() do f
            settings = joinpath(f.d, "driver.json")
            write(settings, JSON3.write(Dict("host" => "fixture", "base" => f.target.base,
                "identity_file" => nothing, "local_dir" => f.d)))
            driver = joinpath(@__DIR__, "..", "helpers", "remote_restart_driver.jl")
            project = dirname(dirname(@__DIR__))
            cmd = `$(Base.julia_cmd()) --startup-file=no --project=$project $driver $settings`
            first_process = run(pipeline(cmd; stdout=devnull, stderr=joinpath(f.d, "first.err")); wait=false)
            try
                _wait_remote(; timeout=180) do
                    process_exited(first_process) && error(read(joinpath(f.d, "first.err"), String))
                    isfile(joinpath(f.d, "stage.log.remote.json")) && isfile(joinpath(_dir(f), "starts"))
                end
                dir = _dir(f)
                kill(first_process, Base.SIGKILL)
                wait(first_process)
                uploads = readlines(joinpath(f.d, "uploads"))
                second_process = run(pipeline(cmd; stdout=devnull, stderr=joinpath(f.d, "second.err")); wait=false)
                try
                    sleep(0.3)
                    @test readlines(joinpath(dir, "starts")) == ["start"]
                    touch(joinpath(dir, "release"))
                    _wait_remote(() -> process_exited(second_process); timeout=180)
                    @test success(second_process)
                    @test readlines(joinpath(f.d, "uploads")) == uploads
                    @test read(joinpath(f.d, "stage.out"), String) == "known input\n"
                finally
                    process_running(second_process) && kill(second_process, Base.SIGKILL)
                end
            finally
                process_running(first_process) && kill(first_process, Base.SIGKILL)
            end
        end
    end

    @testset "a lost launch reply cannot launch a second command" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && echo start >> starts && while ! test -f release; do sleep 0.05; done; cat input > output"
            f.drop_start_reply[] = true
            @test_throws ProcessFailedException f.invoke(command)
            dir = _dir(f)
            _wait_remote(() -> isfile(joinpath(dir, "starts")))
            uploads = f.uploads[]
            resumed = @async f.invoke(command)
            sleep(0.1)
            @test readlines(joinpath(dir, "starts")) == ["start"]
            touch(joinpath(dir, "release"))
            fetch(resumed)
            @test f.uploads[] == uploads
        end
    end

    @testset "completed remote work survives a failed partial download" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && cat input > output && echo second > second"
            f.fail_download[] = true
            other = ["second" => joinpath(f.d, "second.out")]
            @test_throws Exception f.invoke(command; extra_outputs=other)
            dir = _dir(f)
            @test read(joinpath(dir, ".metamanifold", "exit"), String) == "0\n"
            @test !isfile(joinpath(f.d, "stage.out"))
            @test !isfile(joinpath(f.d, "second.out"))
            @test !any(occursin(".part-", p) for p in readdir(f.d))
            uploads = f.uploads[]
            f.fail_download[] = false
            f.invoke(command; extra_outputs=other)
            @test f.uploads[] == uploads
            @test read(joinpath(f.d, "second.out"), String) == "second\n"
        end
    end

    @testset "cancel kills the process group and removes staging" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && sleep 300 & wait"
            t = @async f.invoke(command)
            _wait_remote(() -> isfile(joinpath(f.d, "stage.log.remote.json")) &&
                isfile(joinpath(_dir(f), ".metamanifold", "pid")))
            dir = _dir(f)
            pid = parse(Int, first(split(read(joinpath(dir, ".metamanifold", "pid"), String))))
            schedule(t, InterruptException(); error=true)
            @test !isnothing(_task_error(t))
            @test !isdir(dir)
            @test !isfile(joinpath(f.d, "stage.log.remote.json"))
            @test !success(`bash -c $("kill -0 -- -$pid 2>/dev/null")`)
        end
    end

    @testset "cancel during recovery preparation stops the existing attempt" begin
        for phase in (:identity, :command)
            _remote_fixture() do f
                pause = Ref(false)
                waiting = Ref(false)
                command = function (d)
                    if pause[] && ((d == "/METAMANIFOLD_STAGE") == (phase == :identity))
                        waiting[] = true
                        sleep(300)
                    end
                    "cd '$d' && echo started > started && sleep 300 && cat input > output"
                end
                initial = @async f.invoke(command)
                _wait_remote(() -> isfile(joinpath(f.d, "stage.log.remote.json")) &&
                    isfile(joinpath(_dir(f), "started")))
                dir = _dir(f)
                pid = parse(Int, first(split(read(joinpath(dir, ".metamanifold", "pid"), String))))
                f.drop_status[] = true
                @test !isnothing(_task_error(initial))
                f.drop_status[] = false
                pause[] = true
                recovery = @async f.invoke(command)
                _wait_remote(() -> waiting[])
                schedule(recovery, InterruptException(); error=true)
                @test !isnothing(_task_error(recovery))
                @test !isdir(dir)
                @test !isfile(joinpath(f.d, "stage.log.remote.json"))
                @test !success(`bash -c $("kill -0 -- -$pid 2>/dev/null")`)
            end
        end
    end

    @testset "unacknowledged cancellation is durable and retried" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && echo started > started && sleep 1 && cat input > output"
            t = @async f.invoke(command)
            _wait_remote(() -> isfile(joinpath(f.d, "stage.log.remote.json")) &&
                isfile(joinpath(_dir(f), "started")))
            dir = _dir(f)
            f.fail_cancel[] = true
            schedule(t, InterruptException(); error=true)
            @test !isnothing(_task_error(t))
            @test _journal(f)["cancel_requested"]
            @test isdir(dir)
            f.fail_cancel[] = false
            f.invoke(d -> command(d) * " && echo changed-setting")
            @test !isdir(dir)
            @test !_journal(f)["cancel_requested"]
        end
    end

    @testset "equal-content concurrent workflows never share a command" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && echo started > started && while ! test -f release; do sleep 0.05; done; cat input > output"
            a = @async f.invoke(command; log="a")
            b = @async f.invoke(command; log="b")
            _wait_remote(() -> all(isfile(joinpath(f.d, "$n.log.remote.json")) &&
                isfile(joinpath(_dir(f, n), "started")) for n in ("a", "b")))
            da, db = _dir(f, "a"), _dir(f, "b")
            @test da != db
            @test_throws ErrorException f.invoke(command; log="a")
            touch(joinpath(da, "release"))
            touch(joinpath(db, "release"))
            fetch(a); fetch(b)
            @test read(joinpath(f.d, "a.out")) == read(joinpath(f.d, "b.out"))
        end
    end

    @testset "changed bytes, settings and external data cannot reuse an attempt" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && cat input > output"
            f.fail_download[] = true
            @test_throws Exception f.invoke(command)
            original = _dir(f)
            write(f.input, "other input\n") # same size: mtime/size is not the identity
            @test_throws Exception f.invoke(command)
            changed = _dir(f)
            @test changed != original
            @test_throws Exception f.invoke(d -> command(d) * " && echo setting")
            @test _dir(f) != changed
            settings = _dir(f)
            @test_throws Exception f.invoke(d -> command(d) * " && echo setting"; context=["database-checksum"])
            @test _dir(f) != settings
        end
    end

    @testset "a retry retains the failed tool's own checkpoint" begin
        _remote_fixture() do f
            command = d -> "cd '$d'; if ! test -f checkpoint; then touch checkpoint; exit 7; fi; cat input > output"
            @test_throws ErrorException f.invoke(command)
            dir = _dir(f)
            @test _journal(f)["failed"]
            @test isfile(joinpath(dir, "checkpoint"))
            # Simulate a local death before the nonzero exit was journalled.
            journal = _journal(f)
            journal["failed"] = false
            _RX._save_journal(joinpath(f.d, "stage.log.remote.json"), journal)
            uploads = f.uploads[]
            f.invoke(command)
            @test f.uploads[] == uploads
            @test !isdir(dir)
        end
    end

    @testset "a dead remote wrapper leaves a retryable checkpoint" begin
        _remote_fixture() do f
            command = d -> "cd '$d'; if ! test -f checkpoint; then touch checkpoint; sleep 300; fi; cat input > output"
            t = @async f.invoke(command)
            _wait_remote(() -> isfile(joinpath(f.d, "stage.log.remote.json")) &&
                isfile(joinpath(_dir(f), "checkpoint")))
            dir = _dir(f)
            pid = parse(Int, first(split(read(joinpath(dir, ".metamanifold", "pid"), String))))
            run(`kill -KILL $pid`)
            _wait_remote(() -> istaskdone(t))
            @test !isnothing(_task_error(t))
            @test _RX._remote_state(f.ssh, f.target, dir) == "interrupted"
            uploads = f.uploads[]
            f.invoke(command)
            @test f.uploads[] == uploads
            @test !isdir(dir)
        end
    end

    @testset "retry waits for a terminal worker to release its lock" begin
        _remote_fixture() do f
            command = d -> "cd '$d'; if ! test -f checkpoint; then touch checkpoint; exit 7; fi; cat input > output"
            @test_throws ErrorException f.invoke(command)
            dir = _dir(f)
            # Hold the real execution lock after the failure marker exists,
            # reproducing the wrapper/watchdog's terminal handoff window.
            held, release = joinpath(dir, "held"), joinpath(dir, "release")
            lockfile = joinpath(dir, ".metamanifold", "run.lock")
            holder = run(`bash -c 'exec 9>"$1"; flock -x 9; touch "$2"; while ! test -f "$3"; do sleep 0.02; done' _ $lockfile $held $release`; wait=false)
            retry = nothing
            try
                _wait_remote(() -> isfile(held))
                retry = @async f.invoke(command)
                _wait_remote(() -> f.retries[] > 0)
                @test timedwait(() -> istaskdone(retry), 1.0; pollint=0.02) == :timed_out
            finally
                touch(release)
                wait(holder)
            end
            isnothing(retry) || fetch(retry)
            @test read(joinpath(f.d, "stage.out"), String) == "known input\n"
            @test !isdir(dir)
        end
    end

    @testset "missing outputs never publish a checkpoint or retain stale optional files" begin
        _remote_fixture() do f
            @test_throws ErrorException f.invoke(d -> "cd '$d' && true")
            @test !isfile(joinpath(f.d, "stage.out"))
            optional = joinpath(f.d, "optional.out")
            write(optional, "stale data")
            command = d -> "cd '$d' && cat input > output"
            f.invoke(command; optional_outputs=["optional" => optional])
            @test !isfile(optional)
            uploads = f.uploads[]
            f.invoke(command; optional_outputs=["optional" => optional])
            @test f.uploads[] == uploads
            # A stale optional file appearing after cleanup invalidates the receipt.
            write(optional, "stale again")
            f.invoke(command; optional_outputs=["optional" => optional])
            @test !isfile(optional)
            @test f.uploads[] > uploads
        end
    end

    @testset "retention removes only terminal old attempts" begin
        _remote_fixture() do f
            f.fail_download[] = true
            @test_throws Exception f.invoke(d -> "cd '$d' && cat input > output")
            old = _dir(f)
            run(`touch -d '9 days ago' $(joinpath(old, ".metamanifold", "exit"))`)
            f.fail_download[] = false
            command = d -> "cd '$d' && echo started > started && while ! test -f release; do sleep 0.05; done; cat input > output"
            task = @async f.invoke(command; log="live")
            _wait_remote(() -> isfile(joinpath(f.d, "live.log.remote.json")) &&
                isfile(joinpath(_dir(f, "live"), "started")))
            live = _dir(f, "live")
            run(`touch -d '9 days ago' $live`)
            _RX._sweep(f.ssh, f.target)
            @test !isdir(old)
            @test isdir(live)
            @test _RX._remote_state(f.ssh, f.target, live) == "running"
            touch(joinpath(live, "release"))
            fetch(task)
        end
    end

    @testset "a completed attempt is retrieved again instead of wedging the stage" begin
        # A required output can be missing from a zero-exit attempt (the stage
        # is then retried), and a fetch can be interrupted after the attempt
        # finished. Both leave a terminal `exit 0` marker with no receipt.
        _remote_fixture() do f
            command = d -> "cd '$d' && cat input > output"
            @test_throws ErrorException f.invoke(d -> "cd '$d' && true")
            dir = _dir(f)
            @test _RX._remote_state(f.ssh, f.target, dir) == "exit 0"
            @test !haskey(_journal(f), "retrieved")
            uploads = f.uploads[]
            # Running the same stage again retrieves without uploading or
            # launching, and the checkpoint becomes verified.
            f.invoke(command)
            @test f.uploads[] == uploads
            @test read(joinpath(f.d, "stage.out"), String) == "known input\n"
            @test haskey(_journal(f), "retrieved")
            @test !isdir(dir)
        end
    end

    @testset "a superseded preparation is reclaimed instead of stranded" begin
        _remote_fixture() do f
            # Fail the control-script upload after the inputs are in place, so
            # the attempt has no worker and no launch marker: exactly the state
            # a `run_remote` sweep skips.
            f.fail_control[] = :worker
            @test_throws ProcessFailedException f.invoke(d -> "cd '$d' && cat input > output")
            old = _dir(f)
            @test isfile(joinpath(old, "input"))
            # A different identity supersedes it. Best-effort reclamation runs
            # even though the first attempt never installed a worker.
            f.fail_control[] = :none
            command = d -> "cd '$d' && cat input > output"
            f.invoke(command; log="second")
            @test !isdir(old)
            @test read(joinpath(f.d, "second.out"), String) == "known input\n"
            @test isempty(readdir(f.target.base))
        end
    end

    @testset "a declined reclamation does not block the new attempt" begin
        _remote_fixture() do f
            f.fail_control[] = :worker
            @test_throws ProcessFailedException f.invoke(d -> "cd '$d' && cat input > output")
            old = _dir(f)
            f.fail_control[] = :none
            # A host that refuses the reclaim (here: every control command for
            # the superseded directory) must not fail the fresh attempt.
            ssh = function (host, command)
                occursin(old, command) && return `bash -c 'exit 255'`
                f.ssh(host, command)
            end
            target = f.target
            _RX.run_remote(_ -> nothing, target, d -> "cd '$d' && cat input > output";
                inputs=[f.input => "input"], outputs=["output" => joinpath(f.d, "stage.out")],
                log_path=joinpath(f.d, "stage.log"), stage="denoise", context=["changed"],
                _connection=(ssh, f.scp), _poll_interval=0.02)
            @test read(joinpath(f.d, "stage.out"), String) == "known input\n"
        end
    end

    @testset "a finished attempt keeps its remote diagnostics in the stage log" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && echo stage-diagnostics && cat input > output"
            @test_throws ErrorException f.invoke(d -> "cd '$d' && true")
            # A missing local log is exactly what the stage entry points leave
            # behind: they truncate it before calling the runner.
            log = joinpath(f.d, "stage.log")
            @test isfile(log)
            rm(log)
            f.invoke(command)
            text = read(log, String)
            @test occursin("stage-diagnostics", text) || occursin("cmd:", text)
        end
    end

    @testset "an interrupt that a callback swallows still cancels the attempt" begin
        _remote_fixture() do f
            command = d -> "cd '$d' && cat input > output"
            # Leave a completed attempt with no receipt: the next run reattaches
            # without uploading and calls `before_run`, which is where a stage
            # probes a remote database or tool.
            f.fail_download[] = true
            @test_throws Exception f.invoke(command)
            f.fail_download[] = false
            dir = _dir(f)
            @test _RX._remote_state(f.ssh, f.target, dir) == "exit 0"
            entered = Ref(false)
            before_run = function (_, _)
                entered[] = true
                # A production probe that catches every exception.
                try
                    sleep(300)
                catch
                end
            end
            task = @async f.invoke(command; before_run=before_run)
            _wait_remote(() -> entered[])
            schedule(task, InterruptException(); error=true)
            _wait_remote(() -> istaskdone(task))
            @test _task_error(task) isa InterruptException
            # Without the interrupt escaping the callback, no stop would be
            # requested: the journal and the live attempt would both remain.
            @test !isfile(joinpath(f.d, "stage.log.remote.json"))
            @test !isdir(dir)
        end
    end

    @testset "the interrupt guard passes values and every other failure through" begin
        @test _RX._interruptible(() -> 42) == 42
        @test _RX._interruptible(path -> path * "!", "x") == "x!"
        @test_throws InterruptException _RX._interruptible(() -> throw(InterruptException()))
        @test_throws ErrorException _RX._interruptible(() -> error("probe failed"))
        @test_throws InterruptException _RX._interruptible() do
            throw(CompositeException([ErrorException("probe failed"), InterruptException()]))
        end
    end

    @testset "path traversal and reserved control files are refused" begin
        for path in ("/absolute", "../outside", "a/../../b", "a//b", ".metamanifold/pid", "a\nb", "a\\b")
            @test_throws ErrorException _RX._relative(path)
        end
        root_target = (; base="/")
        journal = _RX._new_journal("root fixture", root_target)
        @test dirname(_RX._check_directory(journal, root_target)) == "/"
    end
end
