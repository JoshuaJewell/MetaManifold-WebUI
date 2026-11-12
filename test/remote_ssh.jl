# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

# Dedicated loopback sshd, keys and staging: no user's SSH configuration changes.
# Run on Linux with OpenSSH server/client installed and /run/sshd available:
# julia --project=. -t2 test/remote_ssh.jl
using Test, JSON3, Sockets, UUIDs
include("helpers/remote_protocol_module.jl")
const RX = RemoteProtocol.RemoteExec

function wait_for(predicate; timeout=180)
    timedwait(predicate, timeout; pollint=0.05) == :ok || error("SSH fixture timed out")
end

@testset "Real SSH recovery" begin
    Sys.islinux() || error("The SSH fixture requires Linux")
    sshd = something(Sys.which("sshd"), isfile("/usr/sbin/sshd") ? "/usr/sbin/sshd" : nothing)
    mktempdir() do root
        hostkey, userkey = joinpath(root, "hostkey"), joinpath(root, "userkey")
        run(`ssh-keygen -q -t ed25519 -N "" -f $hostkey`)
        run(`ssh-keygen -q -t ed25519 -N "" -f $userkey`)
        socket = listen(ip"127.0.0.1", 0)
        port = getsockname(socket)[2]
        close(socket)
        server_config = joinpath(root, "sshd_config")
        client_config = joinpath(root, "ssh_config")
        known_hosts = joinpath(root, "known_hosts")
        write(known_hosts, "mm66-fixture " * read(hostkey * ".pub", String))
        write(server_config, """
        ListenAddress 127.0.0.1
        Port $port
        HostKey $hostkey
        PidFile $(joinpath(root, "sshd.pid"))
        AuthorizedKeysFile $userkey.pub
        StrictModes no
        UsePAM no
        PasswordAuthentication no
        KbdInteractiveAuthentication no
        PubkeyAuthentication yes
        Subsystem sftp internal-sftp
        """)
        write(client_config, """
        Host mm66-fixture
            HostName 127.0.0.1
            Port $port
            User $(strip(read(`id -un`, String)))
            IdentityFile $userkey
            IdentitiesOnly yes
            BatchMode yes
            StrictHostKeyChecking yes
            HostKeyAlias mm66-fixture
            UserKnownHostsFile $known_hosts
        """)
        server = run(pipeline(`$sshd -D -e -f $server_config`; stdout=devnull,
            stderr=joinpath(root, "sshd.log")); wait=false)
        controls = String[]
        remote_dirs = String[]
        ssh = (args...) -> `ssh -F $client_config $args`
        try
            wait_for() do
                process_exited(server) && error(read(joinpath(root, "sshd.log"), String))
                try
                    close(connect(ip"127.0.0.1", port)); true
                catch
                    false
                end
            end
            run(ssh("mm66-fixture", "true"))
            for complete_offline in (false, true)
                @testset "reconnect $(complete_offline ? "completed" : "live") command" begin
                    local_dir = mkpath(joinpath(root, string(complete_offline)))
                    base = mkpath(joinpath(local_dir, "remote"))
                    push!(remote_dirs, base)
                    write(joinpath(local_dir, "input"), "verified over SSH\n")
                    settings = joinpath(local_dir, "driver.json")
                    write(settings, JSON3.write(Dict("host" => "mm66-fixture", "base" => base,
                        "identity_file" => userkey, "ssh_config" => client_config, "local_dir" => local_dir)))
                    driver = joinpath(@__DIR__, "helpers", "remote_restart_driver.jl")
                    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(@__DIR__)) $driver $settings`
                    err = joinpath(local_dir, "driver.err")
                    journal = joinpath(local_dir, "stage.log.remote.json")
                    first_process = run(pipeline(cmd; stdout=devnull, stderr=err); wait=false)
                    dir = ""
                    try
                        wait_for() do
                            process_exited(first_process) && error(read(err, String))
                            isfile(journal) || return false
                            dir = JSON3.read(read(journal, String))["directory"]
                            isfile(joinpath(dir, "starts"))
                        end
                        kill(first_process, Base.SIGKILL)
                        wait(first_process)
                        ctl = read(joinpath(local_dir, "control_path"), String)
                        push!(controls, ctl)
                        run(pipeline(ignorestatus(ssh("-S", ctl, "-O", "exit", "mm66-fixture")); stderr=devnull))
                        uploads = readlines(joinpath(local_dir, "uploads"))
                        @test RX._remote_state(ssh, (; host="mm66-fixture"), dir) == "running"
                        if complete_offline
                            touch(joinpath(dir, "release"))
                            wait_for(() -> isfile(joinpath(dir, ".metamanifold", "exit")))
                            @test read(joinpath(dir, ".metamanifold", "exit"), String) == "0\n"
                        end
                        out = joinpath(local_dir, "resumed.log")
                        second_process = run(pipeline(cmd; stdout=out, stderr=err); wait=false)
                        try
                            wait_for(() -> process_exited(second_process) || occursin("Reattaching", read(out, String)))
                            @test readlines(joinpath(local_dir, "uploads")) == uploads
                            if !complete_offline
                                @test readlines(joinpath(dir, "starts")) == ["start"]
                                touch(joinpath(dir, "release"))
                            end
                            wait_for(() -> process_exited(second_process))
                            success(second_process) || error(read(err, String))
                            @test read(joinpath(local_dir, "stage.out"), String) == "verified over SSH\n"
                            @test readlines(joinpath(local_dir, "uploads")) == uploads
                            @test isempty(readdir(base))
                        finally
                            process_running(second_process) && kill(second_process, Base.SIGKILL)
                            append!(controls, readlines(joinpath(local_dir, "control_paths")))
                        end
                    finally
                        process_running(first_process) && kill(first_process, Base.SIGKILL)
                    end
                end
            end
            @testset "cancel over SSH" begin
                base = mkpath(joinpath(root, "cancel"))
                push!(remote_dirs, base)
                target = (; host="mm66-fixture", base, identity_file=userkey, threads=2)
                ctl = joinpath(tempdir(), "mm_$(uuid4())")
                push!(controls, ctl)
                raw_ssh, raw_scp = RX._transport(target, ctl)
                configured = cmd -> Cmd([first(cmd.exec), "-F", client_config, cmd.exec[2:end]...])
                connection = ((args...) -> configured(raw_ssh(args...)), (args...) -> configured(raw_scp(args...)))
                log = joinpath(root, "cancel.log")
                task = @async RX.run_remote(_ -> nothing, target, _ -> "sleep 300";
                    inputs=Pair{String,String}[], outputs=Pair{String,String}[], log_path=log,
                    stage="denoise", _connection=connection, _poll_interval=0.05)
                wait_for() do
                    istaskdone(task) && fetch(task)
                    isfile(log * ".remote.json") || return false
                    dir = JSON3.read(read(log * ".remote.json", String))["directory"]
                    isfile(joinpath(dir, ".metamanifold", "pid"))
                end
                schedule(task, InterruptException(); error=true)
                @test_throws TaskFailedException fetch(task)
                @test isempty(readdir(base))
                @test !isfile(log * ".remote.json")
            end
        finally
            for base in remote_dirs, dir in readdir(base; join=true)
                isfile(joinpath(dir, ".metamanifold", "worker.sh")) || continue
                try RX._cancel_remote(ssh, (; host="mm66-fixture"), dir) catch end
            end
            for ctl in unique(controls)
                run(pipeline(ignorestatus(ssh("-S", ctl, "-O", "exit", "mm66-fixture")); stderr=devnull))
            end
            process_running(server) && kill(server)
            wait(server)
        end
    end
end
