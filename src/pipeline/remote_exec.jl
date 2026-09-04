module RemoteExec

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

# One command on the bioserver, in a staging directory of its own. The inputs
# are uploaded and checked, the command runs, the outputs come back, and the
# directory is removed. DADA2's remote stages and the phylogeny steps both run
# through here.

export run_remote, verify_remote_sizes, checked_target

    using Logging
    using ..PipelineLog
    using ..Tools: _run_killable
    using ..Validation: is_shell_safe_arg

    """
        checked_target(rc, stage; threads) -> NamedTuple or nothing

    The server named by the remote block `rc`, or nothing when it names none.
    The host and staging root become words of the ssh command string, so each
    is refused here if a shell would split or expand it.
    """
    function checked_target(rc::AbstractDict, stage::AbstractString; threads)
        host = get(rc, "host", nothing)
        isnothing(host) && return nothing

        base = get(rc, "staging_dir", nothing)
        isnothing(base) &&
            error("remote.staging_dir must be set explicitly to run '$stage' on $host")
        base = string(base)
        is_shell_safe_arg(base) ||
            error("remote.staging_dir may not contain shell metacharacters or spaces (got: '$base')")
        startswith(base, "/") ||
            error("remote.staging_dir must be an absolute path on the server (got: '$base')")

        # A leading '-' would be read by ssh as an option (e.g. -oProxyCommand).
        (is_shell_safe_arg(string(host)) && !startswith(string(host), "-")) ||
            error("remote.host may not start with '-' or contain shell metacharacters or spaces (got: '$host')")

        idf = get(rc, "identity_file", nothing)
        (; host = string(host),
           base,
           identity_file = isnothing(idf) ? nothing : string(idf),
           threads)
    end

    # scp reports its own failures, but a transfer that ends early - a dropped
    # connection, a full disk on the server - can leave a short file that the
    # remote side will read the front of without complaint. Compare byte counts
    # before anything runs against them.
    function verify_remote_sizes(ssh, target, staging_dir, uploads, emit)
        isempty(uploads) && return nothing
        remote_paths = ["$staging_dir/$rel" for (_, rel) in uploads]
        quoted = join(["'" * p * "'" for p in remote_paths], " ")
        out = try
            read(ssh(target.host, "wc -c $quoted"), String)
        catch err
            error("Remote transfer check: could not stat uploaded files on " *
                  "$(target.host): $(sprint(showerror, err))")
        end

        sizes = Dict{String,Int}()
        for line in eachline(IOBuffer(out))
            parts = split(strip(line); limit=2)
            length(parts) == 2 || continue
            parts[2] == "total" && continue
            n = tryparse(Int, parts[1])
            isnothing(n) || (sizes[String(parts[2])] = n)
        end

        for ((local_path, rel), remote_path) in zip(uploads, remote_paths)
            want = filesize(local_path)
            got  = get(sizes, remote_path, nothing)
            isnothing(got) &&
                error("Remote transfer check: $remote_path is missing on $(target.host) " *
                      "after upload")
            got == want ||
                error("Remote transfer check: $remote_path is $got bytes on " *
                      "$(target.host) but $want bytes locally - the upload was " *
                      "truncated. Re-run the stage.")
        end
        emit("  Verified $(length(uploads)) uploaded file(s) on $(target.host)")
        return nothing
    end

    """
        run_remote(emit, target, command_for; inputs, outputs, optional_outputs,
                   before_run, log_path, stage)

    Run `command_for(staging_dir)` on `target.host` and bring its outputs back.

    `inputs` are `local_path => remote_relative_path` pairs uploaded before the
    command runs; `outputs` are `remote_relative_path => local_path` pairs
    retrieved after it. Every remote path is formed under a staging directory
    created for this call alone and removed again afterwards, so two stages - or
    two runs of the same stage - never share state on the server.

    `optional_outputs` take the same form but are retrieved only if the command
    produced them. `before_run(ssh, staging_dir)` runs once the inputs are
    verified, over the same connection; checks of files already on the server
    go there.

    The command string reaches a remote login shell, so the caller builds it
    only from values that shell cannot word-split or expand.
    """
    function run_remote(emit, target, command_for::Function;
                        inputs::AbstractVector{<:Pair},
                        outputs::AbstractVector{<:Pair},
                        optional_outputs::AbstractVector{<:Pair}=Pair{String,String}[],
                        before_run::Union{Function,Nothing}=nothing,
                        log_path::AbstractString,
                        stage::AbstractString)
        # Unique per call: the study route denoises runs on parallel threads, so
        # a seconds-resolution stamp alone would collide and two runs would write
        # into one staging directory.
        run_id      = "$(stage)_$(floor(Int, time()))_$(rand(UInt32))"
        staging_dir = "$(target.base)/$run_id"

        # ControlMaster reuses a single SSH connection for every ssh/scp call, so
        # a key-less setup prompts for the password once rather than per transfer.
        ctl      = "/tmp/ssh_mux_$run_id"
        id_opt   = isnothing(target.identity_file) ? `` : `-i $(target.identity_file)`
        ssh_opts = `$id_opt -o ControlMaster=auto -o ControlPath=$ctl -o ControlPersist=yes -o ConnectTimeout=15 -o NumberOfPasswordPrompts=1 -o ServerAliveInterval=30 -o ServerAliveCountMax=3`
        ssh = (args...) -> `ssh $ssh_opts $args`
        scp = (args...) -> `scp $ssh_opts $args`

        if isnothing(target.identity_file)
            emit("  Connecting to $(target.host) (enter SSH password if prompted)...")
        else
            emit("  Connecting to $(target.host) (key: $(target.identity_file))...")
        end

        try
            # Every directory a file is written into must exist before it is
            # used: scp will not create one for an upload, and the remote command
            # creates only the output directories it knows about.
            sub_dirs = unique(filter(!isempty, [
                [dirname(string(last(i)))  for i in inputs];
                [dirname(string(first(o))) for o in vcat(collect(outputs),
                                                         collect(optional_outputs))]
            ]))
            mkdirs   = join(["$staging_dir/$d" for d in sub_dirs], " ")
            emit("  Setting up staging directory on $(target.host)")
            run(ssh(target.host, "mkdir -p $staging_dir $mkdirs"))

            # A crash hard enough to kill the trap below - SIGKILL, or the server
            # rebooting - leaves a staging directory nothing will ever collect.
            # Sweep the ones old enough that no stage could still be using them;
            # a week is far longer than any single stage takes and short enough
            # that a full read set is not left lying around indefinitely.
            try
                run(ssh(target.host,
                        "find $(target.base) -maxdepth 1 -type d -name '*_*_*' " *
                        "-mtime +7 -exec rm -rf {} + 2>/dev/null || true"))
            catch e
                @debug "Remote: Could not sweep stale staging directories" exception=e
            end

            emit("  Transferring $(length(inputs)) input file(s) to $(target.host)")
            for (local_path, remote_rel) in inputs
                isfile(local_path) || error("Remote stage input not found: $local_path")
                run(scp(local_path, "$(target.host):$staging_dir/$remote_rel"))
            end
            verify_remote_sizes(ssh, target, staging_dir, inputs, emit)

            isnothing(before_run) || before_run(ssh, staging_dir)

            # If this machine goes down mid-stage, sshd SIGHUPs the session and
            # takes the command with it - but the `finally` below never runs, so
            # nothing would remove the staging directory, and for a read-consuming
            # stage that is the whole filtered read set left on the server. The
            # trap makes the remote side clean up after itself in that case. It
            # deliberately does not catch EXIT: on a normal finish the outputs are
            # still sitting in the staging directory waiting to be fetched.
            command    = command_for(staging_dir)
            remote_cmd = "trap 'rm -rf $staging_dir' HUP INT TERM PIPE; $command"

            emit("  Running $stage on $(target.host):$staging_dir")
            # The resolved command joins the tool log behind the same marker the
            # local stages use, so a remote run is as reproducible from the log
            # as a local one.
            log_command(remote_cmd, log_path)
            open(log_path, "a") do io
                _run_killable(pipeline(ssh(target.host, remote_cmd); stdout=io, stderr=io))
            end

            emit("  Retrieving $(length(outputs)) output file(s) from $(target.host)")
            fetch = function (remote_rel, local_path)
                mkpath(dirname(local_path))
                # Land the file under a temporary name and rename it into place.
                # A checkpoint is judged current by its mtime, so a transfer cut
                # short by a dropped connection would otherwise leave a truncated
                # file that looks newer than its inputs - and the stage would skip
                # over it on the next run rather than redo it.
                part = local_path * ".part"
                try
                    run(scp("$(target.host):$staging_dir/$remote_rel", part))
                    mv(part, local_path; force=true)
                finally
                    isfile(part) && rm(part; force=true)
                end
            end
            for (remote_rel, local_path) in outputs
                fetch(remote_rel, local_path)
            end
            for (remote_rel, local_path) in optional_outputs
                # `test -f` over the established master connection is cheaper
                # than a failed scp, and distinguishes "the command chose not to
                # write this" from "the transfer broke".
                success(ssh(target.host, "test -f $staging_dir/$remote_rel")) || continue
                fetch(remote_rel, local_path)
            end
        finally
            # Cleanup runs whether or not the stage succeeded: the remote log has
            # already been streamed into log_path, so a retained staging
            # directory would hold nothing a failure diagnosis needs while
            # accumulating a full read set per failed attempt.
            parts = filter(!isempty, split(staging_dir, '/'))
            if length(parts) >= 3
                # Retried once. A transfer killed mid-write leaves the server's
                # end of that scp briefly holding the partial file, and on the
                # network filesystems these staging roots usually live on that
                # turns into a silly-renamed .nfsXXXX entry - so the first
                # `rm -rf` fails with "Directory not empty" and a moment later
                # the same command succeeds. Leaving it is not cheap: a failed
                # read-consuming stage strands its whole upload.
                cleaned = false
                for attempt in 1:2
                    attempt > 1 && sleep(2)
                    try
                        run(ssh(target.host, "rm -rf $staging_dir"))
                        cleaned = true
                        break
                    catch e
                        attempt == 2 && @warn "Remote: cleanup failed - remove it by hand" host=target.host staging_dir exception=e
                    end
                end
                cleaned || emit("  Left $(target.host):$staging_dir behind - remove it by hand")
            else
                @warn "Remote: skipping cleanup - staging path looks too shallow to delete safely" staging_dir
            end
            try
                run(`ssh -o ControlPath=$ctl -O exit $(target.host)`)
            catch
                # The master exits with the last channel; nothing is leaked if it
                # has already gone, and failing here would mask the real error.
            end
        end
        return nothing
    end
end
