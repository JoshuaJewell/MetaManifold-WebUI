# SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

module RemoteExec

# Shared target validation and durable remote execution for DADA2 and phylogeny.

export run_remote, verify_remote_sizes, checked_target

    using Logging
    using ..PipelineLog
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
        quoted = join(_sq.(remote_paths), " ")
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

    include("remote_session.jl")
end
