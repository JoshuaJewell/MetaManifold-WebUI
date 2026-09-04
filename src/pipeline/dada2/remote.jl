# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Remote stage execution (bioserver offload)
#
# Four DADA2 stages can run on a compute server instead of this machine. Each
# one ships its own inputs and collects its own outputs, so the stages are
# independent: a run denoised locally can have its chimeras removed on the
# server, and a run whose error model was learnt on the server can be denoised
# here. Nothing is left behind on the server between stages, and no stage
# assumes the stage before it ran in the same place.
#
# That independence is the whole point, and it costs a transfer. learn_errors
# and denoise read the filtered FASTQ set, so offloading either uploads those
# reads; chimera_removal and assign_taxonomy read checkpoints, which are far
# smaller. Choose per stage in `remote.stages` accordingly.

    # The stages `remote.stages` may name; defined in Validation so a config
    # naming anything else is refused at the write gate rather than silently
    # running that stage locally, which looks identical to a working offload.
    using ..Validation: REMOTE_STAGES

    # Everything this module interpolates into the ssh command string becomes an
    # unquoted word there, so it is checked with Validation's stricter gate: a
    # space cannot make the remote shell run anything, but it does split one
    # path into two arguments. Read paths and sample names never reach that
    # string at all - they travel in an uploaded manifest - so what is checked
    # here is the staging root and the stage parameters.
    const _remote_shell_safe = Validation.is_shell_safe_arg

    _as_dict(x) = x isa AbstractDict ? x : Dict{String,Any}()

    """
        _r_threads(full_cfg) -> Union{Bool,Int}

    The thread count every DADA2 R stage runs with, read from the top-level
    `r_threads`.

    `dada2.taxonomy.multithread` is honoured as a deprecated override for
    assign_taxonomy alone, so a project config written before the key moved
    keeps the thread count it asked for. It is not consulted for any other
    stage: a key that used to mean "threads for assignTaxonomy" would be a
    surprising thing to find governing denoising.
    """
    function _r_threads(full_cfg::AbstractDict; stage::AbstractString="")
        if stage == "assign_taxonomy"
            legacy = get(_as_dict(get(_as_dict(get(full_cfg, "dada2", nothing)), "taxonomy", nothing)),
                         "multithread", nothing)
            isnothing(legacy) || return legacy
        end
        get(full_cfg, "r_threads", 4)
    end

    # DADA2 reads `multithread` as a union of Bool and integer; render it as R
    # source text so the remote script can restore the distinction. A string
    # would be rejected by DADA2 and silently downgraded to a single thread.
    _mt_str(x) = x isa Bool ? (x ? "TRUE" : "FALSE") : string(x)

    """
        _remote_target(full_cfg, stage) -> NamedTuple or nothing

    Resolve where `stage` runs. Returns `nothing` when it runs locally.

    The top-level `remote` block names one server and lists the stages to send
    there. `dada2.taxonomy.remote` is honoured as a deprecated per-stage block
    for assign_taxonomy: when it names a host it wins outright, which is what
    keeps a config written before the block moved behaving exactly as it did.
    """
    function _remote_target(full_cfg::AbstractDict, stage::AbstractString)
        legacy = stage == "assign_taxonomy" ?
                 _as_dict(get(_as_dict(get(_as_dict(get(full_cfg, "dada2", nothing)),
                                           "taxonomy", nothing)), "remote", nothing)) :
                 Dict{String,Any}()
        globl  = _as_dict(get(full_cfg, "remote", nothing))

        # A legacy block that names a host enables the stage on its own, exactly
        # as it did before `remote.stages` existed. Otherwise the stage must be
        # listed, so that configuring a server does not silently move every
        # offloadable stage onto it.
        rc = !isnothing(get(legacy, "host", nothing)) ? legacy :
             (stage in string.(get(globl, "stages", String[])) ? globl : return nothing)

        # `remote.threads` lets the server run wider than this machine, which is
        # usually the reason for offloading in the first place.
        threads = let t = get(rc, "threads", nothing)
            isnothing(t) ? _r_threads(full_cfg; stage) : t
        end
        target = checked_target(rc, stage; threads)
        isnothing(target) && return nothing

        rscript = string(get(rc, "rscript", "Rscript"))
        _remote_shell_safe(rscript) ||
            error("remote.rscript may not contain shell metacharacters or spaces (got: '$rscript')")
        (; target..., rscript)
    end

    # Write a newline-delimited manifest for the remote side. Read paths and
    # sample names go over as a file rather than as command-line arguments so
    # that nothing derived from a filename is ever parsed by a login shell.
    function _write_manifest(dir::String, name::String, entries)
        path = joinpath(dir, name)
        open(path, "w") do io
            for e in entries
                println(io, e)
            end
        end
        return path
    end

    # Reference databases reach the server either by upload or as a pre-existing
    # `remote_path`, and a truncated one there is exactly as damaging as locally:
    # R reads the readable prefix and classifies against a partial reference set.
    # Nothing else on the remote side would notice, so test the stream there.
    function _verify_remote_gzip(ssh, target, paths, emit)
        isempty(paths) && return nothing
        for p in paths
            ok = try
                success(ssh(target.host, "gzip -t '" * p * "'"))
            catch err
                error("Remote database check: could not test $p on $(target.host): " *
                      "$(sprint(showerror, err))")
            end
            ok || error("Remote database check: $p on $(target.host) is a truncated " *
                        "or corrupt gzip stream. Replace it on the server (or clear " *
                        "databases.yml remote_path so the verified local copy is sent).")
        end
        emit("  Verified $(length(paths)) reference file(s) on $(target.host)")
        return nothing
    end

    function _verify_remote_sha256(ssh, target, expected::AbstractDict, emit)
        for (p, want) in expected
            _remote_shell_safe(string(p)) || error("Remote database path may not contain shell metacharacters: $p")
            out = try
                readchomp(ssh(target.host, "sha256sum '" * string(p) * "'"))
            catch err
                error("Remote database check: could not hash $p on $(target.host): $(sprint(showerror, err))")
            end
            got = lowercase(first(split(out)))
            got == lowercase(string(want)) ||
                error("Remote database check: $p on $(target.host) has sha256 $got, but databases.yml " *
                      "expects $want. The server holds a different release; replace it or clear remote_path.")
        end
        isempty(expected) || emit("  Checksum matches for $(length(expected)) reference file(s) on $(target.host)")
        nothing
    end

    """
        _run_remote_stage(emit, target, script, params; inputs, outputs, log_path, stage)

    Run one stage on the server and bring its outputs back.

    `inputs` are `local_path => remote_relative_path` pairs uploaded before the
    script runs; `outputs` are `remote_relative_path => local_path` pairs
    retrieved after it. `params` become `key=value` arguments to the script.
    Every remote path is formed under a staging directory created for this call
    alone and removed again afterwards, so two stages - or two runs of the same
    stage - never share state on the server.

    `optional_outputs` take the same form but are retrieved only if the script
    produced them. The diagnostic plots are conditional on there being data left
    to plot, and a local stage that skips one simply leaves no file; demanding it
    here would turn an empty sequence table from a reported result into a failed
    transfer.
    """
    function _run_remote_stage(emit, target, script::AbstractString,
                               params::AbstractVector{<:Pair};
                               inputs::AbstractVector{<:Pair},
                               outputs::AbstractVector{<:Pair},
                               optional_outputs::AbstractVector{<:Pair}=Pair{String,String}[],
                               verify_remote_gzip::AbstractVector{<:AbstractString}=String[],
                               verify_remote_sha256::AbstractDict=Dict{String,String}(),
                               log_path::AbstractString,
                               stage::AbstractString)
        scripts_dir = @__DIR__
        functions_r = joinpath(scripts_dir, "dada2_functions.r")
        script_r    = joinpath(scripts_dir, script)
        isfile(script_r) || error("Remote stage script not found: $script_r")
        uploads = Pair{String,String}[functions_r => "dada2_functions.r", script_r => script,
                                      (string(first(i)) => string(last(i)) for i in inputs)...]

        # A parameter may need to name a file uploaded into the staging
        # directory, whose path only the runner knows. Expanding a placeholder
        # keeps that path out of the callers' hands, where it could be used to
        # build one that escapes the directory.
        expand(v, staging_dir) = replace(string(v), "REMOTE_STAGING/" => "$staging_dir/")

        before_run = function (ssh, staging_dir)
            _verify_remote_gzip(ssh, target, [expand(p, staging_dir) for p in verify_remote_gzip], emit)
            _verify_remote_sha256(ssh, target, verify_remote_sha256, emit)
        end

        command_for = function (staging_dir)
            expanded = [k => expand(v, staging_dir) for (k, v) in params]
            # Parameters are interpolated into a command string that the remote
            # sshd runs through a login shell. Everything derived from a filename
            # travels in a manifest instead, so what remains here is config
            # values from a constrained vocabulary; refuse anything a shell would
            # word-split or expand rather than let it reach that interpolation.
            for (k, v) in expanded
                _remote_shell_safe(string(v)) ||
                    error("Remote stage parameter '$k' may not contain shell " *
                          "metacharacters or spaces (got: '$v')")
            end
            arg_str = join(["$k=$v" for (k, v) in expanded], " ")
            "$(target.rscript) $staging_dir/$script " *
                "functions=$staging_dir/dada2_functions.r " *
                "staging=$staging_dir " *
                arg_str
        end

        run_remote(emit, target, command_for;
                   inputs=uploads, outputs, optional_outputs, before_run, log_path, stage)
    end
