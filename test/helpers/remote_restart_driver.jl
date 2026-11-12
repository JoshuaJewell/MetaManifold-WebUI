# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

include("remote_protocol_module.jl")
using JSON3
using UUIDs

settings = JSON3.read(read(ARGS[1], String))
target = (; host=String(settings.host), base=String(settings.base),
           identity_file=isnothing(settings.identity_file) ? nothing : String(settings.identity_file),
           threads=2)
connection = if haskey(settings, :ssh_config)
    ctl = joinpath(tempdir(), "mm_$(uuid4())")
    write(joinpath(settings.local_dir, "control_path"), ctl)
    open(io -> println(io, ctl), joinpath(settings.local_dir, "control_paths"), "a")
    ssh, scp = RemoteProtocol.RemoteExec._transport(target, ctl)
    configured = cmd -> Cmd([first(cmd.exec), "-F", String(settings.ssh_config), cmd.exec[2:end]...])
    upload = function (src, dst)
        if !startswith(src, target.host * ":")
            open(io -> println(io, src), joinpath(settings.local_dir, "uploads"), "a")
        end
        configured(scp(src, dst))
    end
    ((args...) -> configured(ssh(args...)), upload)
elseif target.host == "fixture"
    ssh = (_, command) -> `bash -c $command`
    scp = function (src, dst)
        if startswith(src, "fixture:")
            `cp -- $(src[9:end]) $dst`
        else
            open(io -> println(io, src), joinpath(settings.local_dir, "uploads"), "a")
            `cp -- $src $(dst[9:end])`
        end
    end
    (ssh, scp)
else
    nothing
end
command = d -> "cd '$d' && echo start >> starts && " *
    "while ! test -f release; do sleep 0.05; done; cat input > output"
emit = message -> (println(message); flush(stdout))
RemoteProtocol.RemoteExec.run_remote(emit, target, command;
    inputs=[joinpath(settings.local_dir, "input") => "input"],
    outputs=["output" => joinpath(settings.local_dir, "stage.out")],
    log_path=joinpath(settings.local_dir, "stage.log"), stage="denoise",
    _connection=connection, _poll_interval=0.02)
