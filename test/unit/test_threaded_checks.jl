# SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
# SPDX-License-Identifier: MPL-2.0

# threaded_checks must count every check made on a thread, and a failure there must
# fail the run. A failing check cannot sit in this suite, so each case runs in a Julia
# of its own and is judged by what a user or CI sees: the exit code and the summary.

@testset "threaded_checks" begin
    helper = abspath(joinpath(@__DIR__, "..", "integration", "threaded_checks.jl"))
    script = """
    include($(repr(helper)))
    mode = only(ARGS)
    @testset "outer" begin
        @test true
        started = Threads.Atomic{Int}(0)
        threaded_checks(["a", "b"]) do i, item
            Threads.atomic_add!(started, 1)
            # Each item waits for the other, so this holds only if both run at once.
            @test timedwait(() -> started[] == 2, 60.0; pollint=0.01) === :ok
            mode == "error" && item == "a" && error("deliberate")
            @test !(mode == "fail" && item == "a")
            @test item in ("a", "b")
        end
        @test true
        @testset "after the threads" begin
            @test true
        end
    end
    """

    # Test is the only package needed, so the suite's flags are overridden rather than
    # compiling it again under --compiled-modules=no.
    function outcome(mode)
        mktempdir() do dir
            path = joinpath(dir, "script.jl")
            write(path, script)
            cmd = `$(Base.julia_cmd()) --startup-file=no --compiled-modules=yes
                   --code-coverage=none --color=no -t 2 $path $mode`
            out = IOBuffer()
            p = run(pipeline(ignorestatus(cmd); stdout=out, stderr=out))
            p.exitcode, String(take!(out))
        end
    end

    @testset "every check is counted" begin
        code, out = outcome("pass")
        @test code == 0
        @test occursin(r"outer\s+\|\s+9\s+9\b", out)
    end

    @testset "a failing check fails the run and stops nothing else" begin
        code, out = outcome("fail")
        @test code == 1
        @test occursin(r"outer\s+\|\s+8\s+1\s+9\b", out)
        @test occursin("after the threads", out)
    end

    @testset "an exception in one item is an error, and the other item finishes" begin
        code, out = outcome("error")
        @test code == 1
        @test occursin(r"outer\s+\|\s+7\s+1\s+8\b", out)
        @test occursin("deliberate", out)
    end
end
