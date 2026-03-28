@testset "PipelineLog" begin

    @testset "reset_log creates log file" begin
        dir = mktempdir()
        PipelineLog.reset_log(dir)
        log_path = joinpath(dir, "pipeline.log")
        @test isfile(log_path)
        content = read(log_path, String)
        @test startswith(content, "Pipeline log started")
        rm(dir; recursive=true)
    end

    @testset "pipeline_log appends timestamped messages" begin
        dir = mktempdir()
        PipelineLog.reset_log(dir)
        PipelineLog.pipeline_log(dir, "Test message 1")
        PipelineLog.pipeline_log(dir, "Test message 2")
        content = read(joinpath(dir, "pipeline.log"), String)
        @test occursin("Test message 1", content)
        @test occursin("Test message 2", content)
        lines = split(strip(content), '\n')
        @test length(lines) == 3  # header + 2 messages
        rm(dir; recursive=true)
    end

    @testset "pipeline_log with ProjectCtx" begin
        dir = mktempdir()
        ctx = ProjectCtx(dir, "", "", "", "")
        PipelineLog.reset_log(ctx)
        PipelineLog.pipeline_log(ctx, "via context")
        content = read(joinpath(dir, "pipeline.log"), String)
        @test occursin("via context", content)
        rm(dir; recursive=true)
    end

    @testset "log_written records SHA-256" begin
        dir = mktempdir()
        PipelineLog.reset_log(dir)
        test_file = joinpath(dir, "test_output.csv")
        write(test_file, "hello world\n")
        PipelineLog.log_written(dir, test_file)
        content = read(joinpath(dir, "pipeline.log"), String)
        @test occursin("Written:", content)
        @test occursin("sha256:", content)
        # SHA-256 hex is 64 chars
        m = match(r"sha256:([0-9a-f]+)", content)
        @test !isnothing(m)
        @test length(m[1]) == 64
        rm(dir; recursive=true)
    end

    @testset "finalise_log appends tool logs" begin
        dir = mktempdir()
        PipelineLog.reset_log(dir)
        PipelineLog.pipeline_log(dir, "Stage complete")

        # Create a fake tool log
        mkpath(joinpath(dir, "cutadapt", "logs"))
        write(joinpath(dir, "cutadapt", "logs", "cutadapt_primer_trimming_stats.txt"),
              "Trimmed: 100 reads\nRetained: 95 reads\n")

        ctx = ProjectCtx(dir, "", "", "", "")
        PipelineLog.finalise_log(ctx)

        content = read(joinpath(dir, "pipeline.log"), String)
        @test occursin("Tool logs appended", content)
        @test occursin("Trimmed: 100 reads", content)
        @test occursin("End of tool logs", content)
        rm(dir; recursive=true)
    end

    @testset "write_combined_log merges run logs" begin
        study_dir = mktempdir()
        run1_dir = joinpath(study_dir, "run1")
        run2_dir = joinpath(study_dir, "run2")
        mkpath(run1_dir)
        mkpath(run2_dir)

        PipelineLog.reset_log(run1_dir)
        PipelineLog.pipeline_log(run1_dir, "Run 1 done")
        PipelineLog.reset_log(run2_dir)
        PipelineLog.pipeline_log(run2_dir, "Run 2 done")

        ctx1 = ProjectCtx(run1_dir, "", "", study_dir, "")
        ctx2 = ProjectCtx(run2_dir, "", "", study_dir, "")

        path = PipelineLog.write_combined_log([ctx1, ctx2]; study_dir)
        @test isfile(path)
        content = read(path, String)
        @test occursin("Combined pipeline log", content)
        @test occursin("Run 1 done", content)
        @test occursin("Run 2 done", content)
        @test occursin("RUN: run1", content)
        @test occursin("RUN: run2", content)
        rm(study_dir; recursive=true)
    end

end
