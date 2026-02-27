@testset "Tools - get_primer_args mode dispatch" begin
    # Use the real primers.yml — always present in the repo
    primers_path = joinpath(@__DIR__, "..", "..", "config", "primers.yml")

    @testset "paired mode emits -g and -G" begin
        args = Tools.get_primer_args(["EMP"], primers_path; mode="paired")
        @test occursin("-g ", args)
        @test occursin("-G ", args)
    end

    @testset "forward mode emits only -g, no -G" begin
        args = Tools.get_primer_args(["EMP"], primers_path; mode="forward")
        @test occursin("-g ", args)
        @test !occursin("-G", args)
    end

    @testset "reverse mode emits only -g for reverse primer, not forward" begin
        data = YAML.load_file(primers_path)
        rev_seq = data["Reverse"]["EMP806R"]
        args = Tools.get_primer_args(["EMP"], primers_path; mode="reverse")
        @test occursin("-g $rev_seq", args)
        @test !occursin("-G", args)
        fwd_seq = data["Forward"]["EMP515F"]
        @test !occursin(fwd_seq, args)
    end

    @testset "default mode is paired" begin
        args_default = Tools.get_primer_args(["EMP"], primers_path)
        args_paired  = Tools.get_primer_args(["EMP"], primers_path; mode="paired")
        @test args_default == args_paired
    end
end
