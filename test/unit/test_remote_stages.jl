# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Bioserver offload: where a stage runs, how many threads it runs with, and what
# reaches the remote command string. These are the decisions made before any SSH
# connection exists, so they are testable without a server; what happens on the
# far side of the connection belongs to the integration suite.

using CodecZlib

const _R = MetaManifold.DADA2

@testset "Remote stage offload" begin

    @testset "_r_threads reads the global setting" begin
        @test _R._r_threads(Dict{String,Any}()) == 4          # factory fallback
        @test _R._r_threads(Dict("r_threads" => 24)) == 24
        @test _R._r_threads(Dict("r_threads" => true)) === true
    end

    @testset "_r_threads honours the deprecated taxonomy override for that stage alone" begin
        cfg = Dict{String,Any}(
            "r_threads" => 8,
            "dada2" => Dict("taxonomy" => Dict("multithread" => 24)),
        )
        # A config written before the key moved keeps the count it asked for...
        @test _R._r_threads(cfg; stage="assign_taxonomy") == 24
        # ...but a key that meant "threads for assignTaxonomy" does not silently
        # take over denoising.
        @test _R._r_threads(cfg; stage="denoise") == 8
        @test _R._r_threads(cfg; stage="learn_errors") == 8
        @test _R._r_threads(cfg; stage="chimera_removal") == 8
    end

    @testset "_mt_str renders DADA2's union type" begin
        # A string reaches DADA2 as an invalid multithread and is silently
        # downgraded to one thread, so Bool and Integer must stay distinct.
        @test _R._mt_str(true)  == "TRUE"
        @test _R._mt_str(false) == "FALSE"
        @test _R._mt_str(24)    == "24"
    end

    @testset "_remote_target returns nothing without a configured host" begin
        for cfg in (Dict{String,Any}(),
                    Dict{String,Any}("remote" => Dict("host" => nothing)),
                    Dict{String,Any}("remote" => Dict("stages" => ["denoise"])))
            for stage in _R.REMOTE_STAGES
                @test isnothing(_R._remote_target(cfg, stage))
            end
        end
    end

    @testset "_remote_target offloads only the listed stages" begin
        cfg = Dict{String,Any}("remote" => Dict(
            "host" => "user@server", "staging_dir" => "/srv/staging",
            "stages" => ["learn_errors", "chimera_removal"]))

        # Naming a server does not move every stage onto it; each is opted in.
        @test !isnothing(_R._remote_target(cfg, "learn_errors"))
        @test !isnothing(_R._remote_target(cfg, "chimera_removal"))
        @test isnothing(_R._remote_target(cfg, "denoise"))
        @test isnothing(_R._remote_target(cfg, "assign_taxonomy"))

        t = _R._remote_target(cfg, "learn_errors")
        @test t.host == "user@server"
        @test t.base == "/srv/staging"
        @test t.rscript == "Rscript"
        @test isnothing(t.identity_file)
    end

    @testset "_remote_target is modular - each stage stands alone" begin
        # The point of the design: chimeras can be removed on the server for a
        # run denoised here, with nothing carried between the two.
        cfg = Dict{String,Any}("remote" => Dict(
            "host" => "user@server", "staging_dir" => "/srv/staging",
            "stages" => ["chimera_removal"]))
        @test isnothing(_R._remote_target(cfg, "denoise"))
        @test !isnothing(_R._remote_target(cfg, "chimera_removal"))
    end

    @testset "_remote_target threads fall back to r_threads" begin
        base = Dict("host" => "user@server", "staging_dir" => "/srv/staging",
                    "stages" => ["denoise"])
        cfg = Dict{String,Any}("r_threads" => 4, "remote" => base)
        @test _R._remote_target(cfg, "denoise").threads == 4

        # The server is usually wider than this machine, which is the reason to
        # offload at all, so it may be told to run wider.
        cfg2 = Dict{String,Any}("r_threads" => 4,
                                "remote" => merge(base, Dict("threads" => 64)))
        @test _R._remote_target(cfg2, "denoise").threads == 64
    end

    @testset "_remote_target keeps the deprecated taxonomy block working" begin
        # A config written before the block moved must behave exactly as it did,
        # including without any top-level remote block at all.
        cfg = Dict{String,Any}("dada2" => Dict("taxonomy" => Dict(
            "remote" => Dict("host" => "user@old", "staging_dir" => "/srv/old"))))
        t = _R._remote_target(cfg, "assign_taxonomy")
        @test !isnothing(t)
        @test t.host == "user@old"
        # And it stays a taxonomy-only setting: it never moves another stage.
        @test isnothing(_R._remote_target(cfg, "denoise"))
    end

    @testset "_remote_target refuses a staging_dir a shell would mangle" begin
        for bad in ("relative/path", "/srv/sta ging", "/srv/\$(whoami)", "/srv/a;rm -rf /")
            cfg = Dict{String,Any}("remote" => Dict(
                "host" => "user@server", "staging_dir" => bad, "stages" => ["denoise"]))
            @test_throws ErrorException _R._remote_target(cfg, "denoise")
        end
    end

    @testset "_remote_target requires a staging_dir once a host is named" begin
        cfg = Dict{String,Any}("remote" => Dict(
            "host" => "user@server", "stages" => ["denoise"]))
        @test_throws ErrorException _R._remote_target(cfg, "denoise")
    end

    @testset "_write_manifest writes one entry per line" begin
        # Read paths and sample names travel this way rather than as command-line
        # arguments, so that no filename is ever parsed by a remote login shell.
        dir = mktempdir()
        path = _R._write_manifest(dir, "fwd.txt", ["a_R1_filt.fastq.gz", "b_R1_filt.fastq.gz"])
        @test readlines(path) == ["a_R1_filt.fastq.gz", "b_R1_filt.fastq.gz"]
        # A name a shell would split survives intact, because nothing splits it.
        path2 = _R._write_manifest(dir, "s.txt", ["sample one", "sample;two"])
        @test readlines(path2) == ["sample one", "sample;two"]
        rm(dir; recursive=true)
    end

    @testset "_fastq_bases counts sequence bases exactly" begin
        # derepFastq's abundance-weighted total is the sum of the read lengths,
        # so a plain count of the sequence lines has to agree with it.
        dir = mktempdir()
        path = joinpath(dir, "s_R1_filt.fastq.gz")
        open(GzipCompressorStream, path, "w") do io
            for (i, seq) in enumerate(("ACGTACGTAC", "TTTT", "GGGGGGGGGGGGGGG"))
                println(io, "@read$i\n$seq\n+\n", repeat("I", length(seq)))
            end
        end
        @test _R._fastq_bases(path) == 10 + 4 + 15
        rm(dir; recursive=true)
    end

    @testset "_learn_errors_prefix matches dada2's own budget loop" begin
        # learnErrors dereplicates in order and breaks once the cumulative base
        # count EXCEEDS nbases, keeping drps[1:i] - so the file that crosses the
        # budget is included, and nothing after it is read. Sending more than
        # this prefix is wasted transfer; sending less would change the model.
        dir = mktempdir()
        paths = String[]
        for i in 1:4
            p = joinpath(dir, "s$(i)_R1_filt.fastq.gz")
            open(GzipCompressorStream, p, "w") do io
                # 100 bases per file, in one read.
                seq = repeat("ACGT", 25)
                println(io, "@r$i\n$seq\n+\n", repeat("I", length(seq)))
            end
            push!(paths, p)
        end

        # Budget met partway: the crossing file is included, the rest are not.
        @test _R._learn_errors_prefix(paths, 150) == paths[1:2]
        # Exactly on the boundary is not "exceeded", so one more file is taken.
        @test _R._learn_errors_prefix(paths, 200) == paths[1:3]
        # A budget the reads never reach uses every file, as dada2 would.
        @test _R._learn_errors_prefix(paths, 10_000) == paths
        @test _R._learn_errors_prefix(String[], 100) == String[]
        rm(dir; recursive=true)
    end

    @testset "thread count and remote host never mark a stage stale" begin
        # r_threads and the remote block change no result, so they are kept out
        # of every hashed section. If one leaks in, changing the thread count
        # would invite a re-run of hours of work to reproduce identical output.
        hashed = join([join(v, ",") for v in values(Config._STAGE_SECTIONS)], ",")
        for key in ("r_threads", "remote")
            @test !occursin(key, hashed)
        end

        dir = mktempdir()
        cfg_path = joinpath(dir, "run_config.yml")
        hash_file = joinpath(dir, "assign_taxonomy.hash")
        write(cfg_path, """
        r_threads: 4
        remote:
          host: ~
          stages: []
        dada2:
          taxonomy:
            database: pr2
            min_boot: 0
            enabled: true
          output:
            taxa_prefix: taxonomy
        """)
        section = Config.stage_sections(:dada2_assign_taxonomy)
        Config._write_section_hash(cfg_path, section, hash_file)
        @test !Config._section_stale(cfg_path, section, hash_file)

        # Raise the thread count and point at a server: the stage stays current.
        write(cfg_path, """
        r_threads: 24
        remote:
          host: "user@server"
          staging_dir: /srv/staging
          stages:
            - assign_taxonomy
        dada2:
          taxonomy:
            database: pr2
            min_boot: 0
            enabled: true
          output:
            taxa_prefix: taxonomy
        """)
        @test !Config._section_stale(cfg_path, section, hash_file)

        # A key that does change the assignment still invalidates it.
        write(cfg_path, """
        r_threads: 24
        dada2:
          taxonomy:
            database: silva
            min_boot: 0
            enabled: true
          output:
            taxa_prefix: taxonomy
        """)
        @test Config._section_stale(cfg_path, section, hash_file)
        rm(dir; recursive=true)
    end

    @testset "validation accepts the new global settings" begin
        cfg = Dict{String,Any}(
            "cutadapt"  => Dict("primer_pairs" => ["EMP"], "min_length" => 200),
            "r_threads" => 24,
            "remote"    => Dict("host" => "user@server",
                                "staging_dir" => "/srv/staging",
                                "threads" => 64,
                                "stages" => ["learn_errors", "assign_taxonomy"]),
        )
        errors = Validation.ValidationError[]
        Validation._validate_pipeline_cfg(errors, cfg, "test")
        @test isempty(errors)
    end

    @testset "validation rejects an unrunnable remote stage" begin
        # A stage name that is not offloadable would otherwise run locally, which
        # looks from the outside exactly like a working offload.
        cfg = Dict{String,Any}(
            "cutadapt" => Dict("primer_pairs" => ["EMP"], "min_length" => 200),
            "remote"   => Dict("host" => "user@server", "staging_dir" => "/srv/staging",
                               "stages" => ["filter_trim"]),
        )
        errors = Validation.ValidationError[]
        Validation._validate_pipeline_cfg(errors, cfg, "test")
        @test any(e -> occursin("remote.stages", e.message), errors)
    end

    @testset "validation rejects bad thread counts" begin
        for v in ("24", 4.5, 0, -1)
            cfg = Dict{String,Any}(
                "cutadapt"  => Dict("primer_pairs" => ["EMP"], "min_length" => 200),
                "r_threads" => v)
            errors = Validation.ValidationError[]
            Validation._validate_pipeline_cfg(errors, cfg, "test")
            @test any(e -> occursin("r_threads", e.message), errors)
        end
    end

    @testset "validation gates the staging_dir that reaches the ssh command" begin
        for bad in ("relative/path", "/srv/sta ging", "/srv/\$(whoami)")
            cfg = Dict{String,Any}(
                "cutadapt" => Dict("primer_pairs" => ["EMP"], "min_length" => 200),
                "remote"   => Dict("host" => "user@server", "staging_dir" => bad,
                                   "stages" => ["denoise"]))
            errors = Validation.ValidationError[]
            Validation._validate_pipeline_cfg(errors, cfg, "test")
            @test any(e -> occursin("staging_dir", e.message), errors)
        end
    end

    @testset "validation leaves the factory placeholder alone" begin
        # The shipped default names no host, so its placeholder staging_dir must
        # not make every unconfigured project invalid.
        factory = YAML.load_file(joinpath(@__DIR__, "..", "..", "config",
                                          "defaults", "pipeline.yml"))
        errors = Validation.ValidationError[]
        Validation._validate_pipeline_cfg(errors, factory, "factory")
        @test isempty(errors)
    end

    @testset "every offloadable stage has a remote script" begin
        # REMOTE_STAGES is what a config may name; a stage listed there with no
        # script behind it would fail only once someone offloaded it.
        scripts = Dict("learn_errors"    => "learn_errors_remote.r",
                       "denoise"         => "denoise_remote.r",
                       "chimera_removal" => "chimera_remote.r",
                       "assign_taxonomy" => "taxonomy_remote.r")
        script_dir = joinpath(@__DIR__, "..", "..", "src", "pipeline", "dada2")
        @test Set(keys(scripts)) == Set(Validation.DADA2_REMOTE_STAGES)
        for name in values(scripts)
            @test isfile(joinpath(script_dir, name))
        end
    end
end
