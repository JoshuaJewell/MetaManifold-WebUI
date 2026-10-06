# Unit tests for the write gates that keep request values out of file paths:
# config values used as file-name parts or remote command words, and the
# pipeline tables the results routes may not delete. Requests run in-process
# through Oxygen's router.
using MetaManifold
SV = MetaManifold.Server
using JSON3, YAML, DuckDB, DBInterface

_wg_request(method, path, body="") =
    SV.Oxygen.internalrequest(SV.HTTP.Request(method, path, ["Content-Type" => "application/json"], body))

@testset "Config write gate" begin
    root = mktempdir()
    old_root = SV.ServerState._root[]
    SV.ServerState.set_root!(root)
    mkpath(joinpath(root, "data"))
    mkpath(joinpath(root, "config"))
    user_cfg = joinpath(root, "config", "pipeline.yml")
    write(user_cfg, "vsearch:\n  identity: 0.97\n")
    try
        bad = [
            "dada2.output.fasta_prefix"     => "../../tmp/pwn",
            "dada2.output.fasta_prefix"     => "a/b",
            "dada2.output.fasta_prefix"     => "",
            "dada2.output.seq_table_prefix" => "../seqtab",
            "dada2.output.taxa_prefix"      => "/tmp/taxonomy",
            "dada2.output.taxa_prefix"      => "tax onomy",
            "cutadapt.r1_suffix"            => "_R1\\E.*",
            "cutadapt.r2_suffix"            => "\\Q_R2",
            "cutadapt.r2_suffix"            => "_R2/x",
            "remote.staging_dir"            => "/srv/mm; rm -rf ~",
            "remote.staging_dir"            => "relative/dir",
            "remote.staging_dir"            => "/srv/m m",
            "remote.stages"                 => ["denoise", "run_anything"],
            "remote.threads"                => 0,
            "remote.threads"                => "4",
            "remote.tools.mafft"            => "mafft;id",
            "remote.tools.iqtree"           => nothing,
        ]
        before = read(user_cfg, String)
        for (k, v) in bad
            r = _wg_request("PATCH", "/api/v1/config", JSON3.write(Dict(k => v)))
            @test r.status == 400
            @test read(user_cfg, String) == before
        end

        good = [
            "dada2.output.fasta_prefix"     => "asvs_v2",
            "dada2.output.seq_table_prefix" => "seqtab.nochim",
            "dada2.output.taxa_prefix"      => "taxonomy-silva",
            "cutadapt.r1_suffix"            => "_1",
            "cutadapt.r2_suffix"            => "_2",
            "remote.host"                   => "user@server",
            "remote.staging_dir"            => "/scratch/metamanifold",
            "remote.stages"                 => ["denoise", "phylogeny_tree"],
            "remote.threads"                => 16,
            "remote.tools.mafft"            => "/opt/mafft/bin/mafft",
        ]
        for (k, v) in good
            r = _wg_request("PATCH", "/api/v1/config", JSON3.write(Dict(k => v)))
            @test r.status == 200
        end
        saved = YAML.load_file(user_cfg)
        @test saved["dada2"]["output"]["fasta_prefix"] == "asvs_v2"
        @test saved["cutadapt"]["r2_suffix"] == "_2"
        @test saved["remote"]["staging_dir"] == "/scratch/metamanifold"
        @test saved["remote"]["stages"] == ["denoise", "phylogeny_tree"]
        @test saved["remote"]["threads"] == 16
        @test saved["remote"]["tools"]["mafft"] == "/opt/mafft/bin/mafft"
        @test saved["vsearch"]["identity"] == 0.97

        # The gate's remote checks are validation's own: what it stored passes them.
        errs = SV.Validation.ValidationError[]
        SV.Validation._validate_remote!(errs, saved["remote"], "test")
        @test isempty(errs)

        # Naming a host without a staging_dir, or with one that is not a string,
        # leaves the remote stages nowhere to stage, and validation says so.
        for remote in (Dict{String,Any}("host" => "user@server"),
                       Dict{String,Any}("host" => "user@server", "staging_dir" => 42))
            errs = SV.Validation.ValidationError[]
            SV.Validation._validate_remote!(errs, remote, "test")
            @test any(e -> occursin("remote.staging_dir must be set", e.message), errs)
        end
    finally
        SV.ServerState._root[] = old_root
        rm(root; recursive=true, force=true)
    end
end

@testset "Pipeline tables cannot be deleted" begin
    root = mktempdir()
    old_root = SV.ServerState._root[]
    SV.ServerState.set_root!(root)
    data_run = joinpath(root, "data", "studyW", "runA")
    mkpath(data_run)
    touch(joinpath(data_run, "s1_R1.fastq.gz"))
    merge_dir = joinpath(root, "projects", "studyW", "runA", "merged")
    mkpath(merge_dir)
    db_path = joinpath(merge_dir, "results.duckdb")
    db = DuckDB.DB(db_path)
    con = DBInterface.connect(db)
    for t in (SV._PIPELINE_TABLES..., "my_view")
        DBInterface.execute(con, "CREATE TABLE \"$t\" (\"SeqName\" VARCHAR, \"s1\" BIGINT)")
        DBInterface.execute(con, "INSERT INTO \"$t\" VALUES ('seq1', 3)")
        write(joinpath(merge_dir, "$t.csv"), "SeqName,s1\nseq1,3\n")
    end
    DBInterface.close!(con); close(db)
    tables() = let db = DuckDB.DB(db_path), con = DBInterface.connect(db)
        try
            Set(string.(DBInterface.execute(con, "SHOW TABLES") |> DataFrame).name)
        finally
            DBInterface.close!(con); close(db)
        end
    end
    try
        @test Set(SV._PIPELINE_TABLES) == Set(["merged", "merged_otu", "merged_cdhit", "cluster_membership"])
        for t in SV._PIPELINE_TABLES
            r = _wg_request("DELETE", "/api/v1/studies/studyW/runs/runA/results/tables/$t")
            @test r.status == 400
            @test JSON3.read(String(r.body)).error == "reserved_name"
            @test isfile(joinpath(merge_dir, "$t.csv"))
            @test t in tables()
            # The save route refuses the same names.
            r = _wg_request("POST", "/api/v1/studies/studyW/runs/runA/results/tables/my_view/save",
                            JSON3.write(Dict("name" => t)))
            @test r.status == 400
        end
        # A user's own table still deletes.
        r = _wg_request("DELETE", "/api/v1/studies/studyW/runs/runA/results/tables/my_view")
        @test r.status == 200
        @test !isfile(joinpath(merge_dir, "my_view.csv"))
        @test !("my_view" in tables())
    finally
        SV.ServerState._root[] = old_root
        rm(root; recursive=true, force=true)
    end
end
