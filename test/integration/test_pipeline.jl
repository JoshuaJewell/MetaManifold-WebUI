# Integration test: full pipeline run on the DADA2 MiSeq SOP dataset.
#
# Two runs committed under data/MiSeq_SOP/:
#   run_A  - F3D0 (day 0) + F3D1 (day 1)  - 2 samples
#   run_B  - F3D2 (day 2)                  - 1 sample
#
# Pipeline settings are in data/MiSeq_SOP/pipeline.yml (study-level).
# Primers are already removed; cutadapt passes reads through unchanged.
#
# Coverage:
#   - All pipeline stages through merge_taxa (both runs)
#   - DuckDB results loading and query helpers
#   - The differential abundance route on the two runs
#
# Requirements:
#   - External tools installed (cutadapt, vsearch, swarm, cd-hit-est)
#   - At least one taxonomy database downloaded (pr2 or silva)
#   - Julia started with multiple threads (e.g. julia -t4)
#
# The test is skipped gracefully if tools or databases are unavailable.

@testset "Integration: MiSeq SOP pipeline" begin

    PROJECT_ROOT = joinpath(@__DIR__, "..", "..")
    PROJECT_NAME = "MiSeq_SOP"

    ## Check prerequisites; fall back to committed CI configs if real ones absent
    ci_config_dir = joinpath(PROJECT_ROOT, "config", "ci")
    tools_config = let p = joinpath(PROJECT_ROOT, "config", "tools.yml")
        isfile(p) ? p : joinpath(ci_config_dir, "tools.yml")
    end
    db_config = let p = joinpath(PROJECT_ROOT, "config", "databases.yml")
        isfile(p) ? p : joinpath(ci_config_dir, "databases.yml")
    end

    if !isfile(tools_config) || !isfile(db_config)
        @warn "Integration test skipped: tools.yml or databases.yml not found"
        return
    end

    tcfg = YAML.load_file(tools_config)
    vsearch_path = let e = get(tcfg, "vsearch", nothing)
        e isa Dict ? get(e, "path", "vsearch") : something(e, "vsearch")
    end
    if isnothing(Sys.which(string(vsearch_path))) && !isfile(string(vsearch_path))
        @warn "Integration test skipped: vsearch not found at '$vsearch_path'"
        return
    end

    ## Clean previous test outputs
    test_study_dir = joinpath(PROJECT_ROOT, "projects", PROJECT_NAME)
    if isdir(test_study_dir)
        for entry in readdir(test_study_dir; join=true)
            isdir(entry) || continue
            for subdir in ("QC", "cutadapt", "dada2", "swarm", "vsearch", "merged")
                d = joinpath(entry, subdir)
                isdir(d) && rm(d; recursive=true)
            end
        end
    end

    local projects
    try
        projects = new_project(PROJECT_NAME;
                               data_dir     = joinpath(PROJECT_ROOT, "data"),
                               projects_dir = joinpath(PROJECT_ROOT, "projects"),
                               config_dir   = joinpath(PROJECT_ROOT, "config"))
    catch e
        @warn "Integration test skipped: could not create project - $e"
        return
    end

    @test length(projects) == 2

    db_name    = "pr2"
    dbs_loaded = ensure_databases(db_config)
    dada2_key  = "$(db_name)_dada2"
    vsearch_key = "$(db_name)_vsearch"

    if !haskey(dbs_loaded, dada2_key)
        @warn "Integration test: DADA2 stage skipped - database '$dada2_key' not available"
        return
    end
    if !haskey(dbs_loaded, vsearch_key)
        @warn "Integration test: vsearch stage skipped - '$vsearch_key' not available"
        return
    end

    db_meta       = make_db_meta(db_config, db_name)
    r_lock_int    = ReentrantLock()
    skip_taxonomy = haskey(ENV, "CI_SKIP_TAXONOMY")

    merged_results = Vector{Union{MergedTables,Nothing}}(undef, length(projects))
    asvs_results   = Vector{Union{ASVResult,Nothing}}(undef, length(projects))

    ## Run pipeline for each project
    Threads.@threads for (i, project) in collect(enumerate(projects))
        n_samples = length(filter(f -> endswith(f, ".fastq.gz"), readdir(project.data_dir))) ÷ 2

        errors = validate_project(project, db_config)
        if !isempty(errors)
            @warn "Integration test: run $(basename(project.dir)) skipped - validation failed"
            continue
        end

        # Stage 0: FastQC/MultiQC over the raw reads. No later stage reads its output,
        # so it is checked here.
        multiqc(project)
        qc_dir = joinpath(project.dir, "QC")
        @test isfile(joinpath(qc_dir, "multiqc_report.html"))
        for logname in ("fastqc.log", "multiqc.log")
            p = joinpath(qc_dir, "logs", logname)
            @test isfile(p)
            # A stage that dies after creating its log leaves it empty; the command
            # marker is the provenance record.
            @test filesize(p) > 0
            @test any(l -> startswith(l, "[MetaManifold] cmd: "), readlines(p))
        end

        # Stage 1: cutadapt
        trimmed = cutadapt(project)
        @test isdir(trimmed.dir)
        trimmed_files = filter(f -> endswith(f, "_trimmed.fastq.gz"), readdir(trimmed.dir))
        @test length(trimmed_files) == n_samples * 2

        # Stage 2: DADA2 (then cd-hit), followed by SWARM (needs ASVs as input)
        asvs = lock(r_lock_int) do
            dada2(project, trimmed, taxonomy_db=dbs_loaded[dada2_key],
                  skip_taxonomy=skip_taxonomy)
        end
        asvs = cdhit(project, asvs)
        otus = swarm(project, asvs)

        @test isfile(asvs.fasta) && filesize(asvs.fasta) > 0
        @test isfile(asvs.count_table) && filesize(asvs.count_table) > 0

        # Provenance: each DADA2 R stage records the exact call it issued, behind
        # the same marker the shell tools use, recoverable by a single grep.
        dada2_logs = joinpath(project.dir, "dada2", "Logs")
        for stage_log in ("filter_trim.log", "learn_errors.log", "denoise.log",
                          "chimera_removal.log")
            p = joinpath(dada2_logs, stage_log)
            @test isfile(p)
            @test any(l -> startswith(l, "[MetaManifold] cmd: R> "), readlines(p))
        end

        # Stage 3: vsearch taxonomy (ASV + OTU in parallel)
        asv_tax_task = Threads.@spawn vsearch(project, asvs, dbs_loaded[vsearch_key])
        otu_tax_task = isnothing(otus) ? nothing :
                       Threads.@spawn vsearch(project, otus, dbs_loaded[vsearch_key])

        asv_tax = fetch(asv_tax_task)
        @test isfile(asv_tax.tsv) && filesize(asv_tax.tsv) > 0

        # Stage 4: merge_taxa
        asv_merged = merge_taxa(project, asvs, asv_tax, db_meta)
        merged = if !isnothing(otus)
            otu_tax    = fetch(otu_tax_task)
            otu_merged = merge_taxa_otu(project, otus, otu_tax, db_meta)
            MergedTables(merge(asv_merged.tables, otu_merged.tables))
        else
            asv_merged
        end

        @test haskey(merged.tables, "merged")
        merged_csv = merged.tables["merged"]
        @test isfile(merged_csv) && filesize(merged_csv) > 0

        df = CSV.read(merged_csv, DataFrames.DataFrame)
        @test nrow(df) > 0
        @test "SeqName" in DataFrames.names(df)

        # Verify sample columns exist by checking non-taxonomy numeric columns
        non_count = Set(["SeqName", "Pident", "Accession", "rRNA", "Organellum", "specimen",
                         "sequence", "OTU", "ASV",
                         "Domain", "Supergroup", "Division", "Subdivision",
                         "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])
        sample_cols = [c for c in DataFrames.names(df)
                       if !(c in non_count) && !endswith(c, "_dada2") &&
                          !endswith(c, "_boot") && !endswith(c, "_vsearch") &&
                          eltype(df[!, c]) <: Union{Missing, Number}]
        @test length(sample_cols) == n_samples
        total = sum(sum(skipmissing(df[!, c])) for c in sample_cols)
        @test total > 0

        merged_results[i] = merged
        asvs_results[i]   = asvs

        @info "Pipeline [$(basename(project.dir))]: $(nrow(df)) ASVs, $total reads across $n_samples samples"

        ## Stage 5: DuckDB results loading
        merge_dir = joinpath(project.dir, "merged")
        swarm_dir = joinpath(project.dir, "swarm")
        db_path = load_results_db(merge_dir;
            swarm_dir = isdir(swarm_dir) ? swarm_dir : nothing)
        @test isfile(db_path)

        # Verify DuckDB query helpers work
        with_results_db(merge_dir) do con
            scols = Analysis.sample_columns(con, "merged")
            @test length(scols) == n_samples

            levels = Analysis.taxonomy_levels(con, "merged")
            @test !isempty(levels)
            @test "Domain" in levels
        end
    end

    ## Differential abundance route on the two runs. run_A has 2 samples and
    # run_B 1, so this checks only that the route answers with a well-formed
    # result or an explicit error, never that anything is significant.
    if all(i -> isassigned(merged_results, i) && !isnothing(merged_results[i]), eachindex(merged_results))
        @testset "Differential abundance route (run_A vs run_B)" begin
            SV = MetaManifold.Server
            old_root = SV.ServerState._root[]
            try
                SV.ServerState.set_root!(PROJECT_ROOT)
                body = JSON3.read(JSON3.write(Dict(
                    "runs" => [Dict("run" => "run_A"), Dict("run" => "run_B")],
                    "table" => "merged", "rank" => "Genus")))
                r = SV._differential(PROJECT_NAME, body)
                out = JSON3.read(String(r.body))
                if r.status == 200
                    @test haskey(out, :rows) && haskey(out, :figure) && haskey(out, :diagnostics)
                    @test out.n_samples.run_A == 2 && out.n_samples.run_B == 1
                    @test all(x -> x.status in ("ok", "boundary", "failed", "filtered"), out.rows)
                else
                    @test r.status in (400, 404, 422, 503)
                    @test haskey(out, :error) && haskey(out, :message)
                    @info "Differential route on MiSeq_SOP answered $(r.status) $(out.error): $(out.message)"
                end
            finally
                SV.ServerState._root[] = old_root
            end
        end
    end

    ## The server's own entry points, on the run directories the loop above built.
    # Nothing up to here has called them: the loop drives the stages by hand, which
    # is not what pressing a stage in the UI does. The full pipeline runs the
    # vsearch stage's two searches concurrently and a single-stage re-run runs them
    # one after the other (_vsearch_searches) - the two shapes whose peaks the
    # occupancy probe exists to tell apart - and these calls are the only callers of
    # that helper outside test/occupancy_disagreement_test.jl.
    #
    # Each vsearch output is removed first. A stage whose outputs are current skips
    # before doing any work, and a stage that skips records nothing: a re-run that
    # skipped would exercise the skip guard rather than the stage. The stages before
    # vsearch still skip on their own guards, so neither call redoes the DADA2 work.
    #
    # The probe is switched on for both, at the narrowest level of the config
    # cascade (the run's own pipeline.yml) - the way the factory defaults file says
    # to switch it on, since the config editor offers only keys that have a factory
    # default. It sits in no stage's hashed section, so turning it on marks nothing
    # stale, and it is on only here: every other call in the suite records exactly
    # the stage keys it always did.
    @testset "Pipeline and stage routes (vsearch re-run)" begin
        SV = MetaManifold.Server
        PV = MetaManifold.Provenance
        MP = MetaManifold.MemoryProbe
        run_dir(run) = joinpath(PROJECT_ROOT, "projects", PROJECT_NAME, run)
        probe_cfg = Dict("occupancy_probe" => Dict("enabled" => true, "interval_ms" => 10))
        old_root  = SV.ServerState._root[]
        original_configs = Dict{String,Union{Nothing,Vector{UInt8}}}()
        try
            SV.ServerState.set_root!(PROJECT_ROOT)

            # run_A: one stage on its own - the serial shape.
            rm(joinpath(run_dir("run_A"), "vsearch"); recursive=true, force=true)
            rm(joinpath(run_dir("run_A"), "swarm", "vsearch"); recursive=true, force=true)
            yaml_a = joinpath(PROJECT_ROOT, "data", PROJECT_NAME, "run_A", "pipeline.yml")
            original_configs[yaml_a] = isfile(yaml_a) ? read(yaml_a) : nothing
            YAML.write_file(yaml_a, probe_cfg)
            response = SV.Oxygen.internalrequest(SV.HTTP.Request(
                "POST", "/api/v1/studies/$PROJECT_NAME/runs/run_A/stages/vsearch"))
            @test response.status == 200
            job = SV.JobQueue.get_job(String(JSON3.read(String(response.body)).id))
            @test !isnothing(job)
            task = lock(SV.JobQueue._jobs_lock) do
                job.task
            end
            isnothing(task) || wait(task)
            @test job.status == SV.JobQueue.complete
            @test isfile(joinpath(run_dir("run_A"), "vsearch", "taxonomy.tsv"))
            @test isfile(joinpath(run_dir("run_A"), "swarm", "vsearch", "taxonomy.tsv"))

            # run_B: the whole pipeline - the concurrent shape.
            rm(joinpath(run_dir("run_B"), "vsearch"); recursive=true, force=true)
            rm(joinpath(run_dir("run_B"), "swarm", "vsearch"); recursive=true, force=true)
            yaml_b = joinpath(PROJECT_ROOT, "data", PROJECT_NAME, "run_B", "pipeline.yml")
            original_configs[yaml_b] = isfile(yaml_b) ? read(yaml_b) : nothing
            YAML.write_file(yaml_b, probe_cfg)
            merged = SV._run_full_pipeline(PROJECT_NAME, "run_B", db_config, dbs_loaded)
            @test haskey(merged.tables, "merged")
            @test isfile(joinpath(run_dir("run_B"), "vsearch", "taxonomy.tsv"))
            @test isfile(joinpath(run_dir("run_B"), "swarm", "vsearch", "taxonomy.tsv"))

            # What each run recorded: its vsearch stage, and - the probe being on -
            # that stage's occupancy, one labelled task per search.
            for run in ("run_A", "run_B")
                att   = PV.read_attestation(joinpath(run_dir(run), "attestation.yml"))
                stage = only([s for s in att.doc["stages"] if s["name"] == "vsearch"])
                @test haskey(stage, "occupancy")
                occ = get(stage, "occupancy", nothing)
                isnothing(occ) && continue
                @test occ["probe"] == MP.PROBE
                @test occ["interval_ms"] == 10
                labels = sort([String(t["label"]) for t in occ["tasks"]])
                @test labels == ["vsearch:asv", "vsearch:otu"]
            end
        finally
            for (path, contents) in original_configs
                isnothing(contents) ? rm(path; force=true) : write(path, contents)
            end
            SV.ServerState._root[] = old_root
        end
    end

end
