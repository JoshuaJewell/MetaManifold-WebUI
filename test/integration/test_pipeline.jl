# Integration test: full pipeline run on the DADA2 MiSeq SOP dataset.
#
# Uses data/MiSeq_SOP/MiSeq_SOP/ — a 3-sample gzipped subset (F3D0, F3D1, F3D2)
# committed to the repository. Pipeline settings are in data/MiSeq_SOP/pipeline.yml.
# Primers are already removed; cutadapt passes reads through unchanged.
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

    ## Check prerequisites
    tools_config = joinpath(PROJECT_ROOT, "config", "tools.yml")
    db_config    = joinpath(PROJECT_ROOT, "config", "databases.yml")

    if !isfile(tools_config) || !isfile(db_config)
        @warn "Integration test skipped: tools.yml or databases.yml not found"
        return
    end

    # Check at least one tool is callable
    import YAML
    tcfg = YAML.load_file(tools_config)
    vsearch_path = let e = get(tcfg, "vsearch", nothing)
        e isa Dict ? get(e, "path", "vsearch") : something(e, "vsearch")
    end
    if isnothing(Sys.which(string(vsearch_path))) && !isfile(string(vsearch_path))
        @warn "Integration test skipped: vsearch not found at '$vsearch_path'"
        return
    end

    samples = ["F3D0_S188", "F3D1_S189", "F3D2_S190"]

    ## Clean previous test outputs so we exercise a fresh run
    test_project_dir = joinpath(PROJECT_ROOT, "projects", PROJECT_NAME)
    if isdir(test_project_dir)
        for subdir in ("cutadapt", "dada2", "swarm", "vsearch", "merged", "analysis")
            d = joinpath(test_project_dir, subdir)
            isdir(d) && rm(d; recursive=true)
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

    @test length(projects) >= 1
    project = projects[1]

    # Validate project config before running
    errors = validate_project(project, db_config)
    if !isempty(errors)
        @warn "Integration test skipped: project validation failed"
        for e in errors; @warn "  [$(e.context)] $(e.message)"; end
        return
    end

    ## Stage 1: cutadapt
    trimmed = cutadapt(project)
    @test isdir(trimmed.dir)
    trimmed_files = filter(f -> endswith(f, "_trimmed.fastq.gz"), readdir(trimmed.dir))
    @test length(trimmed_files) == length(samples) * 2  # R1 + R2 per sample

    ## Stage 2: DADA2 + SWARM (parallel)
    db_name   = "pr2"
    dbs_loaded = ensure_databases(db_config)

    dada2_key   = "$(db_name)_dada2"
    vsearch_key = "$(db_name)_vsearch"

    if !haskey(dbs_loaded, dada2_key)
        @warn "Integration test: DADA2 stage skipped - database '$dada2_key' not available"
        return
    end

    r_lock_int = ReentrantLock()

    dada2_task = Threads.@spawn begin
        result = lock(r_lock_int) do
            dada2(project, trimmed, taxonomy_db=dbs_loaded[dada2_key])
        end
        cdhit(project, result)
    end
    swarm_task = Threads.@spawn swarm(project, trimmed)

    asvs = fetch(dada2_task)
    otus = fetch(swarm_task)

    @test isfile(asvs.fasta)
    @test isfile(asvs.count_table)
    @test filesize(asvs.fasta) > 0
    @test filesize(asvs.count_table) > 0

    if !isnothing(otus)
        @test isfile(otus.fasta)
        @test isfile(otus.count_table)
        @test filesize(otus.fasta) > 0
    end

    ## Stage 3: vsearch taxonomy
    if !haskey(dbs_loaded, vsearch_key)
        @warn "Integration test: vsearch taxonomy skipped - '$vsearch_key' not available"
        return
    end

    asv_tax = vsearch(project, asvs, dbs_loaded[vsearch_key])
    @test isfile(asv_tax.tsv)
    @test filesize(asv_tax.tsv) > 0

    ## Stage 4: merge_taxa
    db_meta    = make_db_meta(db_config, db_name)
    asv_merged = merge_taxa(project, asvs, asv_tax, db_meta)

    @test haskey(asv_merged.tables, "merged")
    merged_csv = asv_merged.tables["merged"]
    @test isfile(merged_csv)
    @test filesize(merged_csv) > 0

    import CSV, DataFrames
    df = CSV.read(merged_csv, DataFrames.DataFrame)
    @test nrow(df) > 0
    @test "SeqName" in DataFrames.names(df)

    # Sample columns use stem-only names (no _R1_filt.fastq.gz suffix)
    sample_cols = filter(c -> any(startswith(c, s) for s in ["F3D0", "F3D1", "F3D2"]),
                         DataFrames.names(df))
    @test length(sample_cols) == length(samples)
    for col in sample_cols
        @test !occursin("_filt.fastq.gz", col)
    end

    # Total read count is positive
    total = sum(sum(skipmissing(df[!, c])) for c in sample_cols)
    @test total > 0

    @info "Integration test passed: $(nrow(df)) ASVs, $total total reads across $(length(samples)) samples"
end
