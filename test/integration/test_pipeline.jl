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
#   - analyse_run  for each run
#   - analyse_study across both runs
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
            for subdir in ("cutadapt", "dada2", "swarm", "vsearch", "merged", "analysis")
                d = joinpath(entry, subdir)
                isdir(d) && rm(d; recursive=true)
            end
        end
        study_analysis = joinpath(test_study_dir, "analysis")
        isdir(study_analysis) && rm(study_analysis; recursive=true)
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

    db_meta = make_db_meta(db_config, db_name)
    r_lock_int = ReentrantLock()
    plot_lock  = ReentrantLock()

    merged_results = Vector{Union{MergedTables,Nothing}}(undef, length(projects))
    asvs_results   = Vector{Union{ASVResult,Nothing}}(undef, length(projects))

    ## Run pipeline for each project
    Threads.@threads for (i, project) in collect(enumerate(projects))
        n_samples = length(filter(f -> endswith(f, ".fastq.gz"), readdir(project.data_dir))) ÷ 2

        errors = validate_project(project, db_config)
        if !isempty(errors)
            @warn "Integration test: run $run_name skipped - validation failed"
            continue
        end

        # Stage 1: cutadapt
        trimmed = cutadapt(project)
        @test isdir(trimmed.dir)
        trimmed_files = filter(f -> endswith(f, "_trimmed.fastq.gz"), readdir(trimmed.dir))
        @test length(trimmed_files) == n_samples * 2

        # Stage 2: DADA2 + SWARM in parallel
        dada2_task = Threads.@spawn begin
            result = lock(r_lock_int) do
                dada2(project, trimmed, taxonomy_db=dbs_loaded[dada2_key])
            end
            cdhit(project, result)
        end
        swarm_task = Threads.@spawn swarm(project, trimmed)

        asvs = fetch(dada2_task)
        otus = fetch(swarm_task)

        @test isfile(asvs.fasta) && filesize(asvs.fasta) > 0
        @test isfile(asvs.count_table) && filesize(asvs.count_table) > 0

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
            MergedTables(merge(asv_merged.tables, otu_merged.tables),
                         asv_merged.filter_order, asv_merged.filter_colours)
        else
            asv_merged
        end

        @test haskey(merged.tables, "merged")
        merged_csv = merged.tables["merged"]
        @test isfile(merged_csv) && filesize(merged_csv) > 0

        df = CSV.read(merged_csv, DataFrames.DataFrame)
        @test nrow(df) > 0
        @test "SeqName" in DataFrames.names(df)

        sample_cols = Analysis._sample_cols(df, db_meta)
        @test length(sample_cols) == n_samples
        for col in sample_cols
            @test !occursin("_filt.fastq.gz", col)
        end
        total = sum(sum(skipmissing(df[!, c])) for c in sample_cols)
        @test total > 0

        merged_results[i] = merged
        asvs_results[i]   = asvs

        @info "Pipeline [$(basename(project.dir))]: $(nrow(df)) ASVs, $total reads across $n_samples samples"

        ## Stage 5: analyse_run
        analyse_run(project, merged, asvs, db_meta; plot_lock)

        analysis_dir = joinpath(project.dir, "analysis")
        @test isdir(analysis_dir)
        @test isfile(joinpath(analysis_dir, "pipeline_summary.csv"))
        @test isfile(joinpath(analysis_dir, "Figures", "alpha_diversity.pdf"))
    end

    ## Stage 6: analyse_study (both runs together)
    valid = count(!isnothing, merged_results)
    if valid >= 2
        analyse_study(projects, merged_results, fill(db_meta, length(projects)); plot_lock)

        study_analysis_dir = joinpath(test_study_dir, "analysis")
        @test isdir(study_analysis_dir)
        @test isfile(joinpath(study_analysis_dir, "Figures", "nmds.pdf")) ||
              isfile(joinpath(study_analysis_dir, "Figures", "alpha_comparison.pdf"))
        @info "Study analysis complete"
    end

end
