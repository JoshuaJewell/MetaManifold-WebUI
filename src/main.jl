#!/usr/bin/env julia

# Activate the project environment (Project.toml one directory up from src/).
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

include("types.jl")
include("log.jl")
include("config.jl")
include("databases.jl")
include("call_tools.jl")
include("dada2.jl")
include("merge_and_filter_taxa.jl")
include("project.jl")
include("diversity.jl")
include("plots.jl")
include("analysis.jl")

using CSV
using YAML
using .PipelineTypes, .PipelineLog, .Config, .Databases, .Tools, .TaxonomyTableTools, .DADA2, .ProjectSetup
using .DiversityMetrics, .PipelinePlots, .Analysis

## Instantiate parameters
config_dir       = "./config"
databases_config = joinpath(config_dir, "databases.yml")

dbs = ensure_databases(databases_config)

## Main
const r_lock = ReentrantLock()

projects = new_project("PR2_v_SILVA")

merged_results = Vector{Union{MergedTables, Nothing}}(undef, length(projects))
asvs_results   = Vector{Union{ASVResult, Nothing}}(undef, length(projects))
db_metas       = Vector{DatabaseMeta}(undef, length(projects))

Threads.@threads for (i, project) in collect(enumerate(projects))
    reset_log(project)
    multiqc(project.data_dir, joinpath(project.dir, "QC"))

    # Resolve per-project database from its config cascade.
    run_cfg_path = write_run_config(project)
    db_name = string(get(get(get(YAML.load_file(run_cfg_path), "dada2", Dict()), "taxonomy", Dict()), "database", "pr2"))
    db_metas[i] = make_db_meta(databases_config, db_name)

    trimmed = cutadapt(project)

    asvs = lock(r_lock) do
        dada2(project, trimmed, taxonomy_db = dbs["$(db_name)_dada2"])
    end
    asvs = lock(r_lock) do; cdhit(project, asvs; optional_args = "-c 1"); end

    if haskey(dbs, "$(db_name)_vsearch")
        tax    = vsearch(project, asvs, dbs["$(db_name)_vsearch"])
        merged = merge_taxa(project, asvs, tax, db_metas[i])
    else
        @warn "No vsearch database for '$db_name' - using DADA2 tax_counts as merged for $(basename(project.dir))"
        tax_counts = joinpath(project.dir, "dada2", "Tables", "tax_counts.csv")
        merged = isfile(tax_counts) && filesize(tax_counts) > 0 ?
                 MergedTables(Dict("merged" => tax_counts), String[], Dict{String,String}()) :
                 nothing
    end
    merged_results[i] = merged
    asvs_results[i]   = asvs
end

# Analysis: data prep runs in parallel, CairoMakie calls serialized via plot_lock.
const plot_lock = ReentrantLock()

Threads.@threads for (i, project) in collect(enumerate(projects))
    merged = merged_results[i]
    asvs   = asvs_results[i]
    (isnothing(merged) || isnothing(asvs)) && continue
    analyse_run(project, merged, asvs, db_metas[i]; plot_lock)
end

if any(!isnothing, merged_results)
    analyse_study(projects, merged_results, db_metas; plot_lock)
end

write_combined_log(projects)
