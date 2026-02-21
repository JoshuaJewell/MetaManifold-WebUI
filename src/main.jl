#!/usr/bin/env julia

include("types.jl")
include("config.jl")
include("databases.jl")
include("call_tools.jl")
include("dada2.jl")
include("merge_and_filter_taxa.jl")
include("project.jl")

using CSV
using YAML
using .PipelineTypes, .Config, .Databases, .Tools, .TaxonomyTableTools, .DADA2, .ProjectSetup

## Instantiate parameters
config_dir       = "./config"
databases_config = joinpath(config_dir, "databases.yml")

dbs = ensure_databases(databases_config)

## Main
const r_lock = ReentrantLock()

projects = new_project("Multi_v_Vespa")

Threads.@threads for project in projects
    multiqc(project.data_dir, joinpath(project.dir, "QC"))

    trimmed = cutadapt(project)

    asvs = lock(r_lock) do
        dada2(project, trimmed, taxonomy_db = dbs["pr2_dada2"])
    end
    #asvs = lock(r_lock) do; cdhit(project, asvs); end  # optional, uncomment to enable
    tax    = vsearch(project, asvs, dbs["pr2_vsearch"])
    merged = merge_taxa(project, asvs, tax)
end
