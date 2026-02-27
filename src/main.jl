#!/usr/bin/env julia

# Activate the project environment (Project.toml one directory up from src/).
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

import Logging, Dates
struct _TimestampLogger <: Logging.AbstractLogger
    inner::Logging.AbstractLogger
end
Logging.min_enabled_level(l::_TimestampLogger) = Logging.min_enabled_level(l.inner)
Logging.shouldlog(l::_TimestampLogger, args...) = Logging.shouldlog(l.inner, args...)
Logging.catch_exceptions(l::_TimestampLogger) = Logging.catch_exceptions(l.inner)
function Logging.handle_message(l::_TimestampLogger, level, message, _module, group, id, file, line; kwargs...)
    ts = Dates.format(Dates.now(), "HH:MM:SS")
    Logging.handle_message(l.inner, level, "[$ts] $message", _module, group, id, file, line; kwargs...)
end
Logging.global_logger(_TimestampLogger(Logging.global_logger()))

include("core/types.jl")
include("core/graph.jl")
include("core/log.jl")
include("core/config.jl")
include("core/databases.jl")
include("core/validate.jl")
include("pipeline/tools.jl")
include("pipeline/dada2.jl")
include("pipeline/swarm.jl")
include("pipeline/merge_taxa.jl")
include("core/project.jl")
include("analysis/diversity.jl")
include("analysis/plots.jl")
include("analysis/analysis.jl")

using CSV
using YAML
using .PipelineTypes, .PipelineGraph, .PipelineLog, .Config, .Databases, .Validation, .Tools, .TaxonomyTableTools, .DADA2, .OTUPipeline, .ProjectSetup
using .DiversityMetrics, .PipelinePlots, .Analysis

## Instantiate parameters
config_dir       = "./config"
databases_config = joinpath(config_dir, "databases.yml")

dbs = ensure_databases(databases_config)

## Main
const r_lock = ReentrantLock()

projects = new_project("Multi_v_Vespa")

tools_config = joinpath(config_dir, "tools.yml")
n_errors = validate_environment(projects, databases_config, tools_config)
n_errors > 0 && error("Aborting: fix the $n_errors validation error(s) above before running the pipeline.")

merged_results = Vector{Union{MergedTables, Nothing}}(undef, length(projects))
asvs_results   = Vector{Union{ASVResult, Nothing}}(undef, length(projects))
db_metas       = Vector{DatabaseMeta}(undef, length(projects))

Threads.@threads for (i, project) in collect(enumerate(projects))
    try
        reset_log(project)
        multiqc(project)

        # Resolve per-project database from its config cascade.
        run_cfg_path = write_run_config(project)
        db_name = string(get(get(get(YAML.load_file(run_cfg_path), "dada2", Dict()), "taxonomy", Dict()), "database", "pr2"))
        db_metas[i] = make_db_meta(databases_config, db_name)

        trimmed = cutadapt(project)

        # DADA2+cdhit and SWARM are independent - run in parallel.
        dada2_task = Threads.@spawn begin
            result = lock(r_lock) do
                dada2(project, trimmed, taxonomy_db = dbs["$(db_name)_dada2"])
            end
            cdhit(project, result)
        end
        swarm_task = Threads.@spawn swarm(project, trimmed)

        asvs = fetch(dada2_task)
        otus = fetch(swarm_task)

        if haskey(dbs, "$(db_name)_vsearch")
            # vsearch for ASVs and OTUs are independent - run in parallel.
            asv_tax_task = Threads.@spawn vsearch(project, asvs, dbs["$(db_name)_vsearch"])
            otu_tax_task = isnothing(otus) ? nothing :
                           Threads.@spawn vsearch(project, otus, dbs["$(db_name)_vsearch"])

            asv_tax    = fetch(asv_tax_task)
            asv_merged = merge_taxa(project, asvs, asv_tax, db_metas[i])

            if !isnothing(otus)
                otu_tax    = fetch(otu_tax_task)
                otu_merged = merge_taxa_otu(project, otus, otu_tax, db_metas[i])
                # Combine ASV and OTU tables into a single MergedTables.
                # OTU keys carry _otu suffix so there is no collision with ASV keys.
                # filter_order and filter_colours come from the ASV side (used for the
                # priority filter composition chart which is ASV-based).
                merged = MergedTables(
                    merge(asv_merged.tables, otu_merged.tables),
                    asv_merged.filter_order,
                    asv_merged.filter_colours
                )
            else
                merged = asv_merged
            end
        else
            @warn "No vsearch database for '$db_name' - using DADA2 tax_counts as merged for $(basename(project.dir))"
            tax_counts = joinpath(project.dir, "dada2", "Tables", "tax_counts.csv")
            merged = isfile(tax_counts) && filesize(tax_counts) > 0 ?
                     MergedTables(Dict("merged" => tax_counts), String[], Dict{String,String}()) :
                     nothing
        end
        merged_results[i] = merged
        asvs_results[i]   = asvs
    catch e
        @error "Project $(basename(project.dir)) failed" exception=(e, catch_backtrace())
        merged_results[i] = nothing
        asvs_results[i]   = nothing
    end
end

# Analysis: data prep runs in parallel, CairoMakie calls serialized via plot_lock.
const plot_lock = ReentrantLock()

Threads.@threads for (i, project) in collect(enumerate(projects))
    merged = merged_results[i]
    asvs   = asvs_results[i]
    (isnothing(merged) || isnothing(asvs)) && continue
    isassigned(db_metas, i) || continue
    try
        analyse_run(project, merged, asvs, db_metas[i]; plot_lock)
    catch e
        @error "analyse_run failed for $(basename(project.dir))" exception=(e, catch_backtrace())
    end
end

if any(!isnothing, merged_results)
    analyse_study(projects, merged_results, db_metas; plot_lock)
end

write_combined_log(projects)
