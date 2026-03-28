#!/usr/bin/env julia
# MetabarcodingPipeline test suite
#
# Run with:
#   julia --project=. test/runtests.jl
#
# For integration tests (requires tools + databases):
#   julia --project=. -t4 test/runtests.jl --integration

using Test

const RUN_INTEGRATION = "--integration" in ARGS
const RUN_SERVER      = "--server"      in ARGS

## Bootstrap: load the pipeline modules without running server.jl
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

using CSV, DataFrames, JSON3, Logging, YAML, DuckDB, DBInterface

# Load only what tests need
include(joinpath(@__DIR__, "..", "src", "core", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "log.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "config.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "databases.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "duckdb_store.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "validate.jl"))
include(joinpath(@__DIR__, "..", "src", "pipeline", "tools.jl"))
include(joinpath(@__DIR__, "..", "src", "pipeline", "merge_taxa.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "project.jl"))
include(joinpath(@__DIR__, "..", "src", "analysis", "diversity.jl"))
include(joinpath(@__DIR__, "..", "src", "analysis", "analysis.jl"))

using .PipelineTypes, .PipelineLog, .Config, .Databases, .DuckDBStore, .Validation
using .Tools, .TaxonomyTableTools, .ProjectSetup
using .DiversityMetrics, .Analysis

## Unit tests (always run)
@testset "MetabarcodingPipeline" begin

    include("unit/test_diversity.jl")
    include("unit/test_merge_taxa.jl")
    include("unit/test_config.jl")
    include("unit/test_validation.jl")
    include("unit/test_tools.jl")
    include("unit/test_analysis.jl")
    include("unit/test_duckdb_store.jl")
    include("unit/test_analysis_duckdb.jl")
    include("unit/test_config_hashing.jl")
    include("unit/test_project.jl")
    include("unit/test_log.jl")
    include("unit/test_databases.jl")
    include("unit/test_merge_taxa_mappings.jl")

    ## Integration tests (opt-in)
    if RUN_INTEGRATION
        include(joinpath(@__DIR__, "..", "src", "pipeline", "dada2.jl"))
        include(joinpath(@__DIR__, "..", "src", "pipeline", "swarm.jl"))
        using .DADA2, .OTUPipeline
        include("integration/test_pipeline.jl")
    else
        @info "Skipping integration tests (pass --integration to enable)"
    end

    ## Server smoke tests (opt-in, starts a Julia subprocess - slow on first run)
    if RUN_SERVER
        using HTTP
        include("integration/test_server.jl")
    else
        @info "Skipping server smoke tests (pass --server to enable)"
    end

end
