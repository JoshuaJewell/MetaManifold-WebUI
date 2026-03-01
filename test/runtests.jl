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

## Bootstrap: load the pipeline modules without running main.jl
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

# Load only what tests need - avoids CairoMakie startup cost for unit tests
include(joinpath(@__DIR__, "..", "src", "core", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "graph.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "log.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "config.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "databases.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "validate.jl"))
include(joinpath(@__DIR__, "..", "src", "pipeline", "tools.jl"))
include(joinpath(@__DIR__, "..", "src", "pipeline", "merge_taxa.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "project.jl"))
include(joinpath(@__DIR__, "..", "src", "analysis", "diversity.jl"))
include(joinpath(@__DIR__, "..", "src", "analysis", "plotly.jl"))

using CSV, DataFrames, JSON3, Logging, YAML
using .PipelineTypes, .PipelineGraph, .PipelineLog, .Config, .Databases, .Validation
using .Tools, .TaxonomyTableTools, .ProjectSetup
using .DiversityMetrics, .PipelinePlotsPlotly

## Unit tests (always run)
@testset "MetabarcodingPipeline" begin

    include("unit/test_diversity.jl")
    include("unit/test_merge_taxa.jl")
    include("unit/test_config.jl")
    include("unit/test_validation.jl")
    include("unit/test_tools.jl")
    include("unit/test_plotly.jl")

    ## Integration tests (opt-in)
    if RUN_INTEGRATION
        include(joinpath(@__DIR__, "..", "src", "pipeline", "dada2.jl"))
        include(joinpath(@__DIR__, "..", "src", "pipeline", "swarm.jl"))
        include(joinpath(@__DIR__, "..", "src", "analysis", "plots.jl"))
        include(joinpath(@__DIR__, "..", "src", "analysis", "plotly.jl"))
        include(joinpath(@__DIR__, "..", "src", "analysis", "analysis.jl"))
        using JSON3
        using .DADA2, .OTUPipeline, .Analysis
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
