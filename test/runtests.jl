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

using CSV, DataFrames, YAML
using .PipelineTypes, .PipelineGraph, .PipelineLog, .Config, .Databases, .Validation
using .Tools, .TaxonomyTableTools, .ProjectSetup
using .DiversityMetrics

## Unit tests (always run)
@testset "MetabarcodingPipeline" begin

    include("unit/test_diversity.jl")
    include("unit/test_merge_taxa.jl")
    include("unit/test_config.jl")
    include("unit/test_validation.jl")

    ## Integration tests (opt-in)
    if RUN_INTEGRATION
        include(joinpath(@__DIR__, "..", "src", "pipeline", "dada2.jl"))
        include(joinpath(@__DIR__, "..", "src", "pipeline", "swarm.jl"))
        include(joinpath(@__DIR__, "..", "src", "analysis", "plots.jl"))
        include(joinpath(@__DIR__, "..", "src", "analysis", "analysis.jl"))
        using .DADA2, .OTUPipeline
        include("integration/test_pipeline.jl")
    else
        @info "Skipping integration tests (pass --integration to enable)"
    end

end
