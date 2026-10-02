#!/usr/bin/env julia
# MetaManifold test suite
#
# Run with:
#   julia --project=. test/runtests.jl
#
# For integration tests (requires tools + databases):
#   julia --project=. -t4 test/runtests.jl --integration
#
# Name unit test files to run only those, e.g.
#   julia --project=. test/runtests.jl test_routes.jl trees
using Test

const RUN_INTEGRATION = "--integration" in ARGS
const RUN_SERVER      = "--server"      in ARGS

const UNIT_FILES = [
    "test_diversity.jl",
    "test_merge_taxa.jl",
    "test_config.jl",
    "test_validation.jl",
    "test_tools.jl",
    "test_analysis.jl",
    "test_duckdb_store.jl",
    "test_analysis_duckdb.jl",
    "test_config_hashing.jl",
    "test_project.jl",
    "test_log.jl",
    "test_databases.jl",
    "test_merge_taxa_mappings.jl",
    "test_routes.jl",
    "test_publication_tables.jl",
    "test_heatmap.jl",
    "test_trees.jl",
    "test_report_funnel.jl",
    "test_phylogeny.jl",
    "test_composition.jl",
    "test_composition_library.jl",
    "test_primers_library.jl",
    "test_databases_library.jl",
    "test_categories.jl",
    "test_read_conservation.jl",
    "test_r_runtime.jl",
    "test_dada2_commands.jl",
    "test_remote_stages.jl",
    "test_jobs.jl",
    "test_provenance.jl",
    "test_install_pins.jl",
    "test_migrate_composition.jl",
    "test_sample_reads.jl",
    "test_determinism.jl",
    "test_differential.jl",
]

# The files named on the command line ("test_routes.jl", "routes" or
# "unit/test_routes.jl"), or every unit file when none is named.
function _selected_files()
    named = filter(a -> !startswith(a, "--"), ARGS)
    isempty(named) && return UNIT_FILES
    files = map(named) do a
        b = basename(a)
        b = startswith(b, "test_") ? b : "test_" * b
        endswith(b, ".jl") ? b : b * ".jl"
    end
    unknown = filter(f -> !isfile(joinpath(@__DIR__, "unit", f)), files)
    isempty(unknown) || error("No such unit test file: $(join(unknown, ", "))")
    files
end

using MetaManifold
using CSV, DataFrames, JSON3, Logging, YAML, DuckDB, DBInterface, Dates

using MetaManifold.PipelineTypes, MetaManifold.PipelineLog, MetaManifold.Config
using MetaManifold.Databases, MetaManifold.DuckDBStore, MetaManifold.Validation
using MetaManifold.Tools, MetaManifold.TaxonomyTableTools, MetaManifold.ProjectSetup
using MetaManifold.DiversityMetrics, MetaManifold.Analysis
using MetaManifold.Categories, MetaManifold.CompositionLibrary

## Unit tests (always run)
@testset "MetaManifold" begin

    for f in _selected_files()
        include(joinpath(@__DIR__, "unit", f))
    end

    ## Integration tests (opt-in)
    if RUN_INTEGRATION
        using MetaManifold.DADA2, MetaManifold.OTUPipeline
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
