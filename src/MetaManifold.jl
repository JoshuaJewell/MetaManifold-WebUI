# SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

module MetaManifold

# Core
include("core/types.jl")
include("core/r_runtime.jl")
# Provenance probes R through the shared runtime lock, so it follows r_runtime.jl.
include("core/provenance.jl")
# Tools registers the children it starts with the probe, so it precedes pipeline/.
include("core/memory_probe.jl")
include("core/log.jl")
include("core/config.jl")
include("core/databases.jl")
include("core/duckdb_store.jl")
include("core/validate.jl")
include("core/project.jl")
include("core/categories.jl")
include("core/composition_library.jl")
include("core/primers_library.jl")
include("core/databases_library.jl")

# Annotation

# Pipeline
include("pipeline/tools.jl")
include("pipeline/remote_exec.jl")
include("pipeline/merge_taxa.jl")
include("pipeline/dada2.jl")
include("pipeline/swarm.jl")
include("pipeline/phylogeny.jl")

# Analysis
include("canonical.jl")
include("analysis/numeric_policy.jl")
include("analysis/exact_summaries.jl")
include("analysis/high_precision.jl")
include("analysis/derived_analyses.jl")
include("analysis/method_catalogue.jl")
include("analysis/diversity.jl")
include("analysis/analysis.jl")
include("analysis/differential.jl")

# Web server
include("server/server.jl")

end
