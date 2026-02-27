module Analysis

# Post-pipeline analysis: per-run, per-group, and study-level outputs.
#
# Entry points (called from main.jl):
#   analyse_run(project, merged, asvs, db_meta; plot_lock)
#   analyse_study(projects, merged_results, db_metas; plot_lock)
#   load_metadata(start_dir, study_dir)
#
# Internal structure (src/analysis/):
#   helpers.jl     - column helpers, taxonomy view, rank/source key utilities
#   stats.jl       - alpha diversity, count matrices, NMDS/PERMANOVA (R), metadata loading
#   composition.jl - pipeline summary CSV, priority filter composition
#   report.jl      - chart generation, text report writing, table formatters
#   run.jl         - analyse_run   (level 1: per-run charts + reports)
#   group.jl       - _analyse_group (level 2: multi-run comparison)
#   study.jl       - _analyse_study_level + analyse_study (level 3: cross-group)
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

export analyse_run, analyse_study, load_metadata

    using CSV, DataFrames, Dates, Logging, Statistics, YAML, RCall
    using ..PipelineTypes, ..PipelineLog, ..Config, ..DiversityMetrics, ..PipelinePlots

    include("helpers.jl")
    include("stats.jl")
    include("composition.jl")
    include("report.jl")
    include("run.jl")
    include("group.jl")
    include("study.jl")

end
