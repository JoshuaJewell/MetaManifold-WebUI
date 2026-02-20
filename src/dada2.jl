module DADA2

# DADA2 amplicon sequencing pipeline - Julia orchestrator
#
# Stages and their Web UI page groupings:
#
#   QC page:        prefilter_qc, filter_trim
#   Denoising page: learn_errors, denoise, filter_length, chimera_removal
#   Taxonomy page:  assign_taxonomy
#
# Each stage is independently callable and saves a checkpoint so the pipeline
# can be resumed or individual stages re-run across Julia sessions.
#
# Call dada2() to run all stages in sequence without stopping.
#
# Notice:
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).
#
# This work is based on the DADA2 tutorial by Benjamin J. Callahan, et al.,
# available at https://benjjneb.github.io/dada2/tutorial.html, with modification
# into a single module. The original material is licensed under the Creative
# Commons Attribution 4.0 International License (CC BY 4.0):
# https://creativecommons.org/licenses/by/4.0/.

export dada2, prefilter_qc, filter_trim, learn_errors, denoise,
       filter_length, chimera_removal, assign_taxonomy

    import Downloads
    using Logging, RCall, YAML

    include("dada2/context.jl")   # shared helpers and _pipeline_context
    include("dada2/qc.jl")        # prefilter_qc, filter_trim
    include("dada2/denoise.jl")   # learn_errors, denoise, filter_length
    include("dada2/chimera.jl")   # chimera_removal
    include("dada2/taxonomy.jl")  # assign_taxonomy

    """
        dada2(config_path; progress)

    Run the complete DADA2 pipeline from raw reads to a taxonomy-annotated ASV
    table, calling all stages in sequence:

        prefilter_qc -> filter_trim -> learn_errors ->
        denoise -> filter_length -> chimera_removal -> assign_taxonomy

    For interactive use - reviewing intermediate outputs or adjusting parameters
    between steps - call the individual stage functions directly instead.

    ## Outputs
    Written to `workspace.root/Tables/`:
    - `seqtab_nochim.csv`              - chimera-free ASV count table
    - `asvs.fasta` / `asvs.csv`       - ASV sequences with short identifiers
    - `taxonomy.csv`                   - taxonomy assignments
    - `taxonomy_bootstraps.csv`        - bootstrap confidence values
    - `taxonomy_combined.csv`          - taxonomy + bootstraps combined
    - `tax_counts.csv`                 - taxonomy + per-sample counts
    - `asv_counts.csv`                 - ASV sequences + per-sample counts (no taxonomy)
    - `pipeline_stats.csv`             - read counts at each pipeline stage

    Written to `workspace.root/Checkpoints/`:
    - `ckpt_filter.RData`  - filter_stats
    - `ckpt_errors.RData`  - fwd_errors, rev_errors
    - `ckpt_denoise.RData` - dada objects, merged, unfiltered seq_table
    - `ckpt_length.RData`  - dada objects, merged, length-filtered seq_table
    - `ckpt_chimera.RData` - seq_table_nochim, index
    - `checkpoint.RData`   - final R environment snapshot
    """
    function dada2(config_path::String; progress=nothing, input_dir=nothing, workspace_root=nothing, taxonomy_db=nothing)
        #prefilter_qc(config_path; progress, input_dir, workspace_root)
        #filter_trim(config_path; progress, input_dir, workspace_root); R"gc()"
        #learn_errors(config_path; progress, input_dir, workspace_root); R"gc()"
        #denoise(config_path; progress, input_dir, workspace_root); R"gc()"
        #filter_length(config_path; progress, input_dir, workspace_root); R"gc()"
        #chimera_removal(config_path; progress, input_dir, workspace_root); R"gc()"
        assign_taxonomy(config_path; progress, input_dir, workspace_root, taxonomy_db)
    end

end
