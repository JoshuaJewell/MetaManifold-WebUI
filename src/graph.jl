module PipelineGraph

# Declarative pipeline graph.
#
# PIPELINE_STAGES lists every stage in dependency order. Each StageNode records:
#   - what wire types it consumes (inputs) and produces (output)
#   - which pipeline.yml config sections affect it (for cache invalidation and WebUI)
#
# Stages that consume ProjectCtx do so implicitly — it is omitted from `inputs`
# because every high-level overload takes it. The graph is linear for Illumina
# paired-end runs; branching (e.g. dada2 vs swarm, cdhit optional) is expressed
# by optional/alternative nodes that share the same output type.
#
# This module does not drive execution; main.jl still calls stages explicitly.
# The declaration exists so tools (print_pipeline, the WebUI, future schedulers)
# can inspect the graph without parsing source code.
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under AGPLv3.

export PIPELINE_STAGES, print_pipeline, stage_by_name

    using ..PipelineTypes

    const PIPELINE_STAGES = StageNode[
        StageNode(
            :fastqc_multiqc,
            "Quality Control",
            DataType[],           # inputs: raw FASTQs are found via ProjectCtx.data_dir
            Nothing,
            ["fastqc", "multiqc"]
        ),
        StageNode(
            :cutadapt,
            "Primer Trimming",
            DataType[],           # inputs: raw FASTQs via ProjectCtx.data_dir
            TrimmedReads,
            ["cutadapt"]
        ),
        StageNode(
            :dada2,
            "Denoising & ASV Table",
            DataType[TrimmedReads],
            ASVResult,
            ["dada2"]
        ),
        StageNode(
            :cdhit,
            "ASV Dereplication (optional)",
            DataType[ASVResult],
            ASVResult,
            ["cdhit"]
        ),
        StageNode(
            :vsearch,
            "Taxonomy Search",
            DataType[ASVResult],  # dispatches on HasFasta — also works for OTUResult
            TaxonomyHits,
            ["vsearch"]
        ),
        StageNode(
            :merge_taxa,
            "Taxonomy Table Merge & Filter",
            DataType[ASVResult, TaxonomyHits],
            MergedTables,
            ["merge_taxa"]
        ),
        StageNode(
            :analyse_run,
            "Per-Run Analysis",
            DataType[MergedTables, ASVResult],
            Nothing,
            ["analysis"]
        ),
        StageNode(
            :analyse_group,
            "Per-Group Analysis",
            DataType[MergedTables],  # vector of MergedTables across runs in a group
            Nothing,
            ["analysis"]
        ),
        StageNode(
            :analyse_study,
            "Study-Wide Analysis",
            DataType[MergedTables],  # vector across all groups
            Nothing,
            ["analysis"]
        ),
    ]

    """
        stage_by_name(name::Symbol) -> StageNode

    Look up a stage by its machine name. Throws if not found.
    """
    function stage_by_name(name::Symbol)
        idx = findfirst(s -> s.name == name, PIPELINE_STAGES)
        isnothing(idx) && error("No pipeline stage named :$name")
        PIPELINE_STAGES[idx]
    end

    """
        print_pipeline([io])

    Print a human-readable summary of all declared pipeline stages.
    """
    function print_pipeline(io::IO = stdout)
        println(io, "Declared pipeline stages:")
        println(io, "─"^60)
        for (i, s) in enumerate(PIPELINE_STAGES)
            ins  = isempty(s.inputs) ? "raw FASTQs (via ProjectCtx)" :
                   join(string.(s.inputs), ", ")
            out  = s.output === Nothing ? "— (side effects only)" : string(s.output)
            cfgs = join(s.config_sections, ", ")
            println(io, "  $i. $(s.label)  [:$(s.name)]")
            println(io, "     inputs:  $ins")
            println(io, "     output:  $out")
            println(io, "     config:  $cfgs")
            i < length(PIPELINE_STAGES) && println(io)
        end
        println(io, "─"^60)
    end

end
