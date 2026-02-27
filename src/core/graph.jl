module PipelineGraph

# Declarative pipeline graph - single source of truth for stage identity.
#
# PIPELINE_STAGES lists every stage in dependency order. Each StageNode records:
#   - what wire types it consumes (inputs) and produces (output)
#   - which pipeline.yml config sections affect it
#
# config_sections is AUTHORITATIVE: skip guards in each stage function call
# stage_sections(:name) rather than hardcoding section strings. Any change to
# which config keys invalidate a stage is made here and nowhere else.
#
# Stages that consume ProjectCtx do so implicitly - it is omitted from `inputs`
# because every high-level overload takes it.
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under AGPLv3.

export PIPELINE_STAGES, stage_by_name, stage_sections, print_pipeline

    using ..PipelineTypes

    const PIPELINE_STAGES = StageNode[

        ## Quality Control
        StageNode(
            :fastqc_multiqc,
            "Quality Control",
            DataType[],
            Nothing,
            ["fastqc", "multiqc"]
        ),

        ## Primer Trimming
        StageNode(
            :cutadapt,
            "Primer Trimming",
            DataType[],
            TrimmedReads,
            ["cutadapt"]
        ),

        ## DADA2 sub-stages (individually resumable)
        StageNode(
            :dada2_filter_trim,
            "Filter & Trim",
            DataType[TrimmedReads],
            Nothing,            # produces ckpt_filter.RData; wire type is DenoisedASVs at end
            ["dada2.filter_trim"]
        ),
        StageNode(
            :dada2_learn_errors,
            "Learn Error Rates",
            DataType[TrimmedReads],
            Nothing,
            ["dada2.dada"]
        ),
        StageNode(
            :dada2_denoise,
            "Denoise & Merge",
            DataType[TrimmedReads],
            Nothing,
            ["dada2.dada", "dada2.merge"]
        ),
        StageNode(
            :dada2_filter_length,
            "ASV Length Filter",
            DataType[TrimmedReads],
            Nothing,
            ["dada2.asv"]
        ),
        StageNode(
            :dada2_chimera_removal,
            "Chimera Removal",
            DataType[TrimmedReads],
            DenoisedASVs,
            ["dada2.asv", "dada2.output"]
        ),
        StageNode(
            :dada2_assign_taxonomy,
            "Taxonomy Assignment (DADA2)",
            DataType[DenoisedASVs],
            ASVResult,
            ["dada2.taxonomy", "dada2.output"]
        ),

        ## Post-DADA2
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
            DataType[ASVResult],   # dispatches on HasFasta - also works for OTUResult
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

        ## OTU pipeline (parallel to DADA2)
        StageNode(
            :swarm,
            "OTU Clustering (SWARM)",
            DataType[TrimmedReads],
            OTUResult,
            ["swarm"]
        ),

        ## Analysis
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
            DataType[MergedTables],
            Nothing,
            ["analysis"]
        ),
        StageNode(
            :analyse_study,
            "Study-Wide Analysis",
            DataType[MergedTables],
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
        stage_sections(name::Symbol) -> String

    Return the comma-joined config section string for a stage, in the format
    expected by `_section_stale` / `_write_section_hash`.

    This is the single point of truth: skip guards call this rather than
    hardcoding section strings, so changes to which config keys invalidate a
    stage are made here and propagate automatically.

        _section_stale(config_path, stage_sections(:cutadapt), hash_file)
    """
    function stage_sections(name::Symbol)::String
        join(stage_by_name(name).config_sections, ",")
    end

    """
        print_pipeline([io])

    Print a human-readable summary of all declared pipeline stages.
    """
    function print_pipeline(io::IO = stdout)
        println(io, "Declared pipeline stages:")
        println(io, "-"^60)
        for (i, s) in enumerate(PIPELINE_STAGES)
            ins  = isempty(s.inputs) ? "raw FASTQs (via ProjectCtx)" :
                   join(string.(s.inputs), ", ")
            out  = s.output === Nothing ? "- (side effects only)" : string(s.output)
            cfgs = join(s.config_sections, ", ")
            println(io, "  $i. $(s.label)  [:$(s.name)]")
            println(io, "     inputs:  $ins")
            println(io, "     output:  $out")
            println(io, "     config:  $cfgs")
            i < length(PIPELINE_STAGES) && println(io)
        end
        println(io, "-"^60)
    end

end
