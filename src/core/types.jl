module PipelineTypes

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

    export HasFasta, ProjectCtx, TrimmedReads, ASVResult, OTUResult, DenoisedASVs, TaxonomyHits, MergedTables, DatabaseMeta,
           StageNode

    abstract type HasFasta end

    struct ProjectCtx
        dir::String           # projects/{.../run}/ - stage outputs and run_config.yml live here
        config_dir::String    # config/ - global defaults, tools.yml, primers.yml
        data_dir::String      # data/{.../run}/ - FASTQ input files
        study_dir::String     # projects/{study}/ - output root (logs, analysis)
        data_study_dir::String # data/{study}/ - top of the per-run config cascade
    end

    struct TrimmedReads
        dir::String         # .../cutadapt/
    end

    struct ASVResult <: HasFasta
        fasta::String        # .../dada2/Tables/asvs.fasta  (or cdhit/asvs.fasta after cdhit)
        count_table::String  # .../dada2/Tables/seqtab_nochim.csv
        taxonomy::String     # .../dada2/Tables/taxonomy.csv
    end

    struct OTUResult <: HasFasta
        fasta::String        # .../swarm/seeds.fasta (OTU representative sequences)
        count_table::String  # .../swarm/otu_table.csv (per-OTU per-sample counts)
    end

    struct DenoisedASVs
        chimera_ckpt::String    # absolute path to ckpt_chimera.RData
        config_path::String     # path to run_config.yml (merged cascade)
        workspace_root::String  # path to dada2/ output directory
        input_dir::String       # path to trimmed reads directory
    end

    struct TaxonomyHits
        tsv::String          # .../vsearch/taxonomy.tsv
    end

    struct DatabaseMeta
        name::String
        levels::Vector{String}
        vsearch_format::String
        corrections::Vector{Dict{String,Any}}
        noncounts::Set{String}
    end

    struct MergedTables
        tables::Dict{String,String}  # name => CSV path; always includes "merged" (unfiltered)
        filter_order::Vector{String}  # filter stems in pipeline.yml order (for priority)
        filter_colours::Dict{String,String}  # filter stem => hex colour override (from YAML "colour" key)
    end

    """
        StageNode(name, label, inputs, output, config_sections)

    Declarative description of one pipeline stage.

    - `name`:            machine identifier (Symbol), e.g. `:cutadapt`
    - `label`:           human-readable stage name
    - `inputs`:          wire types this stage requires (excluding `ProjectCtx`, which every
                         ProjectCtx-aware stage implicitly receives)
    - `output`:          wire type produced; `Nothing` for side-effect-only stages
    - `config_sections`: pipeline.yml section keys whose content controls this stage
                         (used by skip guards and the WebUI to know what to re-hash)
    """
    struct StageNode
        name::Symbol
        label::String
        inputs::Vector{DataType}
        output::DataType
        config_sections::Vector{String}
    end

end
