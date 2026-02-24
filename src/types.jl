module PipelineTypes

    export HasFasta, ProjectCtx, TrimmedReads, ASVResult, DenoisedASVs, TaxonomyHits, MergedTables, DatabaseMeta

    abstract type HasFasta end

    struct ProjectCtx
        dir::String         # projects/{.../run}/ - configs and stage output dirs live here
        config_dir::String  # config/
        data_dir::String    # data/{.../run}/ - FASTQ input files
        study_dir::String   # projects/{study}/ - top of the per-project config cascade
    end

    struct TrimmedReads
        dir::String         # .../cutadapt/
    end

    struct ASVResult <: HasFasta
        fasta::String        # .../dada2/Tables/asvs.fasta  (or cdhit/asvs.fasta after cdhit)
        count_table::String  # .../dada2/Tables/seqtab_nochim.csv
        taxonomy::String     # .../dada2/Tables/taxonomy.csv
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

end
