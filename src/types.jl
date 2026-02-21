module PipelineTypes

    export HasFasta, ProjectCtx, TrimmedReads, ASVResult, TaxonomyHits, MergedTables

    abstract type HasFasta end

    struct ProjectCtx
        dir::String         # projects/{.../run}/ - configs and stage output dirs live here
        config_dir::String  # config/
        data_dir::String    # data/{.../run}/ - FASTQ input files
    end

    struct TrimmedReads
        dir::String         # .../cutadapt/
    end

    struct ASVResult <: HasFasta
        fasta::String        # .../dada2/Tables/asvs.fasta  (or cdhit/asvs.fasta after cdhit)
        count_table::String  # .../dada2/Tables/seqtab_nochim.csv
        taxonomy::String     # .../dada2/Tables/taxonomy.csv
    end

    struct TaxonomyHits
        tsv::String          # .../vsearch/taxonomy.tsv
    end

    struct MergedTables
        tables::Dict{String,String}  # name => CSV path; always includes "merged" (unfiltered)
    end

end
