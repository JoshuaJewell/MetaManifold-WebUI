module PipelineTypes

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

    export HasFasta, ProjectCtx, TrimmedReads, ASVResult, OTUResult, DenoisedASVs, TaxonomyHits, MergedTables, DatabaseMeta,
           find_fastqs, FastqEntry

    abstract type HasFasta end

    """
        ProjectCtx(dir, config_dir, data_dir, study_dir, data_study_dir[, data_dirs])

    Execution context for a single pipeline run.

    When `pool_children: true` is set in a group-level `pipeline.yml`, the group
    directory becomes the leaf run. `data_dir` points to the group directory (used
    for config cascade resolution) while `data_dirs` lists the child directories
    that actually contain FASTQ files. For non-pooled runs `data_dirs == [data_dir]`.
    """
    struct ProjectCtx
        dir::String            # projects/{.../run}/ - stage outputs and run_config.yml live here
        config_dir::String     # config/ - global defaults, tools.yml, primers.yml
        data_dir::String       # data/{.../run}/ - primary data dir (config cascade anchor)
        study_dir::String      # projects/{study}/ - output root (logs, analysis)
        data_study_dir::String # data/{study}/ - top of the per-run config cascade
        data_dirs::Vector{String} # all FASTQ source dirs; [data_dir] when not pooled
    end

    # Convenience constructor kept for backwards-compatible 5-arg call sites.
    ProjectCtx(dir, config_dir, data_dir, study_dir, data_study_dir) =
        ProjectCtx(dir, config_dir, data_dir, study_dir, data_study_dir, [data_dir])

    struct FastqEntry
        path::String    # absolute path to the .fastq.gz file
        name::String    # output-safe sample filename (prefixed when pooled)
    end

    """
        find_fastqs(project) -> Vector{FastqEntry}

    Collect all `.fastq.gz` files across `project.data_dirs`. When the project
    pools children (`length(data_dirs) > 1`), each file is prefixed with a
    sanitised form of its subdirectory name (spaces to underscores) to prevent
    collisions and allow downstream filtering by sub-group.
    """
    function find_fastqs(project::ProjectCtx)
        entries = FastqEntry[]
        pooled  = length(project.data_dirs) > 1
        for d in project.data_dirs
            prefix = pooled ? replace(basename(d), " " => "_") * "_" : ""
            for f in readdir(d)
                endswith(f, ".fastq.gz") || continue
                push!(entries, FastqEntry(joinpath(d, f), prefix * f))
            end
        end
        return entries
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
    end

end
