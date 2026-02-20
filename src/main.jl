#!/usr/bin/env julia

include("databases.jl")
include("call_tools.jl")
include("dada2.jl")
include("merge_and_filter_taxa.jl")

using CSV
using YAML
using .Databases, .Tools, .TaxonomyTableTools, .DADA2

## Instantiate filesystem
# Root directories
data_dir = "./data"
config_dir = "./config"
output_dir = "./output"

# Project name - all stage outputs live under output/{project_name}/
# Re-running with the same project_name overwrites previous results.
project_name = "Multiplex_pool"
project_dir  = joinpath(output_dir, project_name)
fastq_input_dir = joinpath(data_dir, project_name)

# QC paths
fastqc_dir = joinpath(project_dir, "FastQC")

# Cutadapt paths
trimmed_dir = joinpath(project_dir, "cutadapt")
primers_config = joinpath(config_dir, "primers.yml")

# DADA2 paths
dada2_dir = joinpath(project_dir, "dada2")

# VSEARCH paths
vsearch_dir = joinpath(project_dir, "vsearch")
fasta_outfile = joinpath(dada2_dir, "Tables/asvs.fasta")

# Merge/filter paths (derived from project_dir)
merged_dir = joinpath(project_dir, "merged")
merged_outfile_multi = joinpath(merged_dir, "merged_multi.csv")
filtered_outfile_multi = joinpath(merged_dir, "protist_filtered_multi.csv")

## Instantiate parameters
# Cutadapt parameters
primer_pairs = ["TarEuk", "Meta2"]
cutadapt_optional_args = "-m 200 --discard-untrimmed"

# DADA2 parameters
dada2_config_path = joinpath(config_dir, "dada2.yml")

# VSEARCH parameters
vsearch_optional_args = "--id 0.75 --query_cov 0.8"

# Merge and filter (DADA2-VSEARCH) parameters
mkpath(merged_dir)
multiv = joinpath(vsearch_dir, "taxonomy.tsv")
multid = joinpath(dada2_dir, "Tables/taxonomy.csv")

protist_filter = joinpath(config_dir, "protist_filter.yml")

# Download/use database
dbs = ensure_databases(dada2_config_path)

## Main
#fastqc_all(fastq_input_dir, fastqc_dir)

#cutadapt(primer_pairs, primers_config, fastq_input_dir, trimmed_dir, optional_args = cutadapt_optional_args)

dada2(dada2_config_path, input_dir = trimmed_dir, workspace_root = dada2_dir, taxonomy_db = dbs["pr2_dada2"])

#vsearch(fasta_outfile, dbs["pr2_vsearch"], vsearch_dir, optional_args = vsearch_optional_args)

CSV.write(merged_outfile_multi, merge_taxonomy_counts(multiv, multid));     @info "Written: $merged_outfile_multi"
CSV.write(filtered_outfile_multi, filter_table(merge_taxonomy_counts(multiv, multid), protist_filter)); @info "Written: $filtered_outfile_multi"
