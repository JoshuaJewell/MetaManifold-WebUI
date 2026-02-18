include("run_cutadapt.jl")
include("run_dada2.jl")
include("merge_and_filter_taxa.jl")

using CSV
using .Cutadapt, .TaxonomyTableTools, .DADA2

## Instantiate filesystem
# Root directories
data_dir = "./data"
config_dir = "./config"
output_dir = "./output"

# Cutadapt paths
fastq_input_dir = joinpath(data_dir, "fastq")
cutadapt_dir = joinpath(output_dir, "cutadapt")
primers_config = joinpath(config_dir, "primers.yml")

## Instantiate parameters
# Cutadapt parameters
primer_pairs = ["TarEuk", "Meta2"]
optional_args = "-m 200 --discard-untrimmed"

# DADA2 parameters
dada2_config_dir = joinpath(config_dir, "dada2_config.yml")

# Merge and filter (DADA2-VSEARCH) parameters
multiv = joinpath(output_dir, "vsearch/taxonomy_multi_pool.tsv")
vespav = joinpath(output_dir, "vsearch/taxonomy_vespa_pool_fwdonly.tsv")
multid = joinpath(output_dir, "dada2/tax_counts_fasta_multi_pool.csv")
vespad = joinpath(output_dir, "dada2/tax_counts_fasta_vespa_pool_fwdonly.csv")

protist_filter = joinpath(config_dir, "protist_filter.yml")

merged_outfile_multi = joinpath(output_dir, "merged_multi.csv")
merged_outfile_vespa = joinpath(output_dir, "merged_vespa.csv")

filtered_outfile_multi = joinpath(output_dir, "protist_filtered_vespa.csv")
filtered_outfile_vespa = joinpath(output_dir, "protist_filtered_vespa.csv")

## Main

#cutadapt(primer_pairs, primers_config, fastq_input_dir, cutadapt_dir, optional_args = optional_args)

#dada2(dada_config_dir)

#CSV.write(merged_outfile_multi, merge_taxonomy_counts(multiv, multid))
#CSV.write(merged_outfile_vespa, merge_taxonomy_counts(vespav, vespad))

#CSV.write(filtered_outfile_multi, filter_table(merge_taxonomy_counts(vespav, vespad), protist_filter))
#CSV.write(filtered_outfile_vespa, filter_table(merge_taxonomy_counts(vespav, vespad), protist_filter))