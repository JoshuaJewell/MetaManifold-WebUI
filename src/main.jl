include("run_cutadapt.jl")
include("merge_and_filter_taxa.jl")

using CSV
using .Cutadapt, .TaxonomyTableTools

## Instantiate
primers_path = "./inputs/primers.yml"
primer_pairs = ["TarEuk", "Meta2"]
fastq_input_dir = "./inputs/fastq/"
cutadapt_dir = "./cutadapt/"
optional_args = "-m 200 --discard-untrimmed"

## Main

#cutadapt(primer_pairs, primers_path, fastq_input_dir, cutadapt_dir, optional_args = optional_args)

CSV.write("protist_filtered.csv", filter_table(merge_taxonomy_counts(),"./inputs/protist_filter.yml"))