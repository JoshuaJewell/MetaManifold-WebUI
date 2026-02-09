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

multiv = "./vsearch/taxonomy_multi_pool.tsv"
vespav = "./vsearch/taxonomy_vespa_pool_fwdonly.tsv"
multid = "./DADA2/tax_counts_fasta_multi_pool.csv"
vespad = "./DADA2/tax_counts_fasta_vespa_pool_fwdonly.csv"

## Main

#cutadapt(primer_pairs, primers_path, fastq_input_dir, cutadapt_dir, optional_args = optional_args)

CSV.write("protist_filtered_multi_pool.csv", filter_table(merge_taxonomy_counts(multiv, multid),"./inputs/protist_filter_multiplex.yml"))
