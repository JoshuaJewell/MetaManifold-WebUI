include("run_cutadapt.jl")

using .Cutadapt

## Instantiate
primers_path = "./inputs/primers.yml"
primer_pairs = ["TarEuk", "Meta2"]
fastq_input_dir = "./inputs/fastq/"
cutadapt_dir = "./cutadapt/"
optional_args = "-m 200 --discard-untrimmed"

## Main

cutadapt(primer_pairs, primers_path, fastq_input_dir, cutadapt_dir, optional_args = optional_args)