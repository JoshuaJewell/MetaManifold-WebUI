#!/usr/bin/env julia

include("run_cutadapt.jl")
include("dada2.jl")
include("merge_and_filter_taxa.jl")
include("run_vsearch.jl")

using CSV
using YAML
using .Cutadapt, .TaxonomyTableTools, .DADA2, .VSEARCH

## Tools loading (to move to new module at some point)
"""
    load_tools(config_path) -> Dict{String,String}

Read config/tools.yml and return a Dict of tool name => resolved path.
If the file does not exist, returns an empty Dict so tools fall back to PATH.
Paths with `@` are SSH remote paths (user@host:/path), the calling module is
responsible for routing those calls via SSH.
"""
function load_tools(config_path = joinpath(@__DIR__, "..", "config", "tools.yml"))
    isfile(config_path) || return Dict{String,String}()
    data = YAML.load_file(config_path)
    tools = Dict{String,String}()
    for (name, info) in data
        path = info isa Dict ? get(info, "path", nothing) : nothing
        path !== nothing && (tools[name] = string(path))
    end
    tools
end

# Return the configured binary path for `tool_key`, or `default` if not set.
tool_bin(tools, tool_key, default = tool_key) = get(tools, tool_key, default)

# Load resolved tool paths from config/tools.yml (empty Dict if not present)
tools = load_tools()

## Instantiate filesystem
# Root directories
data_dir = "./data"
config_dir = "./config"
output_dir = "./output"

# Cutadapt paths
fastq_input_dir = joinpath(data_dir, "fastq")
cutadapt_dir = joinpath(output_dir, "cutadapt")
primers_config = joinpath(config_dir, "primers.yml")

# DADA2/VSEARCH paths
fasta_outfile = joinpath(output_dir, "dada2/Tables/asvs.fasta")
reference_database = "./databases/pr2_version_5.0.0_SSU_dada2.fasta.gz"
vsearch_dir = joinpath(output_dir, "vsearch")

## Instantiate parameters
# Cutadapt parameters
primer_pairs = ["TarEuk", "Meta2"]
cutadapt_optional_args = "-m 200 --discard-untrimmed"

# DADA2 parameters
dada2_config_dir = joinpath(config_dir, "dada2.yml")

# VSEARCH parameters
vsearch_optional_args = "--id 0.75 --query_cov 0.8"

# Merge and filter (DADA2-VSEARCH) parameters
multiv = joinpath(output_dir, "vsearch/taxonomy.tsv")
multid = joinpath(output_dir, "dada2/Tables/taxonomy.csv")

protist_filter = joinpath(config_dir, "protist_filter.yml")

merged_outfile_multi = joinpath(output_dir, "merged/merged_multi.csv")
filtered_outfile_multi = joinpath(output_dir, "merged/protist_filtered_multi.csv")

## Main

cutadapt(primer_pairs, primers_config, fastq_input_dir, cutadapt_dir, optional_args = cutadapt_optional_args, cutadapt_bin = tool_bin(tools, "cutadapt"))

dada2(dada2_config_dir)

vsearch(fasta_outfile, reference_database, vsearch_dir, optional_args = vsearch_optional_args, vsearch_bin = tool_bin(tools, "vsearch"))

mkpath(dirname(merged_outfile_multi))

CSV.write(merged_outfile_multi, merge_taxonomy_counts(multiv, multid))
CSV.write(filtered_outfile_multi, filter_table(merge_taxonomy_counts(multiv, multid), protist_filter))