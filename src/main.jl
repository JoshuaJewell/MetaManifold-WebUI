#!/usr/bin/env julia

include("run_cutadapt.jl")
include("dada2.jl")
include("merge_and_filter_taxa.jl")

using CSV
using YAML
using .Cutadapt, .TaxonomyTableTools, .DADA2

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

## Instantiate parameters
# Cutadapt parameters
primer_pairs = ["TarEuk", "Meta2"]
optional_args = "-m 200 --discard-untrimmed"

# DADA2 parameters
dada2_config_dir = joinpath(config_dir, "dada2.yml")

# Merge and filter (DADA2-VSEARCH) parameters
multiv = joinpath(output_dir, "vsearch/taxonomy_multi_pool.tsv")
vespav = joinpath(output_dir, "vsearch/taxonomy_vespa_pool_fwdonly.tsv")
multid = joinpath(output_dir, "dada2/tax_counts_fasta_multi_pool.csv")
vespad = joinpath(output_dir, "dada2/tax_counts_fasta_vespa_pool_fwdonly.csv")

protist_filter = joinpath(config_dir, "protist_filter.yml")

merged_outfile_multi = joinpath(output_dir, "merged_multi.csv")
merged_outfile_vespa = joinpath(output_dir, "merged_vespa.csv")

filtered_outfile_multi = joinpath(output_dir, "protist_filtered_multi.csv")
filtered_outfile_vespa = joinpath(output_dir, "protist_filtered_vespa.csv")

## Main

#cutadapt(primer_pairs, primers_config, fastq_input_dir, cutadapt_dir,
#    optional_args = optional_args, cutadapt_bin = tool_bin(tools, "cutadapt"))

#dada2(dada2_config_dir)

#CSV.write(merged_outfile_multi, merge_taxonomy_counts(multiv, multid))
#CSV.write(merged_outfile_vespa, merge_taxonomy_counts(vespav, vespad))

#CSV.write(filtered_outfile_multi, filter_table(merge_taxonomy_counts(vespav, vespad), protist_filter))
#CSV.write(filtered_outfile_vespa, filter_table(merge_taxonomy_counts(vespav, vespad), protist_filter))