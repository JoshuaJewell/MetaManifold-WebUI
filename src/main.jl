#!/usr/bin/env julia

include("databases.jl")
include("run_cutadapt.jl")
include("dada2.jl")
include("merge_and_filter_taxa.jl")
include("run_vsearch.jl")

using CSV
using YAML
using .Databases, .Cutadapt, .TaxonomyTableTools, .DADA2, .VSEARCH

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

# Project name - all stage outputs live under output/{project_name}/
# Re-running with the same project_name overwrites previous results.
project_name = "project"
project_dir  = joinpath(output_dir, project_name)

# Cutadapt paths
fastq_input_dir = joinpath(data_dir, "fastq")
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
multiv = joinpath(vsearch_dir, "taxonomy.tsv")
multid = joinpath(dada2_dir, "Tables/taxonomy.csv")

protist_filter = joinpath(config_dir, "protist_filter.yml")

# Download/use database
dbs = ensure_databases(dada2_config_path)

## Main
#cutadapt(primer_pairs, primers_config, fastq_input_dir, trimmed_dir, optional_args = cutadapt_optional_args, cutadapt_bin  = tool_bin(tools, "cutadapt"))

#dada2(dada2_config_path, input_dir = trimmed_dir, workspace_root = dada2_dir, taxonomy_db = dbs["pr2_dada2"])

vsearch(fasta_outfile, dbs["pr2_vsearch"], vsearch_dir, optional_args = vsearch_optional_args, vsearch_bin = tool_bin(tools, "vsearch"))

CSV.write(merged_outfile_multi, merge_taxonomy_counts(multiv, multid))
CSV.write(filtered_outfile_multi, filter_table(merge_taxonomy_counts(multiv, multid), protist_filter))