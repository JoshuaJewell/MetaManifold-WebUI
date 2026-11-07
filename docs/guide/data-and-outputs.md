<!-- SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Input data and outputs

## Input data

Place paired-end FASTQ files under `data/{project_name}/` following Illumina naming:

```
data/MyProject/SampleName_*_L001_R1_001.fastq.gz
data/MyProject/SampleName_*_L001_R2_001.fastq.gz
```

For multi-run projects, nest runs in subdirectories. The server detects any directory containing `.fastq.gz` files as a leaf run and creates a matching project directory under `projects/{project_name}/`.

## Output structure

All outputs for a given run live under `projects/{project_name}/{run}/`:

```
projects/{project_name}/{run}/
|-- cutadapt/                    # Trimmed FASTQ pairs and logs
|   `-- logs/
|-- QC/
|   |-- fastqc/                  # Per-file FastQC HTML reports
|   |-- multiqc_report.html      # MultiQC summary across all samples
|   `-- logs/
|-- dada2/
|   |-- Tables/
|   |   |-- seqtab_nochim.csv    # ASV count table
|   |   |-- asvs.fasta           # ASV sequences
|   |   |-- asvs.csv             # ASV sequence index
|   |   |-- taxonomy.csv         # Taxonomy assignments
|   |   |-- taxonomy_bootstraps.csv
|   |   |-- taxonomy_combined.csv
|   |   |-- tax_counts.csv       # Taxonomy + per-sample counts
|   |   |-- asv_counts.csv       # Sequences + per-sample counts
|   |   `-- pipeline_stats.csv
|   |-- Figures/                 # Quality profile and error rate PDFs
|   |-- Checkpoints/             # RData checkpoints for stage resumption
|   `-- Logs/                    # Per-stage R logs
|-- cdhit/
|   |-- asvs.fasta               # Clustered ASV sequences
|   `-- asvs.fasta.clstr         # Cluster membership file
|-- swarm/
|   |-- otus.fasta               # OTU representative sequences
|   |-- otus.count_table.csv     # OTU count table (samples x OTUs)
|   `-- logs/
|-- vsearch/
|   |-- taxonomy.tsv             # Top-hit taxonomy assignments (ASV or OTU)
|   `-- logs/
`-- merged/
    |-- merged.csv               # Merged taxonomy + counts (all taxa)
    |-- protist_filter.csv       # Filtered subset (one per merge_taxa.filters entry)
    `-- results.duckdb           # DuckDB database for API queries
```

Phylogeny outputs sit beside the runs and in the shared library:

```
projects/{project_name}/
|-- trees/                       # Newick and jplace files for the tree viewer
`-- phylogeny/{id}/              # One placement
    |-- placement.json           # Name, reference tree, queries, settings
    |-- queries.fasta
    |-- align/combined.aln.fasta # References plus queries (MAFFT --addfragments)
    |-- trim/                    # combined.trim.fasta and the kept columns
    |-- place/placement.jplace   # RAxML EPA, with RAxML's own files
    |-- accumulate/accumulated.jplace
    `-- qc/, logs/, status.json, attestation.yml

reference_trees/{id}/            # One reference tree
|-- reference.json               # Name, description, settings
|-- references.fasta
|-- align/reference.aln.fasta
|-- trim/                        # reference.trim.fasta and the kept columns
|-- tree/reference.treefile      # With IQ-TREE's report, log and consensus tree
`-- qc/, logs/, status.json, attestation.yml
```
