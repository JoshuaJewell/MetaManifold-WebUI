# MetaManifold

A Julia pipeline for amplicon metabarcoding from raw paired-end Illumina reads to filtered, taxonomy-annotated ASV tables.

## Overview

This project aims to wrap a standard amplicon sequencing workflow into a single, configurable project. It handles multiplex primer trimming amplicon denoising, taxonomy assignment, and taxonomic filtering.

**Pipeline stages**

| Step | Tool | Output |
|------|------|--------|
| Primer trimming | cutadapt | Trimmed FASTQ pairs per sample |
| Quality assessment | FastQC / MultiQC | HTML QC reports |
| Amplicon denoising | DADA2 (R) | ASV count table, FASTA, taxonomy |
| Taxonomy assignment | vsearch | Per-ASV taxonomy TSV |
| Clustering | cd-hit-est | Clustered ASV representatives |
| Merge + filter | Julia | Combined taxonomy/count table, protist-filtered output |

## Prerequisites

- **Julia** ≥ 1.0 - installed automatically by `install.sh` if missing
- **R** ≥ 4.0 - must be installed before running `install.sh`
  - Ubuntu/Debian: `sudo apt install r-base`
  - macOS: `brew install r` or [CRAN package](https://cran.r-project.org/bin/macosx/)

## Installation

```bash
git clone https://github.com/JoshuaJewell/MetaManifold.git
cd MetaManifold
bash install.sh
```

`install.sh` will check for Julia and R, install Julia and R dependencies, locate or download each external tool (cutadapt, FastQC, MultiQC, vsearch, cd-hit-est). Tool paths can be configured manually in `config/tools.yml`, including by SSH if you wish to use a server-hosted binary.

To update:
```bash
bash install.sh --update
```

### Tool paths

`install.sh` generates `config/tools.yml`. You can edit it manually, for example, to point a tool at a remote server:

```yaml
vsearch:
  path: "user@bioserver:/home/user/software/vsearch"
```

When a remote SSH path is set, the pipeline routes that tool's invocations through
`ssh`. `config/tools.example.yml` contains full config format.

## Configuration

All configuration lives in `config/`:

| File | Purpose |
|------|---------|
| `primers.yml` | Primer sequences and pair definitions |
| `dada2.yml` | DADA2 pipeline parameters (truncation, error model, taxonomy DB, etc.) |
| `protist_filter.yml` | Taxonomic filtering rules for protist output |
| `tools.yml` | Tool paths (cutadapt, FastQC, MultiQC, vsearch, cd-hit-est) |

### Defining primer pairs `primers.yml`

 Maps primer names to sequences and defines which forward/reverse sequences constitute a pair:

```yaml
Forward:
  PrimerF: "CCAGCASCYGCGGTAATTCC"

Reverse:
  Primer1R: "ACTTTCGTTCTTGATYRA"
  Primer2R: "DCTKTCGTYCTTGATYRA"

Pairs:
  - PrimerPair1:
      - PrimerF
      - Primer1R
  - PrimerPair2:
      - PrimerF
      - Primer2R
```

This allows you to referenece one or multiple primer pairs within the pipeline by name. This allows you to store all primer pairs you often work on in one place and then refer to them in whichever combinations you need as and when. Shared primers across pairs (same forward in two pairs, as the above example) are automatically deduplicated in the cutadapt invocation since otherwise it complains a bit.

### Configuring DADA2 (`dada2.yml`)

> **Performance tip:** For large datasets and reference databases, consider installing the [optimised DADA2 fork](https://github.com/JoshuaJewell/dada2) in place of the standard Bioconductor package. It provides acceleration for CPU and Nvidia CUDA GPUs for taxonomy assignment, with no changes to the API or configuration required. This fork is experimental, so if you encounter unexpected results, the standard Bioconductor release should be considered the reference implementation. It was too experimental for me, so opted to offload assignTaxonomy() to server in dada2.yml...

Copy `config/dada2.example.yml` to `config/dada2.yml` and edit it. It contains your local paths and optionally remote server credentials and should never be shared.

```yaml
workspace:
  root: "./output/dada2/"      # output directories are created here
  input_dir: "./output/cutadapt/" # trimmed FASTQ input (cutadapt output)

file_patterns:
  forward: "_R1_trimmed.fastq.gz"
  reverse: "_R2_trimmed.fastq.gz"
  sample_name_split: "_"       # character to split filenames on
  sample_name_index: 1         # which element is the sample name (1-based)
  mode: "paired"               # paired | forward | reverse

# Filter and trim - DADA2's filterAndTrim():
filter_trim:
  trunc_q: 2
  trunc_len: [220, 220]        # [forward, reverse]; first value used for single-end
  max_ee: [3, 3]               # maximum expected errors in F and R reads
  min_len: 175
  max_n: 0
  match_ids: true
  rm_phix: true

# Denoising - learnErrors() and dada():
dada:
  seed: 123
  nbases: 200000000
  max_consist: 15
  pool_method: "pseudo"        # none | pseudo | true

# Merging - mergePairs(), paired mode only:
merge:
  min_overlap: 20
  max_mismatch: 0
  trim_overhang: true

# ASV length filtering and chimera removal:
asv:
  band_size_min: 200           # null to skip length filtering
  band_size_max: 430
  denovo_method: "consensus"   # consensus | pooled | per-sample

# Reference databases - downloaded and cached on first use:
databases:
  dir: "./databases"
  pr2:
    dada2:
      uri: "https://..."       # DADA2-format FASTA
      local: ~                 # set to a local path to skip download
    vsearch:
      uri: "https://..."       # vsearch-format FASTA
      local: ~

# Taxonomy - assignTaxonomy() against the configured database:
taxonomy:
  database: pr2                # key into databases: section above
  multithread: 4               # threads for assignTaxonomy(); higher = more RAM
  min_boot: 0                  # minimum bootstrap confidence to retain (0-100)
  levels:
    - "Domain"
    - "Supergroup"
    - "Division"
    - "Subdivision"
    - "Class"
    - "Order"
    - "Family"
    - "Genus"
    - "Species"

  # Optional: offload the memory-intensive assignTaxonomy() step to a remote
  # server via SSH. Omit or set host to null to run locally.
  # DISCLAIMER: You are solely responsible for ensuring you have authorisation
  # to use the configured host. See dada2.example.yml for the full disclaimer.
  remote:
    host: ~                    # user@hostname
    rscript: "Rscript"         # path to Rscript on the server
    staging_dir: "/absolute/path/on/server"
    db_path: ~                 # absolute path to database on server (null = transfer local copy)

# Output filenames (all written to workspace.root/Tables/):
output:
  seq_table_prefix: "seqtab_nochim"
  fasta_prefix: "asvs"
  taxa_prefix: "taxonomy"
  combined_filename: "tax_counts.csv"
  asv_filename: "asv_counts.csv"
```

**Outputs written to `workspace.root/Tables/`:**

| File | Contents |
|------|----------|
| `seqtab_nochim.csv` | Chimera-free ASV count table (samples × ASVs) |
| `asvs.fasta` / `asvs.csv` | ASV sequences with short identifiers (seq1, seq2, …) |
| `taxonomy.csv` | Taxonomy assignments per ASV |
| `taxonomy_bootstraps.csv` | Bootstrap confidence values per rank |
| `taxonomy_combined.csv` | Taxonomy + bootstrap columns combined |
| `tax_counts.csv` | Taxonomy + per-sample counts |
| `asv_counts.csv` | ASV sequences + per-sample counts (no taxonomy) |
| `pipeline_stats.csv` | Read counts retained at each pipeline stage |

### Configuring taxonomic filtering (`protist_filter.yml`)

The filter file controls the `filter_table()` step, which removes non-target taxa
and remaps Supergroup labels for consistency with PR2 division names.

> **Database compatibility:** The default config is tuned specifically for [PR2](https://pr2-database.org/). Column names (`Supergroup`, `Division`, `Subdivision`) and pattern strings (`Eukaryota:plas`, `TSAR:chro`, etc.) reflect PR2's rank structure and nomenclature. If you use a different reference database, you will need to update both the column names here and the `levels` list in `dada2.yml` to match that database's ranks.

**Division → Supergroup remapping**

PR2 assigns a `Supergroup` label that can be coarser than `Division` for some lineages. The `mappings` block overrides `Supergroup` with the `Division` value where they should be treated as equivalent:

```yaml
mappings:
  Rhizaria: Rhizaria
  Alveolata: Alveolata
  Stramenopiles: Stramenopiles
  # add further Division: Supergroup pairs as needed
```

**Exclusion filters**

Each entry specifies a column and a substring. Rows where that column contains the substring (partial match) are removed:

```yaml
filters:
  - column: Domain
    pattern: Bacteria
  - column: Domain
    pattern: Eukaryota:plas
  - column: Domain
    pattern: Eukaryota:mito
  - column: Subdivision
    pattern: Metazoa
  - column: Subdivision
    pattern: Fungi
  # add further column/pattern pairs to exclude additional lineages
```

## Usage

Edit `src/main.jl` to set your  parameters, then uncomment or add pipeline steps you want to run, e.g.:

```julia
# Remove primers from Illumina reads
cutadapt(primer_pairs, primers_config, fastq_input_dir,
         cutadapt_dir, optional_args = optional_args)

# Run dada2 denoising, chimera removal, etc. in accordance with config
dada2(dada2_config)

# Use vsearch local alignment to create a taxonomy table
vsearch(fasta)

# Merge vsearch and idtax tables by ASV
merged = merge_taxonomy_counts(vsearch_tsv, dada2_csv)

# Output table of all ID'd reads and table of reads for protists 
CSV.write(merged_outfile_name, merged)
CSV.write(filtered_outfile_name, filter_table(merged, protist_filter_path))
```

Run from the project root:

```bash
julia --project=. src/main.jl
```

### Input data

Place paired-end FASTQ files in `data/fastq/` following Illumina naming:

```
SampleName_*_L001_R1_001.fastq.gz
SampleName_*_L001_R2_001.fastq.gz
```

### Output structure

```
output/{project_name}/
├── cutadapt/                    # Trimmed FASTQ pairs and logs
├── FastQC/                      # FastQC HTML reports
├── dada2/
│   ├── Tables/
│   │   ├── seqtab_nochim.csv    # ASV count table
│   │   ├── asvs.fasta           # ASV sequences
│   │   ├── asvs.csv             # ASV sequence index
│   │   ├── taxonomy.csv         # Taxonomy assignments
│   │   ├── taxonomy_bootstraps.csv
│   │   ├── taxonomy_combined.csv
│   │   ├── tax_counts.csv       # Taxonomy + per-sample counts
│   │   ├── asv_counts.csv       # Sequences + per-sample counts
│   │   └── pipeline_stats.csv
│   ├── Figures/                 # Quality profile and error rate PDFs
│   └── Checkpoints/             # RData checkpoints for stage resumption
├── vsearch/
│   └── taxonomy.tsv
└── merged/
    ├── merged_multi.csv         # Merged DADA2 + vsearch taxonomy + counts
    └── protist_filtered_multi.csv
```

## Third-party tools

This project orchestrates the following tools. Each is fetched from its upstream source by `install.sh` and is subject to its own license - no third-party binaries are included in this repository.

| Tool | License | Source |
|------|---------|--------|
| [cutadapt](https://github.com/marcelm/cutadapt) | MIT | PyPI |
| [FastQC](https://github.com/s-andrews/FastQC) | GPL v3 | Babraham Bioinformatics |
| [MultiQC](https://github.com/MultiQC/MultiQC) | GPL v3 | PyPI |
| [vsearch](https://github.com/torognes/vsearch) | GPL v3 | GitHub Releases |
| [cd-hit](https://github.com/weizhongli/cdhit) | GPL v2+ | GitHub Releases / apt |
| [DADA2](https://benjjneb.github.io/dada2/) ([optimised fork](https://github.com/JoshuaJewell/dada2)) | LGPL v3 | Bioconductor |

## Acknowledgements

This pipeline draws on the following prior work:

- **Frédéric Mahé** - [Fred's metabarcoding pipeline](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline) informed the overall workflow architecture: the sequencing of primer trimming, vsearch-based taxonomy assignment, and final table merge/filter stages. Fred's pipeline uses swarm + OTUs where MetaManifold uses DADA2 + ASVs.
- **Benjamin J. Callahan _et al._** - [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html), used under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), which the `dada2.r` module is based on.

The following colleagues at the **Department of Parasitology, Charles University** (Faculty of Science, BIOCEV, Vestec, Czech Republic) contributed to this work:

- **Mgr. Jiří Novák** (supervisor) - scripts from which several modules and configurations were adapted:
  - `run_cutadapt.jl`
  - `merge_and_filter_taxa.jl`
  - `dada2.yml`
  - `protist_filter.yml`
- **Bc. Ekaterina Kandaurova** - designed the primer pairs in `primers.yml`.
- **doc. Mgr. Vladimír Hampl** - provided laboratory access and resources.

## License

Copyright © 2026 Joshua Benjamin Jewell.

Source code is licensed under the [GNU Affero General Public License v3.0](LICENSE).

This documentation (README.md) is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).