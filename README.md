# MetaManifold

[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Julia ≥ 1.0](https://img.shields.io/badge/Julia-%E2%89%A51.0-9558B2?logo=julia)](https://julialang.org)
[![R ≥ 4.0](https://img.shields.io/badge/R-%E2%89%A54.0-276DC3?logo=r)](https://www.r-project.org)

A Julia pipeline for amplicon metabarcoding from raw paired-end Illumina reads to filtered, taxonomy-annotated ASV tables.

## Overview

This project aims to wrap a standard amplicon sequencing workflow into a single, configurable project. It handles multiplex primer trimming, amplicon denoising, taxonomy assignment, and taxonomic filtering.

**Pipeline stages**

| Step | Tool | Output |
|------|------|--------|
| Primer trimming | cutadapt | Trimmed FASTQ pairs per sample |
| Quality assessment | FastQC / MultiQC | HTML QC reports |
| Amplicon denoising | DADA2 (R) | ASV count table, FASTA, taxonomy |
| Taxonomy assignment | vsearch | Per-ASV taxonomy TSV |
| Clustering | cd-hit-est | Clustered ASV representatives |
| Merge + filter | Julia | Combined taxonomy/count table, filtered output |

## Prerequisites

- **Julia** ≥ 1.0 - installed automatically by `install.sh` if missing
- **R** ≥ 4.0 - required for the DADA2 stage
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
`ssh`. `config/defaults/tools.yml` contains the full config format.

### Creating a new project

```julia
projects = new_project("MyProject")
# Bootstraps projects/MyProject/ mirroring data/MyProject/ subdirectory structure - ensure "MyProject" name matches folder in data/.
# Creates dada2.yml, cutadapt.yml, vsearch.yml, cdhit.yml, merge_taxa.yml at each level.
```

Re-running `new_project` never overwrites existing configs. To reset a level to its parent's
settings, delete that config file and re-run.

Or manually:
```bash
mkdir -p projects/MyProject
cp config/defaults/dada2.yml projects/MyProject/dada2.yml
cp config/defaults/cutadapt.yml projects/MyProject/cutadapt.yml
# etc.
```

## Configuration

Global configuration lives in `config/`; per-project pipeline configs live in `projects/{name}/`. When per-run configs are missing, they are automatically populated with per-project configs, which are populated by global config if missing:

| File | Purpose |
|------|---------|
| `config/databases.yml` | Database URIs and optional local paths |
| `config/primers.yml` | Primer sequences and pair definitions |
| `config/filters/` | Directory of taxonomic filter configs |
| `config/tools.yml` | Tool binary paths (cutadapt, FastQC, MultiQC, vsearch, cd-hit-est) |
| `projects/{name}/cutadapt.yml` | Primer pair selection and cutadapt optional args |
| `projects/{name}/dada2.yml` | DADA2 pipeline parameters |
| `projects/{name}/vsearch.yml` | vsearch alignment thresholds |
| `projects/{name}/cdhit.yml` | cd-hit-est clustering threshold (optional stage) |
| `projects/{name}/merge_taxa.yml` | Which filter configs to apply in the merge step |

### Configuring databases (`config/databases.yml`)

This is the single place to manage DB URI's shared across all projects.

```yaml
databases:
  dir: "./databases"
  pr2:
    dada2:
      uri: "https://..."       # DADA2-format FASTA (downloaded on first use)
      local: ~                 # set to a local path to skip download
    vsearch:
      uri: "https://..."       # vsearch-format FASTA
      local: ~
```

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

Store all primer pairs in here and reference whichever combinations you need per project. Shared primers across pairs (same forward primer in two pairs) are automatically deduplicated in the `cutadapt` invocation since otherwise it complains a bit. If you need duplicates, you must create the same sequence under a different name.

### Configuring cutadapt (`projects/{name}/cutadapt.yml`)

Selects which primer pairs to apply and passes additional arguments to cutadapt. Inherits from `config/cutadapt.yml` if no per-project file exists.

```yaml
# Names must match keys in the Pairs section of config/primers.yml.
primer_pairs:
  - PrimerPair1
  - PrimerPair2

optional_args: "-m 200 --discard-untrimmed"
```

`optional_args` is passed verbatim to cutadapt after the primer arguments. Here, the `-m` flag sets the minimum read length after trimming; `--discard-untrimmed` drops reads with no primer match.

### Configuring DADA2 (`projects/{name}/dada2.yml`)

> **Performance tip:** For large datasets and reference databases, consider installing my [optimised DADA2 fork](https://github.com/JoshuaJewell/dada2) in place of the standard Bioconductor package. It provides CPU and Nvidia CUDA GPU acceleration for taxonomy assignment with no API or config changes required. This fork is experimental - if you encounter unexpected results, the standard Bioconductor release should be considered the reference implementation. Ultimately, I had to offload `assignTaxonomy()` to a remote server (`taxonomy.remote` setting).

```yaml
file_patterns:
  forward: "_R1_trimmed.fastq.gz"
  reverse: "_R2_trimmed.fastq.gz"
  sample_name_split: "_"       # character to split filenames on
  sample_name_index: 1         # which element is the sample name (1-based)
  mode: "paired"               # paired | forward | reverse

# Filter and trim - DADA2's filterAndTrim():
filter_trim:
  trunc_q: 2
  trunc_len: [220, 220]        # [forward, reverse]; first value used for single-end mode
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

# Taxonomy - assignTaxonomy() against the configured database:
taxonomy:
  database: pr2                # key into config/databases.yml
  multithread: true            # threads for assignTaxonomy()
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
  # to use the configured host. See config/defaults/dada2.yml for the full disclaimer.
  remote:
    host: ~                    # user@hostname
    rscript: "Rscript"         # path to Rscript on the server
    staging_dir: "/absolute/path/on/server"
    db_path: ~                 # absolute path to database on server (null = transfer local copy)

# Output filenames (all written to workspace_root/Tables/):
output:
  seq_table_prefix: "seqtab_nochim"
  fasta_prefix: "asvs"
  taxa_prefix: "taxonomy"
  combined_filename: "tax_counts.csv"
  asv_filename: "asv_counts.csv"
```

**Outputs written to `projects/{name}/{run}/dada2/Tables/`:**

| File | Contents |
|------|----------|
| `seqtab_nochim.csv` | Chimera-free ASV count table (samples x ASVs) |
| `asvs.fasta` / `asvs.csv` | ASV sequences with short identifiers (seq1, seq2, ...) |
| `taxonomy.csv` | Taxonomy assignments per ASV |
| `taxonomy_bootstraps.csv` | Bootstrap confidence values per rank |
| `taxonomy_combined.csv` | Taxonomy + bootstrap columns combined |
| `tax_counts.csv` | Taxonomy + per-sample counts |
| `asv_counts.csv` | ASV sequences + per-sample counts (no taxonomy) |
| `pipeline_stats.csv` | Read counts retained at each pipeline stage |

### Configuring vsearch (`projects/{name}/vsearch.yml`)

Controls the alignment thresholds used when assigning taxonomy against the reference database.

```yaml
optional_args: "--id 0.75 --query_cov 0.8"
```

`optional_args` is passed verbatim to `vsearch --usearch_global`. Key thresholds:
- `--id` - minimum sequence identity (0-1); lower values recover more hits at the cost of specificity
- `--query_cov` - minimum fraction of the query that must be aligned; filters partial matches

### Configuring cd-hit-est (`projects/{name}/cdhit.yml`)

Clustering step that collapses near-identical ASVs before taxonomy assignment.

```yaml
optional_args: "-c 0.9"
```

`optional_args` is passed verbatim to `cd-hit-est`. `-c` sets the sequence identity threshold for clustering (default 0.9 = 90%).

### Configuring merge_taxa (`projects/{name}/merge_taxa.yml`)

Controls which filter configs are applied when merging taxonomy and count tables. `merged.csv` (unfiltered) is always written; each entry in `filters` produces an additional filtered CSV.

```yaml
filters:
  - "protist_filter.yml"   # -> merged/protist_filter.csv
```

Each entry is a filename relative to `config/filters/`. Remove all entries (or set `filters: []`) to produce only the unfiltered `merged.csv`.

### Configuring taxonomic filtering (`config/filters/protist_filter.yml`)

The filter file controls the `filter_table()` step, which removes non-target taxa and remaps Supergroup labels for consistency with PR2 division names.

Place filter configs in `config/filters/` and reference them by filename in `merge_taxa.yml`. The default `protist_filter.yml` targets eukaryotic protists from PR2-annotated data:

```yaml
# Division -> Supergroup remapping.
# This block overrides Supergroup with the Division value where they should be equivalent.
mappings:
  Rhizaria: Rhizaria
  Alveolata: Alveolata
  Stramenopiles: Stramenopiles
  Hemimastigophora: Hemimastigophora
  Discoba: Discoba
  Metamonada: Metamonada
  Telonemia: Telonemia
  Ancyromonadida: Ancyromonadida

# Exclusion filters: rows where the named column contains the pattern are removed.
filters:
  - column: Domain
    pattern: Bacteria
  - column: Domain
    pattern: Archaea
  - column: Domain
    pattern: Eukaryota:plas
  - column: Domain
    pattern: Eukaryota:mito
  - column: Supergroup
    pattern: TSAR:chro
  - column: Subdivision
    pattern: Metazoa
  - column: Subdivision
    pattern: Fungi
  - column: Division
    pattern: Rhodophyta
  - column: Class
    pattern: Embryophyceae

# Remove rows with an empty or unassigned Domain field.
remove_empty_domain: true
```

> **Database compatibility:** Column names and patterns above are tuned for [PR2](https://pr2-database.org/). If you use a different reference database, update the column names to match that database's rank structure and adjust the `levels` list in `dada2.yml` accordingly.

## Usage

Place FASTQ files in `data/MyProject`. If you want to run multiple sets of FASTQ files, you can use subdirectories like `data/MyProject/Primer1`, `data/MyProject/Primer2`. `new_project()` will automatically generate a matching project directory structure, nothing in `data/` is ever overwritten. Edit `src/main.jl` to set your project name and run:

```julia
dbs      = ensure_databases(databases_config)
projects = new_project("MyProject")

const r_lock = ReentrantLock()

Threads.@threads for project in projects
    multiqc(project.data_dir, joinpath(project.dir, "QC"))

    trimmed = cutadapt(project)

    asvs = lock(r_lock) do
        dada2(project, trimmed, taxonomy_db = dbs["pr2_dada2"])
    end
    # asvs = lock(r_lock) do; cdhit(project, asvs); end  # optional clustering

    tax    = vsearch(project, asvs, dbs["pr2_vsearch"])
    merged = merge_taxa(project, asvs, tax)
end
```

Each stage returns a wire type (`TrimmedReads`, `ASVResult`, `TaxonomyHits`, `MergedTables`) and skips automatically if outputs are already up to date relative to their inputs (mtime-based). This means that rerunning pipeline after config change will only perform the minimum necessary steps to output.

Run from the project root:

```bash
julia --project=. src/main.jl
```

Or with threads:

```bash
julia -t 2 --project=. src/main.jl
```

### Input data

Place paired-end FASTQ files under `data/{project_name}/` following Illumina naming:

```
data/MyProject/SampleName_*_L001_R1_001.fastq.gz
data/MyProject/SampleName_*_L001_R2_001.fastq.gz
```

For multi-run projects, nest runs in subdirectories - `new_project` will detect any directory containing `.fastq.gz` files as a leaf run and create a matching project directory under `projects/{project_name}/`.

### Output structure

All outputs for a given run live under `projects/{project_name}/{run}/`:

```
projects/{project_name}/{run}/
├── cutadapt/                    # Trimmed FASTQ pairs and logs
│   └── logs/
│       ├── cutadapt_primer_trimming_stats.txt
│       └── cutadapt_trimmed_percentage.txt
├── QC/
│   ├── fastqc/                  # Per-file FastQC HTML reports
│   ├── multiqc_report.html      # MultiQC summary across all samples
│   └── logs/
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
│   ├── Checkpoints/             # RData checkpoints for stage resumption
│   └── Logs/                    # Per-stage R logs
├── cdhit/
│   ├── asvs.fasta               # Clustered ASV sequences
│   └── asvs.fasta.clstr         # Cluster membership file
├── vsearch/
│   ├── taxonomy.tsv
│   └── logs/
└── merged/
    ├── merged.csv               # Merged vsearch taxonomy + ASV counts (all taxa)
    └── protist_filter.csv       # Filtered subset (one file per entry in merge_taxa.yml)
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