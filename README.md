# MetaManifold

[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Julia 1.12](https://img.shields.io/badge/Julia-1.12-9558B2?logo=julia)](https://julialang.org)
[![R >= 4.0](https://img.shields.io/badge/R-%E2%89%A54.0-276DC3?logo=r)](https://www.r-project.org)
[![CI](https://github.com/JoshuaJewell/MetaManifold-WebUI/actions/workflows/ci.yml/badge.svg)](https://github.com/JoshuaJewell/MetaManifold-WebUI/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JoshuaJewell/MetaManifold-WebUI/graph/badge.svg?token=20F1VLF590)](https://codecov.io/gh/JoshuaJewell/MetaManifold-WebUI)

MetaManifold wraps standard amplicon sequencing workflows into a single configurable Julia orchestrator: from raw paired-end Next Generation Sequencing reads through denoising, taxonomy assignment and taxonomic filtering, with interactive configuration and analysis in the browser.

<p align="center">
  <img src=".github/screenshots/hero.png" width="850" alt="A run page showing per-stage status, launch controls, run-level configuration and the read funnel">
</p>

## Overview

MetaManifold consists of a Julia backend (pipeline engine + REST API) and a TypeScript/React frontend. The pipeline runs FastQC, MultiQC, cutadapt, DADA2, swarm, vsearch, and cd-hit-est under the hood, and places ASVs on reference trees with MAFFT, trimAl, IQ-TREE, RAxML and gappa; results are stored in per-run DuckDB databases and served to the frontend as interactive Plotly charts and filterable tables. Pipeline configuration is editable in the web UI at every cascade level (see [Configuration](#configuration)).

**Pipeline stages**

```
Raw FASTQs  (data/{study}/[{group}/]{run}/*.fastq.gz)
      |
   cutadapt, primer trimming
      |
      +----------------------------+
      |                            |
   DADA2*, ASV;               swarm*, OTU;
   filter & trim              merge pairs
   learn error rates          dereplicate
   denoise + merge            chimera filter
   length filter              cluster OTUs
   chimera removal                 |
   taxonomy assign*                |
      |                            |
   cd-hit-est*, demultiplex        |
      |                            |
   vsearch*                   vsearch, global alignment
      |                            |
      +----------------------------+
      |
   merge_taxa;
   join ASV tables*
   join counts-taxonomy
   apply filters
      |
   DuckDB results store          Reference sequences
      |                                |
      |                          reference tree*;
      |                          MAFFT align
      |                          trimAl trim
      |                          IQ-TREE tree
      |                                |
      +---------------+----------------+
                      |
               placement*;
               MAFFT add queries
               trimAl trim
               RAxML EPA place
               gappa accumulate
                      |
               jplace trees
```
*optional

**Analysis**

Once a run completes, analysis is performed on request through the web UI, both per run and across runs within a study:

- Alpha diversity (richness, Shannon, Simpson), per sample or as cross-run comparison boxplots with significance testing and optional lines joining paired samples
- Taxonomic composition bar charts at a chosen rank, relative or absolute
- Organism-composition charts that classify ASVs/OTUs into biological categories (a "Composition" view, e.g. protozoa, helminths, fungi, host)
- Taxon overlap across runs as proportional Euler or UpSet plots
- Pipeline stage read-count summaries
- NMDS ordination (Bray-Curtis, via R/vegan)
- PERMANOVA, with PERMDISP to check for differences in dispersion (via R/vegan)
- Differential abundance between two conditions: a negative-binomial model per taxon (via R/MASS) with Benjamini-Hochberg adjusted p-values (see [Differential abundance](#differential-abundance))

Counts may be normalised before analysis (none, rarefaction to a fixed or auto-resolved depth, or scaling with ranked subsampling (SRS), which reaches the depth by scaling counts and keeps more of the community structure), and contamination-flagged taxa may be included or excluded. Ordinations can also apply a Hellinger transform before the Bray-Curtis dissimilarity. Composition charts can be drawn as a grid faceted on any two of run, group and sub-group. All analysis charts are returned as Plotly JSON and rendered interactively in the browser.

## Documentation

Longer documentation for users, maintainers, developers and theorists is kept as a wiki in
[metadatastician/MetaManifold-Evidence, under `docs/wikis/`](https://github.com/metadatastician/MetaManifold-Evidence/tree/main/docs/wikis).
Start at its [Home page](https://github.com/metadatastician/MetaManifold-Evidence/blob/main/docs/wikis/Home.md).

## Prerequisites

- **Julia** 1.12 (installed by `install.sh` if missing, and pinned to the Manifest version for this directory)
- **R** >= 4.0 (required for the DADA2 stage and NMDS/PERMANOVA analysis)
  - Ubuntu/Debian: `sudo apt install r-base`
  - macOS: `brew install r` or [CRAN package](https://cran.r-project.org/bin/macosx/)
- **bun** builds the frontend. The installer uses a bun on PATH at the version pinned in `config/defaults/tool_versions.yml`, or downloads that release into `bin/`.

## Installation

```bash
git clone https://github.com/JoshuaJewell/MetaManifold-WebUI.git
cd MetaManifold-WebUI
bash install.sh
```

`install.sh` installs Julia (via juliaup) if missing, checks for R, then hands
off to `install.jl`, which installs the Julia dependencies, locates or downloads
each external tool (cutadapt, FastQC, MultiQC, vsearch, cd-hit-est, and the
phylogeny tools MAFFT, trimAl, IQ-TREE, RAxML and gappa), and reproduces the R
environment from `renv.lock`. It ends with a summary listing every dependency as
`OK`, `SKIPPED`, or `ACTION NEEDED`. Tool paths can be set in
`config/tools.yml`. MAFFT, IQ-TREE and RAxML can be skipped on a machine that
sends those stages to a server (see
[Threads and remote execution](#threads-and-remote-execution-top-level-of-pipelineyml)).

Installing R and the system `-dev` headers the R packages compile against needs
root. If R is missing, the installer stops and prints how to install it. If
headers are missing, the summary prints the command that installs them; run it
and rerun `install.sh`.

To update:
```bash
bash install.sh --update
```

### Reproducing the R environment

The R-side dependencies (DADA2, vegan, and their transitive packages) are
pinned with [`renv`](https://rstudio.github.io/renv/); `renv.lock` is the source
of truth and the committed `.Rprofile` activates the project-local library on
any `R`/`Rscript` invocation from the repository root. `install.sh` runs the
restore; to run it by hand:

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
```

To add or upgrade a package, install it inside the project
(`renv::install(...)`) and record the change with `renv::snapshot()`.

### Tool paths

`install.sh` generates `config/tools.yml`, one local path per tool:

```yaml
vsearch:
  path: "/home/user/software/vsearch"
```

`config/defaults/tools.yml` contains the full config format. Stages that run on a
server are configured in the `remote` block of `pipeline.yml`.

## Quick start

### 1. Place paired-end FASTQs under data/
Put `.fastq.gz` files into `data/MyProject/run_A`. 

### 2. Start the server
`bash start.sh`

The installer builds the frontend into `web/dist`, and `install.sh --update` rebuilds it. `start.sh` rebuilds it whenever a frontend source is newer than the build; `BUILD=1 bash start.sh` forces a rebuild.

### 3. Open http://localhost:8080

The web UI lets you create studies, configure pipeline parameters, launch runs, and explore results. All state lives in the filesystem under `data/` and `projects/`.

| Environment variable      | Default           | Description   |
| ------------------------- | ----------------- | ------------- |
| `JULIA_METAMANIFOLD_PORT` | `8080`            | Server port   |
| `JULIA_METAMANIFOLD_ROOT` | working directory | Project root  |
| `JULIA_THREADS`           | `8`               | Julia threads |

Or run the Julia server directly:

```bash
julia --project=. scripts/serve.jl
```

The server is part of the `MetaManifold` package, so the first start after an
install or a source change compiles the package. That compile runs a workload of
the requests a first page visit makes, which takes a few minutes and caches
their native code. While developing, skip the workload with a
`LocalPreferences.toml` (git-ignored) in the repository root:

```toml
[MetaManifold]
precompile_workload = false
```

The server then makes those requests itself in the background as it starts, so the first page visit does not wait on compilation.

## Web interface

The browser interface is the primary way to drive MetaManifold. Beyond creating studies and launching runs, it offers the following:

### Editing configuration

Every pipeline setting can be edited in the UI without touching a YAML file. Configuration is presented as collapsible accordion sections (study design, primer trimming, DADA2 denoising and taxonomy, OTU clustering, analysis) and can be set at any level of the cascade: the instance-wide defaults, a study, a group, or an individual run. Edits at a finer level override coarser ones (see [Configuration](#configuration) for the cascade rules). When a setting changes, the affected pipeline stages are flagged as stale, with a tooltip listing the changed keys and their level.

### Pipeline runs and outputs

A run's page is the working surface for that run. It is where the run-level configuration above is edited, where the full pipeline or any individual stage (including the DADA2 substages) is launched, and where each stage's status is shown. Long-running jobs report progress live through a server-sent event stream, and the jobs panel lets you watch or cancel them. As stages complete, their outputs become available through the views that follow: quality reports, the results table and composition.

### QC

Raw-read QC (FastQC aggregated by MultiQC) is embedded in each run's page with a rerun control. The DADA2 quality, denoising, merging and taxonomy diagnostics sit beside their stage settings.

<p align="center">
  <img src=".github/screenshots/qc.png" width="800" alt="The MultiQC report embedded in a run's QC tab">
</p>

### Results explorer

Each run's merged taxonomy-and-count table, and any derived tables, can be browsed. The table supports per-column filtering (text search, numeric range, include/exclude lists) and a global text filter, column sorting, configurable pagination, and column visibility toggles including taxonomy-source presets (VSEARCH-only, DADA2-only, or all) and a switch for the per-sample count columns. Sequences carry BLAST links, and OTU rows can be expanded to their constituent sequences. Frequently used filters can be saved as named presets and reapplied; filtered tables can be saved back into the run or exported as `.xlsx`.

<p align="center">
  <img src=".github/screenshots/results-explorer.png" width="800" alt="Results explorer table with per-column filters and taxonomy-source column presets">
</p>

### Phylogenetic placement

Reference trees (under SYSTEM) are built once from uploaded reference sequences
and used by any study. A study's Trees tab places its ASVs on one of them: pick
runs and, for pooled runs, subgroups, a results table, a rank and taxa, or upload
a FASTA. Each step shows where it runs, its state and its log, and the align and
trim steps show a QC view of the alignment. A finished placement adds the
reference tree, the placement jplace and the accumulated jplace to the study's
tree list. See
[Configuring phylogenetic placement](#configuring-phylogenetic-placement-phylogeny-in-pipelineyml).

### Tree viewer

Newick and jplace files in a study's tree list open in an interactive viewer:
rectangular or circular layouts, rerooting, ladderizing, collapsing and naming
clades, per-clade colours and fonts, placement symbols sized by count, and
support values carried over from another tree. The view is saved beside the
tree and exports as SVG, PNG or Newick. Reference trees open in the same viewer
from the library page.

<p align="center">
  <img src=".github/screenshots/tree-viewer.png" width="800" alt="The tree viewer showing a placement tree with collapsed clades, bootstrap support and placement symbols">
</p>

### Composition

The composition view classifies each ASV/OTU into a biological category and renders per-sample or pooled stacked bar charts. Category sets live in `config/composition.yml` (the bundled `default` set covers protozoa, helminths, fungi, host, plants, and invertebrates); each category references a named taxonomic filter from the `filters:` library in that same file. Both are editable from the Compositions page under SYSTEM in the sidebar. A category summary precedes the chart, and a quality filter can cap the number of unresolved taxonomic placeholders admitted.

<p align="center">
  <img src=".github/screenshots/composition.png" width="800" alt="Organism-category composition of two runs as stacked bars">
</p>

### Differential abundance

The Differential Abundance tab of a study's analysis workspace tests each taxon
for a difference in abundance between two conditions. Select exactly two runs or
pooled-run sub-groups, a results table and a rank. The first condition is the
reference, so a positive fold change means the taxon is more abundant in the
second. The rank list shows only the ranks the selected runs share. With
aggregation on, a pooled run counts as one condition rather than one per
sub-group.

Counts are summed per taxon at the chosen rank, with rows that have no name at
that rank pooled as `Unclassified`. This happens after the read-count bounds and
any `analysis.exclude_categories` entry that applies to `differential`. Each
taxon is fitted with its own negative-binomial GLM, `MASS::glm.nb(y ~ group +
offset(log(size factor)))`, on the raw integer read counts. There are no
pseudocounts. The group coefficient gets a Wald test, and Benjamini-Hochberg
adjustment runs over the taxa that produced a p-value. Size factors come from
the whole count matrix and are centred to a geometric mean of 1. The code is in
`src/analysis/differential.jl`, and the route in `src/server/routes/analysis.jl`.

Settings live under `analysis.differential` (see
[Configuring analysis](#configuring-analysis-analysis-in-pipelineyml)):

- `offset`: `tss` (default) uses each sample's total reads. `rle` uses the
  median of ratios over the taxa with reads in every sample.
- `min_prevalence`: a fraction from 0 to 1 (default 0). A taxon with reads in a
  smaller fraction of the samples is not fitted.

The result is a volcano plot and a table. The plot draws log2 fold change
against -log10 of the raw p-value, coloured by whether the adjusted p-value is
below 0.05. Taxa without an adjusted p-value are left off the plot but stay in
the table. Each table row has a status:

- `ok`: the model fitted.
- `boundary`: the model fitted, but the dispersion parameter theta reached a
  bound (at least 1e7, so effectively Poisson, or at most 1e-8). These taxa
  keep their p-value and stay in the adjusted family.
- `failed`: no p-value, with the reason. Causes are counts that are constant
  across samples, an error from `glm.nb`, a fit that did not converge, an
  aliased group coefficient, or a non-finite estimate.
- `filtered`: below `min_prevalence`, so not fitted.

A failed or filtered taxon never gets a stand-in p-value, and it is not counted
in the adjustment. The table downloads as `differential_abundance.csv` and can
be added to the report. The CSV starts with `#` lines that record the
comparison, method and settings.

The test refuses, with an explicit error, when:

- the selection is not exactly two conditions, or both are the same;
- the two conditions resolve to different ranks;
- a condition has no results table, no taxonomy columns or no sample columns,
  or the read-count bounds remove all its samples;
- `offset` or `min_prevalence` is invalid;
- a count is not a non-negative integer, a condition has no samples, or there
  are fewer than 3 samples in total;
- `tss` meets a sample with no reads, or `rle` finds no taxon with reads in
  every sample (use `tss`);
- the R package MASS cannot be loaded;
- no taxon could be fitted at all.

### Figures and report

The Figures tab lays analysis charts out on a page (A4, letter or a journal
column width), in lettered groups of panels with shared colours, fonts and
sizes, and exports it as PDF, PNG or TIFF at a chosen resolution. The
Publication Tables tab builds one table at a time: rows of taxa at a rank or of
categories from a Compositions set, columns of runs, sub-groups or samples, and
any of ASVs, reads and percentage per column, with an optional difference
between two sub-groups. A table can be copied for Word or downloaded as `.xlsx` or
`.csv`. Charts, figures, tables and trees can be added to the study's Report tab,
which downloads them together as one ZIP with a captions file.

<p align="center">
  <img src=".github/screenshots/figures.png" width="800" alt="The figure builder with page settings on the left and a two-panel figure of alpha diversity comparisons">
</p>

## Configuration

Pipeline settings use a cascade: each level overrides the one above it, and any key you omit is inherited from the nearest ancestor. The merged result is written to each run's `run_config.yml`, which records the settings the run used.

Settings can be edited in the web UI (per-study, per-group, or per-run) or in the YAML files.

| File | Purpose |
|------|---------|
| `config/defaults/` | Canonical defaults for every setting; do not edit |
| `config/composition.yml` | Composition library: named taxonomic filters and the category sets that reference them |
| `config/presets/` | Saved table-view filter presets, written from the Tables view |
| `config/databases.yml` | Database URIs and optional local paths. Editable from the Databases page under SYSTEM in the sidebar |
| `config/primers.yml` | Primer sequences and pair definitions. Editable from the Primers page under SYSTEM in the sidebar |
| `config/tools.yml` | Tool binary paths (cutadapt, FastQC, MultiQC, vsearch, cd-hit-est, MAFFT, trimAl, IQ-TREE, RAxML, gappa) |
| `config/pipeline.yml` | Machine-level overrides (lowest user-editable precedence) |
| `data/{name}/pipeline.yml` | Study-level overrides |
| `data/{name}/{group}/pipeline.yml` | Group-level overrides (intermediate directories) |
| `data/{name}/{run}/pipeline.yml` | Run-level overrides (highest precedence) |
| `projects/{name}/{run}/run_config.yml` | Generated merged config (provenance); do not edit |

Each `pipeline.yml` stub is created with a comment block explaining that level's role. Write only the keys you want to change; omit the rest.

### Configuring databases (`config/databases.yml`)

This is the single place to manage DB URIs shared across all projects.

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

Edit this on the Databases page under SYSTEM in the sidebar, or in the YAML. The page edits the shared cache directory (`dir`) and, per database, the dada2 and vsearch source URIs, a `local:` override for a file already on disk, `remote_path` (dada2 only) for a file already present on the remote taxonomy host, the ordered taxonomy `levels`, the `vsearch_format` parser (`pr2` or `generic`), and the taxonomy `corrections`. Adding and removing a database is supported.

Removing or renaming a database, or changing its `levels`, reports the studies that use it, including those that inherit it through the cascade.

Both formats of one database should come from the same reference release, because the consensus of the two classifiers compares their labels as text. The editor warns when the version numbers in the two URIs differ.

### Defining primer pairs (`config/primers.yml`)

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

Store all primer pairs in here and reference whichever combinations you need per project. A primer shared by two pairs is passed to cutadapt once; to pass it twice, add the same sequence under a second name.

Edit this on the Primers page under SYSTEM in the sidebar, or in the YAML. The page adds and removes primers and composes pairs from them, checking each sequence against the IUPAC codes as you type. A pair naming a primer that does not exist is rejected on save.

Pair names are referenced by `cutadapt.primer_pairs` in `pipeline.yml`. Removing or renaming a pair reports the studies, groups and runs that use it. Renaming a primer updates the pairs that use it.

### Configuring cutadapt (`cutadapt:` in `pipeline.yml`)

Selects which primer pairs to apply and controls trimming behaviour.

```yaml
cutadapt:
  # Names must match keys in the Pairs section of config/primers.yml.
  primer_pairs:
    - PrimerPair1
    - PrimerPair2
  min_length: 200           # discard reads shorter than this after trimming (-m)
  discard_untrimmed: true   # drop reads where no adapter was found (--discard-untrimmed)
  cores: 0                  # parallel cores; 0 = auto-detect (-j)
  quality_cutoff: ~         # 3' quality trimming cutoff, null to disable (-q)
  error_rate: ~             # max adapter mismatch rate, null = cutadapt default (-e)
  overlap: ~                # min adapter overlap length, null = cutadapt default (-O)
  optional_args: ""         # additional flags passed verbatim to cutadapt
```

### Threads and remote execution (top level of `pipeline.yml`)

Both settings sit outside the stages. Changing either does not mark completed
work stale.

```yaml
# Worker threads for the DADA2 stages learn_errors, denoise, chimera_removal
# and assign_taxonomy: true = every core, false = one, or a positive integer.
r_threads: 4

# DISCLAIMER: you are solely responsible for ensuring you are authorised to use
# the host configured here, and to place sequencing data on it.
remote:
  host: ~                        # user@hostname; ~ keeps every stage local
  rscript: "Rscript"             # Rscript on the server
  staging_dir: "/absolute/path/on/server"
  identity_file: ~               # SSH key; ~ uses the agent, then a password prompt
  threads: ~                     # threads on the server; ~ uses r_threads (DADA2) or phylogeny.threads
  stages: []                     # any of the stages below
  tools:                         # programs the phylogeny stages call on the server
    mafft: "mafft"
    iqtree: "iqtree"
    raxml: "raxmlHPC-PTHREADS-SSE3"
```

Each stage listed in `remote.stages` uploads its own inputs and downloads its
own outputs, so any combination of stages can run on the server:

| Stage | Uploads | Returns |
| --- | --- | --- |
| `learn_errors` | the reads within the `dada2.dada.nbases` budget | `ckpt_errors.RData` and the error-rate plot |
| `denoise` | filtered reads and `ckpt_errors.RData` | `ckpt_denoise.RData` |
| `chimera_removal` | two checkpoints | a checkpoint |
| `assign_taxonomy` | one checkpoint, and the database unless `dada2.remote_path` in `config/databases.yml` names a copy on the server | a checkpoint |
| `phylogeny_align` | the reference FASTA | the reference alignment (MAFFT) |
| `phylogeny_tree` | the trimmed reference alignment | the reference tree (IQ-TREE) |
| `phylogeny_add` | the queries and the trimmed reference alignment | the combined alignment (MAFFT `--addfragments`) |
| `phylogeny_place` | the trimmed combined alignment and the reference tree | the placement jplace (RAxML EPA) |

trimAl and gappa run on this machine. Reference trees belong to no study and
read `remote` from `config/pipeline.yml`; placements read it through their
study's cascade.

`filter_trim` runs on this machine with DADA2's own threading, so `r_threads`
does not apply to it.

Each stage invocation stages its files in its own directory on the server and
removes it when the stage ends. Directories left by an interrupted stage are
removed on a later connection once they are a week old.

**Migrating from `dada2.taxonomy.multithread` and `dada2.taxonomy.remote`:** both
still work for `assign_taxonomy` and are deprecated. After moving them to the
top level, run this once to keep completed taxonomy work current:

```
julia --project=. -e 'include("scripts/migrate_r_threads.jl"); MigrateRThreads.run_migration("projects")'
```

### Configuring DADA2 (`dada2:` in `pipeline.yml`)

```yaml
dada2:
  file_patterns:
    mode: "paired"               # paired | forward | reverse

  # Filter and trim; DADA2's filterAndTrim():
  filter_trim:
    trunc_q: 2
    trunc_len: [220, 220]        # [forward, reverse]; first value used for single-end mode
    max_ee: [3, 3]               # maximum expected errors in F and R reads
    min_len: 175
    max_n: 0
    match_ids: true
    rm_phix: true

  # Denoising; learnErrors() and dada():
  dada:
    seed: 123
    nbases: 200000000
    max_consist: 15
    pool_method: "pseudo"        # none | pseudo | true

  # Merging; mergePairs(), paired mode only:
  merge:
    min_overlap: 20
    max_mismatch: 0
    trim_overhang: true

  # ASV length filtering and chimera removal:
  asv:
    band_size_min: 200           # null to skip length filtering
    band_size_max: 430
    denovo_method: "consensus"   # consensus | pooled | per-sample

  # Taxonomy; assignTaxonomy() against the configured database:
  taxonomy:
    database: pr2                # key into config/databases.yml
    min_boot: 0                  # minimum bootstrap confidence to retain (0-100)
    # Rank names come from `levels:` in config/databases.yml.

  # Output filename prefixes (all written to dada2/Tables/):
  output:
    seq_table_prefix: "seqtab_nochim"
    fasta_prefix: "asvs"
    taxa_prefix: "taxonomy"
```

**Outputs written to `projects/{name}/{run}/dada2/Tables/`:**

| File                      | Contents                                               |
| ---------------------------| --------------------------------------------------------|
| `seqtab_nochim.csv`       | Chimera-free ASV count table (samples x ASVs)          |
| `asvs.fasta` / `asvs.csv` | ASV sequences with short identifiers (seq1, seq2, ...) |
| `taxonomy.csv`            | Taxonomy assignments per ASV                           |
| `taxonomy_bootstraps.csv` | Bootstrap confidence values per rank                   |
| `taxonomy_combined.csv`   | Taxonomy + bootstrap columns combined                  |
| `tax_counts.csv`          | Taxonomy + per-sample counts                           |
| `asv_counts.csv`          | ASV sequences + per-sample counts (no taxonomy)        |
| `pipeline_stats.csv`      | Read counts retained at each pipeline stage            |

### Configuring vsearch (`vsearch:` in `pipeline.yml`)

Controls the alignment thresholds used when assigning taxonomy against the reference database. The same settings apply to ASVs and OTUs.

```yaml
vsearch:
  identity: 0.75        # minimum sequence identity (--id)
  query_cov: 0.8        # minimum fraction of query covered (--query_cov)
  maxaccepts: ~         # stop after this many hits per query, null = vsearch default
  maxrejects: ~         # max rejected candidates, null = vsearch default
  strand: ~             # "plus" or "both"; null = vsearch default
  optional_args: ""     # additional flags passed verbatim to vsearch
```

### Configuring cd-hit-est (`cdhit:` in `pipeline.yml`)

Optional clustering step that collapses near-identical ASVs before vsearch taxonomy assignment. Use it with multiplexed primers, where one template amplified by different primers yields near-identical ASVs.

```yaml
cdhit:
  identity: 1           # sequence identity threshold (-c)
  threads: 0            # worker threads; 0 = all available (-T)
  optional_args: ""     # additional flags passed verbatim to cd-hit-est
```

### Configuring swarm (`swarm:` in `pipeline.yml`)

OTU clustering pipeline run in parallel with DADA2. Produces an OTU count table and FASTA which are carried through vsearch taxonomy assignment and `merge_taxa` alongside the ASV outputs.

```yaml
swarm:
  differences: 1          # -d: max differences between sequences in the same cluster
  threads: 0              # -t: worker threads; 0 = all available
  chimera_check: true     # run vsearch --uchime_denovo before clustering
  min_abundance: 2        # --minsize: discard singleton dereps before clustering
  fastq_minovlen: 20      # min overlap for paired-end merging
  identity: 0.97          # --id: threshold for mapping reads back to OTU seeds
  optional_args: ""       # additional flags passed verbatim to swarm
```

### Configuring phylogenetic placement (`phylogeny:` in `pipeline.yml`)

Controls building reference trees (`reference`) and placing a study's sequences on them (`placement`). A reference tree or placement can override its own keys from its page.

```yaml
phylogeny:
  threads: 4                  # MAFFT, IQ-TREE, RAxML and gappa on this machine
  reference:
    align:                    # MAFFT
      strategy: "localpair"   # auto | localpair | genafpair | globalpair | 6merpair
      maxiterate: 1000        # --maxiterate; 0 omits it
    trim:                     # trimAl
      method: "manual"
      gap_threshold: 0.3
    tree:                     # IQ-TREE
      model: "MFP"            # -m; MFP selects a model with ModelFinder
      bootstrap: "standard"   # standard (-b) | ultrafast (-bb, 1000 replicates or more)
      replicates: 100
  placement:
    align:                    # MAFFT --addfragments onto the trimmed reference alignment
      strategy: "auto"
      maxiterate: 0
    trim:                     # trimAl
      method: "manual"
      gap_threshold: 0.01
    place:                    # RAxML EPA
      model: "GTRCATI"        # -m
      heuristic: 0.2          # -G; ~ tries every branch
    accumulate:               # gappa edit accumulate
      threshold: 0.8          # --threshold, 0.5 to 1
```

Each step also takes `optional_args`, additional flags passed verbatim to its tool.

### Configuring merge_taxa (`merge_taxa:` in `pipeline.yml`)

Controls which filter configs are applied when merging taxonomy and count tables. `merged.csv` (unfiltered) is always written; each entry in `filters` produces an additional filtered CSV.

```yaml
merge_taxa:
  filters:
    - "protist_filter.yml"   # -> merged/protist_filter.csv
```

Each entry names a filter in the `filters:` library of `config/composition.yml`. Remove all entries (or set `filters: []`) to produce only the unfiltered `merged.csv`.

### Configuring analysis (`analysis:` in `pipeline.yml`)

Controls the defaults applied to the analysis charts (alpha diversity, taxa bar, NMDS, differential abundance, etc.). Per-chart choices such as the taxonomic rank and relative or absolute abundance are made in the UI.

```yaml
analysis:
  exclude_categories:            # composition categories to drop from figures; [] to keep all
    - {set: contamination, category: Contaminant, apply_to: [diversity, taxa, venn, differential]}
                                 # apply_to surfaces: diversity | taxa | composition | venn | differential
                                 # (omit apply_to to act on every surface)
  alpha:                         # richness, Shannon and Simpson
    normalisation: none          # none | rarefy | srs (scaling with ranked subsampling)
    depth: 0                     # rarefy/srs reads per sample; 0 = smallest non-empty sample
    rarefaction_iterations: 100  # rarefy draws averaged per value
    show_points: true            # overlay individual sample points on boxplots
    annotate_significance: false # annotate pairwise significance on grouped alpha
    pairwise_brackets: false     # draw significance brackets between groups
    paired_samples: false        # treat samples as paired in the significance test
    paired_lines: false          # join each sample to itself across groups with a line
    significance_test: "kruskal_wallis"  # test used for group comparison
  beta:                          # Bray-Curtis for NMDS and PERMANOVA
    normalisation: hellinger     # hellinger | none | rarefy | srs
    depth: 0                     # rarefy/srs reads per sample; 0 = smallest non-empty sample
  nmds:
    max_stress: 0.2              # warn if NMDS stress exceeds this value
  differential:                  # negative-binomial differential abundance (MASS::glm.nb)
    offset: tss                  # tss | rle (rle is refused when no taxon has reads in every sample)
    min_prevalence: 0.0          # fraction of samples, 0 to 1, a taxon needs reads in to be tested
```

### Configuring taxonomic filtering (`filters:` in `config/composition.yml`)

Each named filter in the `filters:` library of `config/composition.yml` defines one biological group to extract from the merged table. A category set references these filters by name, and the same filters back the `merge_taxa.filters` stage, which produces one additional CSV per entry. Edit them on the Compositions page under SYSTEM in the sidebar, or in the YAML.

Saved table-view presets are a separate concern and live in `config/presets/`; the Tables view reads and writes them.

#### Database-specific filters

Each filter carries a `databases:` key so that it is only applied when the active database matches. The following filters ship in the library:

| Category | PR2 match |
|----------|-----------|
| `bacteria_archaea` | `Domain` = Bacteria\|Archaea |
| `environmental_protozoa` | `Subdivision` = Cercozoa\|Gyrista\|Ciliophora\|Chrompodellids |
| `fungi` | `Subdivision` = Fungi |
| `helminths` | `Class` = Nematoda (excl. *Miculenchus*) |
| `parasitic_protozoa` | `Subdivision` = Apicomplexa\|Parabasalia\|Fornicata\|Bigyra |
| `plants_invertebrates` | Exclusion-based (PR2 ranks) |
| `protist` | Exclusion-based (PR2 ranks) |
| `vertebrates` | `Class` = Craniata |

Example:

```yaml
# fungi.pr2.yml
databases: [pr2]

filters:
  - column: Subdivision
    pattern: Fungi
    action: keep        # keep rows matching the pattern (default action is exclude)

remove_empty:
  - Subdivision
```

#### Filter file format

```yaml
databases: [pr2]          # omit to apply regardless of active database

mappings:                 # optional column remapping applied before filters
  - source_column: Division
    target_column: Supergroup
    values: { Rhizaria: Rhizaria, Alveolata: Alveolata }

filters:
  - column: Domain
    pattern: "Bacteria|Archaea"
    regex: true           # false (default) = substring match
    action: exclude       # exclude (default) | keep

remove_empty:             # remove rows where this column is blank or "NA"
  - Domain
```

## Deployment 

### Local (single machine)

```bash
bash start.sh
```

Open `http://localhost:8080`. The backend serves the frontend.

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

## REST API

The server exposes a REST API under `/api/v1/`. Key endpoint groups:

| Group          | Endpoints                                                                                                     | Description                                                                              |
| ----------------| ---------------------------------------------------------------------------------------------------------------| ------------------------------------------------------------------------------------------|
| Studies        | `GET/POST/DELETE /studies`, `POST .../rename`                                                                 | List, create, rename, delete studies                                                     |
| Groups         | `POST/DELETE /studies/{study}/groups`, `POST .../rename`                                                      | Create, rename, delete groups                                                            |
| Runs           | `GET/POST/DELETE /studies/{study}/runs`, `POST .../rename`                                                    | List, create, rename, delete runs                                                        |
| Config         | `GET/PATCH/DELETE .../config`, `GET .../config/overrides`                                                     | Read and edit config at any cascade level; list downstream overrides                     |
| Primers        | `GET /primers`, `GET /primers/document`, `PUT /primers`                                                       | List pair names; read and replace the whole primers document (validated before writing)  |
| Pipeline       | `POST .../pipeline`, `POST .../stages/{stage}`                                                                | Launch full-study, single-run, or individual-stage jobs                                  |
| Jobs           | `GET /jobs`, `GET/DELETE /jobs/{id}`                                                                          | Monitor and cancel running pipeline jobs                                                 |
| Events         | `GET /events`                                                                                                 | Server-sent event stream of real-time job and stage updates                              |
| Results        | `GET/POST/DELETE .../results/tables/...`                                                                      | List, query, filter, save, export (`.xlsx`), and delete tables; OTU member drill-down    |
| QC             | `GET .../results/qc`, `GET .../results/dada2`                                                                 | MultiQC report metadata and DADA2 figures, logs, stats                                   |
| Analysis       | `POST .../runs/{run}/analysis/{alpha,chart}`, `GET .../runs/{run}/analysis/{pipeline-stats,ranks}`            | Per-run charts and rank discovery                                                        |
| Cross-run      | `POST /studies/{study}/analysis/{alpha,chart,chart-facet,nmds,permanova,venn,publication-tables,differential}` | Comparison, faceted charts, NMDS, PERMANOVA, taxon overlap, publication tables and differential abundance across runs |
| Composition    | `POST .../runs/{run}/composition/summary`, `POST .../composition/{source}/query`, `POST .../composition/{source}/distinct/{column}` | Organism-composition summaries and tables for a run                             |
| Composition library | `GET /composition`, `POST/DELETE /composition/{filters,sets}/{name}`, `GET /category-sets`, `POST/DELETE /category-sets/{name}` | Edit the filters and category sets in `config/composition.yml`               |
| Reference trees | `GET/POST /reference-trees`, `GET/PUT/DELETE /reference-trees/{id}`, `GET/PUT .../fasta`, `POST .../run`, `GET .../log/{step}`, `GET .../qc/{align,trim}`, `GET .../alignment/{raw,trimmed}`, `GET .../treefile`, `POST .../trim-preview` | Build and inspect the shared reference trees |
| Placements     | `GET/POST /studies/{study}/placements`, `GET/PUT/DELETE .../placements/{id}`, `POST .../placements/preview`, `GET/PUT .../fasta/queries`, `POST .../run`, `GET .../log/{step}`, `GET .../qc/{align,trim}`, `GET .../alignment/{raw,trimmed}`, `POST .../trim-preview` | Place a study's sequences on a reference tree |
| Trees          | `GET/POST /studies/{study}/trees`, `GET/DELETE /studies/{study}/trees/{file}`, `PUT /studies/{study}/trees/{file}/view` | Upload, read and delete Newick and jplace files; save each tree's view state          |
| Filter presets | `GET/POST/DELETE /filter-presets`                                                                             | Save, list, delete reusable table filters                                                |
| Databases      | `GET /databases`, `GET /databases/document`, `PUT /databases`, `POST /databases/{key}/download`               | List and download taxonomy databases; read and replace the whole databases document (validated before writing, returns advisory warnings) |
| System         | `POST /init`, `GET /capabilities`                                                                             | Initialise project directories; report server capabilities (e.g. R availability)         |

All responses are JSON. Analysis endpoints return Plotly chart specifications.

## Running tests

```bash
julia --project=. test/runtests.jl                  # unit tests
julia --project=. test/runtests.jl --integration    # adds the MiSeq SOP pipeline run
julia --project=. test/runtests.jl --server         # adds the HTTP server smoke tests
```

Set `CI_SKIP_TAXONOMY=1` to skip taxonomy assignment in the integration run. The full suite loads R and every package and needs several GB of memory; on a small machine run a few files at a time. `cd frontend && bun run test` runs the frontend unit tests.

Name unit test files to run only those: `julia --project=. test/runtests.jl routes trees`. `dev/dev.sh` runs the development server, rebuilds and tests one job at a time, stopping and restarting the server around builds and tests (`dev/dev.sh` alone lists its commands). Deferred work is listed in `dev/TODO.md`.

## Architecture

```
frontend/           TypeScript + React + Vite (SPA)
src/
  core/             Types, config cascade, validation, DuckDB store, logging
  pipeline/         Pipeline stages (cutadapt, dada2, swarm, vsearch, cd-hit-est, merge_taxa, phylogeny)
  analysis/         Diversity metrics + Plotly chart builders
  server/           Oxygen.jl HTTP server (MetaManifold.Server) and its precompile workload
    routes/         REST API route handlers
scripts/            serve.jl (server entry point) and one-off maintenance scripts
config/             Default configs, filters, CI fixtures
data/               Input FASTQs (user-managed)
projects/           Pipeline outputs (generated)
```

Each pipeline stage returns a typed result (`TrimmedReads`, `ASVResult`, `OTUResult`, `TaxonomyHits`, `MergedTables`) and is skipped when its outputs are up to date (mtime-based for files, content-hash-based for configuration). Rerunning after a config change only re-executes the minimum necessary stages.

## Third-party tools

This project orchestrates the following tools. Each is fetched from its upstream source by `install.sh` and is subject to its own licence; no third-party binaries are included in this repository.

| Tool                                            | Licence | Source                  |
| -------------------------------------------------| ---------| -------------------------|
| [cutadapt](https://github.com/marcelm/cutadapt) | MIT     | PyPI                    |
| [FastQC](https://github.com/s-andrews/FastQC)   | GPL v3  | Babraham Bioinformatics |
| [MultiQC](https://github.com/MultiQC/MultiQC)   | GPL v3  | PyPI                    |
| [DADA2](https://benjjneb.github.io/dada2/)      | LGPL v3 | Bioconductor            |
| [swarm](https://github.com/frederic-mahe/swarm) | GPL v3  | GitHub Releases         |
| [vsearch](https://github.com/torognes/vsearch)  | GPL v3  | GitHub Releases         |
| [cd-hit](https://github.com/weizhongli/cdhit)   | GPL v2+ | GitHub Releases / apt   |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | BSD  | Source / apt            |
| [trimAl](https://github.com/inab/trimal)        | GPL v3  | GitHub Releases         |
| [IQ-TREE](https://github.com/iqtree/iqtree3)    | GPL v2  | GitHub Releases         |
| [RAxML](https://github.com/stamatak/standard-RAxML) | GPL v3 | Source / apt         |
| [gappa](https://github.com/lczech/gappa)        | GPL v3  | Source                  |

## Acknowledgements

This pipeline draws on the following prior work:

- **Frédéric Mahé**: [Fred's metabarcoding pipeline](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline) informed the overall workflow architecture, namely the sequencing of primer trimming, `swarm.jl`, vsearch-based taxonomy assignment, and the final table merge/filter stages.
- **Benjamin J. Callahan _et al._**: [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html), used under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), on which `dada2.jl` and its modules are based.

The following colleagues at the **Department of Parasitology, Charles University** (Faculty of Science, BIOCEV, Vestec, Czech Republic) contributed to this work:

- **Mgr. Jiří Novák** (supervisor): scripts from which several modules and configurations were adapted.
- **doc. Mgr. Vladimír Hampl**: provided laboratory access and resources.
- **Mgr. Paulína Pristašová**: <3.

## Licence

Copyright © 2026 Joshua Benjamin Jewell.

Source code is licensed under the [GNU Affero General Public License v3.0](LICENSE).

This documentation (README.md) is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).
