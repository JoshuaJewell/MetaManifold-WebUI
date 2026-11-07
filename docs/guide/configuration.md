<!-- SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Configuration

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

### Profiling stage memory (`occupancy_probe:` in `pipeline.yml`)

A diagnostic for sizing a machine, off by default and not offered in the config
editor. Add it to a run's or study's `pipeline.yml` by hand:

```yaml
occupancy_probe:
  enabled: true
  interval_ms: 10    # sampling interval in milliseconds (at least 1)
```

Each stage then records an `occupancy` entry in the run's `attestation.yml`:

- `child_peak_rss_kib`: the largest total resident memory of the stage's tool
  processes (vsearch, swarm, cd-hit and so on) at any one sample. Tools run
  concurrently can approach the sum of their individual peaks.
- under `tasks`, each task's `peak_rss_kib` (largest sampled resident memory)
  and `hwm_kib` (the kernel's high-water mark), with `samples: 0` for a process
  too short-lived to be sampled.
- `julia_allocated_bytes`, the Julia allocation over the stage. It is
  process-wide, so it includes other runs' pipelines when several run at once.

Linux only, since it reads `/proc`. For stages run on a remote server, the
figures describe the local ssh client, not the work. With the probe off, the
records are unchanged.

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

Controls the defaults applied to the analysis charts (alpha diversity, taxa bar, NMDS, etc.). Per-chart choices such as the taxonomic rank and relative or absolute abundance are made in the UI.

```yaml
analysis:
  exclude_categories:            # composition categories to drop from figures; [] to keep all
    - {set: contamination, category: Contaminant, apply_to: [diversity, taxa, venn]}
                                 # apply_to surfaces: diversity | taxa | composition | venn
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
