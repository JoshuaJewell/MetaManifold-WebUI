<!-- SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Pipeline and analysis

MetaManifold consists of a Julia backend (pipeline engine + REST API) and a TypeScript/React frontend. The pipeline runs FastQC, MultiQC, cutadapt, DADA2, swarm, vsearch, and cd-hit-est under the hood, and places ASVs on reference trees with MAFFT, trimAl, IQ-TREE, RAxML and gappa; results are stored in per-run DuckDB databases and served to the frontend as interactive Plotly charts and filterable tables. Pipeline configuration is editable in the web UI at every cascade level (see [Configuration](configuration.md)).

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

Counts may be normalised before analysis (none, rarefaction to a fixed or auto-resolved depth, or scaling with ranked subsampling (SRS), which reaches the depth by scaling counts and keeps more of the community structure), and contamination-flagged taxa may be included or excluded. Ordinations can also apply a Hellinger transform before the Bray-Curtis dissimilarity. Composition charts can be drawn as a grid faceted on any two of run, group and sub-group. All analysis charts are returned as Plotly JSON and rendered interactively in the browser.
