# DADA2 amplicon sequencing pipeline
#
# Processes paired- or single-end Illumina reads into ASVs and assigns
# taxonomy using a reference database (default: PR2). Options for output
# of bootstraps or counts only, as well as parameters for filterAndTrim(),
# dada(), mergePairs(), and assignTaxonomy().
#
# Usage:
#   Rscript dada2.r <config.yaml>
#
# Outputs (in workspace.root/Tables/):
#   seqtab_nochim.csv        - chimera-free ASV count table (sequences as rows)
#   asvs.fasta / asvs.csv    - ASV sequences with short identifiers (seq1, seq2 ...)
#   taxonomy.csv             - taxonomy assignments (optionally with bootstrap confidence values)
#   tax_counts.xlsx          - combined taxonomy + per-sample counts
#   pipeline_stats.csv       - read counts at each pipeline stage
#
# Notice:
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).
#
# This work is based on the DADA2 tutorial by Benjamin J. Callahan, et al.,
# available at https://benjjneb.github.io/dada2/tutorial.html, with modification
# into a single module. The original material is licensed under the Creative
# Commons Attribution 4.0 International License (CC BY 4.0):
# https://creativecommons.org/licenses/by/4.0/.

source("./src/dada2_functions.r")

## Pipeline
main <- function() {
  # Load args
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    stop("Usage: Rscript dada2_newer.r <config.yaml>")
  }
  yaml_path <- args[1]

  # Load config
  message("Loading and validating config: ", yaml_path)
  cfg     <- load_config(yaml_path)
  cfg     <- validate_config(cfg)
  verbose <- cfg$verbose %||% TRUE
  mode    <- cfg$file_patterns$mode %||% "paired"

  message("Setting up workspace: ", cfg$workspace$root)
  paths <- setup_workspace(cfg$workspace$root)

  # Find input files
  message("Discovering FASTQ files (mode: ", mode, ")")
  fastq <- find_fastq_files(
    cfg$workspace$input_dir,
    cfg$file_patterns$forward,
    cfg$file_patterns$reverse,
    mode
  )
  validate_sample_files(fastq$forward, fastq$reverse, mode)

  # In forward or reverse-only mode, sample names come from whichever side
  # is present; the split pattern in config determines what counts as a name.
  primary_files <- fastq$forward %||% fastq$reverse
  sample_names  <- extract_sample_names(
    primary_files,
    cfg$file_patterns$sample_name_split,
    cfg$file_patterns$sample_name_index
  )
  message("  Found ", length(sample_names), " samples")

  # Single-sample fallback:
  # dada() returns a bare dada object (not a list) when given a single file.
  # This breaks makeSequenceTable(), mergePairs(), and sapply() downstream.
  # Workaround: duplicate the file paths so the pipeline sees 2 samples, then
  # drop the duplicate row from all results before writing outputs.
  # ASV calls are unaffected because dada() processes each file independently.
  single_sample <- length(sample_names) == 1
  if (single_sample) {
    warning("Only 1 sample found. Duplicating it to work around dada() returning ",
            "a bare object for single-file input. The duplicate will be dropped ",
            "from all outputs.")
    fastq$forward <- rep(fastq$forward, 2)
    if (!is.null(fastq$reverse)) fastq$reverse <- rep(fastq$reverse, 2)
    sample_names <- c(sample_names, paste0(sample_names, "_dup"))
  }

  # Quality assessment (pre-filter)
  # Inspect quality_unfiltered.pdf to choose truncLen and maxEE values.
  message("Plotting unfiltered quality profiles")
  plot_quality_profiles(
    fastq$forward,
    fastq$reverse,
    file.path(paths$Figures, "quality_unfiltered.pdf")
  )

  # Filter and trim
  # Removes low-quality bases and reads. Reads that pass filtering are written
  # to Filtered/; empty samples (0 reads) are silently skipped by DADA2 later.
  message("Filtering and trimming reads")
  filtered <- make_filtered_paths(sample_names, paths$Filtered, mode)

  # For single-end modes, route the relevant files into the fwd slots and
  # pass NULL for rev so run_filter_trim() calls the correct filterAndTrim form.
  in_fwd      <- if (mode != "reverse") fastq$forward else fastq$reverse
  out_fwd     <- if (mode != "reverse") filtered$forward else filtered$reverse
  in_rev_arg  <- if (mode == "paired") fastq$reverse else NULL
  out_rev_arg <- if (mode == "paired") filtered$reverse else NULL

  filter_stats <- run_filter_trim(in_fwd, in_rev_arg, out_fwd, out_rev_arg,
                                  cfg$filter_trim, verbose)

  message("Plotting filtered quality profiles")
  plot_quality_profiles(
    filtered$forward,
    filtered$reverse,
    file.path(paths$Figures, "quality_filtered.pdf")
  )

  # Error model
  # set.seed() should be called before learnErrors()
  message("Learning error rates")
  set.seed(cfg$dada$seed %||% 123L)

  fwd_errors <- if (mode != "reverse") {
    learnErrors(filtered$forward,
                nbases = cfg$dada$nbases, MAX_CONSIST = cfg$dada$max_consist,
                verbose = verbose)
  } else NULL

  rev_errors <- if (mode != "forward") {
    learnErrors(filtered$reverse,
                nbases = cfg$dada$nbases, MAX_CONSIST = cfg$dada$max_consist,
                verbose = verbose)
  } else NULL

  # Inspect error_rates.pdf: the fitted line should follow the observed points.
  plot_error_rates(fwd_errors, rev_errors,
                   file.path(paths$Figures, "error_rates.pdf"))

  # Denoising
  # dada() applies the error model to resolve ASV's from sequencing errors.
  # pool = "pseudo" shares information across samples to improve detection of
  # rare variants.
  message("Denoising reads")
  dada_fwd <- if (!is.null(fwd_errors)) {
    dada(filtered$forward, err = fwd_errors, pool = cfg$dada$pool_method,
         verbose = verbose)
  } else NULL

  dada_rev <- if (!is.null(rev_errors)) {
    dada(filtered$reverse, err = rev_errors, pool = cfg$dada$pool_method,
         verbose = verbose)
  } else NULL

  # Merge and sequence table
  # Paired reads are joined at the overlapping region. Longer overlap and
  # lower maxMismatch = higher confidence merges but fewer total merges.
  message("Building sequence table")
  if (mode == "paired") {
    merged <- mergePairs(
      dada_fwd, filtered$forward,
      dada_rev, filtered$reverse,
      minOverlap   = cfg$merge$min_overlap,
      maxMismatch  = cfg$merge$max_mismatch,
      trimOverhang = cfg$merge$trim_overhang,
      verbose      = verbose
    )
    seq_table <- makeSequenceTable(merged)
  } else {
    merged    <- NULL
    seq_table <- makeSequenceTable(dada_fwd %||% dada_rev)
  }

  # Plot the full length distribution before any length filtering so that
  # off-target bands are visible when choosing band_size_min / band_size_max.
  plot_length_distribution(seq_table,
                           file.path(paths$Figures, "length_distribution.pdf"))

  # Length filtering: retains only ASVs matching the expected amplicon size.
  # Set band_size_min / band_size_max to null in config to skip this step.
  band_min <- cfg$asv$band_size_min
  band_max <- cfg$asv$band_size_max
  if (!is.null(band_min) && !is.null(band_max)) {
    message("  Filtering by length: ", band_min, "-", band_max, " bp")
    seq_table <- filter_by_length(seq_table, band_min, band_max)
    # Second plot confirms off-target lengths were successfully removed.
    plot_length_distribution(seq_table,
                             file.path(paths$Figures, "length_distribution_filtered.pdf"))
  }

  # Chimera removal
  # "consensus" removes sequences identified as chimeric in the majority
  # of samples in which they appear.
  message("Removing chimeras")
  seq_table_nochim <- removeBimeraDenovo(seq_table, method = cfg$asv$denovo_method,
                                         verbose = verbose)

  nochim_pct <- sum(seq_table_nochim) / sum(seq_table) * 100
  message("  Chimeric reads removed: ", round(100 - nochim_pct, 2),
          "% | Retained: ", round(nochim_pct, 2), "%")

  # Drop duplicate row (if using single-sample fallback)
  if (single_sample) {
    filter_stats     <- filter_stats[1, , drop = FALSE]
    if (!is.null(dada_fwd)) dada_fwd <- dada_fwd[1]
    if (!is.null(dada_rev)) dada_rev <- dada_rev[1]
    if (!is.null(merged))   merged   <- merged[1]
    seq_table_nochim <- seq_table_nochim[1, , drop = FALSE]
    sample_names     <- sample_names[1]
  }

  # Pipeline stats
  message("Computing pipeline stats")
  stats <- compute_pipeline_stats(filter_stats, dada_fwd, dada_rev, merged,
                                  seq_table_nochim, sample_names, mode)
  write.csv(stats, file.path(paths$Tables, "pipeline_stats.csv"), quote = FALSE)
  if (verbose) print(stats)

  # Write core outputs
  write_seq_table(seq_table_nochim, paths$Tables,
                  cfg$output$seq_table_prefix %||% "seqtab_nochim")

  message("Writing output tables")
  # write_fasta() returns index (SeqName - Sequence), which is the join key
  # used in write_taxa_table() and write_combined_table() below.
  index <- write_fasta(seq_table_nochim, paths$Tables,
                       cfg$output$fasta_prefix %||% "asvs")

  # combined_input() defaults to index (SeqName + Sequence + counts only),
  # which is the alternative output mode. It is replaced with the full
  # taxa_df when taxonomy is run and combined_mode = "regular".
  combined_input <- index

  # Taxonomy
  if (!(cfg$taxonomy$skip %||% FALSE)) {
    message("Assigning taxonomy")
    db_path     <- fetch_taxonomy_db(cfg$taxonomy$uri, paths$Taxonomy)
    taxa_result <- run_assign_taxonomy(seq_table_nochim, db_path,
                                       cfg$taxonomy, verbose)
    taxa_df <- write_taxa_table(
      taxa_result$tax, taxa_result$boot, index,
      paths$Tables,
      cfg$output$taxa_prefix %||% "taxonomy",
      cfg$output$bootstraps %||% "combined"
    )
    if ((cfg$output$combined_mode %||% "regular") == "regular") {
      combined_input <- taxa_df
    }
  } else {
    message("  Skipping taxonomy (taxonomy.skip = true)")
  }

  # Checkpoint saved after taxonomy so the full environment (including taxa)
  # can be restored with load() without re-running the pipeline.
  checkpoint <- file.path(paths$Analysis, "checkpoint.RData")
  save.image(checkpoint)
  message("  Checkpoint saved: ", checkpoint)

  write_combined_table(combined_input, seq_table_nochim, paths$Tables,
                       cfg$output$combined_filename %||% "tax_counts.xlsx")

  message("Pipeline complete.")
  message("Outputs:")
  message("  ", file.path(paths$Tables, paste0(cfg$output$seq_table_prefix %||% "seqtab_nochim", ".csv")))
  message("  ", file.path(paths$Tables, paste0(cfg$output$fasta_prefix %||% "asvs", ".fasta")))
  message("  ", file.path(paths$Tables, paste0(cfg$output$fasta_prefix %||% "asvs", ".csv")))
  message("  ", file.path(paths$Tables, paste0(cfg$output$taxa_prefix %||% "taxonomy", ".csv")))
  message("  ", file.path(paths$Tables, cfg$output$combined_filename %||% "tax_counts.xlsx"))
  message("  ", file.path(paths$Tables, "pipeline_stats.csv"))
  message("  ", checkpoint)
}

# Run pipeline
if (sys.nframe() == 0L) main()