# DADA2 amplicon sequencing functions
#
# R wrappers and functions for dada2.jl main workflow
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

library(dada2)
library(tidyverse)
# library(yaml)

# Config loading/validation, workspace setup, and file discovery are now
# handled by dada2.jl, passing resolved paths and values directly to R.

## Plotting

# Saves per-base quality score profiles for up to 3 forward and 3 reverse
# samples. Inspect before and after filtering to guide truncLen / maxEE choices.
plot_quality_profiles <- function(fwd_files, rev_files, output_pdf) {
  pdf(output_pdf, width = 8, height = 6)
  if (!is.null(fwd_files) && length(fwd_files) > 0) {
    print(plotQualityProfile(fwd_files[seq_len(min(3L, length(fwd_files)))]))
  }
  if (!is.null(rev_files) && length(rev_files) > 0) {
    print(plotQualityProfile(rev_files[seq_len(min(3L, length(rev_files)))]))
  }
  dev.off()
  invisible(NULL)
}

# Plots the learned substitution error rates against quality scores.
# The fitted line should track the observed points closely; poor fit suggests
# the error model did not converge. If so, try increasing nbases or max_consist.
plot_error_rates <- function(fwd_errors, rev_errors, output_pdf) {
  pdf(output_pdf, width = 8, height = 6)
  if (!is.null(fwd_errors)) print(plotErrors(fwd_errors, nominalQ = TRUE))
  if (!is.null(rev_errors)) print(plotErrors(rev_errors, nominalQ = TRUE))
  dev.off()
  invisible(NULL)
}

# Plots the distribution of merged ASV lengths. Inspect to confirm the
# dominant peak corresponds to the expected amplicon size and to set
# band_size_min / band_size_max for off-target removal.
plot_length_distribution <- function(seq_table, output_pdf) {
  pdf(output_pdf, width = 8, height = 6)
  plot(table(nchar(getSequences(seq_table))),
       xlab = "Merged read length (bp)",
       ylab = "Count",
       main = "ASV length distribution")
  dev.off()
  invisible(NULL)
}

## Filter and trim

# filterAndTrim() is now called directly from dada2.jl, which
# handles modes and parameters without this intermediate wrapper.
# I might delete it when I feel more destructive...

# run_filter_trim <- function(fwd_in, rev_in, fwd_out, rev_out, params, verbose) {
#   trunc_len <- unlist(params$trunc_len)
#   max_ee    <- unlist(params$max_ee)
#
#   if (!is.null(rev_in)) {
#     filterAndTrim(
#       fwd_in, fwd_out,
#       rev_in, rev_out,
#       truncQ   = params$trunc_q,
#       truncLen = trunc_len,
#       maxEE    = max_ee,
#       minLen   = params$min_len,
#       maxN     = params$max_n,
#       matchIDs = params$match_ids,
#       rm.phix  = params$rm_phix,
#       verbose  = verbose
#     )
#   } else {
#     # Single-end: use only the first element of trunc_len / maxEE
#     filterAndTrim(
#       fwd_in, fwd_out,
#       truncQ   = params$trunc_q,
#       truncLen = trunc_len[1],
#       maxEE    = max_ee[1],
#       minLen   = params$min_len,
#       maxN     = params$max_n,
#       rm.phix  = params$rm_phix,
#       verbose  = verbose
#     )
#   }
# }

## Sequence table

# Retains only ASVs whose length falls within [band_min, band_max].
# Amplicons outside this range are typically non-specific or artefactual
# (e.g. primer dimers, chimeras that survived removal, off-target loci).
filter_by_length <- function(seq_table, band_min, band_max) {
  lengths <- nchar(colnames(seq_table))
  seq_table[, lengths >= band_min & lengths <= band_max, drop = FALSE]
}

## Pipeline stats

# Builds a per-sample read-count table showing reads retained at each step:
# input -> filtered -> denoised (F/R) -> merged -> nochim.
# Large drops at any stage can indicate a problem:
#   filtered  - overly strict truncLen / maxEE
#   denoised  - poor error model fit
#   merged    - insufficient overlap, mismatched truncation lengths
#   nochim    - high chimera rate
compute_pipeline_stats <- function(filter_stats, dada_fwd, dada_rev, merged,
                                   seq_table_nochim, sample_names, mode) {
  get_n <- function(x) sum(getUniques(x))

  track <- as.data.frame(filter_stats)
  colnames(track) <- c("input", "filtered")

  if (mode %in% c("paired", "forward") && !is.null(dada_fwd)) {
    track$denoisedF <- sapply(dada_fwd, get_n)
  }
  if (mode %in% c("paired", "reverse") && !is.null(dada_rev)) {
    track$denoisedR <- sapply(dada_rev, get_n)
  }
  if (mode == "paired" && !is.null(merged)) {
    track$merged <- sapply(merged, get_n)
  }
  track$nochim <- rowSums(seq_table_nochim)

  rownames(track) <- sample_names
  track
}

## Taxonomy

# Downloads the reference database if not already present.
# The file is stored in the project's Taxonomy/ folder so it can be reused
# across runs without re-downloading.
fetch_taxonomy_db <- function(uri, local_dir) {
  local_path <- file.path(local_dir, basename(uri))
  if (!file.exists(local_path)) {
    message("  Downloading taxonomy database: ", uri)
    download.file(uri, local_path, mode = "wb")
  }
  local_path
}

# Assigns taxonomy to ASVs using a naive Bayesian classifier.
# outputBootstraps is always TRUE so bootstrap confidence values (0-100 per
# rank) are available for downstream filtering regardless of minBoot.
# minBoot sets the threshold at which assignments are returned; 0 returns all
# assignments.
run_assign_taxonomy <- function(seq_table, db_path, params, verbose) {
  result <- assignTaxonomy(
    seq_table,
    db_path,
    multithread      = params$multithread,
    minBoot          = params$min_boot,
    outputBootstraps = TRUE,
    verbose          = verbose,
    taxLevels        = params$levels
  )
  list(tax = result$tax, boot = result$boot)
}

## Output

# Writes the chimera-free ASV count table as a CSV with sequences as row names
# and samples as columns.
write_seq_table <- function(seq_table, tables_dir, prefix) {
  output_path <- file.path(tables_dir, paste0(prefix, ".csv"))
  write.csv(t(seq_table), output_path, quote = FALSE)
  invisible(NULL)
}

# Writes ASV sequences to a FASTA file and a companion CSV mapping short
# identifiers (seq1, seq2, ...) to full sequences. Returns the index
# data frame, which is used as the join key in later tables.
write_fasta <- function(seq_table, tables_dir, prefix) {
  sequences <- colnames(seq_table)
  n         <- length(sequences)
  seq_names <- sprintf("seq%d", seq_len(n))

  fasta_lines <- c(rbind(paste0(">", seq_names), sequences))
  writeLines(fasta_lines, file.path(tables_dir, paste0(prefix, ".fasta")))

  index <- data.frame(SeqName = seq_names, Sequence = sequences,
                      stringsAsFactors = FALSE)
  write.csv(index, file.path(tables_dir, paste0(prefix, ".csv")),
            quote = FALSE, row.names = FALSE)

  index
}

# Writes taxonomy assignments to CSV. Bootstrap confidence values (0-100 per
# taxonomic rank) are handled based on bootstrap_mode:
#   none     - taxonomy columns only
#   combined - taxonomy and bootstrap columns in one file (suffix: _boot)
#   separate - taxonomy and bootstraps written to two separate files
# Returns taxa_df (SeqName + Sequence + taxonomy), used by write_combined_table().
write_taxa_table <- function(tax_matrix, boot_matrix, index, tables_dir,
                             prefix, bootstrap_mode) {
  taxa_df <- as.data.frame(tax_matrix, stringsAsFactors = FALSE)
  taxa_df$Sequence <- rownames(taxa_df)
  taxa_df <- dplyr::left_join(index, taxa_df, by = "Sequence")

  if (bootstrap_mode == "combined") {
    boot_df <- as.data.frame(boot_matrix, stringsAsFactors = FALSE)
    colnames(boot_df) <- paste0(colnames(boot_df), "_boot")
    boot_df$Sequence <- rownames(boot_matrix)
    combined_df <- dplyr::left_join(taxa_df, boot_df, by = "Sequence")
    write.csv(combined_df, file.path(tables_dir, paste0(prefix, ".csv")),
              quote = FALSE, row.names = FALSE)
  } else {
    write.csv(taxa_df, file.path(tables_dir, paste0(prefix, ".csv")),
              quote = FALSE, row.names = FALSE)
    if (bootstrap_mode == "separate") {
      boot_df <- as.data.frame(boot_matrix, stringsAsFactors = FALSE)
      boot_df$Sequence <- rownames(boot_matrix)
      boot_df <- dplyr::left_join(index, boot_df, by = "Sequence")
      write.csv(boot_df,
                file.path(tables_dir, paste0(prefix, "_bootstraps.csv")),
                quote = FALSE, row.names = FALSE)
    }
  }

  taxa_df
}

# Joins taxonomy (or index for alternative mode) with per-sample counts and
# writes an Excel workbook.
write_combined_table <- function(taxa_df, seq_table, tables_dir, filename) {
  seq_t <- as.data.frame(t(seq_table), stringsAsFactors = FALSE)
  seq_t <- tibble::rownames_to_column(seq_t, "Sequence")
  combined <- dplyr::left_join(taxa_df, seq_t, by = "Sequence")
  output_path <- file.path(tables_dir, filename)
  write.csv(combined, output_path, quote = FALSE, row.names = FALSE)
  invisible(NULL)
}
