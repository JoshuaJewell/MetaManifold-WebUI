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
library(dplyr)
library(tibble)
# library(yaml)

# BiocParallel defaults to MulticoreParam, which forks. ShortRead::qa(), reached
# from plotQualityProfile(), dispatches through bplapply() and so forks children
# out of the embedded R session. That R session lives inside a multithreaded
# Julia process, and fork() carries over only the calling thread: any lock a
# Julia thread happened to hold is inherited already locked, with no owner left
# to release it. The children wedge, and the parent blocks forever in select()
# waiting to collect them, holding the R lock and the whole pipeline with it.
# SerialParam does the same work in-process, where there is nothing to deadlock.
BiocParallel::register(BiocParallel::SerialParam())

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
#
# `kept` indexes the samples that still had reads after filtering; the dada,
# merge and nochim columns cover only those, and the others are recorded as 0.
compute_pipeline_stats <- function(filter_stats, dada_fwd, dada_rev, merged,
                                   seq_table_nochim, sample_names, mode,
                                   kept = seq_along(sample_names)) {
  get_n <- function(x) sum(getUniques(x))
  # dada() and mergePairs() return a bare object for a single input.
  as_list <- function(x) if (inherits(x, "dada") || is.data.frame(x)) list(x) else x
  spread <- function(values) {
    out <- numeric(length(sample_names))
    out[kept] <- values
    out
  }

  track <- as.data.frame(filter_stats)
  colnames(track) <- c("input", "filtered")

  if (mode %in% c("paired", "forward") && !is.null(dada_fwd)) {
    track$denoisedF <- spread(sapply(as_list(dada_fwd), get_n))
  }
  if (mode %in% c("paired", "reverse") && !is.null(dada_rev)) {
    track$denoisedR <- spread(sapply(as_list(dada_rev), get_n))
  }
  if (mode == "paired" && !is.null(merged)) {
    track$merged <- spread(sapply(as_list(merged), get_n))
  }
  track$nochim <- if (nrow(seq_table_nochim) == length(kept)) spread(rowSums(seq_table_nochim)) else 0

  rownames(track) <- sample_names
  track
}

## Taxonomy

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
# taxonomic rank)
# Returns taxa_df (SeqName + Sequence + taxonomy), used by write_combined_table().
write_taxa_table <- function(tax_matrix, boot_matrix, index, tables_dir, prefix) {
  taxa_df <- as.data.frame(tax_matrix, stringsAsFactors = FALSE)
  taxa_df$Sequence <- rownames(taxa_df)
  taxa_df <- dplyr::left_join(index, taxa_df, by = "Sequence")

  write.csv(taxa_df, file.path(tables_dir, paste0(prefix, ".csv")),
            quote = FALSE, row.names = FALSE)

  boot_df <- as.data.frame(boot_matrix, stringsAsFactors = FALSE)
  colnames(boot_df) <- paste0(colnames(boot_df), "_boot")
  boot_df$Sequence <- rownames(boot_matrix)
  boot_df <- dplyr::left_join(index, boot_df, by = "Sequence")
  write.csv(boot_df, file.path(tables_dir, paste0(prefix, "_bootstraps.csv")),
            quote = FALSE, row.names = FALSE)

  combined_df <- dplyr::left_join(taxa_df, boot_df, by = c("SeqName", "Sequence"))
  write.csv(combined_df, file.path(tables_dir, paste0(prefix, "_combined.csv")),
            quote = FALSE, row.names = FALSE)

  taxa_df
}

# Joins taxonomy with per-sample counts
write_combined_table <- function(taxa_df, index, seq_table, tables_dir, tax_filename, asv_filename) {
  seq_t <- as.data.frame(t(seq_table), stringsAsFactors = FALSE)
  seq_t <- tibble::rownames_to_column(seq_t, "Sequence")
  write.csv(dplyr::left_join(taxa_df, seq_t, by = "Sequence"),
            file.path(tables_dir, tax_filename), quote = FALSE, row.names = FALSE)
  write.csv(dplyr::left_join(index, seq_t, by = "Sequence"),
            file.path(tables_dir, asv_filename), quote = FALSE, row.names = FALSE)
  invisible(NULL)
}

## Remote stage helpers
# Shared by the *_remote.r scripts, which Julia uploads alongside this file and
# runs through Rscript on the bioserver. They are defined here rather than
# repeated in each script so the four remote stages parse their arguments and
# coerce their thread count identically.

# Parses `key=value` command-line arguments into a named list, splitting on the
# FIRST `=` only so a value may itself contain one. Stops when a required key is
# absent, which is a great deal easier to read than the NULL dereference that
# would otherwise surface several lines later.
parse_remote_args <- function(args, required = character(0)) {
  p <- list()
  for (a in args) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) >= 2) p[[kv[1]]] <- paste(kv[-1], collapse = "=")
  }
  missing_args <- setdiff(required, names(p))
  if (length(missing_args) > 0)
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
  p
}

# DADA2's `multithread` is a union type: TRUE/FALSE select "every core" and
# "one core", an integer pins the count. The value arrives over the command line
# as text, so restore the distinction rather than passing a string, which DADA2
# would reject as invalid and silently downgrade to a single thread.
remote_multithread <- function(x) {
  if (tolower(x) %in% c("true", "false")) tolower(x) == "true" else as.integer(x)
}

# Reads a newline-delimited manifest written by the Julia side. Sample names and
# read paths travel this way instead of being interpolated into the ssh command
# string, which a remote login shell would word-split and glob.
read_manifest <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(character(0))
  lines <- readLines(path, warn = FALSE)
  lines[nzchar(lines)]
}
