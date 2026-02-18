# DADA2 amplicon sequencing functions
#
# Functions for dada2.r main workflow
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
library(openxlsx)
library(tidyverse)
library(yaml)

# Null-coalescing operator: returns x if not NULL, otherwise y.
# Used to apply config defaults.
`%||%` <- function(x, y) if (!is.null(x)) x else y

## Configuration
load_config <- function(yaml_path) {
  if (!file.exists(yaml_path)) {
    stop("Config file not found: ", yaml_path)
  }
  yaml.load_file(yaml_path)
}

validate_config <- function(cfg) {
  required <- c("workspace", "file_patterns", "filter_trim",
                "dada", "merge", "asv", "taxonomy", "output")
  missing_sections <- setdiff(required, names(cfg))
  if (length(missing_sections) > 0) {
    stop("Missing required config sections: ",
         paste(missing_sections, collapse = ", "))
  }

  if (is.null(cfg$workspace$root)) {
    stop("workspace.root is required")
  }
  if (is.null(cfg$workspace$input_dir)) {
    stop("workspace.input_dir is required")
  }
  if (!dir.exists(cfg$workspace$input_dir)) {
    stop("workspace.input_dir does not exist: ", cfg$workspace$input_dir)
  }

  mode <- cfg$file_patterns$mode %||% "paired"
  if (!mode %in% c("paired", "forward", "reverse")) {
    stop("file_patterns.mode must be one of: paired, forward, reverse")
  }

  boot_mode <- cfg$output$bootstraps %||% "combined"
  if (!boot_mode %in% c("none", "combined", "separate")) {
    stop("output.bootstraps must be one of: none, combined, separate")
  }

  # taxonomy.uri is only required when taxonomy assignment is not skipped
  if (!(cfg$taxonomy$skip %||% FALSE) && is.null(cfg$taxonomy$uri)) {
    stop("taxonomy.uri is required when taxonomy.skip is not true")
  }

  combined_mode <- cfg$output$combined_mode %||% "regular"
  if (!combined_mode %in% c("regular", "alternative")) {
    stop("output.combined_mode must be one of: regular, alternative")
  }

  invisible(cfg)
}

# Creates the standard directory tree under root and returns paths as a named
# list so downstream functions can reference them by:
# directory <- paths$Directory (e.g. paths$Tables).
setup_workspace <- function(root) {
  dirs <- list(
    Tables   = file.path(root, "Tables"),
    Analysis = file.path(root, "Analysis"),
    Taxonomy = file.path(root, "Taxonomy"),
    Figures  = file.path(root, "Figures"),
    Filtered = file.path(root, "Filtered")
  )
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  dirs
}

# Finds FASTQ files matching forward/reverse patterns. Files are sorted so that
# each forward file always corresponds its reverse file by name.
find_fastq_files <- function(input_dir, fwd_pattern, rev_pattern, mode) {
  fwd <- NULL
  rev <- NULL
  if (mode %in% c("paired", "forward")) {
    fwd <- sort(list.files(input_dir, pattern = fwd_pattern, full.names = TRUE))
  }
  if (mode %in% c("paired", "reverse")) {
    rev <- sort(list.files(input_dir, pattern = rev_pattern, full.names = TRUE))
  }
  list(forward = fwd, reverse = rev)
}

# Checks that at least one sample was found and that forward/reverse counts
# match in paired mode. With one sample, dada() would return a bare object
# instead of a list so duplication is performed in main() by the file paths
# before the pipeline runs.
validate_sample_files <- function(fwd, rev, mode) {
  if (mode %in% c("paired", "forward")) {
    if (length(fwd) < 1) {
      stop("No forward FASTQ files found. ",
           "Check workspace.input_dir and file_patterns.forward.")
    }
    missing_fwd <- fwd[!file.exists(fwd)]
    if (length(missing_fwd) > 0) {
      stop("Forward files not found:\n  ", paste(missing_fwd, collapse = "\n  "))
    }
  }
  if (mode %in% c("paired", "reverse")) {
    if (length(rev) < 1) {
      stop("No reverse FASTQ files found. ",
           "Check workspace.input_dir and file_patterns.reverse.")
    }
    missing_rev <- rev[!file.exists(rev)]
    if (length(missing_rev) > 0) {
      stop("Reverse files not found:\n  ", paste(missing_rev, collapse = "\n  "))
    }
  }
  if (mode == "paired" && length(fwd) != length(rev)) {
    stop("Forward file count (", length(fwd), ") does not match ",
         "reverse file count (", length(rev), ").")
  }
  invisible(NULL)
}

# Extracts sample names. (e.g. "Sample1_R1_trimmed.fastq.gz" split by "_"
# at index 1 gives "Sample1").
extract_sample_names <- function(fwd_files, split_char, split_index) {
  sapply(strsplit(basename(fwd_files), split_char), `[`, split_index)
}

# Constructs output paths for filtered reads.
make_filtered_paths <- function(sample_names, filtered_dir, mode = "paired") {
  fwd <- NULL
  rev <- NULL
  if (mode != "reverse") {
    fwd <- file.path(filtered_dir, paste0(sample_names, "_R1_filt.fastq.gz"))
    names(fwd) <- sample_names
  }
  if (mode != "forward") {
    rev <- file.path(filtered_dir, paste0(sample_names, "_R2_filt.fastq.gz"))
    names(rev) <- sample_names
  }
  list(forward = fwd, reverse = rev)
}

## Plotting

# Saves per-base quality score profiles for up to 3 forward and 3 reverse
# samples. Inspect before and after filtering to guide truncLen / maxEE choices.
plot_quality_profiles <- function(fwd_files, rev_files, output_pdf) {
  pdf(output_pdf, width = 8, height = 6)
  if (!is.null(fwd_files) && length(fwd_files) > 0) {
    plotQualityProfile(fwd_files[seq_len(min(3L, length(fwd_files)))])
  }
  if (!is.null(rev_files) && length(rev_files) > 0) {
    plotQualityProfile(rev_files[seq_len(min(3L, length(rev_files)))])
  }
  dev.off()
  invisible(NULL)
}

# Plots the learned substitution error rates against quality scores.
# The fitted line should track the observed points closely; poor fit suggests
# the error model did not converge. If so, try increasing nbases or max_consist.
plot_error_rates <- function(fwd_errors, rev_errors, output_pdf) {
  pdf(output_pdf, width = 8, height = 6)
  if (!is.null(fwd_errors)) plotErrors(fwd_errors, nominalQ = TRUE)
  if (!is.null(rev_errors)) plotErrors(rev_errors, nominalQ = TRUE)
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

# Wrapper around filterAndTrim that handles both paired and single-end modes.
# Key parameters:
#   truncQ   - truncate reads at the first base with Phred quality ≤ this value
#   maxEE    - maximum expected errors per read (calculated from Phred scores);
#              lower = stricter, recommended 2–5 for typical Illumina data
#   matchIDs - (paired only) discard read pairs where one read was filtered out,
#              ensuring F and R files remain perfectly synchronised
run_filter_trim <- function(fwd_in, rev_in, fwd_out, rev_out, params, verbose) {
  trunc_len <- unlist(params$trunc_len)
  max_ee    <- unlist(params$max_ee)

  if (!is.null(rev_in)) {
    filterAndTrim(
      fwd_in, fwd_out,
      rev_in, rev_out,
      truncQ   = params$trunc_q,
      truncLen = trunc_len,
      maxEE    = max_ee,
      minLen   = params$min_len,
      maxN     = params$max_n,
      matchIDs = params$match_ids,
      rm.phix  = params$rm_phix,
      verbose  = verbose
    )
  } else {
    # Single-end: use only the first element of trunc_len / maxEE
    filterAndTrim(
      fwd_in, fwd_out,
      truncQ   = params$trunc_q,
      truncLen = trunc_len[1],
      maxEE    = max_ee[1],
      minLen   = params$min_len,
      maxN     = params$max_n,
      rm.phix  = params$rm_phix,
      verbose  = verbose
    )
  }
}

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
# outputBootstraps is always TRUE so bootstrap confidence values (0–100 per
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

# Writes taxonomy assignments to CSV. Bootstrap confidence values (0–100 per
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
  openxlsx::write.xlsx(combined, output_path, overwrite = TRUE,
                       asTable = FALSE, sheetName = "1.sampling",
                       firstRow = TRUE, zoom = 90, keepNA = TRUE)
  invisible(NULL)
}
