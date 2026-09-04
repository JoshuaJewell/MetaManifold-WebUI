#!/usr/bin/env Rscript
# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Standalone chimera removal and core table writing for remote execution.
# Called by chimera_removal() in chimera.jl via SSH, through _run_remote_stage().
# All paths refer to the remote filesystem.
#
# This is the cheapest stage to offload: it reads two checkpoints rather than
# the read set, so the transfer is small however many samples the run carries.
#
# Arguments (key=value):
#   functions      path to dada2_functions.r
#   staging        staging directory on this host
#   mode           paired|forward|reverse
#   denovo_method  consensus|pooled|per-sample
#   seq_prefix     stem for the ASV count table
#   fasta_prefix   stem for the ASV FASTA and index
#   multithread    TRUE|FALSE|<integer>
#   verbose        true|false

args <- commandArgs(trailingOnly = TRUE)

functions_path <- sub("^functions=", "", grep("^functions=", args, value = TRUE)[1])
if (is.na(functions_path)) stop("Missing required argument: functions")
source(functions_path)

p <- parse_remote_args(args, required = c("functions", "staging", "mode",
                                          "denovo_method",
                                          "seq_prefix", "fasta_prefix", "multithread"))

staging    <- p[["staging"]]
tables_dir <- file.path(staging, "Tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

mode          <- p[["mode"]]
verbose       <- tolower(p[["verbose"]]) == "true"
multithread   <- remote_multithread(p[["multithread"]])
seq_prefix    <- p[["seq_prefix"]]
fasta_prefix  <- p[["fasta_prefix"]]
# Sample names arrive as $staging/samples.txt, in row order; kept.txt indexes
# the samples that still had reads after filtering.
final_names   <- read_manifest(file.path(staging, "samples.txt"))
kept          <- as.integer(read_manifest(file.path(staging, "kept.txt")))

load(file.path(staging, "ckpt_filter.RData"))   # filter_stats
load(file.path(staging, "ckpt_length.RData"))   # dada_fwd, dada_rev, merged, seq_table

has_data <- isTRUE(sum(seq_table, na.rm = TRUE) > 0)
if (has_data) {
    seq_table_nochim <- removeBimeraDenovo(seq_table, method = p[["denovo_method"]],
                                           multithread = multithread, verbose = verbose)
    nochim_pct <- sum(seq_table_nochim) / sum(seq_table) * 100
    message("  Chimeric reads removed: ", round(100 - nochim_pct, 2),
            "% | Retained: ", round(nochim_pct, 2), "%")
} else {
    seq_table_nochim <- seq_table
    message("Chimera removal skipped - seq_table is empty")
}

stats <- compute_pipeline_stats(filter_stats, dada_fwd, dada_rev, merged,
                                seq_table_nochim, final_names, mode, kept)
write.csv(stats, file.path(tables_dir, "pipeline_stats.csv"), quote = FALSE)
if (verbose) print(stats)

if (has_data) {
    write_seq_table(seq_table_nochim, tables_dir, seq_prefix)
    index <- write_fasta(seq_table_nochim, tables_dir, fasta_prefix)
} else {
    # The local stage touches these so downstream stages find the files they
    # expect and fail on content rather than on absence; do the same here so a
    # remote run leaves the run directory in the identical state.
    index <- data.frame(SeqName = character(0), sequence = character(0))
    file.create(file.path(tables_dir, paste0(seq_prefix, ".csv")))
    file.create(file.path(tables_dir, paste0(fasta_prefix, ".fasta")))
    file.create(file.path(tables_dir, paste0(fasta_prefix, ".csv")))
}

save(seq_table_nochim, index, file = file.path(staging, "ckpt_chimera.RData"))
message("Remote chimera_removal complete.")
