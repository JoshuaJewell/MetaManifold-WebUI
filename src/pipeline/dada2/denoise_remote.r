#!/usr/bin/env Rscript
# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Standalone denoising and sequence-table construction for remote execution.
# Called by denoise() in denoise.jl via SSH, through _run_remote_stage().
# All paths refer to the remote filesystem.
#
# The error model is uploaded with the reads rather than assumed to be on the
# server, so this stage runs whether or not learn_errors ran here.
#
# Arguments (key=value):
#   functions     path to dada2_functions.r
#   staging       staging directory on this host
#   mode          paired|forward|reverse
#   pool_method   value for dada(pool=), as written in the config
#   pool_bool     true when pool_method is to be read as a logical, not a string
#   min_overlap   mergePairs minOverlap    (paired mode)
#   max_mismatch  mergePairs maxMismatch   (paired mode)
#   trim_overhang mergePairs trimOverhang  (paired mode)
#   seed          random seed
#   multithread   TRUE|FALSE|<integer>
#   verbose       true|false

args <- commandArgs(trailingOnly = TRUE)

functions_path <- sub("^functions=", "", grep("^functions=", args, value = TRUE)[1])
if (is.na(functions_path)) stop("Missing required argument: functions")
source(functions_path)

p <- parse_remote_args(args, required = c("functions", "staging", "mode",
                                          "pool_method", "pool_bool", "seed",
                                          "multithread"))

filtered_dir <- file.path(p[["staging"]], "Filtered")
# Read basenames arrive as $staging/fwd.txt and $staging/rev.txt, written by the
# Julia side and uploaded with the reads. A mode that carries no reads on one
# side simply has no manifest there.
fwd_names    <- read_manifest(file.path(p[["staging"]], "fwd.txt"))
rev_names    <- read_manifest(file.path(p[["staging"]], "rev.txt"))
fwd_files    <- if (length(fwd_names)) file.path(filtered_dir, fwd_names) else character(0)
rev_files    <- if (length(rev_names)) file.path(filtered_dir, rev_names) else character(0)

mode        <- p[["mode"]]
verbose     <- tolower(p[["verbose"]]) == "true"
multithread <- remote_multithread(p[["multithread"]])

# dada(pool=) distinguishes the logical TRUE from the string "pseudo", and the
# config may spell either. The command line flattens both to text, so the Julia
# side sends the type alongside the value rather than having this script guess:
# guessing would silently turn a config that says pool_method: "true" into a
# different denoising run from the one the local stage performs.
pool_method <- if (tolower(p[["pool_bool"]]) == "true") {
    tolower(p[["pool_method"]]) == "true"
} else {
    p[["pool_method"]]
}

load(file.path(p[["staging"]], "ckpt_errors.RData"))
set.seed(as.integer(p[["seed"]]), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

if (length(fwd_files) > 0) {
    dada_fwd <- dada(fwd_files, err = fwd_errors, pool = pool_method,
                     multithread = multithread, verbose = verbose)
} else {
    dada_fwd <- NULL
}

if (length(rev_files) > 0) {
    dada_rev <- dada(rev_files, err = rev_errors, pool = pool_method,
                     multithread = multithread, verbose = verbose)
} else {
    dada_rev <- NULL
}

# A single input comes back as a bare object; name it so the sequence table
# carries the sample name.
name_one <- function(x, files) if (inherits(x, "dada") || is.data.frame(x)) setNames(list(x), basename(files)) else x
if (!is.null(dada_fwd)) dada_fwd <- name_one(dada_fwd, fwd_files)
if (!is.null(dada_rev)) dada_rev <- name_one(dada_rev, rev_files)

if (mode == "paired") {
    merged <- mergePairs(dada_fwd, fwd_files, dada_rev, rev_files,
                         minOverlap   = as.integer(p[["min_overlap"]]),
                         maxMismatch  = as.integer(p[["max_mismatch"]]),
                         trimOverhang = tolower(p[["trim_overhang"]]) == "true",
                         verbose      = verbose)
    merged    <- name_one(merged, fwd_files)
    seq_table <- makeSequenceTable(merged)
} else {
    merged    <- NULL
    seq_table <- makeSequenceTable(if (!is.null(dada_fwd)) dada_fwd else dada_rev)
}

# Written only when there is something to plot, matching the local stage. The
# Julia side retrieves this file as an optional output for that reason.
if (sum(seq_table) > 0) {
    plot_length_distribution(seq_table, file.path(p[["staging"]], "length_distribution.pdf"))
} else {
    message("Skipping length distribution plot: seq_table is empty after merging")
}

save(dada_fwd, dada_rev, merged, seq_table,
     file = file.path(p[["staging"]], "ckpt_denoise.RData"))
message("Remote denoise complete.")
