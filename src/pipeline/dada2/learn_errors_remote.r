#!/usr/bin/env Rscript
# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Standalone error-model learning for remote execution.
# Called by learn_errors() in denoise.jl via SSH, through _run_remote_stage().
# All paths refer to the remote filesystem.
#
# The filtered reads this stage learns from are uploaded to $staging/Filtered
# and named by the manifests, so the stage depends on nothing the server may
# have kept from an earlier one.
#
# Arguments (key=value):
#   functions    path to dada2_functions.r
#   staging      staging directory on this host
#   nbases       bases sampled when learning the error model
#   max_consist  self-consistency iterations
#   seed         random seed
#   multithread  TRUE|FALSE|<integer>
#   verbose      true|false

args <- commandArgs(trailingOnly = TRUE)

# dada2_functions.r carries parse_remote_args, so its own path is pulled out by
# hand: the helper that would read it does not exist until it is sourced.
functions_path <- sub("^functions=", "", grep("^functions=", args, value = TRUE)[1])
if (is.na(functions_path)) stop("Missing required argument: functions")
source(functions_path)

p <- parse_remote_args(args, required = c("functions", "staging", "nbases",
                                          "max_consist", "seed", "multithread"))

filtered_dir <- file.path(p[["staging"]], "Filtered")
# Read basenames arrive as $staging/fwd.txt and $staging/rev.txt, written by the
# Julia side and uploaded with the reads. A mode that carries no reads on one
# side simply has no manifest there.
fwd_names    <- read_manifest(file.path(p[["staging"]], "fwd.txt"))
rev_names    <- read_manifest(file.path(p[["staging"]], "rev.txt"))
fwd_files    <- if (length(fwd_names)) file.path(filtered_dir, fwd_names) else character(0)
rev_files    <- if (length(rev_names)) file.path(filtered_dir, rev_names) else character(0)

verbose     <- tolower(p[["verbose"]]) == "true"
multithread <- remote_multithread(p[["multithread"]])
nbases      <- as.numeric(p[["nbases"]])
max_consist <- as.integer(p[["max_consist"]])

set.seed(as.integer(p[["seed"]]), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

# Mirrors the local stage exactly: a mode that carries no reads on one side
# leaves that error model NULL rather than learning from an empty vector.
if (length(fwd_files) > 0) {
    fwd_errors <- learnErrors(fwd_files, nbases = nbases, MAX_CONSIST = max_consist,
                              multithread = multithread, verbose = verbose)
} else {
    fwd_errors <- NULL
}

if (length(rev_files) > 0) {
    rev_errors <- learnErrors(rev_files, nbases = nbases, MAX_CONSIST = max_consist,
                              multithread = multithread, verbose = verbose)
} else {
    rev_errors <- NULL
}

plot_error_rates(fwd_errors, rev_errors, file.path(p[["staging"]], "error_rates.pdf"))
save(fwd_errors, rev_errors, file = file.path(p[["staging"]], "ckpt_errors.RData"))
message("Remote learn_errors complete.")
