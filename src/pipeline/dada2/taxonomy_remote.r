#!/usr/bin/env Rscript
# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Standalone taxonomy assignment for remote execution.
# Called by assign_taxonomy() in taxonomy.jl via SSH, through _run_remote_stage().
# All paths refer to the remote filesystem.
#
# Arguments (key=value):
#   functions   path to dada2_functions.r
#   staging     staging directory on this host
#   db          path to the taxonomy database; either uploaded under staging or
#               a pre-existing path named by databases.yml dada2.remote_path
#   prefix      taxonomy output prefix (e.g. "taxonomy")
#   multithread TRUE|FALSE|<integer>
#   min_boot    minimum bootstrap threshold
#   levels      comma-separated taxonomy level names
#   seed        master RNG seed
#   verbose     true|false

args <- commandArgs(trailingOnly = TRUE)

functions_path <- sub("^functions=", "", grep("^functions=", args, value = TRUE)[1])
if (is.na(functions_path)) stop("Missing required argument: functions")
source(functions_path)

p <- parse_remote_args(args, required = c("functions", "staging", "db", "prefix",
                                          "multithread", "min_boot", "levels", "seed"))

staging    <- p[["staging"]]
tables_dir <- file.path(staging, "Tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(staging, "ckpt_chimera.RData"))

verbose     <- tolower(p[["verbose"]]) == "true"
multithread <- remote_multithread(p[["multithread"]])
min_boot    <- as.integer(p[["min_boot"]])
levels      <- strsplit(p[["levels"]], ",", fixed = TRUE)[[1]]

# assignTaxonomy bootstraps by randomly subsampling kmers, so without a seed
# both the bootstrap values and the winning genus vary between runs on identical
# input. Seed here exactly as learn_errors and denoise do.
set.seed(as.integer(p[["seed"]]), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

taxa_result <- run_assign_taxonomy(
    seq_table_nochim, p[["db"]],
    list(multithread = multithread, min_boot = min_boot, levels = levels),
    verbose
)
taxa_df <- write_taxa_table(taxa_result$tax, taxa_result$boot, index,
                             tables_dir, p[["prefix"]])

save(seq_table_nochim, index, taxa_df, file = file.path(staging, "checkpoint.RData"))
message("Remote taxonomy assignment complete.")
