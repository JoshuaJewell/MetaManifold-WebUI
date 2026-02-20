#!/usr/bin/env Rscript
# Standalone taxonomy assignment for remote execution.
# Called by _assign_taxonomy_remote() in taxonomy.jl via SSH.
# All paths refer to the remote filesystem.
#
# Arguments (key=value):
#   functions   path to dada2_functions.r
#   ckpt        path to ckpt_chimera.RData
#   db          path to taxonomy database
#   tables      output directory for taxonomy CSV files
#   save        path to write checkpoint.RData
#   prefix      taxonomy output prefix (e.g. "taxonomy")
#   multithread number of threads
#   min_boot    minimum bootstrap threshold
#   levels      comma-separated taxonomy level names
#   verbose     true|false

args <- commandArgs(trailingOnly = TRUE)
p <- list()
for (a in args) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) >= 2) p[[kv[1]]] <- paste(kv[-1], collapse = "=")
}

required <- c("functions", "ckpt", "db", "tables", "save",
               "prefix", "multithread", "min_boot", "levels")
missing_args <- setdiff(required, names(p))
if (length(missing_args) > 0)
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))

source(p[["functions"]])

load(p[["ckpt"]])
dir.create(p[["tables"]], recursive = TRUE, showWarnings = FALSE)

verbose     <- tolower(p[["verbose"]]) == "true"
multithread <- as.integer(p[["multithread"]])
min_boot    <- as.integer(p[["min_boot"]])
levels      <- strsplit(p[["levels"]], ",", fixed = TRUE)[[1]]

taxa_result <- run_assign_taxonomy(
    seq_table_nochim, p[["db"]],
    list(multithread = multithread, min_boot = min_boot, levels = levels),
    verbose
)
taxa_df <- write_taxa_table(taxa_result$tax, taxa_result$boot, index,
                             p[["tables"]], p[["prefix"]])

save(seq_table_nochim, index, taxa_df, file = p[["save"]])
message("Remote taxonomy assignment complete.")
