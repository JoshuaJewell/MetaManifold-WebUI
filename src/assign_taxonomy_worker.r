# Worker script for assignTaxonomy().
#
# Invoked as a subprocess by assign_taxonomy() in dada2.jl. On Unix this
# process is launched with `nice -n 10` so that the CPU-intensive naive
# Bayesian classification does not monopolise all cores and lock up the
# system. Memory for the reference database is also freed when this process
# exits rather than accumulating in the Julia/RCall session.
#
# Usage (internal - do not call directly):
#   Rscript assign_taxonomy_worker.r \
#       <functions_r> <chimera_ckpt> <db_path> <taxa_ckpt> \
#       <multithread> <min_boot> <levels_csv> <verbose>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8L) {
  stop("assign_taxonomy_worker.r: expected 8 arguments")
}

functions_r  <- args[[1L]]
chimera_ckpt <- args[[2L]]
db_path      <- args[[3L]]
taxa_ckpt    <- args[[4L]]
multithread  <- as.integer(args[[5L]])
min_boot     <- as.numeric(args[[6L]])
levels       <- strsplit(args[[7L]], ",", fixed = TRUE)[[1L]]
verbose      <- as.logical(args[[8L]])

source(functions_r)
load(chimera_ckpt)   # provides seq_table_nochim

taxa_result <- run_assign_taxonomy(seq_table_nochim, db_path,
                                   list(multithread = multithread,
                                        min_boot = min_boot,
                                        levels = levels),
                                   verbose)

save(taxa_result, file = taxa_ckpt)
message("Taxonomy assignment complete. Saved to ", taxa_ckpt)
