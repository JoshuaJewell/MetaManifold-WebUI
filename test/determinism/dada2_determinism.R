# DADA2 determinism experiment: run on the compute server, not in CI.
#
# Answers, on real reads, the questions the Julia unit suite cannot:
#
#   repeat      the same stage twice with the same seed and threads is identical
#   threads     1 thread versus N threads is identical
#   rngkind     a session whose RNG kind was changed first is identical, given
#               the pinned set.seed() call the pipeline now makes
#   seed        assignTaxonomy's bootstraps genuinely depend on the seed (so a
#               fixed seed is load-bearing, and the repeat check is not vacuous)
#   production  re-running assign_taxonomy with the production run's seed and
#               thread count reproduces the stored taxonomy and bootstraps
#
# Usage:
#   Rscript dada2_determinism.R <reads_dir> <n_samples> <threads> <tax_db> \
#           <chimera_ckpt> <tables_dir> <prod_threads> <out_dir>
#
# Writes <out_dir>/report.tsv (check, stage, status, detail) and exits non-zero
# if any check fails.

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 8)
reads_dir    <- args[1]
n_samples    <- as.integer(args[2])
threads      <- as.integer(args[3])
tax_db       <- args[4]
chimera_ckpt <- args[5]
tables_dir   <- args[6]
prod_threads <- as.integer(args[7])
out_dir      <- args[8]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SEED <- 123L
# Exactly the call the pipeline makes (src/pipeline/dada2/*.jl and *_remote.r).
pinned_seed <- function(seed) set.seed(seed, kind = "Mersenne-Twister",
                                       normal.kind = "Inversion", sample.kind = "Rejection")
perturb_rng <- function() {
  RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rounding")
  invisible(runif(12345))
}

report <- data.frame(check = character(), stage = character(), status = character(),
                     detail = character(), stringsAsFactors = FALSE)
record <- function(check, stage, ok, detail = "") {
  status <- if (isTRUE(ok)) "PASS" else "FAIL"
  report[nrow(report) + 1, ] <<- list(check, stage, status, detail)
  cat(sprintf("[%s] %-10s %-20s %s\n", status, check, stage, detail))
  flush.console()
}
timed <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- force(expr)
  cat(sprintf("    %s: %.1fs\n", label, proc.time()[["elapsed"]] - t0))
  flush.console()
  v
}
# identical() alone gives no detail on failure.
same <- function(a, b) {
  ok <- identical(a, b)
  detail <- if (ok) "identical" else
    paste(head(capture.output(print(all.equal(a, b))), 3), collapse = " | ")
  list(ok = ok, detail = detail)
}

cat(sprintf("R %s, dada2 %s, RcppParallel %s, host %s, %d threads\n",
            getRversion(), packageVersion("dada2"), packageVersion("RcppParallel"),
            Sys.info()[["nodename"]], threads))

## Reads: the n_samples largest forward/reverse pairs, for a realistic workload.
fwd_all <- list.files(reads_dir, pattern = "_R1_filt\\.fastq\\.gz$", full.names = TRUE)
fwd <- sort(head(fwd_all[order(-file.size(fwd_all))], n_samples))
rev <- sub("_R1_filt", "_R2_filt", fwd)
stopifnot(length(fwd) > 0, all(file.exists(rev)))
cat(sprintf("Using %d samples, %.1f MB of reads\n", length(fwd),
            sum(file.size(c(fwd, rev))) / 1e6))

## learnErrors: repeat, and 1 vs N threads.
NBASES <- 2e7
learn <- function(th) {
  pinned_seed(SEED)
  list(f = learnErrors(fwd, nbases = NBASES, MAX_CONSIST = 15, multithread = th, verbose = FALSE),
       r = learnErrors(rev, nbases = NBASES, MAX_CONSIST = 15, multithread = th, verbose = FALSE))
}
errN  <- timed("learnErrors N",        learn(threads))
errN2 <- timed("learnErrors N repeat", learn(threads))
err1  <- timed("learnErrors 1",        learn(1L))
s <- same(getErrors(errN$f), getErrors(errN2$f)); record("repeat",  "learnErrors", s$ok, s$detail)
s <- same(getErrors(errN$f), getErrors(err1$f));  record("threads", "learnErrors.fwd", s$ok, s$detail)
s <- same(getErrors(errN$r), getErrors(err1$r));  record("threads", "learnErrors.rev", s$ok, s$detail)

## dada (pseudo-pooled, as production) + mergePairs + chimera removal.
denoise <- function(th, err) {
  pinned_seed(SEED)
  df <- dada(fwd, err = err$f, pool = "pseudo", multithread = th, verbose = FALSE)
  dr <- dada(rev, err = err$r, pool = "pseudo", multithread = th, verbose = FALSE)
  mg <- mergePairs(df, fwd, dr, rev, verbose = FALSE)
  st <- makeSequenceTable(mg)
  list(st = st,
       nochim = removeBimeraDenovo(st, method = "consensus", multithread = th, verbose = FALSE))
}
dnN  <- timed("denoise N",        denoise(threads, errN))
dnN2 <- timed("denoise N repeat", denoise(threads, errN))
dn1  <- timed("denoise 1",        denoise(1L, errN))
s <- same(dnN$st, dnN2$st);        record("repeat",  "dada+merge", s$ok, s$detail)
s <- same(dnN$st, dn1$st);         record("threads", "dada+merge", s$ok, s$detail)
s <- same(dnN$nochim, dn1$nochim); record("threads", "removeBimera", s$ok,
                                          sprintf("%s (%d ASVs)", s$detail, ncol(dnN$nochim)))

perturb_rng()
dnK <- timed("denoise N after RNG perturbation", denoise(threads, errN))
s <- same(dnN$nochim, dnK$nochim); record("rngkind", "denoise", s$ok, s$detail)

## assignTaxonomy on the production run's chimera-free table.
load(chimera_ckpt)   # seq_table_nochim, index
seqs <- colnames(seq_table_nochim)
cat(sprintf("assignTaxonomy on %d production ASVs against %s\n", length(seqs), basename(tax_db)))
levels <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order",
            "Family", "Genus", "Species")
tax <- function(sq, th, seed = SEED, perturb = FALSE) {
  if (perturb) perturb_rng()
  pinned_seed(seed)
  assignTaxonomy(sq, tax_db, multithread = th, minBoot = 0, outputBootstraps = TRUE,
                 taxLevels = levels, verbose = FALSE)
}
tN  <- timed("assignTaxonomy N",        tax(seqs, threads))
tN2 <- timed("assignTaxonomy N repeat", tax(seqs, threads))
s <- same(tN, tN2); record("repeat", "assignTaxonomy", s$ok, s$detail)

tK <- timed("assignTaxonomy N perturbed RNG", tax(seqs, threads, perturb = TRUE))
s <- same(tN, tK); record("rngkind", "assignTaxonomy", s$ok, s$detail)

tS <- timed("assignTaxonomy N seed 124", tax(seqs, threads, seed = 124L))
nboot <- sum(tS$boot != tN$boot)
record("seed", "assignTaxonomy", nboot > 0,
       sprintf("%d bootstrap cells differ between seed 123 and 124", nboot))

# Single-threaded assignTaxonomy over the whole table is slow; a 500-ASV slice
# answers the thread question just as well.
sub <- head(seqs, 500)
tsN <- timed("assignTaxonomy slice N", tax(sub, threads))
ts1 <- timed("assignTaxonomy slice 1", tax(sub, 1L))
s <- same(tsN, ts1); record("threads", "assignTaxonomy", s$ok, s$detail)
tsP <- timed("assignTaxonomy slice prod threads", tax(sub, prod_threads))
s <- same(tsN, tsP); record("threads", "assignTaxonomy.prod", s$ok,
                            sprintf("%s (%d vs %d threads)", s$detail, threads, prod_threads))

## Production reproduction against the tables the real run wrote.
prod_tax  <- read.csv(file.path(tables_dir, "taxonomy.csv"), stringsAsFactors = FALSE,
                      check.names = FALSE)
prod_boot <- read.csv(file.path(tables_dir, "taxonomy_bootstraps.csv"), stringsAsFactors = FALSE,
                      check.names = FALSE)
cat("Stored taxonomy.csv columns:  ", paste(head(colnames(prod_tax), 12), collapse = ", "), "\n")
cat("Stored bootstraps.csv columns:", paste(head(colnames(prod_boot), 12), collapse = ", "), "\n")
tP <- timed("assignTaxonomy prod threads", tax(seqs, prod_threads))

rank_cols <- intersect(levels, colnames(prod_tax))
boot_cols <- colnames(prod_boot)[sub("_boot$", "", colnames(prod_boot)) %in% levels]
if (length(rank_cols) > 0 && nrow(prod_tax) == length(seqs)) {
  mine <- tP$tax[, rank_cols, drop = FALSE]; mine[is.na(mine)] <- ""
  theirs <- as.matrix(prod_tax[, rank_cols, drop = FALSE]); theirs[is.na(theirs)] <- ""
  diff_cells <- sum(mine != theirs)
  record("production", "taxonomy.csv", diff_cells == 0,
         sprintf("%d of %d rank cells differ", diff_cells, length(mine)))
} else {
  record("production", "taxonomy.csv", FALSE,
         sprintf("layout not comparable: %d stored rows, %d ASVs, rank cols: %s",
                 nrow(prod_tax), length(seqs), paste(rank_cols, collapse = ",")))
}
if (length(boot_cols) > 0 && nrow(prod_boot) == length(seqs)) {
  mine_b <- tP$boot[, sub("_boot$", "", boot_cols), drop = FALSE]
  diff_b <- sum(mine_b != as.matrix(prod_boot[, boot_cols, drop = FALSE]))
  record("production", "bootstraps.csv", diff_b == 0,
         sprintf("%d of %d bootstrap cells differ", diff_b, length(mine_b)))
} else {
  record("production", "bootstraps.csv", FALSE,
         sprintf("layout not comparable: %d stored rows, boot cols: %s",
                 nrow(prod_boot), paste(boot_cols, collapse = ",")))
}

write.table(report, file.path(out_dir, "report.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(errN = errN, dnN = dnN, tN = tN, tP = tP), file.path(out_dir, "results.rds"))
cat(sprintf("\n%d checks, %d failed\n", nrow(report), sum(report$status == "FAIL")))
quit(status = if (any(report$status == "FAIL")) 1L else 0L)
