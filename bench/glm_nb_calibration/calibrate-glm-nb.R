# SPDX-License-Identifier: MPL-2.0
# Fresh reproducer for the small-n question raised during upstream PR 62.
# This is not the missing script that produced the earlier quoted rates.
# Run from a checkout of PR 62: Rscript /path/to/calibrate-glm-nb.R
args <- commandArgs(trailingOnly = TRUE)
source_path <- if (length(args)) args[[1]] else "src/analysis/differential.jl"
text <- paste(readLines(source_path, warn = FALSE), collapse = "\n")
parts <- strsplit(text, 'const _FIT_R = raw"""', fixed = TRUE)[[1]]
stopifnot(length(parts) == 2L)
fit_code <- strsplit(parts[[2]], '"""', fixed = TRUE)[[1]][[1]]
read_constant <- function(name) {
  pattern <- paste0("const ", name, " = ([0-9.eE+-]+)")
  match <- regmatches(text, regexec(pattern, text))[[1]]
  stopifnot(length(match) == 2L)
  as.numeric(match[[2]])
}
stopifnot(requireNamespace("MASS", quietly = TRUE))
# Enable the optional comparison explicitly; never substitute another test.
compare_gp <- identical(Sys.getenv("COMPARE_GLMGAMPOI"), "1")
if (compare_gp) {
  stopifnot(requireNamespace("glmGamPoi", quietly = TRUE))
  stopifnot(packageVersion("glmGamPoi") == "1.22.0")
}
cat(R.version.string, "\nMASS:", as.character(packageVersion("MASS")), "\n")
if (compare_gp) cat("glmGamPoi:", as.character(packageVersion("glmGamPoi")), "\n")
cat("Source MD5:", unname(tools::md5sum(source_path)), "\n")
reps <- as.integer(Sys.getenv("CALIBRATION_REPS", "10"))
taxa <- as.integer(Sys.getenv("CALIBRATION_TAXA", "600"))
stopifnot(reps >= 2L, taxa >= 10L)
rows <- list()
for (n in c(4L, 12L)) {
  for (replicate_id in seq_len(reps)) {
    seed <- 620000L + 1000L * n + replicate_id
    set.seed(seed)
    # All taxa are null: the same mean and dispersion in both groups.
    mu <- exp(runif(taxa, log(5), log(1000)))
    dispersion <- 0.05 + 1 / sqrt(mu)
    exposure <- exp(rnorm(2L * n, sd = 0.3))
    counts <- matrix(rnbinom(2L * n * taxa,
                            mu = as.vector(outer(exposure, mu)),
                            size = rep(1 / dispersion, each = 2L * n)),
                     nrow = 2L * n)
    library_sizes <- rowSums(counts)
    stopifnot(all(library_sizes > 0))
    offset <- log(library_sizes / exp(mean(log(library_sizes))))
    retained <- colMeans(counts > 0) >= 0.1
    for (method in if (compare_gp) c("glm_nb", "glmGamPoi") else "glm_nb") {
      env <- new.env(parent = globalenv())
      env$da_counts <- counts[, retained, drop = FALSE]
      env$da_group <- rep(c("A", "B"), each = n)
      env$da_levels <- c("A", "B")
      env$da_offset <- offset
      env$da_theta_upper <- read_constant("THETA_UPPER")
      env$da_theta_lower <- read_constant("THETA_LOWER")
      env$da_use_glp <- method == "glmGamPoi"
      eval(parse(text = fit_code), envir = env)
      result <- env$da_result
      valid <- result$status != "failed" & is.finite(result$pvalue)
      rejected <- sum(result$pvalue[valid] < 0.05)
      rows[[length(rows) + 1L]] <- data.frame(
        n_per_group = n, replicate = replicate_id, seed = seed, method = method,
        generated = taxa, retained = sum(retained), valid = sum(valid),
        failed = sum(!valid), boundary = sum(result$status == "boundary"),
        rejected = rejected, rate = rejected / sum(valid))
    }
  }
}
results <- do.call(rbind, rows)
write.csv(results, "calibration-replicates.csv", row.names = FALSE)
for (n in unique(results$n_per_group)) {
  for (method in unique(results$method)) {
    d <- subset(results, n_per_group == n & results$method == method)
    cat(sprintf("n=%d method=%s valid=%d failed=%d boundary=%d rejected=%d pooled_rate=%.6f mean_replicate_rate=%.6f replicate_SE=%.6f\n",
                n, method, sum(d$valid), sum(d$failed), sum(d$boundary),
                sum(d$rejected), sum(d$rejected) / sum(d$valid),
                mean(d$rate), sd(d$rate) / sqrt(nrow(d))))
  }
}
# Replicates are independent; taxa within one matrix share estimated TSS
# offsets (and, optionally, a dispersion trend). Do not use a naive pooled
# binomial standard error as though all taxon tests were independent.
print(sessionInfo())
