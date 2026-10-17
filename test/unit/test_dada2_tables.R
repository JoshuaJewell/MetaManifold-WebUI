# Regression checks for CSV assembly, without loading DADA2 or dplyr.
# Run from the repository root: Rscript --vanilla test/unit/test_dada2_tables.R
expressions <- parse("src/pipeline/dada2/dada2_functions.r")
wanted <- c("base_left_join", "write_taxa_table", "write_combined_table")
loaded <- character()
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) {
    eval(expr)
    loaded <- c(loaded, as.character(expr[[2]]))
  }
}
stopifnot(setequal(loaded, wanted))

# Unsorted keys, repeated matches, unmatched rows, overlapping column names,
# NA keys, and helper-column name collisions must preserve the intended order.
x <- data.frame(value = 1:4, key = c("b", "a", "z", NA),
                .left_join_position = 11:14)
y <- data.frame(key = c("a", "b", "b", NA), value = 21:24,
                .right_join_position = 31:34)
got <- base_left_join(x, y, "key")
stopifnot(identical(got$key, c("b", "b", "a", "z", NA_character_)),
          identical(got$value.x, c(1L, 1L, 2L, 3L, 4L)),
          identical(got$value.y, c(22L, 23L, 21L, NA_integer_, 24L)),
          identical(names(got), c("value.x", "key", ".left_join_position",
                                   "value.y", ".right_join_position")))
stopifnot(nrow(base_left_join(x[FALSE, ], y, "key")) == 0L,
          nrow(base_left_join(x, y[FALSE, ], "key")) == nrow(x))
# Multiple keys: matching only one component must not produce a match.
a <- data.frame(id = c("a", "a"), key = c("x", "y"), value = 1:2)
b <- data.frame(id = c("a", "b"), key = c("y", "x"), other = 3:4)
stopifnot(identical(base_left_join(a, b, c("id", "key"))$other,
                    c(NA_integer_, 3L)))

# Rank names and their suffixed variants are also valid sample IDs.
# Repeated suffixing must preserve both the taxonomy strings and count columns.
tax_collision <- data.frame(Sequence = c("AA", "GG"), Species = c("taxA", "taxG"))
count_collision <- data.frame(Sequence = c("GG", "AA"), Species = c(10L, 20L),
                              Species.x = c(30L, 40L), Species.x.x = c(50L, 60L))
collision <- base_left_join(tax_collision, count_collision, "Sequence")
stopifnot(identical(names(collision),
                    c("Sequence", "Species.x.x.x", "Species.y", "Species.x", "Species.x.x")),
          identical(collision[[2]], c("taxA", "taxG")),
          identical(collision[[3]], c(20L, 10L)),
          identical(collision[[4]], c(40L, 30L)),
          identical(collision[[5]], c(60L, 50L)))

# Fixed expected CSV bytes cover every writer output, including non-syntactic
# sample names and an unmatched taxonomy row. No sequencing dependencies needed.
out <- tempfile("dada2-tables-")
dir.create(out)
index <- data.frame(SeqName = c("ASV3", "ASV1", "ASV2"),
                    Sequence = c("TT", "GG", "AA"))
tax <- matrix(c("G", "A"), ncol = 1,
              dimnames = list(c("GG", "AA"), "Kingdom"))
boot <- matrix(c(90L, 80L), ncol = 1, dimnames = dimnames(tax))
taxa <- write_taxa_table(tax, boot, index, out, "taxonomy")
counts <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), nrow = 2,
                 dimnames = list(c("sample-1", "sample 2"), c("AA", "TT", "GG")))
write_combined_table(taxa, index, counts, out, "tax_counts.csv", "asv_counts.csv")
expected <- list(
  taxonomy.csv = c("SeqName,Sequence,Kingdom", "ASV3,TT,NA", "ASV1,GG,G", "ASV2,AA,A"),
  taxonomy_bootstraps.csv = c("SeqName,Sequence,Kingdom_boot", "ASV3,TT,NA", "ASV1,GG,90", "ASV2,AA,80"),
  taxonomy_combined.csv = c("SeqName,Sequence,Kingdom,Kingdom_boot", "ASV3,TT,NA,NA", "ASV1,GG,G,90", "ASV2,AA,A,80"),
  tax_counts.csv = c("SeqName,Sequence,Kingdom,sample-1,sample 2", "ASV3,TT,NA,3,4", "ASV1,GG,G,5,6", "ASV2,AA,A,1,2"),
  asv_counts.csv = c("SeqName,Sequence,sample-1,sample 2", "ASV3,TT,3,4", "ASV1,GG,5,6", "ASV2,AA,1,2")
)
for (name in names(expected)) {
  path <- file.path(out, name)
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  stopifnot(identical(bytes, charToRaw(paste0(paste(expected[[name]], collapse = "\n"), "\n"))))
}
unlink(out, recursive = TRUE)
cat("DADA2 join and five CSV fixtures passed\n")
