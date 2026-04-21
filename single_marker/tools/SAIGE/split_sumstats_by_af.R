#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript split_sumstats_by_af.R /path/to/sum_stats.txt", call. = FALSE)
}

input_file <- args[1]

if (!file.exists(input_file)) {
  stop(paste("Input file does not exist:", input_file), call. = FALSE)
}

cat("Reading input file:", input_file, "\n")

dat <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

required_cols <- c("SNP", "CHR", "BP", "AF", "PVAL")
missing_cols <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop(
    paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
    call. = FALSE
  )
}

out_dir <- dirname(input_file)

file_af_lt <- file.path(out_dir, "sum_stats_AF_lt_0.01.txt")
file_af_ge <- file.path(out_dir, "sum_stats_AF_ge_0.01.txt")

file_af_lt_min <- file.path(out_dir, "sum_stats_AF_lt_0.01_CHR_BP_SNP_PVAL.txt")
file_af_ge_min <- file.path(out_dir, "sum_stats_AF_ge_0.01_CHR_BP_SNP_PVAL.txt")

dat_af_lt <- dat[dat$AF < 0.01, , drop = FALSE]
dat_af_ge <- dat[dat$AF >= 0.01, , drop = FALSE]

write.table(
  dat_af_lt,
  file = file_af_lt,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  dat_af_ge,
  file = file_af_ge,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cols_to_keep <- c("CHR", "BP", "SNP", "PVAL")

dat_af_lt_min <- dat_af_lt[, cols_to_keep, drop = FALSE]
dat_af_ge_min <- dat_af_ge[, cols_to_keep, drop = FALSE]

write.table(
  dat_af_lt_min,
  file = file_af_lt_min,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  dat_af_ge_min,
  file = file_af_ge_min,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cat("Done.\n")
cat("Files written:\n")
cat("  ", file_af_lt, "\n")
cat("  ", file_af_ge, "\n")
cat("  ", file_af_lt_min, "\n")
cat("  ", file_af_ge_min, "\n")
cat("Rows with AF < 0.01 :", nrow(dat_af_lt), "\n")
cat("Rows with AF >= 0.01:", nrow(dat_af_ge), "\n")
