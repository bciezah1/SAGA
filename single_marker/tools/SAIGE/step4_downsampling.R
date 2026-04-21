#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript step4_downsampling.R <input_frq_file> <output_snp_list>")
}

input_file <- args[1]
output_file <- args[2]

if (!file.exists(input_file)) {
  stop(paste("Input file does not exist:", input_file))
}

maf_data <- read.table(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("SNP", "MAF")
missing_cols <- setdiff(required_cols, colnames(maf_data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in .frq file:", paste(missing_cols, collapse = ", ")))
}

maf_data <- maf_data %>%
  filter(!is.na(SNP), SNP != "", !is.na(MAF)) %>%
  mutate(
    maf_bin = cut(
      MAF,
      breaks = c(0, 0.01, 0.05, 0.10, 0.20, 0.50),
      include.lowest = TRUE,
      right = TRUE
    )
  )

if (nrow(maf_data) == 0) {
  stop("No valid SNPs found after filtering MAF table.")
}

set.seed(123)

downsampled_snps <- maf_data %>%
  group_by(maf_bin) %>%
  group_modify(function(df, key) {
    n_bin <- nrow(df)
    n_keep <- max(1, ceiling(n_bin * 0.01))
    if (n_keep >= n_bin) {
      df
    } else {
      df %>% slice_sample(n = n_keep)
    }
  }) %>%
  ungroup() %>%
  distinct(SNP, .keep_all = TRUE) %>%
  arrange(SNP)

if (nrow(downsampled_snps) == 0) {
  stop("Downsampling produced zero SNPs.")
}

write.table(
  downsampled_snps$SNP,
  file = output_file,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

summary_file <- paste0(output_file, ".summary.tsv")
bin_summary <- maf_data %>%
  count(maf_bin, name = "n_before") %>%
  left_join(
    downsampled_snps %>% count(maf_bin, name = "n_after"),
    by = "maf_bin"
  )

write.table(
  bin_summary,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Downsampling completed successfully.\n")
cat(paste("SNP list written to:", output_file, "\n"))
cat(paste("Summary written to:", summary_file, "\n"))