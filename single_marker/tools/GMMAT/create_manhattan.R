#!/usr/bin/env Rscript

suppressMessages({
  library(CMplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_manhattan.R <input_file> [output_dir]", call. = FALSE)
}

input_file <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], getwd())

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

prefix <- sub("\\.txt$", "", basename(input_file))

cat("Reading data from:", input_file, "\n")
cat("Writing plots to:", output_dir, "\n")

data <- read.table(input_file, header = TRUE, sep = "\t")

required_cols <- c("SNP", "CHR", "BP", "PVAL")
if (!all(required_cols %in% colnames(data))) {
  stop("Input file must contain columns: SNP, CHR, BP, PVAL")
}

chisq <- qchisq(1 - data$PVAL, 1)
lambda_gc <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, 1), 3)
cat("Lambda GC =", lambda_gc, "\n")

highlight_threshold <- 5e-8
snps_to_highlight <- data$SNP[data$PVAL < highlight_threshold]
cat("Number of SNPs to highlight (P < ", highlight_threshold, "): ", length(snps_to_highlight), "\n")

old_wd <- getwd()
setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)

output_file <- paste0(prefix, "_manhattan_highlight.jpg")
cat("Creating Manhattan plot:", file.path(output_dir, output_file), "\n")

thresholds <- c(1e-5, 5e-8)
threshold_colors <- c("black", "black")
threshold_lty <- c(2, 2)

chrom_colors <- c("lightblue", "blue")

CMplot(
  data,
  type = "p",
  plot.type = "m",
  LOG10 = TRUE,
  col = chrom_colors,
  highlight = snps_to_highlight,
  highlight.type = "p",
  highlight.col = "red",
  highlight.cex = 0.7,
  highlight.pch = 20,
  file = "jpg",
  file.name = paste0(prefix, "_manhattan_highlight"),
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  width = 14,
  height = 6,
  band = 0.6,
  cex = 0.45,
  threshold = thresholds,
  threshold.col = threshold_colors,
  threshold.lty = threshold_lty,
  main = paste0("Manhattan Plot (Lambda = ", lambda_gc, ")")
)

cat("Done! Manhattan plot saved to:", file.path(output_dir, output_file), "\n")