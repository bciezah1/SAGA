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

if (!file.exists(input_file)) {
  stop(paste("Input file does not exist:", input_file), call. = FALSE)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

prefix <- sub("\\.[^.]+$", "", basename(input_file))

cat("Reading data from:", input_file, "\n")
cat("Writing plots to:", output_dir, "\n")

data <- read.table(
  input_file,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

colnames_present <- colnames(data)

has_metal_format <- all(c("MarkerName", "P-value") %in% colnames_present)
has_manhattan_format <- all(c("CHR", "BP", "SNP", "PVAL") %in% colnames_present)

if (!has_metal_format && !has_manhattan_format) {
  stop(
    paste(
      "Input file must contain either:",
      "(1) MarkerName and P-value, or",
      "(2) CHR, BP, SNP, and PVAL."
    ),
    call. = FALSE
  )
}

if (has_metal_format) {
  cat("Detected METAL-style input format: MarkerName + P-value\n")
  cat("Parsing chromosome and base-pair position from MarkerName...\n")

  marker_parts <- strsplit(data$MarkerName, ":", fixed = TRUE)

  data$CHR <- suppressWarnings(as.integer(
    vapply(marker_parts, function(x) {
      if (length(x) >= 1) x[1] else NA_character_
    }, character(1))
  ))

  data$BP <- suppressWarnings(as.integer(
    vapply(marker_parts, function(x) {
      if (length(x) >= 2) x[2] else NA_character_
    }, character(1))
  ))

  data$SNP <- data$MarkerName
  data$PVAL <- suppressWarnings(as.numeric(data[["P-value"]]))
} else {
  cat("Detected Manhattan-ready input format: CHR + BP + SNP + PVAL\n")

  data$CHR <- suppressWarnings(as.integer(data$CHR))
  data$BP <- suppressWarnings(as.integer(data$BP))
  data$SNP <- as.character(data$SNP)
  data$PVAL <- suppressWarnings(as.numeric(data$PVAL))
}

data <- data[
  !is.na(data$SNP) &
  !is.na(data$CHR) &
  !is.na(data$BP) &
  !is.na(data$PVAL),
  ,
  drop = FALSE
]

data <- data[data$PVAL > 0 & data$PVAL <= 1, , drop = FALSE]
data <- data[data$CHR >= 1 & data$CHR <= 22, , drop = FALSE]

if (nrow(data) == 0) {
  stop("No valid variants remain after filtering.", call. = FALSE)
}

data <- data[order(data$CHR, data$BP), , drop = FALSE]

cat("Total variants for plotting:", nrow(data), "\n")

chisq <- qchisq(1 - data$PVAL, 1)
lambda_gc <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, 1), 3)
cat("Lambda GC =", lambda_gc, "\n")

highlight_threshold <- 5e-8
snps_to_highlight <- data$SNP[data$PVAL < highlight_threshold]

cat(
  "Number of SNPs to highlight (P < ",
  highlight_threshold,
  "): ",
  length(snps_to_highlight),
  "\n",
  sep = ""
)

if (length(snps_to_highlight) > 0) {
  cat("Highlighted SNPs:\n")
  print(snps_to_highlight)
} else {
  cat("No genome-wide significant SNPs to highlight.\n")
}

old_wd <- getwd()
setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)

output_file_prefix <- prefix
cat("Creating Manhattan plot:", file.path(output_dir, paste0(output_file_prefix, ".jpg")), "\n")

thresholds <- c(1e-5, 5e-8)
threshold_colors <- c("black", "black")
threshold_lty <- c(2, 2)
chrom_colors <- c("lightblue", "blue")

plot_data <- data[, c("SNP", "CHR", "BP", "PVAL"), drop = FALSE]
colnames(plot_data) <- c("SNP", "CHR", "BP", "P")

CMplot(
  plot_data,
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
  file.name = output_file_prefix,
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
  main = paste0(output_file_prefix, " Manhattan Plot (Lambda = ", lambda_gc, ")")
)

cat("Done! Manhattan plot saved to:", file.path(output_dir, paste0(output_file_prefix, ".jpg")), "\n")