#!/usr/bin/env Rscript

# ==============================
# Manhattan Plot with highlighted SNPs and threshold lines
# ==============================

suppressMessages({
  library(CMplot)
})

# ======== Parse command-line arguments ========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_manhattan_highlight.R <input_file>", call. = FALSE)
}

input_file <- args[1]
prefix <- sub("\\.txt$", "", basename(input_file))

cat("Reading data from:", input_file, "\n")

# ======== Load data ========
data <- read.table(input_file, header = TRUE, sep = "\t")

# Check required columns
required_cols <- c("SNP", "CHR", "BP", "PVAL")
if (!all(required_cols %in% colnames(data))) {
  stop("Input file must contain columns: SNP, CHR, BP, PVAL")
}

# ======== Compute lambda GC ========
chisq <- qchisq(1 - data$PVAL, 1)
lambda_gc <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, 1), 3)
cat("Lambda GC =", lambda_gc, "\n")

# ======== Identify SNPs to highlight ========
highlight_threshold <- 1e-4
SNPs_to_highlight <- data$SNP[data$PVAL < highlight_threshold]
cat("Number of SNPs to highlight (P < ", highlight_threshold, "): ", length(SNPs_to_highlight), "\n")

# ======== Manhattan plot ========
output_file <- paste0(prefix, "_manhattan_highlight.jpg")
cat("Creating Manhattan plot:", output_file, "\n")

# Threshold lines
thresholds <- c(1e-5, 5e-8)           # black dotted line at 1e-5, red dotted line at 5e-8
threshold_colors <- c("black", "red")
threshold_lty <- c(2, 2)

CMplot(
  data,
  type = "h",                          # vertical lines
  plot.type = "m",                      # Manhattan plot
  LOG10 = TRUE,                         # -log10(P)
  highlight = SNPs_to_highlight,
  highlight.type = "p",
  highlight.col = NULL,
  highlight.cex = 1.2,
  highlight.pch = 19,
  file = "jpg",
  file.name = paste0(prefix, "_manhattan_highlight"),
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  width = 14,
  height = 6,
  band = 0.6,
  threshold = thresholds,
  threshold.col = threshold_colors,
  threshold.lty = threshold_lty,
  main = paste0("Manhattan Plot (Lambda = ", lambda_gc, ")")
)

cat("Done! Manhattan plot saved to:", output_file, "\n")
