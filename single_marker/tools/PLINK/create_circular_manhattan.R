#!/usr/bin/env Rscript

# ==============================
# Circular Manhattan Plot with lambda
# ==============================

suppressMessages({
  library(CMplot)
})

# ======== Parse command-line arguments ========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_circular_manhattan.R <input_file>", call. = FALSE)
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

# ======== Circular Manhattan plot ========
output_file <- paste0(prefix, "_circular.jpg")
cat("Creating circular Manhattan plot:", output_file, "\n")

CMplot(
  data,
  type = "p",
  plot.type = "c",
  chr.labels = paste("Chr", c(1:22), sep = ""),   # your chromosomes
  r = 0.4,
  cir.chr.h = 1.3,
  chr.den.col = "black",
  file = "jpg",
  file.name = paste0(prefix, "_circular"),
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  width = 10,
  height = 10,
  main = paste0("Circular Manhattan Plot (? = ", lambda_gc, ")")
)

cat("Done! Plot saved to:", output_file, "\n")
