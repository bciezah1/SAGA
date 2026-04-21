#!/usr/bin/env Rscript

suppressMessages({
  library(CMplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_circular_manhattan.R <input_file> [output_dir]", call. = FALSE)
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

old_wd <- getwd()
setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)

output_file <- paste0(prefix, "_circular.jpg")
cat("Creating circular Manhattan plot:", file.path(output_dir, output_file), "\n")
data <- data[, c("SNP", "CHR", "BP", "PVAL")]
colnames(data) <- c("SNP", "CHR", "BP", "P")

CMplot(
  data,
  type = "p",
  plot.type = "c",
  chr.labels = paste("Chr", c(1:22), sep = ""),
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
  main = paste0("Circular Manhattan Plot (Lambda = ", lambda_gc, ")")
)

cat("Done! Plot saved to:", file.path(output_dir, output_file), "\n")