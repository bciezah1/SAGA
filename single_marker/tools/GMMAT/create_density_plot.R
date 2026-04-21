#!/usr/bin/env Rscript

suppressMessages({
  library(CMplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_density_plot.R <input_file> [output_dir]", call. = FALSE)
}

input_file <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], getwd())

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

output_prefix <- tools::file_path_sans_ext(basename(input_file))

cat("Reading data from:", input_file, "\n")
cat("Writing plots to:", output_dir, "\n")

data <- read.table(input_file, header = TRUE, sep = "\t")

required_cols <- c("CHR", "BP")
if (!all(required_cols %in% names(data))) {
  stop("Input file must contain columns: CHR, BP", call. = FALSE)
}

old_wd <- getwd()
setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)

output_file_prefix <- paste0(output_prefix, "_density")
cat("Creating SNP density plot:", file.path(output_dir, paste0(output_file_prefix, ".jpg")), "\n")

CMplot(
  data,
  plot.type = "d",
  bin.size = 1e6,
  chr.den.col = c("darkgreen", "yellow", "red"),
  file = "jpg",
  file.name = output_file_prefix,
  dpi = 300,
  main = paste0("SNP Density Plot - ", output_prefix),
  file.output = TRUE,
  verbose = TRUE,
  width = 9,
  height = 6
)

cat("Done! Plot saved to:", file.path(output_dir, paste0(output_file_prefix, ".jpg")), "\n")