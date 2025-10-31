#!/usr/bin/env Rscript

# ======== Load necessary libraries ========
suppressMessages({
  library(CMplot)
})

# ======== Parse command-line arguments ========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_density_plot.R <input_file> [output_prefix]", call. = FALSE)
}

input_file <- args[1]
output_prefix <- ifelse(length(args) >= 2, args[2], tools::file_path_sans_ext(basename(input_file)))

# ======== Load data ========
cat("Reading data from:", input_file, "\n")
data <- read.table(input_file, header = TRUE, sep = "\t")

# Must contain at least CHR and BP
required_cols <- c("CHR", "BP")
if (!all(required_cols %in% names(data))) {
  stop("Input file must contain columns: CHR, BP", call. = FALSE)
}

# ======== Create SNP density plot ========
output_file_prefix <- paste0(output_prefix, "_density")

cat("Creating SNP density plot:", output_file_prefix, "\n")

windinfo <- CMplot(
  data,
  plot.type = "d",                      # Density plot
  bin.size = 1e6,                       # 1 Mb window
  chr.den.col = c("darkgreen", "yellow", "red"),  # Color gradient by density
  file = "jpg",                         # Output format
  file.name = output_file_prefix,       # Output file name
  dpi = 300,
  main = paste0("SNP Density Plot - ", output_prefix),
  file.output = TRUE,                   # Save to file
  verbose = TRUE,
  width = 9,
  height = 6
)

cat("Done! Plot saved to:", paste0(output_file_prefix, ".jpg"), "\n")
