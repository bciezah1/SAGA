#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)


# Define input and output filenames
input_file <- paste0("temp_snps_dataset.freq.frq")

output_file <- paste0("list_snps_for_grm.txt")

# Load MAF data
maf_data <- read.table(input_file, header=TRUE)

# Define bins for MAF stratification
maf_data <- maf_data %>%
  mutate(maf_bin = cut(MAF, breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.5), include.lowest = TRUE))

# Sample SNPs proportionally from each MAF bin
set.seed(123)  # For reproducibility
downsampled_snps <- maf_data %>%
  group_by(maf_bin) %>%
  sample_frac(0.01)  # Adjust percentage (e.g., 0.05 for 5%)

# Save SNPs to file
write.table(downsampled_snps, output_file, col.names=FALSE, row.names=FALSE, quote=FALSE)

cat(paste("downsampling completed!\n"))
