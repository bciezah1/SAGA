#!/usr/bin/env Rscript

suppressMessages({
  library(CMplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript create_qq.plot.R <input_file> [output_dir]", call. = FALSE)
}

input_file <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], getwd())

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

prefix <- tools::file_path_sans_ext(basename(input_file))

cat("Reading data from:", input_file, "\n")
cat("Writing plots to:", output_dir, "\n")

data <- read.table(input_file, header = TRUE, sep = "\t")

if (!"PVAL" %in% names(data)) {
  stop("Input file must contain column: PVAL", call. = FALSE)
}

data <- data[!is.na(data$PVAL) & data$PVAL > 0 & data$PVAL <= 1, ]

chisq <- qchisq(1 - data$PVAL, 1)
lambda_gc <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, 1), 3)
cat("Lambda GC =", lambda_gc, "\n")

output_file <- file.path(output_dir, paste0(prefix, "_qq_lambda.jpg"))
cat("Creating QQ plot:", output_file, "\n")
data <- data[, c("SNP", "CHR", "BP", "PVAL")]
colnames(data) <- c("SNP", "CHR", "BP", "P")

jpeg(output_file, width = 1500, height = 1500, res = 300)

CMplot(
  data,
  plot.type = "q",
  col = c("darkorchid"),
  box = FALSE,
  file.output = FALSE,
  conf.int = TRUE,
  conf.int.col = NULL,
  threshold.col = "red",
  threshold.lty = 2,
  verbose = TRUE,
  main = "QQ Plot of GWAS Results"
)

usr <- par("usr")
text(
  x = usr[1] + 0.1 * (usr[2] - usr[1]),
  y = usr[4] - 0.1 * (usr[4] - usr[3]),
  labels = bquote(lambda[GC] == .(lambda_gc)),
  cex = 1.4,
  col = "black",
  adj = 0
)

dev.off()

cat("Done! QQ plot with lambda saved to:", output_file, "\n")