#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    paste(
      "Usage:",
      "Rscript merge_pheno_pcs.R",
      "<PHENO_INPUT>",
      "<QCED_ASSOC_FAM>",
      "<QCED_KINSHIP_FAM>",
      "<PCA_EIGENVEC>",
      "<OUT_PHENO_MATCHED>",
      "<OUT_PHENO_WITH_PCS>"
    )
  )
}

pheno_file <- args[1]
assoc_fam_file <- args[2]
kinship_fam_file <- args[3]
pca_file <- args[4]
out_pheno_matched <- args[5]
out_pheno_with_pcs <- args[6]

read_table_safe <- function(path, header = FALSE) {
  tryCatch(
    {
      read.table(
        path,
        header = header,
        stringsAsFactors = FALSE,
        sep = "",
        check.names = FALSE,
        comment.char = "",
        quote = ""
      )
    },
    error = function(e) {
      stop(paste("Failed to read file:", path, "|", conditionMessage(e)))
    }
  )
}

if (!file.exists(pheno_file)) {
  stop(paste("Phenotype file not found:", pheno_file))
}
if (!file.exists(assoc_fam_file)) {
  stop(paste("QCed.assoc.fam not found:", assoc_fam_file))
}
if (!file.exists(kinship_fam_file)) {
  stop(paste("QCed.kinship.fam not found:", kinship_fam_file))
}
if (!file.exists(pca_file)) {
  stop(paste("PCA eigenvec file not found:", pca_file))
}

pheno <- read_table_safe(pheno_file, header = TRUE)

if (ncol(pheno) < 2) {
  stop("Phenotype file must have at least two columns: FID and IID")
}

colnames(pheno)[1:2] <- c("FID", "IID")
pheno$FID <- as.character(pheno$FID)
pheno$IID <- as.character(pheno$IID)

assoc_fam <- read_table_safe(assoc_fam_file, header = FALSE)
kinship_fam <- read_table_safe(kinship_fam_file, header = FALSE)
pca <- read_table_safe(pca_file, header = FALSE)

if (ncol(assoc_fam) < 2) {
  stop("QCed.assoc.fam must have at least 2 columns")
}
if (ncol(kinship_fam) < 2) {
  stop("QCed.kinship.fam must have at least 2 columns")
}
if (ncol(pca) < 12) {
  stop("PCA eigenvec file must have at least 12 columns: FID IID PC1-PC10")
}

assoc_ids <- assoc_fam[, 1:2, drop = FALSE]
colnames(assoc_ids) <- c("FID", "IID")
assoc_ids$FID <- as.character(assoc_ids$FID)
assoc_ids$IID <- as.character(assoc_ids$IID)

kinship_ids <- kinship_fam[, 1:2, drop = FALSE]
colnames(kinship_ids) <- c("FID", "IID")
kinship_ids$FID <- as.character(kinship_ids$FID)
kinship_ids$IID <- as.character(kinship_ids$IID)

pca <- pca[, 1:12, drop = FALSE]
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))
pca$FID <- as.character(pca$FID)
pca$IID <- as.character(pca$IID)

make_key <- function(df) {
  paste(df$FID, df$IID, sep = "\t")
}

assoc_keys <- make_key(assoc_ids)
kinship_keys <- make_key(kinship_ids)
pca_keys <- make_key(pca)
pheno_keys <- make_key(pheno)

common_keys <- Reduce(
  intersect,
  list(assoc_keys, kinship_keys, pca_keys, pheno_keys)
)

if (length(common_keys) == 0) {
  stop("No common participants found across QCed.assoc, QCed.kinship, PCA, and phenotype")
}

# Keep final order equal to QCed.assoc
assoc_ids_common <- assoc_ids[assoc_keys %in% common_keys, , drop = FALSE]
assoc_common_keys <- make_key(assoc_ids_common)

pheno_sub <- pheno[pheno_keys %in% common_keys, , drop = FALSE]
pca_sub <- pca[pca_keys %in% common_keys, , drop = FALSE]

pheno_sub$key <- make_key(pheno_sub)
pca_sub$key <- make_key(pca_sub)

pheno_sub <- pheno_sub[match(assoc_common_keys, pheno_sub$key), , drop = FALSE]
pca_sub <- pca_sub[match(assoc_common_keys, pca_sub$key), , drop = FALSE]

if (any(is.na(pheno_sub$key))) {
  stop("Phenotype alignment failed after matching common participants")
}
if (any(is.na(pca_sub$key))) {
  stop("PCA alignment failed after matching common participants")
}

if (!all(pheno_sub$key == assoc_common_keys)) {
  stop("Phenotype order does not match QCed.assoc after alignment")
}
if (!all(pca_sub$key == assoc_common_keys)) {
  stop("PCA order does not match QCed.assoc after alignment")
}

pheno_sub$key <- NULL
pca_sub$key <- NULL

pheno_matched <- pheno_sub
pheno_with_pcs <- cbind(pheno_sub, pca_sub[, paste0("PC", 1:10), drop = FALSE])

write.table(
  pheno_matched,
  file = out_pheno_matched,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  pheno_with_pcs,
  file = out_pheno_with_pcs,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Merge complete\n")
cat("Participants in phenotype input: ", nrow(pheno), "\n", sep = "")
cat("Participants in QCed.assoc.fam: ", nrow(assoc_ids), "\n", sep = "")
cat("Participants in QCed.kinship.fam: ", nrow(kinship_ids), "\n", sep = "")
cat("Participants in PCA eigenvec: ", nrow(pca), "\n", sep = "")
cat("Common participants kept: ", length(common_keys), "\n", sep = "")
cat("Wrote: ", out_pheno_matched, "\n", sep = "")
cat("Wrote: ", out_pheno_with_pcs, "\n", sep = "")