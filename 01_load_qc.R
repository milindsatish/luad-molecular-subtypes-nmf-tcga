# ============================================================
# Project: LUAD Integrative Genomics
# Purpose: Load TCGA-LUAD expression and survival data,
#          harmonize sample IDs, and save clean datasets
# ============================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(data.table)
library(tidyverse)

# Define paths
raw_data_dir <- "data/raw"
processed_data_dir <- "data/processed"

expr_file <- file.path(raw_data_dir, "TCGA.LUAD.sampleMap-HiSeqV2.gz")
surv_file <- file.path(raw_data_dir, "LUAD_survival.txt")

# ------------------------------------------------------------
# Load expression data
# ------------------------------------------------------------
expr <- fread(expr_file)

cat("Expression data dimensions:\n")
print(dim(expr))

# First column = gene names
expr_genes <- expr[[1]]
expr_sample_ids <- colnames(expr)[-1]

cat("\nExample expression sample IDs:\n")
print(head(expr_sample_ids))

# ------------------------------------------------------------
# Load survival data
# ------------------------------------------------------------
surv <- fread(surv_file)

cat("\nSurvival data dimensions:\n")
print(dim(surv))

cat("\nSurvival data column names:\n")
print(colnames(surv))

# ------------------------------------------------------------
# Harmonize TCGA sample IDs
# ------------------------------------------------------------

expr_ids_short <- substr(expr_sample_ids, 1, 12)
surv_ids_short <- substr(surv$sample, 1, 12)

common_ids <- intersect(expr_ids_short, surv_ids_short)

cat("\nNumber of matched samples:\n")
print(length(common_ids))

# ------------------------------------------------------------
# Create clean expression matrix
# ------------------------------------------------------------
expr_mat <- as.matrix(expr[, -1])
rownames(expr_mat) <- expr_genes
colnames(expr_mat) <- expr_ids_short

expr_mat <- expr_mat[, common_ids]

cat("\nFiltered expression matrix dimensions:\n")
print(dim(expr_mat))

# ------------------------------------------------------------
# Create clean survival table
# ------------------------------------------------------------
surv$sample_short <- substr(surv$sample, 1, 12)
surv_filt <- surv[surv$sample_short %in% common_ids, ]

# Reorder survival rows to match expression columns
surv_filt <- surv_filt[match(common_ids, surv_filt$sample_short), ]

cat("\nFiltered survival data dimensions:\n")
print(dim(surv_filt))

# ------------------------------------------------------------
# Save processed objects
# ------------------------------------------------------------
saveRDS(expr_mat, file.path(processed_data_dir, "LUAD_expression_mat.rds"))
saveRDS(surv_filt, file.path(processed_data_dir, "LUAD_survival_matched.rds"))

cat("\nProcessed data saved successfully.\n")
