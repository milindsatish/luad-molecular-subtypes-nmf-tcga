# ============================================================
# Project: LUAD Integrative Genomics
# Script: 02_feature_filtering.R
# Purpose: Variance-based gene filtering and initial PCA
# ============================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(tidyverse)

# Define paths
processed_data_dir <- "data/processed"

# Load processed expression matrix
expr_mat <- readRDS(file.path(processed_data_dir, "LUAD_expression_mat.rds"))

cat("Expression matrix dimensions (genes x samples):\n")
print(dim(expr_mat))

# ------------------------------------------------------------
# Compute gene-wise variance
# ------------------------------------------------------------
gene_var <- apply(expr_mat, 1, var)

cat("\nSummary of gene-wise variance:\n")
print(summary(gene_var))

# ------------------------------------------------------------
# Variance-based gene filtering
# Keep top 25% most variable genes
# ------------------------------------------------------------
var_cutoff <- quantile(gene_var, 0.75)

expr_filt <- expr_mat[gene_var > var_cutoff, ]

cat("\nFiltered expression matrix dimensions:\n")
print(dim(expr_filt))

# ------------------------------------------------------------
# PCA for sanity check
# ------------------------------------------------------------
pca <- prcomp(t(expr_filt), scale. = TRUE)

cat("\nVariance explained by first 5 PCs:\n")
print(summary(pca)$importance[2, 1:5])

pca_var <- summary(pca)$importance[2, 1:5]

pca_var_df <- data.frame(
  PC = names(pca_var),
  VarianceExplained = as.numeric(pca_var)
)

write.csv(
  pca_var_df,
  "data/processed/LUAD_PCA_variance_explained.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Save filtered expression matrix
# ------------------------------------------------------------
saveRDS(expr_filt, file.path(processed_data_dir, "LUAD_expression_filt.rds"))

cat("\nFiltered expression matrix saved successfully.\n")
