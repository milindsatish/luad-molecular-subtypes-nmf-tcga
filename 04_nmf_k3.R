# ============================================================
# Script: 04_nmf_k3.R
# Purpose: Final NMF factorization at k = 3 and cluster assignment
# ============================================================

# Clean start
rm(list = ls())

# Load libraries
library(NMF)

# Set seed for reproducibility
set.seed(123)

# Paths
processed_data_dir <- "data/processed"

# Load filtered expression matrix
expr_filt <- readRDS(
  file.path(processed_data_dir, "LUAD_expression_filt.rds")
)

# Sanity check: NMF requires non-negative values
stopifnot(all(expr_filt >= 0))

# ------------------------------------------------------------
# Final NMF at k = 3
# ------------------------------------------------------------
nmf_k3 <- nmf(
  expr_filt,
  rank = 3,
  method = "brunet",
  nrun = 30
)

# Save full NMF object
saveRDS(
  nmf_k3,
  file.path(processed_data_dir, "LUAD_nmf_k3.rds")
)

# ------------------------------------------------------------
# Derive hard cluster assignments
# ------------------------------------------------------------
clusters <- apply(coef(nmf_k3), 2, which.max)

# Save clusters
saveRDS(
  clusters,
  file.path(processed_data_dir, "LUAD_clusters_k3.rds")
)

# Quick check
cat("Cluster sizes (k = 3):\n")
print(table(clusters))
