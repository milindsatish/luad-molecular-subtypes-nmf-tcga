# ============================================================
# Project: TCGA-LUAD Molecular Subtyping with NMF
# Purpose:
#   1. Run final NMF factorization at k = 3
#   2. Derive tumor subtype assignments
#   3. Extract gene programs (W matrix)
#   4. Save all outputs for downstream biology & survival
# ============================================================

library(NMF)
library(tidyverse)

set.seed(123)


# -------------------------------
# 1. Load processed input data
# -------------------------------

processed_data_dir <- "data/processed"

expr_filt <- readRDS(
  file.path(processed_data_dir, "LUAD_expression_filt.rds")
)

surv_filt <- readRDS(
  file.path(processed_data_dir, "LUAD_survival_matched.rds")
)

# Sanity checks
stopifnot(all(expr_filt >= 0))
stopifnot(ncol(expr_filt) == nrow(surv_filt))


# -------------------------------
# 2. Final NMF at k = 3
# -------------------------------

# Rationale:
# k = 3 was selected based on stability, silhouette,
# reconstruction error, and interpretability from rank survey

nmf_k3 <- nmf(
  expr_filt,
  rank   = 3,
  method = "brunet",
  nrun   = 30
)

# Save full NMF object
saveRDS(
  nmf_k3,
  file.path(processed_data_dir, "LUAD_nmf_k3.rds")
)


# -------------------------------
# 3. Derive tumor subtype assignments
# -------------------------------

# H matrix: program activity per tumor
H <- coef(nmf_k3)

# Assign each tumor to its dominant program
clusters <- apply(H, 2, which.max)

# Inspect cluster sizes
cat("NMF cluster sizes (k = 3):\n")
print(table(clusters))

# Save cluster assignments
saveRDS(
  clusters,
  file.path(processed_data_dir, "LUAD_clusters_k3.rds")
)


# -------------------------------
# 4. Extract gene programs (W matrix)
# -------------------------------

# W matrix: gene loadings per program
W <- basis(nmf_k3)

# Rank genes by contribution to each program
top_genes <- apply(W, 2, function(x) {
  names(sort(x, decreasing = TRUE))
})

# Keep top genes for interpretation
top20_genes <- lapply(top_genes, head, 20)
top50_genes <- lapply(top_genes, head, 50)

# Save gene programs
saveRDS(
  top20_genes,
  file.path(processed_data_dir, "LUAD_nmf_top20_genes_k3.rds")
)

saveRDS(
  top50_genes,
  file.path(processed_data_dir, "LUAD_nmf_top50_genes_k3.rds")
)


# -------------------------------
# 5. Prepare data for survival analysis
# -------------------------------

# Add NMF cluster labels to survival table
surv_filt$NMF_cluster <- factor(
  clusters,
  levels = 1:3,
  labels = c("Subtype_1", "Subtype_2", "Subtype_3")
)

# Save annotated survival data
saveRDS(
  surv_filt,
  file.path(processed_data_dir, "LUAD_survival_with_NMF_clusters.rds")
)


# -------------------------------
# 6. Summary for downstream scripts
# -------------------------------

cat("\nSummary:\n")
cat("- Final NMF completed at k = 3\n")
cat("- Tumor subtypes saved\n")
cat("- Gene programs extracted (top 20 / 50 genes)\n")
cat("- Survival table annotated with NMF clusters\n")
cat("\nReady for pathway analysis and survival modeling.\n")
