# ============================================================
# Project: TCGA-LUAD Molecular Subtyping with NMF
# 
# Purpose:
#   1. Run final NMF factorization at k = 3
#   2. Assign tumors to dominant NMF programs
#   3. Extract gene programs (W matrix)
#   4. Perform pathway enrichment using MSigDB Hallmark gene sets
#   5. Save all outputs for downstream interpretation
# ============================================================


# -------------------------------
# 0. Clean start and libraries
# -------------------------------

rm(list = ls())

library(NMF)
library(tidyverse)
library(fgsea)

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

# k = 3 selected based on rank survey stability and interpretability
nmf_k3 <- nmf(
  expr_filt,
  rank   = 3,
  method = "brunet",
  nrun   = 3
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
ranked_genes <- apply(W, 2, function(x) {
  sort(x, decreasing = TRUE)
})

# Save top genes for interpretation
top20_genes <- lapply(ranked_genes, function(x) names(x)[1:20])
top50_genes <- lapply(ranked_genes, function(x) names(x)[1:50])

saveRDS(
  top20_genes,
  file.path(processed_data_dir, "LUAD_nmf_top20_genes_k3.rds")
)

saveRDS(
  top50_genes,
  file.path(processed_data_dir, "LUAD_nmf_top50_genes_k3.rds")
)


# -------------------------------
# 5. Pathway enrichment (Hallmark)
# -------------------------------

# MSigDB Hallmark gene sets (stored in project root)
hallmark_gmt <- "h.all.v2025.1.Hs.symbols.gmt.txt"

pathways <- fgsea::gmtPathways(hallmark_gmt)

enrichment_results <- list()

for (k in seq_len(ncol(W))) {
  
  # Extract gene weights for program k
  gene_stats <- W[, k]
  
  # Name by gene symbols
  names(gene_stats) <- rownames(W)
  
  # Rank genes (required for fgsea)
  gene_stats <- sort(gene_stats, decreasing = TRUE)
  
  # Run fgsea
  fgsea_res <- fgsea(
    pathways = pathways,
    stats    = gene_stats
  )
  
  # Order by adjusted p-value
  fgsea_res <- fgsea_res %>%
    arrange(padj)
  
  enrichment_results[[paste0("Program_", k)]] <- fgsea_res
}
# Save enrichment results
saveRDS(
  enrichment_results,
  file.path(processed_data_dir, "LUAD_nmf_hallmark_enrichment_k3.rds")
)


# -------------------------------
# 6. Annotate survival data
# -------------------------------

surv_filt$NMF_cluster <- factor(
  clusters,
  levels = 1:3,
  labels = c("Subtype_1", "Subtype_2", "Subtype_3")
)

saveRDS(
  surv_filt,
  file.path(processed_data_dir, "LUAD_survival_with_NMF_clusters.rds")
)

# ------------------------------------------------------------
# Program Activity Heatmap (H matrix)
# ------------------------------------------------------------

png("figures/fig1_program_activity_heatmap.png", width = 1200, height = 800)

aheatmap(
  coef(nmf_k3),
  scale  = "row",
  labRow = paste0("Program ", 1:3),
  labCol = NA
)

dev.off()

# -------------------------------
# 7. Summary
# -------------------------------

cat("\nPipeline complete:\n")
cat("- Final NMF model saved\n")
cat("- Tumor subtypes assigned\n")
cat("- Gene programs extracted\n")
cat("- Hallmark pathway enrichment completed\n")
cat("- Survival data annotated\n")
