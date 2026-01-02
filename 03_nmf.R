# Project: LUAD Integrative Genomics

# Clear workspace 
rm(list = ls())

# Load libraries
library(NMF)

expr_filt <- readRDS("data/processed/LUAD_expression_filt.rds")

set.seed(123)

nmf_res <- nmf(
  expr_filt,
  rank = 2:6,
  method = "brunet",
  nrun = 5   # reduced to be faster
)

# SAVE
saveRDS(nmf_res, "data/processed/LUAD_nmf_all.rds")

