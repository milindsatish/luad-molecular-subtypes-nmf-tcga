# ============================================================
# Survival Analysis of NMF-defined LUAD Programs
# ============================================================

library(survival)
library(survminer)

# ------------------------------------------------------------
# Create survival object
# ------------------------------------------------------------
surv_obj <- Surv(
  time  = as.numeric(surv$OS.time),
  event = as.numeric(surv$OS)
)

# ============================================================
# 1) Kaplan–Meier: Three NMF subtypes
# ============================================================

km_fit_3group <- survfit(
  surv_obj ~ NMF_cluster,
  data = surv
)

km_plot_3group <- ggsurvplot(
  km_fit_3group,
  data       = surv,
  pval       = TRUE,
  risk.table = TRUE,
  conf.int   = FALSE,
  palette    = "Dark2",
  title      = "Overall Survival by NMF-defined LUAD Subtype",
  xlab       = "Time (days)",
  ylab       = "Overall survival probability"
)

dir.create("figures", showWarnings = FALSE)

ggsave(
  filename = "figures/LUAD_KM_3group_NMF_subtypes.png",
  plot     = km_plot_3group$plot,
  width    = 7,
  height   = 6,
  dpi      = 300
)

print(km_plot_3group)

# ============================================================
# 2) Kaplan–Meier: Program 1 vs Programs 2 + 3
# ============================================================

surv$Program1_vs_Others <- ifelse(
  surv$NMF_cluster == "Subtype_1",
  "Program_1",
  "Program_2_3"
)

surv$Program1_vs_Others <- factor(
  surv$Program1_vs_Others,
  levels = c("Program_2_3", "Program_1")
)

km_fit_2group <- survfit(
  surv_obj ~ Program1_vs_Others,
  data = surv
)

km_plot_2group <- ggsurvplot(
  km_fit_2group,
  data       = surv,
  pval       = TRUE,
  risk.table = TRUE,
  conf.int   = FALSE,
  palette    = c("#2C7BB6", "#D7191C"),
  title      = "Overall Survival: Program 1 vs Programs 2+3",
  xlab       = "Time (days)",
  ylab       = "Overall survival probability"
)

ggsave(
  filename = "figures/LUAD_KM_Program1_vs_Program23.png",
  plot     = km_plot_2group$plot,
  width    = 7,
  height   = 6,
  dpi      = 300
)

print(km_plot_2group)

# ============================================================
# 3) Cox Proportional Hazards Models using NMF H scores
# ============================================================

# Load final NMF model
nmf_k3 <- readRDS("data/processed/LUAD_nmf_k3.rds")

# Extract H matrix (program activity per tumor)
H <- coef(nmf_k3)  # programs x samples

# Ensure sample alignment
common_ids <- intersect(colnames(H), surv$sample_short)
H <- H[, common_ids]
surv_cox <- surv[match(common_ids, surv$sample_short), ]

# ------------------------------------------------------------
# Cox model: Program 1 activity only
# ------------------------------------------------------------
cox_p1 <- coxph(
  Surv(OS.time, OS) ~ H[1, ],
  data = surv_cox
)

summary(cox_p1)

# ------------------------------------------------------------
# Cox model: All programs jointly
# ------------------------------------------------------------
cox_all <- coxph(
  Surv(OS.time, OS) ~ H[1, ] + H[2, ] + H[3, ],
  data = surv_cox
)

summary(cox_all)

