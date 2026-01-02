# luad-molecular-subtypes-nmf-tcga
Identification and biological interpretation of transcriptional programs in TCGA lung adenocarcinoma using non-negative matrix factorization, with pathway enrichment and survival analysis to assess translational relevance.
# TCGA LUAD Transcriptional Programs and Clinical Associations

## Key Takeaways

- Identified **three biologically interpretable transcriptional programs** in TCGA lung adenocarcinoma using unsupervised non-negative matrix factorization (NMF).
- Programs correspond to **invasive/EMT-like**, **proliferative**, and **immune-associated** tumor states.
- Pathway enrichment revealed strong concordance with known LUAD biology, including EMT, KRAS signaling, cell-cycle regulation, and immune response pathways.
- Survival analyses showed **modest, suggestive associations** between program activity and overall survival, consistent with the multifactorial nature of clinical outcomes in TCGA cohorts.
- This project demonstrates a **reproducible, end-to-end translational genomics workflow** suitable for pharma research or PhD-level study.

---

## Overview

This project investigates transcriptional heterogeneity in lung adenocarcinoma (LUAD) using publicly available TCGA RNA-seq and clinical data. The primary goal is to identify latent gene expression programs using Non-negative Matrix Factorization (NMF), interpret their biological meaning via pathway enrichment, and evaluate their potential clinical relevance using survival analyses.

This work is designed as a **pharma- and PhD-ready translational genomics project**, emphasizing reproducibility, interpretability, and honest statistical assessment.

---

## Data Sources

- **Gene expression:** TCGA-LUAD RNA-seq (Illumina HiSeq; UCSC Xena Hub)
- **Clinical survival data:** TCGA-CDR curated survival endpoints
  - Overall Survival (OS)
  - Progression-Free Interval (PFI)
  - Disease-Specific Survival (DSS)
  - Disease-Free Interval (DFI)
- **Pathway annotations:** MSigDB Hallmark gene sets (v2025)

All data were downloaded from UCSC Xena and processed locally in R.

---

## Methods at a Glance

- Downloaded TCGA-LUAD RNA-seq and curated survival data from UCSC Xena
- Matched expression and clinical data by TCGA sample barcodes
- Filtered genes by variance (top 25%) to focus on informative signal
- Applied Non-negative Matrix Factorization (NMF) for unsupervised program discovery
- Interpreted programs using MSigDB Hallmark pathway enrichment (FGSEA)
- Assessed clinical relevance using Kaplan–Meier and Cox proportional hazards models

---

## Preprocessing and Quality Control

- Expression matrix aligned to curated survival data using TCGA sample barcodes
- Final matched cohort: ~516 LUAD tumors
- Gene-wise variance computed across samples
- Top 25% most variable genes retained (≈5,100 genes)
- PCA performed as a sanity check to assess global structure and variance distribution

---

## Unsupervised Program Discovery (NMF)

### Method

Non-negative Matrix Factorization was applied to the filtered expression matrix:

- **W (genes × programs):** defines gene expression programs  
- **H (programs × tumors):** defines program activity per tumor  

NMF was run using the Brunet algorithm across ranks *k* = 2–6 with multiple random initializations. Model selection was guided by:

- Cophenetic correlation
- Consensus clustering stability
- Silhouette scores
- Interpretability of resulting programs

### Selected Model

- **k = 3 programs** provided the best balance of stability and interpretability
- Tumors were assigned to subtypes based on their dominant program activity

---

## Biological Interpretation

### Pathway Enrichment

Genes were ranked by program-specific weights (W matrix) and analyzed using FGSEA against MSigDB Hallmark pathways.

**Key findings:**

- **Program 1:** Enriched for epithelial–mesenchymal transition (EMT), KRAS signaling, coagulation, and invasive pathways
- **Program 2:** Enriched for cell-cycle programs, including E2F targets and G2M checkpoint
- **Program 3:** Enriched for immune-related and inflammatory pathways, including interferon signaling

These programs correspond to biologically interpretable tumor states:

- Invasive / EMT-like  
- Proliferative  
- Immune-influenced  

---

## Survival Analysis

### Figures

- **Figure 1:** Kaplan–Meier survival curves comparing three NMF-defined LUAD subtypes
- **Figure 2:** Kaplan–Meier survival curves comparing Program 1 vs Programs 2+3
- **Figure 3:** Forest plot summarizing Cox proportional hazards model for Program 1 activity

### Kaplan–Meier Analysis

Two Kaplan–Meier analyses were performed using overall survival:

- **Three-group comparison** across NMF-defined subtypes  
  - Trend toward survival differences (log-rank p ≈ 0.08)
  - Separation largely driven by the invasive (Program 1–dominant) subtype

- **Binary comparison:** Program 1 vs Programs 2+3 combined  
  - Trend toward worse survival for Program 1
  - Not statistically significant (p ≈ 0.14)

### Cox Proportional Hazards Models

To preserve continuous information, Cox models were fit using program activity scores (H matrix):

- **Single-program model (Program 1):**
  - Hazard ratio > 1 (directionally consistent with aggressive biology)
  - Not statistically significant

- **Multivariable model (Programs 1–3):**
  - No individual program independently significant
  - Global tests (likelihood ratio, Wald, score) borderline significant (p ≈ 0.06–0.07)
  - Concordance ≈ 0.57, indicating modest prognostic signal

These results suggest that transcriptional programs collectively capture biologically meaningful heterogeneity, though survival associations are modest in this cohort.

---

## Limitations

- Survival effects are limited by censoring and modest effect sizes typical of TCGA cohorts
- Smoking history was explored but not included due to incomplete and inconsistent annotation in publicly available TCGA LUAD datasets
- NMF programs are correlated by design, limiting independent effect estimation in multivariable models

---

## Reproducibility

- All analyses performed in R
- Scripts are modular and numbered by analysis stage
- Intermediate processed data objects saved as RDS files
- Figures generated programmatically and saved to disk

---

## Summary

This project demonstrates a complete unsupervised-to-translational genomics workflow:

- Discovery of latent transcriptional programs
- Biological interpretation via pathway enrichment
- Careful evaluation of clinical relevance using survival analysis
- Transparent reporting of statistical limitations

The results highlight biologically distinct LUAD tumor states and illustrate how unsupervised learning can inform cancer biology even when clinical associations are modest.
