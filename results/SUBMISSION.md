# Reproducibility index

Mapping between each element of the manuscript and the file in this repository that
produced it. Figures are PNG; tables are provided as CSV (machine-readable) and, where
submitted, DOCX (typeset).

---

## 1. Main results

### 1.1 Three-tier signal framework (§3.2)

| Manuscript element | Repository path |
|---|---|
| Figure 1 (three-tier AUC) | `results/fig_three_tier_auc.png` |
| Table 3 summary (A1–A4 / B1–B4) | `results/results_c4_master_summary.csv` |
| Cross-cohort transfer AUC with bootstrap 95% CI (A2, A3, B2, B3) | `results/results_c14_external_ci.csv` |

### 1.2 Feature families and sensitivity analyses (§3.3, Table S8)

| Manuscript element | Repository path |
|---|---|
| Figure 2 (label-coupling control, with vs. without *Lactobacillus*) | `results/fig_label_coupling.png` |
| Figure 3 (feature-family and timepoint contrast) | `results/fig_feature_family.png` |
| Figure 4 (richness forest plot) | `results/fig_richness_forest.png` |
| Cross-family / CLR / feature-frequency sensitivity | `results/results_c20_*.csv` |
| Two-way single-cohort transfer matrix (12 cells; not leave-one-cohort-out) | `results/results_c19_loco_matrix.csv` |
| Cross-cohort directional consistency (ω) | `results/results_c8_direction_consistency.csv` |

### 1.3 Clinical landmarks and power (§3.4)

| Manuscript element | Repository path |
|---|---|
| Figure 5 (nine taxa × three cohorts, outcome correlation heat map) | `results/fig_c18_outcome_heatmap.png` |
| D7 dynamics and transient-clearance statistics (Table S5) | `results/table_S5_d7_dynamic.csv`, `results/fig_S5_d7_dynamic.png` |
| PPV for recovery at D30 among D7-responders (0.350, 7/20) | `results/table_S5_d7_dynamic.csv` |
| Targeted prediction of key clinical taxa (Table S7) | `results/results_c18_univariate.csv`, `results/results_c18_combo.csv` |
| Cross-cohort pooled ρ (metafor) | `results/c6_metafor_pooled.csv`, `results/fig_c6_meta_forest.png` |
| Minimum detectable AUC (MDE) | `results/results_mde_power.csv` (produced by `scripts/mde_power_analysis.py`) |

### 1.4 Descriptive quantities cited in Methods and Results

| Manuscript element | Repository path |
|---|---|
| Baseline CST IV counts per cohort | `results/table_S2_three_tier_framework.csv` |
| Final metric set (A1–A4) | `results/results_c9_final_metrics.csv` |
| 41 / 65 / 83-genus feature-space matrices | `results/results_c8_richness_targets.csv` |
| Cross-family variance decomposition | `results/results_c8_variance_decomp.csv` |

**Primary outputs and sensitivity outputs.** Each reported quantity has one designated
primary output file, listed in the tables above; any additional file covering the same
quantity is a sensitivity analysis under a different resampling protocol and is labelled
as such here.

- Within-cohort CV AUC: primary = `results_c9_final_metrics.csv` (20-seed repeated 5-fold CV).
- Cross-cohort transfer AUC: primary = `results_c14_external_ci.csv` (33-seed bootstrap,
  with 95% CI); this is the file the manuscript reports.
- All reported AUC values are traceable to one of these two files; no value in the manuscript
  is carried over from an intermediate working table.
- Remaining small differences between files reflect the different resampling protocols
  (20 vs. 33 resamples) rather than different underlying estimates; the direction and
  significance of every reported comparison are identical across protocols.

---

## 2. Supplementary material (S1–S9)

### 2.1 Tables S1–S9

| Element | Repository path |
|---|---|
| Table S1 | `results/table_S1_dada2_params.{csv,docx}` |
| Table S2 | `results/table_S2_three_tier_framework.{csv,docx}` |
| Table S3 | `results/table_S3_feature_list.{csv,docx}` |
| Table S4 | `results/table_S4_resolution_annotation.{csv,docx}` |
| Table S5 | `results/table_S5_d7_dynamic.{csv,docx}` |
| Table S6 | `results/table_S6_alpha_family_matrix.csv` |
| Table S7 | `results/results_c18_univariate.csv`, `results/results_c18_combo.csv` (formatted table in Supplementary Material) |
| Tables S8–S9 | `results/results_c20_model_family.csv`, `results/results_c20_clr.csv`, `results/results_c20_clr_pseudo.csv`, `results/results_c20_feat_freq.csv` (formatted tables in Supplementary Material) |

### 2.2 Figures S1–S5

| Element | Repository path |
|---|---|
| Figure S1 | `results/fig_S1_quality_profiles.png` |
| Figure S2 | `results/fig_S2_three_tier.png` |
| Figure S3 | `results/fig_S3_feature_distribution.png` |
| Figure S4 | `results/fig_S4_resolution_contrast.png` |
| Figure S5 | `results/fig_S5_d7_dynamic.png` |

### 2.3 Code S1–S2

| Element | Repository path |
|---|---|
| Code S1 — DADA2 pipeline (R) | `results/code_S1_dada2_pipeline.R` |
| Code S2 — LASSO pipeline (Python) | `results/code_S2_lasso_pipeline.py` |

---

## 3. Sequencing quality profiles

Read-quality verification of the raw FASTQ files, one file per cohort with public data:

| Cohort | Repository path |
|---|---|
| Xiao (PRJNA302078) | `results/PRJNA302078_quality.pdf`, `results/PRJNA302078_quality_py.png` |
| Plummer / SToPBV (PRJNA592384) | `results/PRJNA592384_quality_py.png` |
| Gustin (PRJNA691964) | `results/PRJNA691964_quality_py.png` |

The Mollin cohort has no raw-read quality profile because its raw data are not public;
it was reconstructed from the published supplementary tables and is used only as a
pipeline-heterogeneity reference.

---

## 4. Software versions

Python 3.12 (scikit-learn, xgboost, pandas, numpy); R 4.6.1 (DADA2 ≥ 1.30, phyloseq,
metafor, pROC); SILVA v138.1 reference database.
