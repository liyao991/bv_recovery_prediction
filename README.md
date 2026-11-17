# Code and results — pretreatment 16S rRNA and BV microbiological recovery

Analysis code, derived result tables and figures for:

> **Limited predictive value of pretreatment 16S rRNA for BV microbiological recovery:
> a four-cohort outcome–feature de-coupling audit**

The study asks whether the pretreatment vaginal microbiota carries information that
predicts *microbiological recovery* after standard antibiotic treatment of bacterial
vaginosis (BV), and whether previously reported predictive signal survives removal of
**label coupling** — the mechanical overlap between the outcome definition
(recovery = *Lactobacillus* relative abundance ≥ 50%) and the features used to predict it
(*Lactobacillus* abundance itself).

The headline finding is a negative one: within the baseline dysbiotic (CST IV) subgroup,
prediction was near chance both within and across cohorts, and the moderate whole-cohort
performance was largely attributable to the coupling between outcome and feature.

---

## 1. Cohorts

Four published longitudinal 16S rRNA cohorts were re-analysed with a harmonised
microbiological outcome.

| Cohort | Accession | Platform | Region | Paired samples | Baseline CST IV | Outcome |
|---|---|---|---|---|---|---|
| Xiao 2016 | PRJNA302078 | 454 GS Junior | V1–V3 SE | 68 (D0→D30) | 54/68 | *Lactobacillus* ≥ 50% |
| Plummer 2023 (SToPBV) | PRJNA592384 | Illumina MiSeq | V3–V4 PE | 49 (M0→M1) | 43/49 | *Lactobacillus* ≥ 50% |
| Mollin 2022 | — (no raw reads) | — | — | 33 | 30/33 | *Lactobacillus* ≥ 50% |
| Gustin 2022 (CONRAD BV trial) | PRJNA691964 | Illumina MiSeq | V4 2×150 | 33 (V1→V2→V3; 101 runs) | 33/33 | *Lactobacillus* ≥ 50% |

Notes:

- **Mollin 2022** raw sequence data were not available; this cohort was reconstructed
  from the original publication's supplementary tables and serves as a
  pipeline-heterogeneity reference only.
- **Excluded cohorts**: StepUp-1 (PRJNA398590) and StepUp-2 (PRJNA735440) were excluded
  because they involve partner treatment / male samples; NuSwab (PRJNA1171405) lacked
  outcome labels; water controls were excluded.
- The Gustin cohort is a dedicated treatment cohort in which all baseline samples are
  CST IV, so its "whole cohort" and its "CST IV subgroup" are the same set — this is used
  as a key control for separating label coupling from genuine within-subgroup signal.

---

## 2. Repository layout

```
scripts/    analysis code (Python, R)
results/    derived result tables (CSV), figures (PNG), supplementary tables (CSV/DOCX)
SUBMISSION.md   index mapping each paper element (main figures, Tables S1–S9, Code S1–S2)
                to the file that produced it
```

Pipeline in order:

1. `scripts/download_sra_cohorts.py` — fetch FASTQ from ENA/SRA for the three cohorts with
   public raw data
2. `scripts/dada2_pipeline_template.R` — DADA2 denoising and SILVA v138.1 assignment
   (`code_S1_dada2_pipeline.R` in `results/` is the copy deposited as Supplementary Code S1;
   see the note below on its status)
3. `scripts/export_genus_matrix.R`, `export_species_matrix.R`, `export_alpha_family.R` —
   build genus/species relative-abundance matrices and α-diversity features
4. `scripts/train_main.py`, `four_cohort_matrix.py`, `c20_sensitivity.py`, `c19_loco_matrix.py`
   — modelling, cross-cohort transfer and sensitivity analyses
5. `scripts/c6_metafor.R`, `scripts/mde_power_analysis.py` — cross-cohort correlation
   pooling and minimum detectable effect (MDE) calculation

### Status of the deposited code

Two levels of reproducibility are distinguished here, deliberately:

| Code | Status | What that means |
|---|---|---|
| `results/code_S2_lasso_pipeline.py` | **Executable** | Self-contained given the derived feature matrix; includes a hard assertion that outcome/label columns never enter the feature set |
| `results/code_S1_dada2_pipeline.R` | **Illustrative** | Documents the denoising parameters reported in Supplementary Table S1, but has not been re-executed against raw reads for this deposit. Raw FASTQ for the public cohorts should be fetched with `scripts/download_sra_cohorts.py` and re-processed; exact per-run truncation lengths for `PRJNA691964` depend on read quality and should be re-validated |
| `scripts/*.py`, `scripts/*.R` | **Analysis scripts as run** | Steps 3–5 of the pipeline. Some were written iteratively during the study; where two implementations of the same quantity exist, the manuscript reports the one named in `SUBMISSION.md` |
| `results/results_c4_master_summary.csv` | **Deposited summary (full precision)** | Assembled from the modelling runs at full precision; `results_c9_final_metrics.csv` records the same nine quantities at three decimals for cross-checking. Re-deriving the fourth decimal requires the saved prediction arrays, which are available from the corresponding author |

We state this explicitly rather than labelling everything "frozen code", because a script that
has not been re-run is not a reproduction guarantee. No derived numbers reported in the paper
are affected by this distinction.

---

## 3. Path placeholders

Scripts use placeholders instead of absolute paths. Replace them with your local paths
before running:

| Placeholder | Meaning |
|---|---|
| `<DATA_DIR>` | root of raw/derived data (FASTQ, SILVA reference, feature matrices) |
| `<OUT_DIR>` | computation output directory (DADA2 products, modelling CSVs) |
| `<REPO_DIR>` | root of this repository |
| `<USER_HOME>` | user home directory |

Raw FASTQ files and large intermediate objects (ASV tables, phyloseq RData) are **not**
included; they are regenerated by steps 1–3 above.

---

## 4. Reproducibility notes

- Random seed 42 was used for permutation tests, rarefaction and bootstrap; primary
  results are reported as the seed-mean over 20 random seeds, with cross-validation
  splits independently stratified per seed.
- λ (the LASSO penalty) was re-selected inside each training split by cross-validation,
  with the test fold excluded — no information from the test fold enters model selection.
- Minimum detectable AUC was computed by solving
  `AUC − 0.5 = (z₁₋α/₂ + z₁₋β) · SE_HM(AUC)` with the Hanley–McNeil standard error
  evaluated at the AUC being tested (α = 0.05 two-sided, power 0.80);
  `scripts/mde_power_analysis.py` reproduces these values (output:
  `results/results_mde_power.csv`).
- Internal cross-validation AUC is the mean ± SD across 20×5 CV folds (a dispersion across
  folds, not a confidence interval); external transfer reports a multi-seed bootstrap 95% CI.
- Small differences in the third decimal between result files reflect different analysis
  protocols (e.g. matrix-protocol versus re-check protocol); `SUBMISSION.md` states which
  file corresponds to which reported value.

---

## 5. Data availability

Raw sequencing data are publicly available from the NCBI Sequence Read Archive under
BioProject accessions **PRJNA302078**, **PRJNA592384** and **PRJNA691964**. The Mollin
cohort was reconstructed from the published supplementary tables of the original study.
Derived result tables and figures supporting the findings are provided in `results/`.

## 6. License

Released under the MIT License — see `LICENSE`.

## 7. Contact

Correspondence: Li Yao, Department of Clinical Laboratory, Fujian Provincial Hospital,
Fuzhou, Fujian, China.
