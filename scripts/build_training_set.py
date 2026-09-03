# -*- coding: utf-8 -*-
"""
build_training_set.py - Join CST labels (sample_meta.csv) with genus-level
feature matrices to build the training X/y for BV recurrence prediction.

  * baseline (D0/M0) genus relative abundances  -> log1p(100*x) transform
  * alpha diversity (shannon / chao1 / n_asv)
  * union of genera across cohorts, missing -> 0
  * outcome = 1 recovered (non-CST4 at follow-up), 0 not recovered

Outputs (ASCII path to avoid R/Windows encoding issues):
  <DATA_DIR>/feat/X_features.csv   (sample x features, log1p + scaled later)
  <DATA_DIR>/feat/y_labels.csv     (cohort, subject, base, outcome)
  <DATA_DIR>/feat/feature_meta.csv (feature source: genus vs alpha)
"""
import os
import re
import numpy as np
import pandas as pd

META = "<DATA_DIR>/meta"
FEAT = "<DATA_DIR>/feat"
OUT_X = os.path.join(FEAT, "X_features.csv")
OUT_Y = os.path.join(FEAT, "y_labels.csv")
OUT_F = os.path.join(FEAT, "feature_meta.csv")

ALPHA_COLS = ["shannon", "chao1", "n_asv"]
BV_GENERA = ["Gardnerella", "Atopobium", "Prevotella", "Sneathia",
             "Megasphaera", "Mobiluncus", "Fannyhessea", "Dialister",
             "Prevotella_7", "Prevotella_9", "Aerococcus", "Fastidiosipila"]

def main():
    meta = pd.read_csv(os.path.join(META, "sample_meta.csv"))
    print(f"pairs: {len(meta)}")

    # ---- load genus + alpha per cohort ----
    genus_dict, alpha_dict = {}, {}
    for proj in ["PRJNA302078", "PRJNA592384"]:
        g = pd.read_csv(os.path.join(FEAT, f"{proj}_genus_rel.csv"), index_col=0)
        a = pd.read_csv(os.path.join(FEAT, f"{proj}_alpha.csv"))
        a = a.set_index(a.columns[0]) if a.columns[0] != "sample" else a.set_index("sample")
        genus_dict[proj] = g
        alpha_dict[proj] = a
        print(f"  {proj}: genus {g.shape}, alpha {a.shape}")

    # ---- assemble X rows: one per pair, using BASELINE features ----
    all_genus_cols = sorted(set().union(*[set(g.columns) for g in genus_dict.values()]))
    # drop non-biological catch-alls kept out of modelling (keep Unassigned out? keep it: it is informative)
    print(f"  union of genera: {len(all_genus_cols)}")

    rows = []
    for r in meta.itertuples():
        proj = r.cohort
        g = genus_dict[proj]
        a = alpha_dict[proj]
        base_s = r.base_sample
        if base_s not in g.index:
            print(f"  [warn] {proj} base sample {base_s} not in genus table")
            continue
        feat = {}
        for col in all_genus_cols:
            feat[col] = g.loc[base_s, col] if col in g.columns else 0.0
        for c in ALPHA_COLS:
            feat[c] = a.loc[base_s, c] if base_s in a.index and c in a.columns else np.nan
        feat["cohort"] = proj
        feat["subject"] = r.subject
        feat["base_sample"] = base_s
        feat["outcome_sample"] = r.outcome_sample
        feat["outcome"] = int(r.outcome)
        feat["base_cst4"] = int(r.base_cst4)
        feat["outcome_cst4"] = int(r.outcome_cst4)
        feat["base_LB"] = float(r.base_LB)
        feat["base_BV"] = float(r.base_BV)
        rows.append(feat)

    X = pd.DataFrame(rows)
    print(f"\nX assembled: {X.shape}")

    # ---- feature transform: log1p(100*x) for abundances ----
    genus_feats = [c for c in all_genus_cols]
    X[genus_feats] = np.log1p(100.0 * X[genus_feats].clip(lower=0))
    # alpha: keep as-is (chao1/n_asv already count-scaled); shannon 0-5 range
    print(f"genus features: {len(genus_feats)}, alpha features: {len(ALPHA_COLS)}")

    # ---- feature meta ----
    feat_meta = pd.DataFrame({
        "feature": genus_feats + ALPHA_COLS,
        "source": ["genus"] * len(genus_feats) + ["alpha"] * len(ALPHA_COLS),
        "is_bv_genera": [c in BV_GENERA for c in genus_feats] + [False] * len(ALPHA_COLS),
    })
    # ---- check for zero-variance / NaN ----
    Xc = X[genus_feats + ALPHA_COLS]
    nan_cols = [c for c in Xc.columns if Xc[c].isna().any()]
    zv_cols = [c for c in Xc.columns if Xc[c].nunique() < 2]
    print(f"  NaN columns: {nan_cols}")
    print(f"  zero-variance columns: {len(zv_cols)}")
    if nan_cols:
        print(f"  -> imputing NaN with column mean")
        Xc = Xc.fillna(Xc.mean())
        X[nan_cols] = Xc[nan_cols]

    # ---- save ----
    X.to_csv(OUT_X, index=False)
    y = X[["cohort", "subject", "base_sample", "outcome_sample",
           "outcome", "base_cst4", "outcome_cst4", "base_LB", "base_BV"]].copy()
    y.to_csv(OUT_Y, index=False)
    feat_meta.to_csv(OUT_F, index=False)

    print(f"\nsaved -> {OUT_X} ({X.shape})")
    print(f"saved -> {OUT_Y} ({y.shape})")
    print(f"saved -> {OUT_F} ({feat_meta.shape})")
    print(f"\noutcome dist: {y['outcome'].value_counts().to_dict()}")
    print(f"cohort x outcome:\n{pd.crosstab(y['cohort'], y['outcome'])}")

    # ---- unresolved alias diagnostics (which 592384 subjects dropped?) ----
    cst = pd.read_csv(os.path.join(META, "sample_cst_all.csv"))
    cst = cst[cst["cohort"] == "PRJNA592384"]
    paired_subjects = set(y.loc[y["cohort"] == "PRJNA592384", "subject"])
    all_subjects = set(cst["subject"])
    dropped = sorted(all_subjects - paired_subjects)
    print(f"\nPRJNA592384 subjects in CST table: {len(all_subjects)}, paired: {len(paired_subjects)}")
    print(f"dropped/unpaired subjects: {dropped}")

if __name__ == "__main__":
    main()
