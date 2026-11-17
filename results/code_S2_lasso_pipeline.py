#!/usr/bin/env python3
# ============================================================
# Supplementary Code S2 (corrected)
# LASSO logistic regression for BV microbiological-recovery prediction
#
# STATUS: illustrative implementation of the modelling steps described in
#         Methods. It documents the estimator, the feature-space guard, and the
#         resampling scheme. It is NOT a turnkey script that reproduces the
#         manuscript tables: that requires the derived feature matrices, which
#         are available from the corresponding author on reasonable request.
#
# Reference implementation for the quantities reported in the manuscript:
#   scripts/four_cohort_matrix.py  (feature selection by whitelist, see §1)
#
# Depends: scikit-learn >= 1.3, pandas, numpy
# ============================================================

import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

# ---- 0. Paths ------------------------------------------------------------
# Set LINEB_DATA to the directory holding feat/ and meta/. No absolute paths
# are hard-coded so that the script is portable.
DATA_DIR = os.environ.get("LINEB_DATA", "lineB_data")
XFEAT = os.path.join(DATA_DIR, "feat", "X_features.csv")
Y_LABELS = os.path.join(DATA_DIR, "feat", "y_labels.csv")
FEAT_META = os.path.join(DATA_DIR, "feat", "feature_meta.csv")
XIAO_REL = os.path.join(DATA_DIR, "feat", "PRJNA302078_genus_rel.csv")

SEED_LIST = list(range(20))   # 20-seed repeated 5-fold cross-validation

# ---- 1. Feature-space guard ---------------------------------------------
# Every column that encodes, or can be used to reconstruct, the outcome must be
# excluded from the feature matrix. Listed explicitly rather than by prefix
# matching, and asserted after selection, because a single string-vs-tuple slip
# (e.g. `c not in ("out")`, which tests for a substring of "out") silently lets
# `outcome`, `base_cst4`, `outcome_cst4` and `base_LB` into the model.
LABEL_COLS = {
    "cohort", "subject", "base_sample", "outcome_sample",
    "outcome", "outcome_cst4", "base_cst4", "base_LB", "base_BV",
}


def feature_cols(X: pd.DataFrame, feat_meta: pd.DataFrame | None = None) -> list:
    """Return the feature columns, preferring an explicit whitelist.

    The whitelist route (feat_meta source == genus|alpha) is the one used for
    the manuscript analyses, because it is closed under newly added metadata
    columns. The blacklist route is kept only for inputs without metadata.
    """
    if feat_meta is not None and {"source", "feature"}.issubset(feat_meta.columns):
        allowed = set(feat_meta.loc[feat_meta["source"].isin(["genus", "alpha"]), "feature"])
        cols = [c for c in X.columns if c in allowed]
    else:
        cols = [c for c in X.columns if c not in LABEL_COLS]
        cols = [c for c in cols if not c.startswith(("outcome", "base_"))]

    leaked = sorted(set(cols) & LABEL_COLS)
    if leaked:
        raise ValueError(f"Label leakage into the feature matrix: {leaked}")
    if not cols:
        raise ValueError("No feature columns selected; check feat_meta.")
    return cols


# ---- 2. Outcome definition ----------------------------------------------
# Outcome = 1 if follow-up Lactobacillus relative abundance >= 50%.
# Computed on the LINEAR 0-1 relative-abundance scale, then assigned to whole
# pairs, so that the label vector is aligned with the sample rows by
# construction rather than by position.
def make_outcome(genus_rel: pd.DataFrame, threshold: float = 0.5) -> pd.Series:
    lb = genus_rel["Lactobacillus"].astype(float)
    return (lb >= threshold).astype(int)


# ---- 3. Cohort x subset helper ------------------------------------------
def subset(X: pd.DataFrame, y: pd.DataFrame, cohort: str, cst4_only: bool = False):
    if cohort not in set(y["cohort"]):
        raise ValueError(f"Cohort {cohort!r} not present in the label table; "
                         f"available: {sorted(set(y['cohort']))}")
    sel = (y["cohort"] == cohort).values
    if cst4_only:
        sel &= (y["base_cst4"] == 1).values
    Xs = X.loc[sel].reset_index(drop=True)
    ys = y.loc[sel].reset_index(drop=True)
    if len(Xs) != len(ys):
        raise ValueError(f"Feature/label length mismatch: {len(Xs)} vs {len(ys)}")
    return Xs, ys


# ---- 4. 20-seed repeated 5-fold CV --------------------------------------
def lasso_cv(X: pd.DataFrame, y: np.ndarray, n_splits: int = 5,
             seeds: list = SEED_LIST) -> dict:
    if len(X) != len(y):
        raise ValueError(f"Feature/label length mismatch: {len(X)} vs {len(y)}")
    if len(np.unique(y)) < 2:
        raise ValueError("Only one outcome class present in this subset.")
    aucs = []
    for seed in seeds:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        oof = np.zeros(len(y))
        for tr, te in cv.split(X, y):
            scaler = StandardScaler().fit(X.iloc[tr])
            clf = LogisticRegressionCV(Cs=10, cv=3, penalty="l1", solver="liblinear",
                                       scoring="roc_auc", max_iter=5000, random_state=seed)
            clf.fit(scaler.transform(X.iloc[tr]), y[tr])
            oof[te] = clf.predict_proba(scaler.transform(X.iloc[te]))[:, 1]
        aucs.append(roc_auc_score(y, oof))
    return {"auc_mean": float(np.mean(aucs)), "auc_sd": float(np.std(aucs)),
            "n": int(len(y)), "n_seeds": len(seeds)}


# ---- 5. Example: Xiao, baseline-CST-IV subgroup, within-cohort CV --------
def run_xiao_cst4():
    X_all = pd.read_csv(XFEAT)
    y_all = pd.read_csv(Y_LABELS)
    feat_meta = pd.read_csv(FEAT_META) if os.path.exists(FEAT_META) else None

    X, y = subset(X_all, y_all, cohort="PRJNA302078", cst4_only=True)
    cols = feature_cols(X, feat_meta)

    # Outcome is taken from the label table where available; the recomputation
    # from the genus matrix is a cross-check, not a second definition.
    y_bin = y["outcome"].to_numpy().astype(int)
    if os.path.exists(XIAO_REL):
        g = pd.read_csv(XIAO_REL, index_col=0)
        lb = g["Lactobacillus"].astype(float) if "Lactobacillus" in g.columns else None
        if lb is not None:
            recomputed = (lb.reindex(y["outcome_sample"]).to_numpy() >= 0.5).astype(int)
            if len(recomputed) == len(y_bin) and not np.array_equal(recomputed, y_bin):
                raise ValueError("Outcome mismatch between y_labels and genus_rel.")
    return lasso_cv(X[cols], y_bin)


if __name__ == "__main__":
    print("Xiao CST-IV within-cohort CV:", run_xiao_cst4())
