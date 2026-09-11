#!/usr/bin/env python3
# ============================================================
# Supplementary Code S2
# LASSO logistic regression modelling for BV recovery prediction
# Cohorts: Xiao (PRJNA302078), Mollin (custom label), Gustin (PRJNA691964)
# Depends: scikit-learn >= 1.3, pandas, numpy
# ============================================================

import os, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

# ---- 1. Load pre-built feature matrices (ASCII paths) ----
XFEAT = "D:/lineB_data/feat/X_features.csv"     # sample x features (genus + alpha)
Y_LABELS = "D:/lineB_data/feat/y_labels.csv"   # cohort, subject, base_sample, outcome_sample, base_cst4, outcome
SEED_LIST = list(range(20))                      # 20-seed internal CV

# ---- 2. Outcome definition ----
# Outcome = 1 if follow-up Lactobacillus >= 50% (linear scale; never use log1p space)
def make_outcome(g_rel: pd.DataFrame) -> np.ndarray:
    return (g_rel["Lactobacillus"] >= 0.5).astype(int).values

# ---- 3. Cohort × subset helper ----
def subset(X: pd.DataFrame, y: pd.DataFrame, cohort: str, cst4_only: bool = False):
    sel = (X["cohort"] == cohort)
    if cst4_only:
        sel &= (y["base_cst4"] == 1)
    return X[sel], y[sel]

# ---- 4. 20-seed internal CV ----
def lasso_cv(X: pd.DataFrame, y: np.ndarray, n_splits: int = 5, seeds: list = SEED_LIST) -> dict:
    aucs = []
    for seed in seeds:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        oof = np.zeros(len(y))
        for tr, te in cv.split(X, y):
            scaler = StandardScaler().fit(X.iloc[tr])
            Xtr = scaler.transform(X.iloc[tr])
            Xte = scaler.transform(X.iloc[te])
            clf = LogisticRegressionCV(Cs=10, cv=3, penalty="l1",
                                       solver="liblinear", scoring="roc_auc",
                                       max_iter=5000, random_state=seed)
            clf.fit(Xtr, y[tr])
            oof[te] = clf.predict_proba(Xte)[:, 1]
        aucs.append(roc_auc_score(y, oof))
    return {"auc_mean": np.mean(aucs), "auc_sd": np.std(aucs), "aucs": aucs}

# ---- 5. Three-tier example: Xiao within-cohort CST-IV ----
def run_xiao_cst4():
    X_all = pd.read_csv(XFEAT)
    y_all = pd.read_csv(Y_LABELS)
    X, y = subset(X_all, y_all, cohort="Xiao", cst4_only=True)
    feat_cols = [c for c in X.columns if c not in ("cohort", "subject", "base_sample", "outcome_sample", "out")]
    Xfeat = X[feat_cols]
    y_bin = make_outcome(pd.read_csv("D:/lineB_data/feat/PRJNA302078_genus_rel.csv"))
    return lasso_cv(Xfeat, y_bin)

# ---- 6. Cross-platform transfer example: Xiao → Mollin ----
def run_xiao_to_mollin():
    X_all = pd.read_csv(XFEAT)
    y_all = pd.read_csv(Y_LABELS)
    Xtr_all, ytr = subset(X_all, y_all, cohort="Xiao", cst4_only=True)
    Xte_all, yte = subset(X_all, y_all, cohort="Mollin", cst4_only=True)
    feat_cols = [c for c in Xtr_all.columns if c not in ("cohort", "subject", "base_sample", "outcome_sample", "out")]
    Xtr, Xte = Xtr_all[feat_cols].align(Xte_all[feat_cols], join="inner", axis=1, fill_value=0)
    aucs = []
    for seed in SEED_LIST:
        scaler = StandardScaler().fit(Xtr)
        clf = LogisticRegressionCV(Cs=10, cv=3, penalty="l1",
                                   solver="liblinear", scoring="roc_auc",
                                   max_iter=5000, random_state=seed)
        clf.fit(scaler.transform(Xtr), ytr)
        y_score = clf.predict_proba(scaler.transform(Xte))[:, 1]
        aucs.append(roc_auc_score(yte, y_score))
    return {"auc_mean": np.mean(aucs), "auc_sd": np.std(aucs)}

# ---- 7. Sensitivity: label-coupling-removal (Section 3.2) ----
def run_no_lacto_xiao():
    """Whole-cohort Xiao with Lactobacillus features removed"""
    X_all = pd.read_csv(XFEAT)
    y_all = pd.read_csv(Y_LABELS)
    X, y = subset(X_all, y_all, cohort="Xiao", cst4_only=False)
    feat_cols = [c for c in X.columns if c not in ("cohort", "subject", "base_sample", "outcome_sample", "out") and "Lactobacillus" not in c]
    Xfeat = X[feat_cols]
    y_bin = make_outcome(pd.read_csv("D:/lineB_data/feat/PRJNA302078_genus_rel.csv"))
    return lasso_cv(Xfeat, y_bin)

if __name__ == "__main__":
    print("Xiao CST-IV within-cohort:", run_xiao_cst4())
    print("Xiao → Mollin CST-IV cross-platform:", run_xiao_to_mollin())
    print("Xiao whole-cohort no-Lactobacillus:", run_no_lacto_xiao())
