# -*- coding: utf-8 -*-
"""
xiao_holdout.py - Repeated holdout validation within Xiao (the only robust signal).
68 Xiao subjects -> stratified 50 train / 18 test, repeated over 50 random seeds.
Models: LASSO-LR / RF / XGB (same grids as train_main.py).
Reports test AUC/PR-AUC/Brier distribution vs within-CV estimate.
"""
import os
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split, GridSearchCV
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from xgboost import XGBClassifier

RNG = 42
FEAT = "<DATA_DIR>/feat"
OUT = "<OUT_DIR>/modeling"
N_SEEDS = 50
TEST_SIZE = 18
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "xiao_holdout_log.txt"), "w", buffering=1)

def log(*args):
    msg = " ".join(str(a) for a in args)
    print(msg, flush=True)
    LOG.write(msg + "\n")

# ---------------- data (Xiao only) ----------------
X = pd.read_csv(os.path.join(FEAT, "X_features.csv"))
y = pd.read_csv(os.path.join(FEAT, "y_labels.csv"))
feat = pd.read_csv(os.path.join(FEAT, "feature_meta.csv"))
genus_cols = feat.loc[feat["source"] == "genus", "feature"].tolist()
alpha_cols = feat.loc[feat["source"] == "alpha", "feature"].tolist()
keep = [c for c in genus_cols if (X[c] > 0).mean() >= 0.05]
feats = keep + alpha_cols

m = (y["cohort"] == "PRJNA302078").values
X_xiao = X.loc[m, feats].reset_index(drop=True)
y_xiao = y.loc[m, "outcome"].values
log(f"Xiao: n={len(X_xiao)}, positive={int(y_xiao.sum())}/{len(y_xiao)}")

# ---------------- models ----------------
def make_models():
    lasso = Pipeline([
        ("sc", StandardScaler()),
        ("clf", LogisticRegression(solver="liblinear", l1_ratio=1.0,
                                   class_weight="balanced", max_iter=5000, random_state=RNG)),
    ])
    lasso_grid = {"clf__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]}
    rf = Pipeline([
        ("sc", StandardScaler()),
        ("clf", RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                       random_state=RNG, n_jobs=1)),
    ])
    rf_grid = {"clf__max_depth": [2, 3, 5],
               "clf__min_samples_leaf": [2, 5],
               "clf__max_features": ["sqrt", 0.3]}
    xgb = Pipeline([
        ("sc", StandardScaler()),
        ("clf", XGBClassifier(eval_metric="logloss", random_state=RNG, n_jobs=1)),
    ])
    xgb_grid = {"clf__max_depth": [2, 3],
                "clf__learning_rate": [0.05, 0.1],
                "clf__subsample": [0.8, 1.0],
                "clf__colsample_bytree": [0.7, 1.0],
                "clf__n_estimators": [150, 300]}
    return [("LASSO-LR", lasso, lasso_grid), ("RF", rf, rf_grid), ("XGB", xgb, xgb_grid)]

# ---------------- repeated holdout ----------------
log(f"\n=== Repeated stratified holdout: 50 train / {TEST_SIZE} test, {N_SEEDS} seeds ===")
rows = []
for s in range(N_SEEDS):
    seed = RNG + s
    tr_idx, te_idx = train_test_split(np.arange(len(X_xiao)), test_size=TEST_SIZE,
                                      stratify=y_xiao, random_state=seed)
    Xtr, Xte = X_xiao.iloc[tr_idx], X_xiao.iloc[te_idx]
    ytr, yte = y_xiao[tr_idx], y_xiao[te_idx]
    # sanity: test must contain both classes
    for name, pipe, grid in make_models():
        inner = StratifiedKFold(3, shuffle=True, random_state=seed)
        gs = GridSearchCV(pipe, grid, cv=inner, scoring="roc_auc", n_jobs=1, refit=True)
        gs.fit(Xtr, ytr)
        proba = gs.predict_proba(Xte)[:, 1]
        try:
            a = roc_auc_score(yte, proba)
        except ValueError:
            a = np.nan
        p = average_precision_score(yte, proba)
        b = brier_score_loss(yte, proba)
        rows.append({"seed": s, "model": name, "AUC": a, "PR_AUC": p, "Brier": b})
    if s % 10 == 0:
        log(f"  seed {s}/{N_SEEDS} done")

res = pd.DataFrame(rows)
summ = res.groupby("model").agg(
    AUC_mean=("AUC", "mean"), AUC_sd=("AUC", "std"), AUC_median=("AUC", "median"),
    AUC_q25=("AUC", lambda x: x.quantile(0.25)), AUC_q75=("AUC", lambda x: x.quantile(0.75)),
    AUC_ge06=("AUC", lambda x: (x >= 0.6).mean()),
    PR_AUC_mean=("PR_AUC", "mean"), PR_AUC_sd=("PR_AUC", "std"),
    Brier_mean=("Brier", "mean")).round(4).reset_index()
log(summ.to_string(index=False))
res.to_csv(os.path.join(OUT, "results_xiao_holdout_raw.csv"), index=False)
summ.to_csv(os.path.join(OUT, "results_xiao_holdout.csv"), index=False)
log("\nDONE holdout.")
LOG.close()
