# -*- coding: utf-8 -*-
"""
cv_within_clr.py - Diagnose batch effects & test compositional (CLR) fix.
Part A: within-cohort repeated 5-fold CV on RAW relative-abundance features.
Part B: CLR-transformed features - within-cohort CV + LOCO.
Models: LASSO-LR / RF / XGB (same grids as train_main.py).
Output: <OUT_DIR>/modeling/results_within_cohort.csv, results_clr.csv
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
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from xgboost import XGBClassifier

RNG = 42
FEAT = "<DATA_DIR>/feat"
OUT = "<OUT_DIR>/modeling"
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "cv_within_clr_log.txt"), "w", buffering=1)

def log(*args):
    msg = " ".join(str(a) for a in args)
    print(msg, flush=True)
    LOG.write(msg + "\n")

# ---------------- data ----------------
X = pd.read_csv(os.path.join(FEAT, "X_features.csv"))
y = pd.read_csv(os.path.join(FEAT, "y_labels.csv"))
feat = pd.read_csv(os.path.join(FEAT, "feature_meta.csv"))
genus_cols = feat.loc[feat["source"] == "genus", "feature"].tolist()
alpha_cols = feat.loc[feat["source"] == "alpha", "feature"].tolist()
keep = [c for c in genus_cols if (X[c] > 0).mean() >= 0.05]
feats = keep + alpha_cols
Xm = X[feats].reset_index(drop=True)
ym = y["outcome"].values
cohort = y["cohort"].values
log(f"X {Xm.shape}, positive={int(ym.sum())}/{len(ym)}")

# ---------------- CLR transform (per-sample, genus only) ----------------
def clr_transform(df):
    eps = 1e-6
    gm = df.values.copy()
    gm[gm == 0] = eps
    lg = np.log(gm)
    clr = lg - lg.mean(axis=1, keepdims=True)
    return pd.DataFrame(clr, columns=df.columns, index=df.index)

X_clr = Xm.copy()
X_clr[keep] = clr_transform(Xm[keep])
log(f"CLR computed on {len(keep)} genera")

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

def nested_oof(Xd, yv, models, outer_cv, inner_cv):
    oofs = {}
    for name, pipe, grid in models:
        oof = np.zeros(len(yv))
        for tr, te in outer_cv.split(Xd, yv):
            gs = GridSearchCV(pipe, grid, cv=inner_cv, scoring="roc_auc", n_jobs=1, refit=True)
            gs.fit(Xd.iloc[tr], yv[tr])
            oof[te] = gs.predict_proba(Xd.iloc[te])[:, 1]
        oofs[name] = oof
    return oofs

def safe_metrics(yv, oof):
    try:
        return roc_auc_score(yv, oof), average_precision_score(yv, oof), brier_score_loss(yv, oof)
    except ValueError:
        return np.nan, np.nan, np.nan

# ---------------- Part A: within-cohort CV (raw) ----------------
log("\n=== Part A: within-cohort repeated 5-fold CV (raw features, 5 repeats) ===")
rows = []
for c in np.unique(cohort):
    m = cohort == c
    Xc, yc = Xm[m].reset_index(drop=True), ym[m]
    if len(np.unique(yc)) < 2:
        log(f"{c}: single-class ({yc.sum()}/{len(yc)}) - skipped")
        continue
    log(f"\n--- {c}: n={len(yc)}, positive={int(yc.sum())}/{len(yc)} ---")
    for rep in range(5):
        skf = StratifiedKFold(5, shuffle=True, random_state=RNG + rep)
        inner = StratifiedKFold(3, shuffle=True, random_state=RNG + rep)
        oofs = nested_oof(Xc, yc, make_models(), skf, inner)
        for name, oof in oofs.items():
            a, p, b = safe_metrics(yc, oof)
            rows.append({"part": "within_raw", "cohort": c, "repeat": rep,
                         "model": name, "AUC": a, "PR_AUC": p, "Brier": b})
within_raw = pd.DataFrame(rows)
log(within_raw.groupby(["cohort", "model"]).agg(
    AUC=("AUC", "mean"), PR_AUC=("PR_AUC", "mean")).round(4).to_string())

# ---------------- Part B: CLR - within-cohort CV ----------------
log("\n=== Part B1: CLR features - within-cohort CV (5 repeats) ===")
rowsB = []
for c in np.unique(cohort):
    m = cohort == c
    Xc, yc = X_clr[m].reset_index(drop=True), ym[m]
    if len(np.unique(yc)) < 2:
        continue
    log(f"\n--- {c}: n={len(yc)}, positive={int(yc.sum())}/{len(yc)} ---")
    for rep in range(5):
        skf = StratifiedKFold(5, shuffle=True, random_state=RNG + rep)
        inner = StratifiedKFold(3, shuffle=True, random_state=RNG + rep)
        oofs = nested_oof(Xc, yc, make_models(), skf, inner)
        for name, oof in oofs.items():
            a, p, b = safe_metrics(yc, oof)
            rowsB.append({"part": "within_clr", "cohort": c, "repeat": rep,
                          "model": name, "AUC": a, "PR_AUC": p, "Brier": b})
within_clr = pd.DataFrame(rowsB)
log(within_clr.groupby(["cohort", "model"]).agg(
    AUC=("AUC", "mean"), PR_AUC=("PR_AUC", "mean")).round(4).to_string())

# ---------------- Part B2: CLR - LOCO ----------------
log("\n=== Part B2: CLR features - LOCO ===")
loco_rows = []
for name, pipe, grid in make_models():
    oof = np.zeros(len(ym))
    for c in np.unique(cohort):
        tr = np.where(cohort != c)[0]
        te = np.where(cohort == c)[0]
        inner = StratifiedKFold(3, shuffle=True, random_state=RNG)
        gs = GridSearchCV(pipe, grid, cv=inner, scoring="roc_auc", n_jobs=1, refit=True)
        gs.fit(X_clr.iloc[tr], ym[tr])
        oof[te] = gs.predict_proba(X_clr.iloc[te])[:, 1]
    a, p, b = safe_metrics(ym, oof)
    row = {"part": "loco_clr", "model": name, "AUC": a, "PR_AUC": p, "Brier": b}
    for c in np.unique(cohort):
        m = cohort == c
        row[f"AUC_{c}"] = safe_metrics(ym[m], oof[m])[0]
    loco_rows.append(row)
    log(f"{name}: AUC_LOCO={a:.4f}  Xiao={row['AUC_PRJNA302078']:.4f}  Plummer={row['AUC_PRJNA592384']:.4f}")

loco_clr = pd.DataFrame(loco_rows)
all_res = pd.concat([within_raw, within_clr, loco_clr], ignore_index=True)
all_res.to_csv(os.path.join(OUT, "results_within_clr.csv"), index=False)
log("\n=== SUMMARY: CLR-LOCO ===")
log(loco_clr.round(4).to_string(index=False))
log("\nDONE within/CLR analysis.")
LOG.close()
