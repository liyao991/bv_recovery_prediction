# -*- coding: utf-8 -*-
"""
train_main.py - MAIN analysis: BV recurrence prediction (CST4 labels).
Merged training + LOCO CV + repeated 5-fold stratified CV, nested tuning.
Models: LASSO-LR (l1_ratio=1.0), RandomForest, XGBoost.
Metrics: AUC, PR-AUC, Brier.  Saves CSVs + ROC/PR figure.
"""
import os
import sys
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             brier_score_loss, roc_curve, precision_recall_curve, auc)
from xgboost import XGBClassifier

RNG = 42
np.random.seed(RNG)
FEAT = "<DATA_DIR>/feat"
OUT = "<OUT_DIR>/modeling"
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "train_main_log.txt"), "w", buffering=1)

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
ym = y["outcome"].reset_index(drop=True)
log(f"features: {len(keep)} genera (presence>=5%) + {len(alpha_cols)} alpha = {len(feats)}")
log(f"n={len(Xm)}, positive={int(ym.sum())}/{len(ym)}")
log(f"cohort split: {y['cohort'].value_counts().to_dict()}")
log(f"cohort x outcome:\n{pd.crosstab(y['cohort'], y['outcome'])}")

# ---------------- models (sklearn>=1.8 API) ----------------
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

def nested_oof(Xm, ym, models, outer_cv, inner_cv):
    oofs = {}
    for name, pipe, grid in models:
        oof = np.zeros(len(ym))
        for tr, te in outer_cv.split(Xm, ym):
            gs = GridSearchCV(pipe, grid, cv=inner_cv, scoring="roc_auc", n_jobs=1, refit=True)
            gs.fit(Xm.iloc[tr], ym.iloc[tr])
            oof[te] = gs.predict_proba(Xm.iloc[te])[:, 1]
        oofs[name] = oof
    return oofs

def evaluate(oof, yv):
    return (roc_auc_score(yv, oof), average_precision_score(yv, oof),
            brier_score_loss(yv, oof))

# ---------------- 1) repeated 5-fold ----------------
log("\n=== Repeated 5-fold stratified CV (10 repeats) ===")
rows = []
for rep in range(10):
    skf = StratifiedKFold(5, shuffle=True, random_state=RNG + rep)
    inner = StratifiedKFold(3, shuffle=True, random_state=RNG + rep)
    oofs = nested_oof(Xm, ym, make_models(), skf, inner)
    for name, oof in oofs.items():
        a, p, b = evaluate(oof, ym)
        rows.append({"repeat": rep, "model": name, "AUC": a, "PR_AUC": p, "Brier": b})
cv_df = pd.DataFrame(rows)
summ = cv_df.groupby("model").agg(
    AUC=("AUC", "mean"), AUC_SD=("AUC", "std"),
    PR_AUC=("PR_AUC", "mean"), PR_AUC_SD=("PR_AUC", "std"),
    Brier=("Brier", "mean")).round(4).reset_index()
log(summ.to_string(index=False))
cv_df.to_csv(os.path.join(OUT, "results_main_cv.csv"), index=False)
summ.to_csv(os.path.join(OUT, "results_main_summary.csv"), index=False)

# ---------------- 2) LOCO ----------------
log("\n=== LOCO (leave-one-cohort-out) ===")
cohort = y["cohort"].values
loco_rows = []
oof_loco = {}
for name, pipe, grid in make_models():
    oof = np.zeros(len(ym))
    for c in np.unique(cohort):
        tr = np.where(cohort != c)[0]
        te = np.where(cohort == c)[0]
        inner = StratifiedKFold(3, shuffle=True, random_state=RNG)
        gs = GridSearchCV(pipe, grid, cv=inner, scoring="roc_auc", n_jobs=1, refit=True)
        gs.fit(Xm.iloc[tr], ym.iloc[tr])
        oof[te] = gs.predict_proba(Xm.iloc[te])[:, 1]
    a, p, b = evaluate(oof, ym)
    oof_loco[name] = oof
    # per-cohort metrics
    row = {"model": name, "AUC_LOCO": a, "PR_AUC_LOCO": p, "Brier_LOCO": b}
    for c in np.unique(cohort):
        m = cohort == c
        row[f"AUC_{c}"] = roc_auc_score(ym[m], oof[m]) if len(np.unique(ym[m])) > 1 else np.nan
    loco_rows.append(row)
loco_df = pd.DataFrame(loco_rows)
log(loco_df.round(4).to_string(index=False))
loco_df.to_csv(os.path.join(OUT, "results_loco.csv"), index=False)

# ROC/PR figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for name, oof in oof_loco.items():
    fpr, tpr, _ = roc_curve(ym, oof)
    axes[0].plot(fpr, tpr, lw=1.8, label=f"{name} (AUC={roc_auc_score(ym, oof):.3f})")
    prec, rec, _ = precision_recall_curve(ym, oof)
    axes[1].plot(rec, prec, lw=1.8, label=f"{name} (PR-AUC={auc(rec, prec):.3f})")
axes[0].plot([0, 1], [0, 1], "k--", lw=0.8)
axes[0].set(xlabel="False positive rate", ylabel="True positive rate", title="ROC — LOCO CV")
axes[1].set(xlabel="Recall", ylabel="Precision", title="Precision–Recall — LOCO CV")
for ax in axes:
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_roc_pr_loco.png"), dpi=150)
plt.close(fig)
log("figure saved -> fig_roc_pr_loco.png")

# ---------------- 3) subgroup: baseline CST4 ----------------
log("\n=== Subgroup: baseline CST4 (dysbiosis at baseline) ===")
idx = y["base_cst4"] == 1
Xs, ys_ = Xm[idx.values].reset_index(drop=True), y.loc[idx, "outcome"].reset_index(drop=True)
log(f"n={len(Xs)}, positive={int(ys_.sum())}/{len(ys_)}")
rows = []
for rep in range(10):
    skf = StratifiedKFold(5, shuffle=True, random_state=RNG + rep)
    inner = StratifiedKFold(3, shuffle=True, random_state=RNG + rep)
    oofs = nested_oof(Xs, ys_, make_models(), skf, inner)
    for name, oof in oofs.items():
        a, p, b = evaluate(oof, ys_)
        rows.append({"repeat": rep, "model": name, "AUC": a, "PR_AUC": p, "Brier": b})
cv_sub = pd.DataFrame(rows)
summ_sub = cv_sub.groupby("model").agg(
    AUC=("AUC", "mean"), AUC_SD=("AUC", "std"),
    PR_AUC=("PR_AUC", "mean"), PR_AUC_SD=("PR_AUC", "std"),
    Brier=("Brier", "mean")).round(4).reset_index()
log(summ_sub.to_string(index=False))
summ_sub.to_csv(os.path.join(OUT, "results_sub_cst4.csv"), index=False)

log("\nDONE main analysis.")
LOG.close()
