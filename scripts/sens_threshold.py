# -*- coding: utf-8 -*-
"""
sens_threshold.py - Sensitivity analysis of CST4 LB threshold on model performance.
Fixes the 117 subject pairs; rebuilds outcome labels at LB thresholds
0.2/0.3/0.4/0.5/0.6/0.7; repeated 5-fold stratified CV (5 repeats) per threshold.
Models: LASSO-LR / RF / XGB (same grids as train_main.py).
Output: <OUT_DIR>/modeling/results_sens_threshold.csv
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
META = "<DATA_DIR>/meta"
OUT = "<OUT_DIR>/modeling"
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "sens_log.txt"), "w", buffering=1)

def log(*args):
    msg = " ".join(str(a) for a in args)
    print(msg, flush=True)
    LOG.write(msg + "\n")

# ---------------- fixed X ----------------
X = pd.read_csv(os.path.join(FEAT, "X_features.csv"))
feat = pd.read_csv(os.path.join(FEAT, "feature_meta.csv"))
genus_cols = feat.loc[feat["source"] == "genus", "feature"].tolist()
alpha_cols = feat.loc[feat["source"] == "alpha", "feature"].tolist()
keep = [c for c in genus_cols if (X[c] > 0).mean() >= 0.05]
feats = keep + alpha_cols
Xm = X[feats].reset_index(drop=True)
log(f"X fixed: {Xm.shape}")

# ---------------- outcome LB map ----------------
meta = pd.read_csv(os.path.join(META, "sample_meta.csv"))
allcst = pd.read_csv(os.path.join(META, "sample_cst_all.csv"))
cst_map = {(r.cohort, r.sample): r.Lactobacillus_pct for r in allcst.itertuples()}
out_lb = {(r.cohort, r.subject): cst_map[(r.cohort, r.outcome_sample)] for r in meta.itertuples()}
y = pd.read_csv(os.path.join(FEAT, "y_labels.csv"))
keys = list(zip(y["cohort"].values, y["subject"].values))

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

def nested_oof(Xm, ym, models, outer_cv, inner_cv):
    oofs = {}
    for name, pipe, grid in models:
        oof = np.zeros(len(ym))
        for tr, te in outer_cv.split(Xm, ym):
            gs = GridSearchCV(pipe, grid, cv=inner_cv, scoring="roc_auc", n_jobs=1, refit=True)
            gs.fit(Xm.iloc[tr], ym[tr])
            oof[te] = gs.predict_proba(Xm.iloc[te])[:, 1]
        oofs[name] = oof
    return oofs

# ---------------- sensitivity loop ----------------
thresholds = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
rows = []
for t in thresholds:
    yt = np.array([1 if out_lb[k] >= t else 0 for k in keys], dtype=int)
    npos = int(yt.sum())
    if npos in (0, len(yt)):
        log(f"threshold {t}: degenerate outcome ({npos}/{len(yt)}) - skipped")
        continue
    log(f"\n=== LB threshold {t}: positive(recovered)={npos}/{len(yt)} ===")
    for rep in range(5):
        skf = StratifiedKFold(5, shuffle=True, random_state=RNG + rep)
        inner = StratifiedKFold(3, shuffle=True, random_state=RNG + rep)
        oofs = nested_oof(Xm, yt, make_models(), skf, inner)
        for name, oof in oofs.items():
            try:
                a = roc_auc_score(yt, oof)
                p = average_precision_score(yt, oof)
                b = brier_score_loss(yt, oof)
            except ValueError:
                a, p, b = np.nan, np.nan, np.nan
            rows.append({"lb_threshold": t, "repeat": rep, "model": name,
                         "AUC": a, "PR_AUC": p, "Brier": b, "n_pos": npos})
    # print running summary for this threshold
    tmp = pd.DataFrame([r for r in rows if r["lb_threshold"] == t])
    log(tmp.groupby("model").agg(AUC=("AUC", "mean"), PR_AUC=("PR_AUC", "mean")).round(4).to_string())

res = pd.DataFrame(rows)
summ = res.groupby(["lb_threshold", "model"]).agg(
    AUC=("AUC", "mean"), AUC_SD=("AUC", "std"),
    PR_AUC=("PR_AUC", "mean"), PR_AUC_SD=("PR_AUC", "std"),
    Brier=("Brier", "mean")).round(4).reset_index()
log("\n=== SUMMARY: AUC/PR-AUC by LB threshold ===")
log(summ.to_string(index=False))
res.to_csv(os.path.join(OUT, "results_sens_threshold_raw.csv"), index=False)
summ.to_csv(os.path.join(OUT, "results_sens_threshold.csv"), index=False)
log("\nDONE sensitivity.")
LOG.close()
