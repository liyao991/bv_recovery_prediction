# -*- coding: utf-8 -*-
"""
confirm_clr_loco.py - Robustness check of LASSO + CLR cross-cohort AUC.
Re-runs LOCO with 5 different inner-CV seeds; also tests variants:
  - LASSO with class_weight=None (raw class prior)
  - LASSO with ElasticNet (L1+L2) for stability
  - LASSO + raw features (baseline check)
"""
import os
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

RNG = 42
FEAT = "<DATA_DIR>/feat"
OUT = "<OUT_DIR>/modeling"

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

def clr_transform(df):
    eps = 1e-6
    gm = df.values.copy(); gm[gm == 0] = eps
    return pd.DataFrame(np.log(gm) - np.log(gm).mean(axis=1, keepdims=True),
                        columns=df.columns, index=df.index)

X_clr = Xm.copy()
X_clr[keep] = clr_transform(Xm[keep])

def make_lasso(class_weight="balanced"):
    return Pipeline([
        ("sc", StandardScaler()),
        ("clf", LogisticRegression(solver="liblinear", l1_ratio=1.0,
                                   class_weight=class_weight,
                                   max_iter=5000, random_state=RNG)),
    ])
GRID = {"clf__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]}

def loco_auc(X, y, cohort, pipe, n_seeds=5):
    aucs = []
    for s in range(n_seeds):
        oof = np.zeros(len(y))
        for c in np.unique(cohort):
            tr = np.where(cohort != c)[0]
            te = np.where(cohort == c)[0]
            inner = StratifiedKFold(3, shuffle=True, random_state=s * 17 + 7)
            gs = GridSearchCV(pipe, GRID, cv=inner, scoring="roc_auc", n_jobs=1, refit=True)
            gs.fit(X.iloc[tr], y[tr])
            oof[te] = gs.predict_proba(X.iloc[te])[:, 1]
        aucs.append(roc_auc_score(y, oof))
    return np.array(aucs)

variants = [
    ("LASSO+CLR, balanced",     X_clr, make_lasso("balanced")),
    ("LASSO+CLR, no_cw",        X_clr, make_lasso(None)),
    ("LASSO+raw, balanced",     Xm,    make_lasso("balanced")),
    ("LASSO+raw, no_cw",        Xm,    make_lasso(None)),
]

print("=== Robustness: LOCO AUC over 5 random inner-CV seeds ===", flush=True)
rows = []
for name, Xd, pipe in variants:
    a = loco_auc(Xd, ym, cohort, pipe, n_seeds=5)
    print(f"{name:30s}  mean={a.mean():.4f}  sd={a.std():.4f}  vals={a.round(3).tolist()}", flush=True)
    rows.append({"variant": name, "mean": a.mean(), "sd": a.std(),
                 "v1": a[0], "v2": a[1], "v3": a[2], "v4": a[3], "v5": a[4]})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "results_clr_robustness.csv"), index=False)
print("DONE.", flush=True)
