# -*- coding: utf-8 -*-
"""
xiao_within_interpret.py - Within-Xiao LASSO coefficients for the only robust signal.
Xiao is the only cohort with usable within-cohort CV AUC (~0.65);
Plummer has only 5 negatives (unusable). Fits LASSO on Xiao with inner-CV tuned C,
extracts coefficients and saves a bar chart.
"""
import os
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
from sklearn.model_selection import StratifiedKFold, GridSearchCV

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

# Xiao only
m = (y["cohort"] == "PRJNA302078").values
X_xiao = X.loc[m, feats].reset_index(drop=True)
y_xiao = y.loc[m, "outcome"].reset_index(drop=True)
print(f"Xiao: n={len(X_xiao)}, positive={int(y_xiao.sum())}/{len(y_xiao)}", flush=True)

# refit on all Xiao with inner-CV tuned C
pipe = Pipeline([
    ("sc", StandardScaler()),
    ("clf", LogisticRegression(solver="liblinear", l1_ratio=1.0,
                               class_weight="balanced", max_iter=5000, random_state=RNG)),
])
inner = StratifiedKFold(3, shuffle=True, random_state=RNG)
gs = GridSearchCV(pipe, {"clf__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]},
                  cv=inner, scoring="roc_auc", n_jobs=1, refit=True)
gs.fit(X_xiao, y_xiao)
print(f"Best C: {gs.best_params_}", flush=True)
estimator = gs.best_estimator_.named_steps["clf"]
coef = estimator.coef_.ravel()
nonzero = int((coef != 0).sum())
print(f"Non-zero coefs: {nonzero}/{len(coef)}", flush=True)

imp = pd.DataFrame({"feature": feats, "coef": coef, "abs_coef": np.abs(coef)})
imp = imp.sort_values("abs_coef", ascending=False).reset_index(drop=True)
imp.to_csv(os.path.join(OUT, "feature_importance_xiao_lasso.csv"), index=False)
print("\nTop 20 features (positive coef -> predicts recovered):", flush=True)
print(imp.head(20).round(4).to_string(index=False), flush=True)

# Bar chart: top 20 (or fewer if non-zero < 20)
top = imp[imp["coef"] != 0].head(20).iloc[::-1]
fig, ax = plt.subplots(figsize=(7, 6))
colors = ["#c0392b" if c > 0 else "#2980b9" for c in top["coef"]]
ax.barh(top["feature"], top["coef"], color=colors)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("LASSO coefficient (standardized)")
ax.set_title("Within-Xiao LASSO: predictors of recovery (n=68, CV AUC=0.65)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_xiao_lasso_coef.png"), dpi=150)
plt.close(fig)
print("\nFigure saved -> fig_xiao_lasso_coef.png", flush=True)
