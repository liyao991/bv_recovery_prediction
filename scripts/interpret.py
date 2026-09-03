# -*- coding: utf-8 -*-
"""
interpret.py - Model interpretation for the BV recurrence prediction.
1. Picks best model by LOCO AUC (most honest, cross-cohort metric).
2. Refits on full data with inner-CV tuned hyperparams.
3. SHAP summary (RF/XGB) or LASSO coefficients + top features.
4. Saves: fig_shap.png, fig_importance.png, feature_importance.csv
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from xgboost import XGBClassifier
import shap

RNG = 42
FEAT = "<DATA_DIR>/feat"
OUT = "<OUT_DIR>/modeling"
os.makedirs(OUT, exist_ok=True)

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
print(f"X {Xm.shape}, positive={int(ym.sum())}/{len(ym)}", flush=True)

# ---------------- model zoo ----------------
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
    return {"LASSO-LR": (lasso, lasso_grid), "RF": (rf, rf_grid), "XGB": (xgb, xgb_grid)}

# ---------------- pick best model by LOCO AUC ----------------
loco = pd.read_csv(os.path.join(OUT, "results_loco.csv")).set_index("model")
best = loco["AUC_LOCO"].idxmax()
print(f"Best model by LOCO AUC: {best} (AUC_LOCO={loco.loc[best,'AUC_LOCO']:.3f})", flush=True)

models = make_models()
pipe, grid = models[best]

# ---------------- refit on full data (inner CV for hyperparams) ----------------
inner = StratifiedKFold(3, shuffle=True, random_state=RNG)
gs = GridSearchCV(pipe, grid, cv=inner, scoring="roc_auc", n_jobs=1, refit=True)
gs.fit(Xm, ym)
print(f"Best params: {gs.best_params_}", flush=True)
clf = gs.best_estimator_

# ---------------- SHAP / coefficients ----------------
estimator = clf.named_steps["clf"]
Xs = clf.named_steps["sc"].transform(Xm)
feature_labels = feats

if best == "LASSO-LR":
    coef = estimator.coef_.ravel()
    imp = pd.DataFrame({"feature": feature_labels, "coef": coef, "abs_coef": np.abs(coef)})
    imp = imp.sort_values("abs_coef", ascending=False).reset_index(drop=True)
    imp.to_csv(os.path.join(OUT, "feature_importance.csv"), index=False)
    print("LASSO non-zero coefs:", (imp["coef"] != 0).sum(), "of", len(imp), flush=True)
    print(imp.head(20).to_string(index=False), flush=True)
    # coefficient bar chart (top 20)
    top = imp.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 6.5))
    colors = ["#c0392b" if c > 0 else "#2980b9" for c in top["coef"]]
    ax.barh(top["feature"], top["coef"], color=colors)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_title("LASSO coefficients (positive -> predicts recovered)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_lasso_coef.png"), dpi=150)
    plt.close(fig)
else:
    if best == "RF":
        explainer = shap.TreeExplainer(estimator)
        shap_values = explainer.shap_values(Xs)[1]  # class 1 = recovered
    else:
        explainer = shap.TreeExplainer(estimator)
        shap_values = explainer.shap_values(Xs)
    # importance
    imp = pd.DataFrame({"feature": feature_labels,
                        "mean_abs_shap": np.abs(shap_values).mean(axis=0)})
    imp = imp.sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    imp.to_csv(os.path.join(OUT, "feature_importance.csv"), index=False)
    print(imp.head(20).to_string(index=False), flush=True)
    # SHAP summary plot
    fig, ax = plt.subplots(figsize=(8, 6.5))
    shap.summary_plot(shap_values, Xs, feature_names=feature_labels,
                      max_display=20, show=False)
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_shap.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("SHAP figure saved -> fig_shap.png", flush=True)

print("\nDONE interpretation.", flush=True)
