# -*- coding: utf-8 -*-
"""
train_models.py - BV recurrence prediction from baseline vaginal microbiota
(CST4-based outcome labels).

Analysis design (user-approved: merged training + LOCO CV):
  * Population:
      - MAIN  : all 117 baseline->outcome pairs
      - SUB   : baseline CST4 subset (n=97)  [dysbiosis-at-baseline -> recovery]
  * Features: baseline genus rel. abund. (log1p(100x)) + alpha diversity
      - filter genera present in >=5% samples
      - StandardScaler inside CV folds (no leakage)
  * Models : LASSO-LR, RandomForest, XGBoost (nested CV: inner GridSearch)
  * CV     : (a) LOCO (leave-one-cohort-out)  (b) repeated 5-fold stratified
  * Metrics: AUC, PR-AUC (average precision), Brier, with SD across repeats
  * Sensitivity: CST threshold 20/30/40/50/70% -> rebuild labels -> 5-fold AUC
  * Interpretability: SHAP (best model) + LASSO coefficients

Outputs (ASCII):
  <OUT_DIR>/modeling/results_main.csv
  <OUT_DIR>/modeling/results_threshold_sweep.csv
  <OUT_DIR>/modeling/lasso_coefs.csv
  <OUT_DIR>/modeling/shap_importance.csv
  <OUT_DIR>/modeling/fig_roc_pr.png
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV, RepeatedStratifiedKFold
from sklearn.metrics import (roc_auc_score, average_precision_score, brier_score_loss,
                             roc_curve, precision_recall_curve, auc)
from xgboost import XGBClassifier

warnings.filterwarnings("ignore")
RNG = 42
np.random.seed(RNG)

FEAT = "<DATA_DIR>/feat"
META = "<DATA_DIR>/meta"
OUT = "<OUT_DIR>/modeling"
os.makedirs(OUT, exist_ok=True)

BV_GENERA = ["Gardnerella", "Atopobium", "Prevotella", "Sneathia",
             "Megasphaera", "Mobiluncus", "Fannyhessea", "Dialister",
             "Prevotella_7", "Prevotella_9", "Aerococcus", "Fastidiosipila"]

# ---------------- data ----------------
def load_data():
    X = pd.read_csv(os.path.join(FEAT, "X_features.csv"))
    y = pd.read_csv(os.path.join(FEAT, "y_labels.csv"))
    feat = pd.read_csv(os.path.join(FEAT, "feature_meta.csv"))
    genus_cols = feat.loc[feat["source"] == "genus", "feature"].tolist()
    alpha_cols = feat.loc[feat["source"] == "alpha", "feature"].tolist()
    return X, y, genus_cols, alpha_cols

def filter_features(X, genus_cols, min_presence=0.05):
    keep = [c for c in genus_cols if (X[c] > 0).mean() >= min_presence]
    drop = [c for c in genus_cols if c not in keep]
    return keep, drop

# ---------------- models ----------------
def make_models():
    lasso = Pipeline([
        ("sc", StandardScaler()),
        ("clf", LogisticRegression(penalty="l1", solver="liblinear",
                                   class_weight="balanced", max_iter=5000, random_state=RNG)),
    ])
    lasso_grid = {"clf__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]}

    rf = Pipeline([
        ("sc", StandardScaler()),
        ("clf", RandomForestClassifier(n_estimators=600, class_weight="balanced",
                                       random_state=RNG, n_jobs=1)),
    ])
    rf_grid = {"clf__max_depth": [2, 3, 5, None],
               "clf__min_samples_leaf": [2, 5, 10],
               "clf__max_features": ["sqrt", 0.3]}

    xgb = Pipeline([
        ("sc", StandardScaler()),
        ("clf", XGBClassifier(eval_metric="logloss", random_state=RNG,
                              n_jobs=1, use_label_encoder=False)),
    ])
    xgb_grid = {"clf__max_depth": [2, 3],
                "clf__learning_rate": [0.02, 0.05, 0.1],
                "clf__subsample": [0.7, 1.0],
                "clf__colsample_bytree": [0.5, 1.0],
                "clf__n_estimators": [150, 300]}

    return [("LASSO-LR", lasso, lasso_grid), ("RF", rf, rf_grid), ("XGB", xgb, xgb_grid)]

def nested_cv(Xm, ym, models, outer_cv, inner_cv, verbose=False):
    """Outer CV with inner GridSearch; returns OOF preds per model."""
    results = {}
    for name, pipe, grid in models:
        oof = np.zeros(len(ym))
        for tr, te in outer_cv.split(Xm, ym):
            gs = GridSearchCV(pipe, grid, cv=inner_cv, scoring="roc_auc",
                              n_jobs=1, refit=True)
            gs.fit(Xm.iloc[tr], ym.iloc[tr])
            oof[te] = gs.predict_proba(Xm.iloc[te])[:, 1]
        results[name] = oof
        if verbose:
            print(f"  {name}: OOF computed")
    return results

def evaluate(oof, ym):
    auc = roc_auc_score(ym, oof)
    prauc = average_precision_score(ym, oof)
    brier = brier_score_loss(ym, oof)
    return auc, prauc, brier

def run_cv_repeated(Xm, ym, models, n_repeats=10):
    """Repeated stratified 5-fold CV -> mean +/- SD."""
    out = []
    for rep in range(n_repeats):
        skf = StratifiedKFold(5, shuffle=True, random_state=RNG + rep)
        res = nested_cv(Xm, ym, models, skf, StratifiedKFold(3, shuffle=True, random_state=RNG + rep))
        for name, oof in res.items():
            a, p, b = evaluate(oof, ym)
            out.append({"repeat": rep, "model": name, "AUC": a, "PR_AUC": p, "Brier": b})
    return pd.DataFrame(out)

def run_loco(Xm, ym, models):
    """Leave-one-cohort-out: cohort labels from y_labels."""
    cohort = ym["cohort"].values if "cohort" in ym.columns else None
    # ym passed with cohort column attached externally
    cohorts = np.unique(cohort)
    rows = []
    oof_all = {}
    for name, pipe, grid in models:
        oof = np.zeros(len(ym))
        for c in cohorts:
            tr = np.where(cohort != c)[0]
            te = np.where(cohort == c)[0]
            gs = GridSearchCV(pipe, grid, cv=StratifiedKFold(3, shuffle=True, random_state=RNG),
                              scoring="roc_auc", n_jobs=1, refit=True)
            gs.fit(Xm.iloc[tr], ym.iloc[tr])
            oof[te] = gs.predict_proba(Xm.iloc[te])[:, 1]
        a, p, b = evaluate(oof, ym)
        oof_all[name] = oof
        rows.append({"model": name, "AUC_LOCO": a, "PR_AUC_LOCO": p, "Brier_LOCO": b})
    return pd.DataFrame(rows), oof_all

# ---------------- sensitivity sweep ----------------
def rebuild_labels(threshold):
    """Rebuild outcome from per-sample LB% at new CST4 threshold."""
    cst = pd.read_csv(os.path.join(META, "sample_cst_all.csv"))
    cst["cst4"] = (cst["Lactobacillus_pct"] < threshold).astype(int)
    pairs = []
    for cohort, bt, ot in [("PRJNA302078", "D0", "D30"), ("PRJNA592384", "M0", "M1")]:
        sub = cst[cst["cohort"] == cohort]
        base = sub[sub["timepoint"] == bt].set_index("subject")
        outc = sub[sub["timepoint"] == ot].set_index("subject")
        for sid in base.index.intersection(outc.index):
            o = outc.loc[sid]
            pairs.append({"cohort": cohort, "subject": sid,
                          "outcome": 1 if int(o["cst4"]) == 0 else 0})
    df = pd.DataFrame(pairs)
    # align to the same 117 subjects used in main analysis
    y_main = pd.read_csv(os.path.join(FEAT, "y_labels.csv"))
    key = y_main["cohort"].astype(str) + "_" + y_main["subject"].astype(str)
    df["key"] = df["cohort"].astype(str) + "_" + df["subject"].astype(str)
    df = df[df["key"].isin(set(key))].reset_index(drop=True)
    return df["outcome"].values

def threshold_sweep(Xm, genus_cols, alpha_cols, models):
    threshs = [0.20, 0.30, 0.40, 0.50, 0.60, 0.70]
    y_main = pd.read_csv(os.path.join(FEAT, "y_labels.csv"))
    rows = []
    for t in threshs:
        yv = rebuild_labels(t)
        # use same feature matrix; drop subjects not in yv? all 117 retained by design
        skf = StratifiedKFold(5, shuffle=True, random_state=RNG)
        aucs, praucs = {}, {}
        for name, pipe, grid in models:
            oof = np.zeros(len(yv))
            for tr, te in skf.split(Xm, yv):
                gs = GridSearchCV(pipe, grid, cv=StratifiedKFold(3, shuffle=True, random_state=RNG),
                                  scoring="roc_auc", n_jobs=1, refit=True)
                gs.fit(Xm.iloc[tr], yv[tr])
                oof[te] = gs.predict_proba(Xm.iloc[te])[:, 1]
            aucs[name] = roc_auc_score(yv, oof)
            praucs[name] = average_precision_score(yv, oof)
        for name in aucs:
            rows.append({"threshold": t, "model": name, "AUC": aucs[name], "PR_AUC": praucs[name]})
    df = pd.DataFrame(rows)
    df["n_positive"] = [None] * len(df)
    for t in threshs:
        yv = rebuild_labels(t)
        df.loc[df["threshold"] == t, "n_positive"] = int(yv.sum())
    return df

# ---------------- plots ----------------
def plot_roc_pr(oof_dict, ym, fname):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for name, oof in oof_dict.items():
        fpr, tpr, _ = roc_curve(ym, oof)
        axes[0].plot(fpr, tpr, lw=1.8, label=f"{name} (AUC={roc_auc_score(ym, oof):.3f})")
        prec, rec, _ = precision_recall_curve(ym, oof)
        axes[1].plot(rec, prec, lw=1.8, label=f"{name} (PR-AUC={auc(rec, prec):.3f})")
    axes[0].plot([0, 1], [0, 1], "k--", lw=0.8)
    axes[0].set(xlabel="FPR", ylabel="TPR", title="ROC (LOCO)")
    axes[1].set(xlabel="Recall", ylabel="Precision", title="PR curve (LOCO)")
    for ax in axes:
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=150)
    plt.close(fig)

def main():
    X, y, genus_cols, alpha_cols = load_data()
    keep, drop = filter_features(X, genus_cols, min_presence=0.05)
    print(f"genus: {len(genus_cols)} -> keep {len(keep)}, drop {len(drop)} (presence<5%)")
    feats = keep + alpha_cols
    Xm = X[feats].reset_index(drop=True)
    ym = y["outcome"].reset_index(drop=True)
    models = make_models()

    print(f"\n[MAIN] n={len(Xm)}, features={len(feats)}, pos={ym.sum()}/{len(ym)}")
    # --- repeated 5-fold ---
    cv_df = run_cv_repeated(Xm, ym, models, n_repeats=10)
    summ = cv_df.groupby("model").agg(
        AUC=("AUC", "mean"), AUC_SD=("AUC", "std"),
        PR_AUC=("PR_AUC", "mean"), PR_AUC_SD=("PR_AUC", "std"),
        Brier=("Brier", "mean")).round(4).reset_index()
    print("\nRepeated 5-fold CV (mean±SD):")
    print(summ.to_string(index=False))

    # --- LOCO ---
    ym2 = y[["outcome", "cohort"]].reset_index(drop=True)
    loco_df, oof_all = run_loco(Xm, ym2, models)
    print("\nLOCO (leave-one-cohort-out):")
    print(loco_df.round(4).to_string(index=False))
    for _, r in loco_df.iterrows():
        print(f"  {r['model']}: LOCO AUC={r['AUC_LOCO']:.3f}, PR-AUC={r['PR_AUC_LOCO']:.3f}")

    # save main results
    cv_df.to_csv(os.path.join(OUT, "results_main_cv.csv"), index=False)
    summ.to_csv(os.path.join(OUT, "results_main_summary.csv"), index=False)
    loco_df.to_csv(os.path.join(OUT, "results_loco.csv"), index=False)
    plot_roc_pr(oof_all, ym.values, "fig_roc_pr_loco.png")

    # --- subgroup: baseline CST4 ---
    sub_idx = y["base_cst4"] == 1
    Xs, ys = Xm[sub_idx.values].reset_index(drop=True), y.loc[sub_idx, "outcome"].reset_index(drop=True)
    print(f"\n[SUB baseline-CST4] n={len(Xs)}, pos={ys.sum()}/{len(ys)}")
    cv_sub = run_cv_repeated(Xs, ys, models, n_repeats=10)
    summ_sub = cv_sub.groupby("model").agg(
        AUC=("AUC", "mean"), AUC_SD=("AUC", "std"),
        PR_AUC=("PR_AUC", "mean"), PR_AUC_SD=("PR_AUC", "std"),
        Brier=("Brier", "mean")).round(4).reset_index()
    print(summ_sub.to_string(index=False))
    summ_sub.to_csv(os.path.join(OUT, "results_sub_cst4.csv"), index=False)

    # --- threshold sensitivity ---
    print("\n[Threshold sensitivity]")
    sw = threshold_sweep(Xm, keep, alpha_cols, models)
    print(sw.round(4).to_string(index=False))
    sw.to_csv(os.path.join(OUT, "results_threshold_sweep.csv"), index=False)

    # --- LASSO coefficients (full-data refit) ---
    from sklearn.preprocessing import StandardScaler as SS
    sc = SS().fit(Xm)
    lr = LogisticRegression(penalty="l1", solver="liblinear", class_weight="balanced",
                            C=0.3, max_iter=5000, random_state=RNG).fit(sc.transform(Xm), ym)
    coefs = pd.DataFrame({"feature": feats, "coef": lr.coef_[0]})
    coefs = coefs.reindex(coefs["coef"].abs().sort_values(ascending=False).index)
    coefs.to_csv(os.path.join(OUT, "lasso_coefs.csv"), index=False)
    print("\nTop 12 LASSO |coef|:")
    print(coefs.head(12).round(4).to_string(index=False))

    # --- SHAP (XGB refit on full data, 30 features) ---
    try:
        import shap
        from sklearn.feature_selection import SelectKBest, f_classif
        k = min(30, len(feats))
        sel = SelectKBest(f_classif, k=k).fit(Xm, ym)
        Xk = sel.transform(Xm)
        cols_k = [feats[i] for i in sel.get_support(indices=True)]
        xgb = XGBClassifier(max_depth=2, learning_rate=0.05, n_estimators=300,
                            subsample=0.8, colsample_bytree=0.8, eval_metric="logloss",
                            random_state=RNG, n_jobs=1)
        xgb.fit(Xk, ym)
        explainer = shap.TreeExplainer(xgb)
        sh = explainer.shap_values(Xk)
        imp = pd.DataFrame({"feature": cols_k,
                            "mean_abs_shap": np.abs(sh).mean(0)}).sort_values("mean_abs_shap", ascending=False)
        imp.to_csv(os.path.join(OUT, "shap_importance.csv"), index=False)
        print("\nTop 10 SHAP importance (XGB, top-30 features):")
        print(imp.head(10).round(4).to_string(index=False))
        fig = plt.figure(figsize=(6, 5))
        shap.summary_plot(sh, Xk, feature_names=cols_k, show=False, max_display=15)
        plt.tight_layout()
        fig.savefig(os.path.join(OUT, "fig_shap.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"SHAP skipped: {e}")

    print(f"\nAll outputs -> {OUT}")

if __name__ == "__main__":
    main()
