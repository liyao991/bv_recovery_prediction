# -*- coding: utf-8 -*-
"""
diag_plummer_within.py - 快速归因 Plummer 队列内 CV 不可解读的原因。
问题：Plummer 队列内 CV LASSO 0.431 / RF 0.093 / XGB 0.089，不可解读。
假设A：类别极端不平衡（44/49 阳性）+ class_weight=balanced 把少数类预测反了 → no_cw 应回到 ~0.5
假设B：任务本身无信号（结局几乎恒定/基线 CST4 可预测性上限）→ no_cw 仍 ~0.5
额外检查：Plummer 基线 CST4 亚组结局分布 + 预测分数方向。
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

FEAT = "<DATA_DIR>/feat"
X = pd.read_csv(FEAT + "/X_features.csv")
y = pd.read_csv(FEAT + "/y_labels.csv")
feat = pd.read_csv(FEAT + "/feature_meta.csv")
genus_cols = feat.loc[feat["source"] == "genus", "feature"].tolist()
alpha_cols = feat.loc[feat["source"] == "alpha", "feature"].tolist()
keep = [c for c in genus_cols if (X[c] > 0).mean() >= 0.05]
feats = keep + alpha_cols
Xm = X[feats].reset_index(drop=True)
ym = y["outcome"].values
cohort = y["cohort"].values
sub = y["subject"].values if "subject" in y.columns else None

def run_lasso(Xd, yv, cw, reps=10):
    pipe = Pipeline([("sc", StandardScaler()),
                     ("clf", LogisticRegression(solver="liblinear", l1_ratio=1.0,
                                                class_weight=cw, max_iter=5000, random_state=42))])
    grid = {"clf__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]}
    aucs, oof_all = [], []
    for rep in range(reps):
        skf = StratifiedKFold(5, shuffle=True, random_state=42 + rep)
        inner = StratifiedKFold(3, shuffle=True, random_state=42 + rep)
        oof = np.zeros(len(yv))
        for tr, te in skf.split(Xd, yv):
            gs = GridSearchCV(pipe, grid, cv=inner, scoring="roc_auc", n_jobs=1, refit=True)
            gs.fit(Xd.iloc[tr], yv[tr])
            oof[te] = gs.predict_proba(Xd.iloc[te])[:, 1]
        aucs.append(roc_auc_score(yv, oof))
        oof_all.append(oof)
    return np.mean(aucs), np.std(aucs), aucs, oof_all

# ---------- 1) 队列内 CV：balanced vs no_cw ----------
for c in np.unique(cohort):
    m = cohort == c
    Xc, yc = Xm[m].reset_index(drop=True), ym[m]
    print(f"\n=== {c}: n={len(yc)}, pos={int(yc.sum())}/{len(yc)} (neg={int((1-yc).sum())}) ===")
    for cw in ["balanced", None]:
        mean_, sd_, vals, oofs = run_lasso(Xc, yc, cw)
        tag = "balanced" if cw else "no_cw"
        print(f"  LASSO {tag:9s}: AUC {mean_:.3f} +/- {sd_:.3f}  per-rep={[round(v,3) for v in vals]}")

# ---------- 2) 预测分数方向检查：阴性样本 OOF 分数位置 ----------
m = cohort == "PRJNA592384"
_, _, _, oofs_bal = run_lasso(Xm[m].reset_index(drop=True), ym[m], "balanced", reps=5)
oof0 = oofs_bal[0]
neg_mask = ym[m] == 0
print(f"\nPlummer rep0: 阴性样本 OOF 分数 = {np.round(oof0[neg_mask], 3)}")
print(f"Plummer rep0: 阳性样本 OOF 分数中位数 = {np.median(oof0[~neg_mask]):.3f}")
print(f"Plummer rep0: 全样本分数范围 = [{oof0.min():.3f}, {oof0.max():.3f}]")

# ---------- 3) Plummer 5 例阴性在 CV 折中的分布 ----------
neg_idx = np.where(ym[m] == 0)[0]
print(f"\nPlummer 阴性样本索引: {neg_idx}")
skf = StratifiedKFold(5, shuffle=True, random_state=42)
for k, (tr, te) in enumerate(skf.split(Xm[m], ym[m])):
    print(f"  fold{k}: test 阴性={int((ym[m][te] == 0).sum())}/5, 训练阴性={int((ym[m][tr] == 0).sum())}/5")
