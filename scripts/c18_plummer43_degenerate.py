# -*- coding: utf-8 -*-
"""c18_plummer43_degenerate.py —— diag_plummer_within 0.431 的空模型折比例实测
协议复刻：10 reps x 5 折外层 + 内层 3 折网格(7C)，OOF 聚合 AUC。
统计：每外层折预测是否常数（空模型）、best_C 分布。"""
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
m = (y["cohort"] == "PRJNA592384").values
Xd = X[feats][m].reset_index(drop=True)
yv = y["outcome"].values[m]
print(f'Plummer 全人群: n={len(yv)}, pos={int(yv.sum())}, neg={int((1 - yv).sum())}')

pipe = Pipeline([("sc", StandardScaler()),
                 ("clf", LogisticRegression(solver="liblinear", l1_ratio=1.0,
                                            class_weight="balanced", max_iter=5000, random_state=42))])
grid = {"clf__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]}

n_const, n_fold, bestcs = 0, 0, []
rep_aucs = []
for rep in range(10):
    skf = StratifiedKFold(5, shuffle=True, random_state=42 + rep)
    inner = StratifiedKFold(3, shuffle=True, random_state=42 + rep)
    oof = np.zeros(len(yv))
    for tr, te in skf.split(Xd, yv):
        gs = GridSearchCV(pipe, grid, cv=inner, scoring="roc_auc", n_jobs=1, refit=True)
        gs.fit(Xd.iloc[tr], yv[tr])
        oof[te] = gs.predict_proba(Xd.iloc[te])[:, 1]
        n_fold += 1
        if len(np.unique(oof[te])) == 1:
            n_const += 1
        bestcs.append(gs.best_params_["clf__C"])
    rep_aucs.append(roc_auc_score(yv, oof))

print(f'OOF AUC (10 reps): {np.mean(rep_aucs):.3f} ± {np.std(rep_aucs):.3f}')
print(f'总外层折数 {n_fold}, 其中常数预测(空模型)折 {n_const} ({100 * n_const / n_fold:.0f}%)')
print(f'best_C 分布: {dict(pd.Series(bestcs).value_counts().sort_index())}')
