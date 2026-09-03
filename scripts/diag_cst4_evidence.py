# -*- coding: utf-8 -*-
"""
diag_cst4_evidence.py v2 - 正式检验"治疗前 16S 属级特征与治疗结局相关性有限"。
v1 的 bootstrap+CV 有泄漏（重采样重复样本跨折）→ 弃用。
v2 用无泄漏方法：
  1) 重复分层 holdout（200 次，60/40，LASSO C=1.0）：Xiao all + 基线 CST4 亚组
     → AUC 分布 mean / 95% CI / P(AUC>0.5)（训练-测试无重叠，无泄漏）
  2) 去 Lactobacillus 特征对照（残余信号是否全来自 LB 耦合）
  3) 基线 LB vs 结局 Spearman（range-restricted 内的直接相关）
输出：<OUT_DIR>/modeling/results_cst4_evidence.csv
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

FEAT = "<DATA_DIR>/feat"
X = pd.read_csv(FEAT + "/X_features.csv")
y = pd.read_csv(FEAT + "/y_labels.csv")
feat = pd.read_csv(FEAT + "/feature_meta.csv")
genus_cols = feat.loc[feat["source"] == "genus", "feature"].tolist()
alpha_cols = feat.loc[feat["source"] == "alpha", "feature"].tolist()
keep = [c for c in genus_cols if (X[c] > 0).mean() >= 0.05]
feats = keep + alpha_cols
LB_LIKE = [c for c in keep if "Lactobacillus" in c or "Lactiplantibacillus" in c or "Limosilactobacillus" in c]
print(f"keep genera: {len(keep)}, LB-like: {LB_LIKE}")

def holdout_auc(Xd, yv, feats_use, n_rep=200, seed=42, test_size=0.4):
    rng = np.random.default_rng(seed)
    aucs = []
    for r in range(n_rep):
        tr, te = train_test_split(np.arange(len(yv)), test_size=test_size, stratify=yv,
                                  random_state=int(rng.integers(1e9)))
        pipe = Pipeline([("sc", StandardScaler()),
                         ("clf", LogisticRegression(solver="liblinear", C=1.0, l1_ratio=1.0,
                                                    class_weight="balanced", max_iter=5000,
                                                    random_state=seed))])
        pipe.fit(Xd.iloc[tr][feats_use], yv[tr])
        try:
            aucs.append(roc_auc_score(yv[te], pipe.predict_proba(Xd.iloc[te][feats_use])[:, 1]))
        except ValueError:
            pass
    aucs = np.array(aucs)
    return aucs.mean(), np.percentile(aucs, 2.5), np.percentile(aucs, 97.5), (aucs > 0.5).mean()

rows = []
m = y["cohort"] == "PRJNA302078"
for subset, mask in [("all", m), ("base_cst4", m & (y["base_cst4"] == 1))]:
    Xd, yv = X.loc[mask].reset_index(drop=True), y.loc[mask, "outcome"].values
    n_pos, n_neg = int(yv.sum()), int((1 - yv).sum())
    if n_pos < 2 or n_neg < 2:
        continue
    full_feats = [c for c in feats if c in Xd.columns]
    no_lb = [c for c in full_feats if c not in LB_LIKE]
    m1, lo1, hi1, p1 = holdout_auc(Xd, yv, full_feats)
    m2, lo2, hi2, p2 = holdout_auc(Xd, yv, no_lb)
    lb = y.loc[mask, "base_LB"].values
    rho, p_rho = spearmanr(lb, yv)
    rows.append({"cohort": "PRJNA302078", "subset": subset, "n": len(yv),
                 "pos": n_pos, "neg": n_neg,
                 "AUC_all_feat": round(m1, 3), "CI95": f"[{lo1:.3f},{hi1:.3f}]", "P_gt_0.5": round(p1, 3),
                 "AUC_no_LB": round(m2, 3), "CI95_noLB": f"[{lo2:.3f},{hi2:.3f}]", "P_gt_0.5_noLB": round(p2, 3),
                 "Spearman_LB": round(rho, 3), "p_spear": round(p_rho, 4)})
    print(f"[{subset}] n={len(yv)} ({n_pos}+/{n_neg}-)  "
          f"AUC(all) {m1:.3f} CI[{lo1:.3f},{hi1:.3f}] P(>0.5)={p1:.3f} | "
          f"AUC(noLB) {m2:.3f} CI[{lo2:.3f},{hi2:.3f}] P(>0.5)={p2:.3f} | "
          f"Spearman(LB,outcome)={rho:.3f} (p={p_rho:.4f})")

out = pd.DataFrame(rows)
out.to_csv("<OUT_DIR>/modeling/results_cst4_evidence.csv", index=False)
print("\nsaved -> <OUT_DIR>/modeling/results_cst4_evidence.csv")
