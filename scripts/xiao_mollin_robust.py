# -*- coding: utf-8 -*-
"""
Xiao->Mollin 外推稳健性验证
====================================
1. 基线 CST4 亚组（剔除"基线健康"样本，压缩 LB 耦合）
2. 多种子重复（验证 0.681 非偶然）
3. 预测分布检查（康复 vs 未康复分层质量）
4. 去掉 Lactobacillus 特征后的 AUC（量化标签耦合程度）
"""
import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.pipeline import Pipeline

# ---------- 数据加载与特征构建（复用 xiao_mollin_cross 逻辑） ----------
X_all = pd.read_csv('<DATA_DIR>/feat/X_features.csv')
y_all = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
feat_meta = pd.read_csv('<DATA_DIR>/feat/feature_meta.csv')
genus_cols_x = feat_meta[feat_meta['source'] == 'genus']['feature'].tolist()

xiao_mask = (y_all['cohort'] == 'PRJNA302078').values
X_xiao = X_all[xiao_mask].reset_index(drop=True)
y_xiao = y_all[xiao_mask].reset_index(drop=True)

M_rel = pd.read_csv('<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
M_pairs = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')


def norm_genus(name):
    s = str(name).strip()
    if s in ('NA_unclassified', 'nan', 'NA', 'NaN'):
        return 'NA_unclassified'
    s = s.replace('[', '').replace(']', '')
    if s.lower().startswith('uncultured') or s.lower().startswith('unclassified') or s == 'Bacteria':
        return 'NA_unclassified'
    first = s.split(' ')[0]
    if first in ('Sub', 'sp', 'group', 'incertae'):
        return 'NA_unclassified'
    return first


xiao_genus_map = {c: norm_genus(c) for c in genus_cols_x}
mollin_genus_map = {c: norm_genus(c) for c in M_rel.columns}

X_rel = np.expm1(X_xiao[genus_cols_x]).clip(lower=0)
X_rel = pd.DataFrame(X_rel, columns=genus_cols_x, index=X_xiao.index)


def agg_rel_to_log1p(rel_df, col_map):
    df = rel_df.copy()
    df.columns = [col_map[c] for c in df.columns]
    agg = df.T.groupby(level=0).sum().T
    return np.log1p(agg)


xiao_agg = agg_rel_to_log1p(X_rel, xiao_genus_map)
mollin_agg = agg_rel_to_log1p(M_rel, mollin_genus_map)
inter = sorted(set(xiao_agg.columns) & set(mollin_agg.columns))
Xtr = xiao_agg[inter].reset_index(drop=True)
Mfeat = mollin_agg.reindex(columns=inter).fillna(0)
M_Xi = Mfeat.set_index(M_rel.index)
mollin_v1 = M_pairs[M_pairs['base_sample'].isin(M_Xi.index)]
M_v1 = M_Xi.loc[mollin_v1['base_sample'].tolist()].reset_index(drop=True)
M_y = mollin_v1['outcome'].reset_index(drop=True).values
ytr = y_xiao['outcome'].values
print(f'交集属: {len(inter)} | Xiao {Xtr.shape} | Mollin 测试 {M_v1.shape}, 阳性率 {M_y.mean():.3f}')

# ---------- 1. 多种子稳健性（全量 + 基线 CST4 亚组） ----------
print('\n=== LASSO 20 种子重复：Xiao 训练 -> Mollin 外推 ===')
for label, tr_mask, te_mask in [
    ('全量', np.ones(len(Xtr), bool), np.ones(len(M_v1), bool)),
    ('基线CST4亚组', y_xiao['base_cst4'].values == 1, mollin_v1['base_cst4'].values == 1),
]:
    Xt, yt = Xtr[tr_mask], ytr[tr_mask]
    Xe, ye = M_v1[te_mask], M_y[te_mask]
    if len(ye) == 0:
        continue
    aucs, prs, briers = [], [], []
    for seed in range(20):
        pipe = Pipeline([('sc', StandardScaler()),
                         ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                    class_weight='balanced', max_iter=5000))])
        skf = StratifiedKFold(5, shuffle=True, random_state=seed)
        gs = GridSearchCV(pipe, {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}, cv=skf,
                          scoring='roc_auc', n_jobs=1)
        gs.fit(Xt, yt)
        p = gs.predict_proba(Xe)[:, 1]
        aucs.append(roc_auc_score(ye, p))
        prs.append(average_precision_score(ye, p))
        briers.append(brier_score_loss(ye, p))
    print(f'  {label}: 训练 n={len(yt)} (阳性 {yt.sum()}) | 测试 n={len(ye)} (阳性 {ye.sum()})')
    print(f'    AUC {np.mean(aucs):.3f}±{np.std(aucs):.3f} [min {np.min(aucs):.3f}, max {np.max(aucs):.3f}]')
    print(f'    PR-AUC {np.mean(prs):.3f}±{np.std(prs):.3f} | Brier {np.mean(briers):.3f}±{np.std(briers):.3f}')

# ---------- 2. 去掉 Lactobacillus 特征（量化标签耦合） ----------
print('\n=== 去掉 Lactobacillus 后的外推（标签耦合检验）===')
inter_nolb = [c for c in inter if c != 'Lactobacillus']
Xtr_nolb = Xtr[inter_nolb]
M_nolb = M_v1[inter_nolb]
aucs = []
for seed in range(20):
    pipe = Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    gs = GridSearchCV(pipe, {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}, cv=skf,
                      scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr_nolb, ytr)
    p = gs.predict_proba(M_nolb)[:, 1]
    aucs.append(roc_auc_score(M_y, p))
print(f'  无 LB: AUC {np.mean(aucs):.3f}±{np.std(aucs):.3f} (vs 含 LB 的 {0.681:.3f})')

# ---------- 3. 预测分布（LASSO 最优 C，全量训练） ----------
print('\n=== LASSO 预测分数分层（康复 vs 未康复，Mollin）===')
pipe = Pipeline([('sc', StandardScaler()),
                 ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                            class_weight='balanced', max_iter=5000))])
gs = GridSearchCV(pipe, {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}, cv=StratifiedKFold(5, shuffle=True, random_state=42),
                  scoring='roc_auc', n_jobs=1)
gs.fit(Xtr, ytr)
p = gs.predict_proba(M_v1)[:, 1]
df = pd.DataFrame({'patient': mollin_v1['patient'].values, 'outcome': M_y,
                   'pred': p, 'meta_class': mollin_v1['meta_class_base'].values})
print(df.groupby('outcome')['pred'].describe().round(3).to_string())
print()
print('按临床硬结局分组的预测分数:')
print(df.groupby('meta_class')['pred'].agg(['count', 'mean', 'median']).round(3).to_string())
print()
# 预测分数与 LB 相关性
print('预测分数 vs 基线 LB 丰度 (Pearson r):')
from scipy.stats import pearsonr
lb_vals = M_Xi.loc[mollin_v1['base_sample'].tolist(), 'Lactobacillus'].values if 'Lactobacillus' in M_Xi.columns else None
if lb_vals is not None:
    r, pv = pearsonr(p, lb_vals)
    print(f'  r={r:.3f}, p={pv:.4f}')
print()
print('top 5 预测特征（LASSO 系数绝对值，标准缩放后）:')
coef = gs.best_estimator_.named_steps['clf'].coef_[0]
sc = gs.best_estimator_.named_steps['sc']
# 用标准化后的系数排序
scaled_coef = coef  # 系数已经是在标准化特征上的
idx = np.argsort(-np.abs(scaled_coef))[:8]
for i in idx:
    print(f'  {inter[i]}: {scaled_coef[i]:+.3f}')
