# -*- coding: utf-8 -*-
"""
c19_gustin_ci_probe.py —— #C19 追加：Gustin 起源 CI 异常窄的机制探查
=================================================================
疑点：Gustin 起源的 4 个 LOCO cell bootstrap CI 全部异常窄
  （[0.485,0.563] / [0.434,0.544] / [0.449,0.540] / [0.436,0.551]）
假设：33 例训练 + 强 L1 收缩 -> 预测概率大量并列（ties） ->
      bootstrap 重采样后 AUC 变化小 -> CI 宽度被低估。
验证：Gustin->Xiao 单 cell，统计 20 种子预测概率的唯一值数量与并列结构。
"""
import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

np.random.seed(42)
SEEDS = 20

# ---- 数据加载（复刻 c19） ----
X_all = pd.read_csv('<DATA_DIR>/feat/X_features.csv')
y_all = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
feat_meta = pd.read_csv('<DATA_DIR>/feat/feature_meta.csv')
genus_cols_x = feat_meta[feat_meta['source'] == 'genus']['feature'].tolist()
xiao_mask = (y_all['cohort'] == 'PRJNA302078').values
y_xiao = y_all[xiao_mask].reset_index(drop=True)
X_xiao_log1p = X_all[xiao_mask].reset_index(drop=True)

M_rel = pd.read_csv('<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
M_pairs = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')
G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [re.sub(r'_F\.fastq\.gz$', '', str(i)) for i in G_rel.index]
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')

def norm_genus(name):
    s = str(name).strip()
    if s in ('NA_unclassified', 'Unassigned', 'nan', 'NA', 'NaN', ''):
        return 'NA_unclassified'
    s = s.replace('[', '').replace(']', '')
    if s.lower().startswith('uncultured') or s.lower().startswith('unclassified') or s == 'Bacteria':
        return 'NA_unclassified'
    first = s.split(' ')[0]
    if first in ('Sub', 'sp', 'group', 'incertae', 'Family', 'Order', 'Class', 'Phylum'):
        return 'NA_unclassified'
    return first

def agg_rel_to_log1p(rel_df):
    df = rel_df.copy()
    df.columns = [norm_genus(c) for c in df.columns]
    return np.log1p(df.T.groupby(level=0).sum().T)

X_rel = np.expm1(X_xiao_log1p[genus_cols_x]).clip(lower=0)
X_rel = pd.DataFrame(X_rel, columns=genus_cols_x, index=X_xiao_log1p.index)
X_log = agg_rel_to_log1p(X_rel)
M_log = agg_rel_to_log1p(M_rel)
G_log = agg_rel_to_log1p(G_rel)
inter = sorted(set(X_log.columns) & set(M_log.columns) & set(G_log.columns))
X_f = X_log[inter].reset_index(drop=True)
M_f = M_log.reindex(columns=inter).fillna(0)
G_f = G_log.reindex(columns=inter).fillna(0)

y_x = y_xiao['outcome'].values.astype(int)
M_v1 = M_f.loc[M_pairs['base_sample'].tolist()].reset_index(drop=True)
y_m = M_pairs['outcome'].values.astype(int)
g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)
G_base = G_f.loc[list(g_base)].reset_index(drop=True)
G_out_df = G_f.loc[list(g_out)].reset_index(drop=True)
g_outcome = (G_out_df['Lactobacillus'] >= np.log1p(0.5)).values.astype(int)

# ---- 模型协议（同 c19） ----
def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])
GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}

def fit_tuned(Xtr, ytr, seed):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    gs = GridSearchCV(make_pipe(), GRID, cv=skf, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_

print('===== Gustin->Xiao：预测概率并列结构探查 =====')
P = np.zeros((SEEDS, len(y_x)))
best_cs, n_coefs = [], []
for s in range(SEEDS):
    mdl = fit_tuned(G_base, g_outcome, seed=s)
    P[s] = mdl.predict_proba(X_f)[:, 1]
    best_cs.append(mdl.named_steps['clf'].C)
    coef = mdl.named_steps['clf'].coef_[0]
    n_coefs.append(int((coef != 0).sum()))

print(f'best_C 分布: {pd.Series(best_cs).value_counts().to_dict()}')
print(f'非零系数个数（20 种子）: {n_coefs}')
print(f'平均非零系数: {np.mean(n_coefs):.1f} / 83 属')
print()
for s in [0, 5, 10, 19]:
    nuniq = len(np.unique(P[s]))
    vals, counts = np.unique(P[s], return_counts=True)
    top = sorted(counts, reverse=True)[:3]
    print(f'种子{s}: 预测唯一值 {nuniq}/{len(y_x)}，最大并列组 {top}（各占样本数）')
print()
mean_uniq = np.mean([len(np.unique(P[s])) for s in range(SEEDS)])
print(f'20 种子平均预测唯一值: {mean_uniq:.1f} / {len(y_x)} 样本')

# 对照：Xiao->Gustin（宽 CI cell）
print('\n===== 对照 Xiao->Gustin（CI 较宽的 cell） =====')
P2 = np.zeros((SEEDS, len(g_outcome)))
for s in range(SEEDS):
    mdl = fit_tuned(X_f, y_x, seed=s)
    P2[s] = mdl.predict_proba(G_base)[:, 1]
mean_uniq2 = np.mean([len(np.unique(P2[s])) for s in range(SEEDS)])
n_coefs2 = [int((fit_tuned(X_f, y_x, seed=s).named_steps['clf'].coef_[0] != 0).sum()) for s in range(SEEDS)]
print(f'20 种子平均预测唯一值: {mean_uniq2:.1f} / {len(g_outcome)} 样本')
print(f'平均非零系数: {np.mean(n_coefs2):.1f} / 83 属')

# 空模型检查：非零系数是否为 0（空模型 -> 预测恒定 -> 唯一值=1）
print()
print('判定：若 Gustin 起源非零系数≈0 且唯一值极少 -> CI 窄系并列/退化所致，')
print('     应在框架脚注说明「CI 宽度低估不确定性」，避免审稿人质疑窄 CI=高精度。')
