# -*- coding: utf-8 -*-
"""
c14_a2_ci.py —— #C14 补 A2 测试集 bootstrap CI + 定外推 AUC 口径
=================================================================
背景：#C13 F1 要求所有外推 AUC 改报「点估计 + 测试集 bootstrap 95%CI」，
但当时只给了 B2/B3/A3，缺 A2（Xiao→Mollin full，报告值 0.528±0.010）。
且 C13 Task A 用的点估计是「20 种子平均概率的 AUC」（mean-prob），
与矩阵协议报告值「20 个各种子 AUC 的均值 ± SD」（seed-mean）口径不同
（B2: 0.496 vs 0.505；A3: 0.537 vs 0.556）。

本脚本对全部 4 条外推路径（A2/A3/B2/B3）同时算两种口径：
  口径1 seed-mean：每种子一个测试集 AUC，20 种子取均值（=矩阵协议报告值）
                   CI = 测试集 bootstrap 内对"种子平均 AUC"这一统计量取分位数
  口径2 mean-prob：20 种子概率先平均，再算一个 AUC（=C13 Task A 口径）
                   CI = 测试集 bootstrap（C13 已给 B2/B3/A3，此处复核 + 补 A2）

口径与 four_cohort_matrix.py / c13_reviewer_checks.py 完全一致。
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
OUT = '<REPO_DIR>/results/results_c14_external_ci.csv'
rows = []

# ================= 数据加载（复刻 four_cohort_matrix 口径） =================
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
cst4_x = y_xiao['base_cst4'].values.astype(bool)
M_v1 = M_f.loc[M_pairs['base_sample'].tolist()].reset_index(drop=True)
y_m = M_pairs['outcome'].values.astype(int)
cst4_m = M_pairs['base_cst4'].values.astype(bool)

g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)
G_base = G_f.loc[list(g_base)].reset_index(drop=True)
G_out = G_f.loc[list(g_out)].reset_index(drop=True)
g_outcome = (G_out['Lactobacillus'] > np.log1p(0.5)).values.astype(int)
g_cst4 = (G_base['Lactobacillus'] < np.log1p(0.5)).values.astype(bool)
print(f'数据: Xiao {len(y_x)}(CST4 {cst4_x.sum()}) Mollin {len(y_m)}(CST4 {cst4_m.sum()}) '
      f'Gustin {len(g_outcome)}(CST4 {g_cst4.sum()}) 交集 {len(inter)} 属')

# ================= 模型 =================
def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])
GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}
def fit_tuned(Xtr, ytr, seed):
    gs = GridSearchCV(make_pipe(), GRID,
                      cv=StratifiedKFold(5, shuffle=True, random_state=seed),
                      scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_

# ================= 双口径外推评估 =================
def external_dual(Xtr, ytr, Xte, yte, tag, seeds=20, nboot=5000):
    Xte = Xte.reset_index(drop=True)
    yte = np.asarray(yte).astype(int)
    n = len(yte)
    P = np.zeros((seeds, n))  # 每种子概率矩阵
    for s in range(seeds):
        mdl = fit_tuned(Xtr, ytr, seed=s)
        P[s] = mdl.predict_proba(Xte)[:, 1]

    # 口径1: seed-mean（每种子 AUC 取均值，=矩阵协议报告值）
    per_seed = np.array([roc_auc_score(yte, P[s]) for s in range(seeds)])
    sm_pt, sm_sd = per_seed.mean(), per_seed.std(ddof=1)
    # 口径2: mean-prob（概率先平均，=C13 Task A 口径）
    mp_pt = roc_auc_score(yte, P.mean(axis=0))

    rng = np.random.RandomState(123)
    sm_boots, mp_boots = [], []
    for _ in range(nboot):
        idx = rng.randint(0, n, n)
        if len(np.unique(yte[idx])) < 2:
            continue
        # seed-mean 统计量的 bootstrap：重采样内每种子 AUC 再平均
        sm_boots.append(np.mean([roc_auc_score(yte[idx], P[s][idx]) for s in range(seeds)]))
        mp_boots.append(roc_auc_score(yte[idx], P.mean(axis=0)[idx]))
    sm_lo, sm_hi = np.percentile(sm_boots, [2.5, 97.5])
    mp_lo, mp_hi = np.percentile(mp_boots, [2.5, 97.5])

    print(f'\n[{tag}] n_test={n}')
    print(f'  口径1 seed-mean : {sm_pt:.3f} ± {sm_sd:.3f}   bootstrap CI [{sm_lo:.3f}, {sm_hi:.3f}]')
    print(f'  口径2 mean-prob : {mp_pt:.3f}            bootstrap CI [{mp_lo:.3f}, {mp_hi:.3f}]')
    rows.append(dict(path=tag, n_test=n,
                     seedmean_auc=round(sm_pt, 4), seedmean_sd=round(sm_sd, 4),
                     seedmean_ci_lo=round(sm_lo, 3), seedmean_ci_hi=round(sm_hi, 3),
                     meanprob_auc=round(mp_pt, 4),
                     meanprob_ci_lo=round(mp_lo, 3), meanprob_ci_hi=round(mp_hi, 3)))

print('\n===== 4 条外推路径 × 双口径 =====')
external_dual(X_f, y_x, M_v1, y_m, 'A2 Xiao->Mollin full (报告 0.528±0.010)')
external_dual(X_f, y_x, G_base, g_outcome, 'A3 Xiao->Gustin full (报告 0.556±0.021)')
external_dual(X_f[cst4_x], y_x[cst4_x], M_v1[cst4_m], y_m[cst4_m], 'B2 Xiao->Mollin CST4 (报告 0.505±0.020)')
external_dual(X_f[cst4_x], y_x[cst4_x], G_base[g_cst4], g_outcome[g_cst4], 'B3 Xiao->Gustin CST4 (报告 0.512)')

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f'\n已保存: {OUT}')
