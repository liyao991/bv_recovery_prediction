# -*- coding: utf-8 -*-
"""
c19_loco_matrix.py —— #C19 leave-one-cohort-out 双向交叉验证全矩阵
=================================================================
派单：ASUS #C19（留言板 2026-09-03 16:50）
目的：证明阴性结论不依赖选哪个队列当训练集（堵"只试了一个训练队列"质疑）。

补 4 个方向 cell（每 cell A/B 两起源）：
  Mollin→Gustin / Mollin→Xiao / Gustin→Mollin / Gustin→Xiao
另同协议重跑 Xiao→Mollin / Xiao→Gustin（A2/A3/B2/B3）作 sanity check。

口径与 four_cohort_matrix.py 完全一致：
  - 特征 = 三队列交集 83 属，相对丰度 log1p
  - LASSO-LR（liblinear + l1_ratio=1.0，#C18 定案 = 真 L1），
    class_weight='balanced'，内层 StratifiedKFold(5) GridSearchCV C∈{0.01,0.05,0.1,0.5,1.0}
  - 20 种子（s=0..19）external AUC 的 seed-mean ± SD
  - CI = 测试集 bootstrap 5000 次，对"种子平均 AUC"统计量取 2.5/97.5 分位（#C14 定论口径）
  - A 系：训练/验证均全人群；B 系：训练/验证均基线 CST4 亚组
  - 结局 = 随访 LB≥50%（线性）；Gustin 用 >= log1p(0.5)（#C17 验证与 y_labels 逐例一致）

输出：<OUT_DIR>/modeling/results_c19_loco_matrix.csv
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
N_BOOT = 5000
OUT = '<OUT_DIR>/modeling/results_c19_loco_matrix.csv'

# ================= 1. 数据加载（复刻 four_cohort_matrix.py） =================
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

# ================= 2. 属名标准化 + log1p =================
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
print(f'三队列属名交集: {len(inter)}')

X_f = X_log[inter].reset_index(drop=True)
M_f = M_log.reindex(columns=inter).fillna(0)
G_f = G_log.reindex(columns=inter).fillna(0)

# ================= 3. 样本与标签 =================
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
G_out_df = G_f.loc[list(g_out)].reset_index(drop=True)
g_outcome = (G_out_df['Lactobacillus'] >= np.log1p(0.5)).values.astype(int)   # 50% 口径（#C17 逐例一致）
g_cst4 = (G_base['Lactobacillus'] < np.log1p(0.5)).values.astype(bool)

print(f'Xiao: n={len(y_x)}, 康复率 {y_x.mean():.3f}, CST4 {cst4_x.sum()}')
print(f'Mollin: n={len(y_m)}, 康复率 {y_m.mean():.3f}, CST4 {cst4_m.sum()}')
print(f'Gustin: n={len(g_outcome)}, 康复率 {g_outcome.mean():.3f}, CST4 {g_cst4.sum()}')

# ================= 4. 队列数据注册表 =================
DATA = {
    'Xiao':   {'full': (X_f, y_x),                    'cst4': (X_f[cst4_x].reset_index(drop=True), y_x[cst4_x])},
    'Mollin': {'full': (M_v1, y_m),                   'cst4': (M_v1[cst4_m].reset_index(drop=True), y_m[cst4_m])},
    'Gustin': {'full': (G_base, g_outcome),           'cst4': (G_base[g_cst4].reset_index(drop=True), g_outcome[g_cst4])},
}

# ================= 5. 模型协议（与主分析逐字一致） =================
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

def external_with_ci(Xtr, ytr, Xte, yte, tag, scope, n_train, n_test):
    """20 种子 seed-mean ± SD + 测试集 bootstrap CI（对种子平均 AUC 统计量）"""
    P = np.zeros((SEEDS, len(yte)))
    seed_aucs = []
    for s in range(SEEDS):
        mdl = fit_tuned(Xtr, ytr, seed=s)
        p = mdl.predict_proba(Xte)[:, 1]
        P[s] = p
        seed_aucs.append(roc_auc_score(yte, p))
    seed_aucs = np.array(seed_aucs)
    sm_pt, sm_sd = seed_aucs.mean(), seed_aucs.std()

    rng = np.random.default_rng(42)
    sm_boots = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(yte), size=len(yte))
        sm_boots.append(np.mean([roc_auc_score(yte[idx], P[s][idx]) for s in range(SEEDS)]))
    lo, hi = np.percentile(sm_boots, [2.5, 97.5])

    row = dict(analysis=tag, train_scope=scope, n_train=n_train, n_test=n_test,
               auc=round(float(sm_pt), 4), auc_sd=round(float(sm_sd), 4),
               ci_lo=round(float(lo), 3), ci_hi=round(float(hi), 3))
    print(f"  {tag:38s} AUC {sm_pt:.3f} ± {sm_sd:.3f}  CI [{lo:.3f}, {hi:.3f}]  (train {n_train} / test {n_test})")
    return row

# ================= 6. 跑全矩阵 =================
rows = []

print('\n===== sanity check：Xiao 起源（对照锚点 A2 0.521 / A3 0.556 / B2 0.501 / B3 0.512）=====')
for scope, sys_tag in [('full', 'A'), ('cst4', 'B')]:
    Xtr, ytr = DATA['Xiao'][scope]
    for test in ['Mollin', 'Gustin']:
        Xte, yte = DATA[test][scope]
        cell = f'{sys_tag}: Xiao->{test} ({scope})'
        rows.append(external_with_ci(Xtr, ytr, Xte, yte, cell, scope, len(ytr), len(yte)))

print('\n===== 新增 cell：Mollin / Gustin 起源 =====')
for train in ['Mollin', 'Gustin']:
    for scope, sys_tag in [('full', 'A'), ('cst4', 'B')]:
        Xtr, ytr = DATA[train][scope]
        for test in ['Xiao', 'Mollin', 'Gustin']:
            if test == train:
                continue
            Xte, yte = DATA[test][scope]
            cell = f'{sys_tag}: {train}->{test} ({scope})'
            rows.append(external_with_ci(Xtr, ytr, Xte, yte, cell, scope, len(ytr), len(yte)))

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print(f'\n已保存: {OUT}')

# ================= 7. 汇总透视 =================
print('\n===== LOCO 矩阵总览（AUC [CI]）=====')
for scope in ['full', 'cst4']:
    sub = res[res['train_scope'] == scope]
    print(f'\n-- 起源口径: {scope} --')
    for _, r in sub.iterrows():
        flag = ' ←新' if not r['analysis'].startswith(('A: Xiao', 'B: Xiao')) else '  (sanity)'
        print(f"  {r['analysis']:38s} {r['auc']:.3f} [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]{flag}")
