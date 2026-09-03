# -*- coding: utf-8 -*-
"""
c7_panel_combo.py —— #C7 靶标精简子集组合检验
=================================================
检验「几个 PCR+CE 靶标物种的属级代理组合」是否优于单物种/全集（预期均 ≈0.5）。

队列：Xiao(PRJNA302078) · Plummer(PRJNA592384) · Gustin(PRJNA691964)
口径（与 four_cohort_matrix.py / species_panel_scan.py 一致）：
  - 基线 CST4 亚组 = base_LB < 0.5（Gustin 全队列即 CST4）
  - 二元结局 = 随访 LB >= 50% (康复=1)；Gustin 用修正口径 LB>=0.5
  - 属级相对丰度 -> log1p（属名归一 norm_genus），三队列取交集 = 65 属
三层对照（LASSO-LR，20 种子 × 5 折分层 CV 内层选 C）：
  ① 单物种（属级代理）：Gardnerella / Atopobium / Megasphaera / Lactobacillus
  ② 组合：上述 4 属；加 Prevotella(聚合) 的 5 属备选
  ③ 全集：三队列属名交集 65 属

输出：<OUT_DIR>/modeling/results_c7_panel_combo.csv
"""
import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

RNG = 42
SEEDS = 20
OUT = '<OUT_DIR>/modeling/results_c7_panel_combo.csv'
np.random.seed(RNG)


def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))


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


# ================= 1. 加载三队列属级矩阵（norm_genus 归一 -> log1p） =================
def load_genus_log1p(pid):
    g = pd.read_csv(f'<DATA_DIR>/feat/{pid}_genus_rel.csv', index_col=0)
    g.index = [norm_sid(i) for i in g.index]
    g.columns = [norm_genus(c) for c in g.columns]
    agg = g.T.groupby(level=0).sum().T
    return np.log1p(agg)


G = {c: load_genus_log1p(p) for c, p in
     [('Xiao', 'PRJNA302078'), ('Plummer', 'PRJNA592384'), ('Gustin', 'PRJNA691964')]}

inter = sorted(set(G['Xiao'].columns) & set(G['Plummer'].columns) & set(G['Gustin'].columns))
print(f'三队列属名交集: {len(inter)} 属')
print('核心靶标属是否在交集:', {g: (g in inter) for g in
      ['Lactobacillus', 'Gardnerella', 'Atopobium', 'Megasphaera']})

# Prevotella 聚合（所有 Prevotella* 列求和）
def prevotella_cols(cols):
    return [c for c in cols if c.startswith('Prevotella')]

# ================= 2. 构建 subject 表（CST4 亚组 + 结局） =================
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
G_map = G_map[G_map['control'].isna()]

# 各队列 base_run -> (X_row, y)
data = {}  # cohort -> (X_df [CST4], y_vec)

# Xiao
xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
xiao = xiao[xiao['base_cst4'] == 1]
X_x = G['Xiao'].reindex(columns=inter).fillna(0).loc[xiao['base'].tolist()]
data['Xiao'] = (X_x.reset_index(drop=True), xiao['outcome'].values.astype(int))

# Plummer
plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
plum = plum[plum['base_cst4'] == 1]
X_p = G['Plummer'].reindex(columns=inter).fillna(0).loc[plum['base'].tolist()]
data['Plummer'] = (X_p.reset_index(drop=True), plum['outcome'].values.astype(int))

# Gustin（全队列即 CST4；v1 基线，结局 LB>=0.5）
pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
gb, go = zip(*pairs)
X_g = G['Gustin'].reindex(columns=inter).fillna(0).loc[list(gb)]
g_outcome = (G['Gustin'].loc[list(go), 'Lactobacillus'] >= np.log1p(0.5)).values.astype(int)
# 基线 CST4 筛选（LB<0.5）
g_cst4 = (G['Gustin'].loc[list(gb), 'Lactobacillus'] < np.log1p(0.5)).values
X_g = X_g[g_cst4].reset_index(drop=True)
g_outcome = g_outcome[g_cst4]
data['Gustin'] = (X_g, g_outcome)

print('\n=== 三队列 CST4 亚组 n / 康复率 ===')
for c in ['Xiao', 'Plummer', 'Gustin']:
    Xd, yv = data[c]
    print(f'  {c}: n={len(yv)}, rec={yv.sum()}, rate={yv.mean():.2f}')

# ================= 3. LASSO-LR internal CV =================
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


def internal_cv(X, yv, seeds=SEEDS):
    aucs = []
    for s in range(seeds):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, yv):
            if len(np.unique(yv[tr])) < 2 or len(np.unique(yv[te])) < 2:
                continue
            mdl = fit_tuned(X.iloc[tr], yv[tr], seed=s)
            p = mdl.predict_proba(X.iloc[te])[:, 1]
            aucs.append(roc_auc_score(yv[te], p))
    return np.mean(aucs), np.std(aucs), len(aucs)


# ================= 4. 三层对照 =================
TARGETS = ['Gardnerella', 'Atopobium', 'Megasphaera', 'Lactobacillus']

rows = []
for cohort in ['Xiao', 'Plummer', 'Gustin']:
    Xd, yv = data[cohort]
    # ① 单物种
    for t in TARGETS:
        mu, sd, nf = internal_cv(Xd[[t]], yv)
        rows.append(dict(cohort=cohort, layer='① 单物种', features=t,
                         auc=round(mu, 4), sd=round(sd, 4), n=len(yv), n_folds=nf))
    # ② 组合 4 / 5（+Prevotella）
    pv = prevotella_cols(Xd.columns)
    X4 = Xd[TARGETS]
    mu, sd, nf = internal_cv(X4, yv)
    rows.append(dict(cohort=cohort, layer='② 组合4', features='+'.join(TARGETS),
                     auc=round(mu, 4), sd=round(sd, 4), n=len(yv), n_folds=nf))
    if pv:
        Xd5 = Xd[TARGETS]
        Xd5 = pd.concat([Xd5, Xd[pv].sum(axis=1).rename('Prevotella_agg')], axis=1)
        mu, sd, nf = internal_cv(Xd5, yv)
        rows.append(dict(cohort=cohort, layer='② 组合5(+Prev)', features='4+Prevotella',
                         auc=round(mu, 4), sd=round(sd, 4), n=len(yv), n_folds=nf))
    # ③ 全集 65 属
    mu, sd, nf = internal_cv(Xd, yv)
    rows.append(dict(cohort=cohort, layer='③ 全集', features=f'{len(inter)}属',
                     auc=round(mu, 4), sd=round(sd, 4), n=len(yv), n_folds=nf))

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print('\n=== #C7 靶标组合检验（LASSO-LR 20种子×5折CV AUC）===')
print(res.to_string(index=False))

# 汇总：三队列均值
print('\n=== 三层对照三队列均值（单物种取 4 属均值）===')
pivot = res.pivot_table(index='features', columns='layer', values='auc', aggfunc='mean')
print(pivot.round(4).to_string())
print(f'\n已保存: {OUT}')
