# -*- coding: utf-8 -*-
"""
alpha_family.py —— α 多样性家族对照（家族 × 分辨率矩阵）
================================================================
队列：Xiao(68) · Plummer(49) · Gustin(33 v1->v3)，均为配对结局（D30 LB>=50% 康复）
α 家族：
  richness : observed, chao1
  evenness : pielou, simpson_even(=invsimpson/observed)
  entropy  : shannon, simpson(1-D), invsimpson
  full     : 全部 7 指标
分辨率：asv / genus / species（三级聚合分别计算 α）
问题：
  1. 哪个 α 家族最有预测信号？（对照 C1 混合 0.689 与论文"richness 高->失败"）
  2. 分辨率提升（ASV->属->种）对 α 信号有何影响？（构成特征上已证"种不优于属"）
  3. 与属级构成特征（LASSO 0.53-0.63）对照：失调程度 vs 构成的预测力对比
方法：20 种子 × 分层 5 折 CV，LogisticRegression（标准化 + 类平衡）
输出：<OUT_DIR>/modeling/results_alpha_family.csv
"""
import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

RNG = 42
SEEDS = 20
OUT = '<OUT_DIR>/modeling/results_alpha_family.csv'
np.random.seed(RNG)

# ================= 1. 加载 α 家族 =================
def norm_sid(s):
    return re.sub(r'_F\.fastq\.gz$', '', str(s))

af = {}
for co in ['PRJNA302078', 'PRJNA592384', 'PRJNA691964']:
    d = pd.read_csv(f'<DATA_DIR>/feat/{co}_alpha_family.csv')
    d['run'] = [norm_sid(s) for s in d['sample']]
    af[co] = d.set_index('run')

# ================= 2. 配对与标签 =================
y_all = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')

# Xiao（68）与 Plummer（49）
pairs = {}
for co, tag in [('PRJNA302078', 'Xiao'), ('PRJNA592384', 'Plummer')]:
    sub = y_all[y_all['cohort'] == co]
    base = [norm_sid(s) for s in sub['base_sample']]
    yv = sub['outcome'].values.astype(int)
    cst4 = sub['base_cst4'].values.astype(bool)
    ok = [b in af[co].index for b in base]
    pairs[tag] = dict(base=[b for b, k in zip(base, ok) if k],
                      y=yv[ok], cst4=cst4[ok], n=int(sum(ok)))
    print(f'{tag}: 配对 {len(sub)} -> α 可匹配 {pairs[tag]["n"]}, 康复率 {pairs[tag]["y"].mean():.3f}, CST4 {pairs[tag]["cst4"].sum()}')

# Gustin（33 v1->v3）
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [norm_sid(i) for i in G_rel.index]
g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)
G_out = G_rel.loc[list(g_out)]
g_outcome = (G_out['Lactobacillus'] > np.log1p(0.5)).values.astype(int)
pairs['Gustin'] = dict(base=list(g_base), y=g_outcome,
                       cst4=np.ones(len(g_outcome), dtype=bool), n=len(g_outcome))
print(f'Gustin: {len(g_pairs)} 对, 康复率 {g_outcome.mean():.3f} (基线全 CST4)')

# ================= 3. 指标与家族 =================
METRICS = ['observed', 'chao1', 'shannon', 'simpson', 'invsimpson', 'pielou']
LEVELS = ['asv', 'genus', 'species']
FAMILIES = {
    'richness': ['observed', 'chao1'],
    'evenness': ['pielou'],          # simpson_even 由 Python 派生
    'entropy':  ['shannon', 'simpson', 'invsimpson'],
    'full':     METRICS + ['simpson_even'],
}

def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])

def cv_auc(X, y):
    """20 种子 × 5 折分层 CV；折内单类返回 nan"""
    aucs = []
    for s in range(SEEDS):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, y):
            if len(np.unique(y[tr])) < 2:
                aucs.append(np.nan); continue
            mdl = make_pipe().fit(X.iloc[tr], y[tr])
            p = mdl.predict_proba(X.iloc[te])[:, 1]
            if len(np.unique(y[te])) < 2:
                aucs.append(np.nan); continue
            aucs.append(roc_auc_score(y[te], p))
    aucs = np.array(aucs); aucs = aucs[~np.isnan(aucs)]
    if len(aucs) == 0:
        return np.nan, np.nan
    return aucs.mean(), aucs.std()

rows = []
def rec(tag, mu, sd, n, note=''):
    rows.append(dict(analysis=tag, auc_mean=round(float(mu), 4) if mu == mu else np.nan,
                     auc_sd=round(float(sd), 4) if sd == sd else np.nan,
                     n_test=int(n), note=note))
    print(f'  {tag}: AUC {mu:.3f} ± {sd:.3f} (n={n}) {note}')

# ================= 4. 家族 × 分辨率矩阵 =================
for co_tag in ['Gustin', 'Xiao', 'Plummer']:
    P = pairs[co_tag]
    d = af[{'Gustin': 'PRJNA691964', 'Xiao': 'PRJNA302078', 'Plummer': 'PRJNA592384'}[co_tag]]
    base_idx = [r for r in P['base'] if r in d.index]
    if len(base_idx) != len(P['base']):
        print(f'!! {co_tag}: {len(P["base"]) - len(base_idx)} 个 base 样本无 α 数据')
    a = d.loc[base_idx].copy()
    yv = P['y']
    # 深度列（测序深度混杂参考，见报告说明）
    a['total_reads'] = 0  # 占位（α 导出未含深度，见报告说明）

    print(f'\n===== {co_tag} α 家族 × 分辨率（全人群 n={len(yv)}）=====')
    for lv in LEVELS:
        cols = {m: f'{lv}_{m}' for m in METRICS + ['simpson_even']}
        X_lv = a[[cols[m] for m in METRICS]].copy()
        X_lv[cols['simpson_even']] = X_lv[cols['invsimpson']] / X_lv[cols['observed']].clip(lower=1)
        # 单指标
        for m in METRICS:
            mu, sd = cv_auc(X_lv[[cols[m]]], yv)
            rec(f'{co_tag} {lv} single:{m}', mu, sd, len(yv))
        # 家族组合
        for fam, ms in FAMILIES.items():
            mu, sd = cv_auc(X_lv[[cols[m] for m in ms]], yv)
            rec(f'{co_tag} {lv} family:{fam}', mu, sd, len(yv))

# ================= 5. Xiao CST4 亚组（专门治疗对象）================
print('\n===== Xiao 基线 CST4 亚组 α 家族 × 分辨率（n=54）=====')
P = pairs['Xiao']
d = af['PRJNA302078']
base_idx = [r for r in P['base'] if r in d.index]
a = d.loc[base_idx].copy()
yv = P['y']; cst = P['cst4']
for lv in LEVELS:
    cols = {m: f'{lv}_{m}' for m in METRICS + ['simpson_even']}
    X_lv = a[[cols[m] for m in METRICS]].copy()
    X_lv[cols['simpson_even']] = X_lv[cols['invsimpson']] / X_lv[cols['observed']].clip(lower=1)
    Xc = X_lv[cst].reset_index(drop=True); yc = yv[cst]
    for m in METRICS:
        mu, sd = cv_auc(Xc[[cols[m]]], yc)
        rec(f'XiaoCST4 {lv} single:{m}', mu, sd, int(cst.sum()))
    for fam, ms in FAMILIES.items():
        mu, sd = cv_auc(Xc[[cols[m] for m in ms]], yc)
        rec(f'XiaoCST4 {lv} family:{fam}', mu, sd, int(cst.sum()))

# ================= 6. 保存 =================
res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print(f'\n保存: {OUT} ({len(res)} 行)')
