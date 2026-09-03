# -*- coding: utf-8 -*-
"""
resolution_matrix.py —— 任务 1：种级 vs 属级分辨率对照（队列内）
================================================
队列：Xiao(68 对) · Plummer(49 对) · Gustin(33 对 v1->v3)
分辨率：属级(log1p rel) vs 种级(log1p rel，未定种退化聚合为 "Genus unclassified")
人群：全人群 + 基线 CST4 亚组
协议：20 种子 LASSO（内层 5 折 CV 选 C）内部 CV —— 与 7.7 矩阵一致
输出：<OUT_DIR>/modeling/results_resolution_matrix.csv
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
OUT = '<OUT_DIR>/modeling/results_resolution_matrix.csv'
np.random.seed(RNG)


def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))


def load_feat(co, level):
    df = pd.read_csv(f'<DATA_DIR>/feat/{co}_{level}_rel.csv', index_col=0)
    df.index = [norm_sid(i) for i in df.index]
    return np.log1p(df.clip(lower=0))


# ================= 加载配对标签 =================
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
cohorts = {}
for co in ['PRJNA302078', 'PRJNA592384']:
    sub = y[y['cohort'] == co].reset_index(drop=True)
    cohorts[co] = dict(
        base=[norm_sid(s) for s in sub['base_sample']],
        y=sub['outcome'].values.astype(int),
        cst4=sub['base_cst4'].values.astype(bool))

# Gustin 33 对（v1->v3）
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
g_pairs = []
for pat, s in G_map.groupby('patient'):
    v1 = s[s['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = s[s['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)
cohorts['PRJNA691964'] = dict(base=list(g_base), y=None, cst4=None)


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
            if len(np.unique(yv[te])) < 2:
                continue
            mdl = fit_tuned(X.iloc[tr], yv[tr], seed=s)
            p = mdl.predict_proba(X.iloc[te])[:, 1]
            aucs.append(roc_auc_score(yv[te], p))
    if not aucs:
        return np.nan, np.nan, len(yv)
    return np.mean(aucs), np.std(aucs), len(yv)


def lb_log1p(F, level):
    """Lactobacillus 属级/种级聚合后的 log1p 相对丰度（从 log1p 域还原加和再 log1p）"""
    if level == 'genus':
        cols = [c for c in F.columns if c == 'Lactobacillus']
    else:
        cols = [c for c in F.columns if c.startswith('Lactobacillus ')]
    if not cols:
        return pd.Series(0.0, index=F.index)
    rel = np.expm1(F[cols]).clip(lower=0).sum(axis=1)
    return np.log1p(rel)


rows = []


def rec(tag, mu, sd, n, note=''):
    rows.append(dict(analysis=tag, auc_mean=round(float(mu), 4), auc_sd=round(float(sd), 4),
                     n=int(n), note=note))
    print(f'  {tag}: AUC {mu:.3f} ± {sd:.3f} (n={n}) {note}')


# ================= 主循环：队列 x 分辨率 x 人群 =================
for co, info in cohorts.items():
    short = {'PRJNA302078': 'Xiao', 'PRJNA592384': 'Plummer', 'PRJNA691964': 'Gustin'}[co]
    for level in ['genus', 'species']:
        F = load_feat(co, level)
        base = F.loc[info['base']].reset_index(drop=True)
        if co == 'PRJNA691964':
            G_out = F.loc[list(g_out)].reset_index(drop=True)
            lb_v1 = lb_log1p(base, level)
            lb_v3 = lb_log1p(G_out, level)
            yv = (lb_v3 > np.log1p(0.5)).values.astype(int)
            cst4 = (lb_v1 < np.log1p(0.5)).values.astype(bool)
        else:
            yv = info['y']
            cst4 = info['cst4']
        mu, sd, n = internal_cv(base, yv)
        rec(f'{short} {level} full', mu, sd, n)
        if cst4.sum() >= 20:
            mu, sd, n = internal_cv(base[cst4], yv[cst4])
            rec(f'{short} {level} CST4', mu, sd, n)
        else:
            print(f'  {short} {level} CST4: n={cst4.sum()} 跳过')

# ================= 跨队列种级外推（Xiao -> Gustin） =================
print('\n===== 跨队列种级外推 =====')
X_sp = load_feat('PRJNA302078', 'species')
G_sp = load_feat('PRJNA691964', 'species')
inter = sorted(set(X_sp.columns) & set(G_sp.columns))
print(f'种级特征交集: {len(inter)}（Xiao {X_sp.shape[1]} / Gustin {G_sp.shape[1]}）')
if len(inter) >= 40:
    x = cohorts['PRJNA302078']
    Xf = X_sp.loc[x['base']][inter].reset_index(drop=True)
    gb = G_sp.loc[list(g_base)][inter].reset_index(drop=True)
    go = G_sp.loc[list(g_out)][inter].reset_index(drop=True)
    yv_g = (lb_log1p(go, 'species') > np.log1p(0.5)).values.astype(int)
    aucs = []
    for s in range(SEEDS):
        mdl = fit_tuned(Xf, x['y'], seed=s)
        p = mdl.predict_proba(gb)[:, 1]
        aucs.append(roc_auc_score(yv_g, p))
    aucs = np.array(aucs)
    rec('Xiao species -> Gustin species external', aucs.mean(), aucs.std(), len(yv_g), f'{len(inter)} 种交集')
else:
    print('  交集过小，跨队列种级外推跳过')

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print(f'\n已保存: {OUT}')
print(res.to_string(index=False))
