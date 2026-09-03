# -*- coding: utf-8 -*-
"""
d7_dynamic.py —— 任务 2：Gustin D7 动态监测（Xiao 原文"D0 弱 / D7 有价值"的延续）
================================================
配对：v1(D0) -> v2(D7) -> v3(D30)；结局 = v3 LB>=50% 康复
特征组（属级 log1p）：
  D0 only / D7 only / D0+D7 拼接 / Δ(D7-D0)
附加：D7 CST 状态转换（基线 CST4 -> D7 LB>=0.5 转好）与结局的 2x2 关联
协议：20 种子 LASSO（内层 5 折 CV 选 C）内部 CV
输出：<OUT_DIR>/modeling/results_d7_dynamic.csv
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
OUT = '<OUT_DIR>/modeling/results_d7_dynamic.csv'
np.random.seed(RNG)


def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))


G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [norm_sid(i) for i in G_rel.index]
G_log = np.log1p(G_rel.clip(lower=0))

G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
tri = {}
for pat, s in G_map.groupby('patient'):
    v = {}
    for i, tag in [(1, '01'), (2, '02'), (3, '03')]:
        r = s[s['visit'].astype(int).astype(str).str.zfill(2) == tag]
        if len(r) == 1:
            v[i] = r.iloc[0]['run']
    if set(v) == {1, 2, 3}:
        tri[pat] = v
print(f'Gustin v1-v2-v3 完整配对: {len(tri)} 例')

D0 = G_log.loc[[v[1] for v in tri.values()]].reset_index(drop=True)
D7 = G_log.loc[[v[2] for v in tri.values()]].reset_index(drop=True)
D30 = G_log.loc[[v[3] for v in tri.values()]].reset_index(drop=True)
lb = 'Lactobacillus'
yv = (D30[lb] > np.log1p(0.5)).values.astype(int)
cst4_d0 = (D0[lb] < np.log1p(0.5)).values.astype(bool)
print(f'结局康复率: {yv.mean():.3f} ({yv.sum()}/{len(yv)}); 基线 CST4: {cst4_d0.sum()}/{len(cst4_d0)}')


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


rows = []


def rec(tag, mu, sd, n, note=''):
    rows.append(dict(analysis=tag, auc_mean=round(float(mu), 4), auc_sd=round(float(sd), 4),
                     n=int(n), note=note))
    print(f'  {tag}: AUC {mu:.3f} ± {sd:.3f} (n={n}) {note}')


print('\n===== 特征组比较（全人群） =====')
rec('D0 only (基线属级)', *internal_cv(D0, yv), note='静态基线')
rec('D7 only (D7属级)', *internal_cv(D7, yv), note='动态中间点')
rec('D0+D7 拼接', *internal_cv(pd.concat([D0.add_suffix('_d0'), D7.add_suffix('_d7')], axis=1), yv), note='静态+动态')
rec('Δ(D7-D0) 变化量', *internal_cv(D7 - D0, yv), note='log1p 域变化')

print('\n===== 基线 CST4 亚组 =====')
rec('D0 only CST4', *internal_cv(D0[cst4_d0], yv[cst4_d0]))
rec('D7 only CST4', *internal_cv(D7[cst4_d0], yv[cst4_d0]))
rec('D0+D7 CST4', *internal_cv(pd.concat([D0.add_suffix('_d0'), D7.add_suffix('_d7')], axis=1)[cst4_d0], yv[cst4_d0]))
rec('Δ(D7-D0) CST4', *internal_cv((D7 - D0)[cst4_d0], yv[cst4_d0]))

print('\n===== D7 CST 状态转换与结局关联 =====')
d7_ok = (D7[lb] > np.log1p(0.5)).values  # D7 时 LB>=50% = 转好
cst4_d0v = cst4_d0
mask = cst4_d0v
conv = d7_ok[mask]
outc = yv[mask]
tab = np.zeros((2, 2))
for c, o in zip(conv, outc):
    tab[int(c), int(o)] += 1
print('行=D7转好(1)/未转好(0), 列=D30康复(1)/未康复(0)')
print(pd.DataFrame(tab, index=['未转好', '转好'], columns=['未康复', '康复']).to_string())
tp, fp = tab[1, 1], tab[1, 0]
fn, tn = tab[0, 1], tab[0, 0]
sens = tp / (tp + fn) if (tp + fn) else np.nan
spec = tn / (tn + fp) if (tn + fp) else np.nan
ppv = tp / (tp + fp) if (tp + fp) else np.nan
print(f'  灵敏度(转好->康复): {sens:.3f} ({int(tp)}/{int(tp+fn)})')
print(f'  特异度(未转好->未康复): {spec:.3f} ({int(tn)}/{int(tn+fp)})')
print(f'  阳性预测值: {ppv:.3f}')
try:
    from scipy.stats import fisher_exact
    orr, pv = fisher_exact([[tp, fp], [fn, tn]])
    print(f'  Fisher OR: {orr:.3f} (p={pv:.3f})')
except Exception as e:
    print('  fisher_exact 失败:', e)

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print(f'\n已保存: {OUT}')
print(res.to_string(index=False))
