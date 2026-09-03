# -*- coding: utf-8 -*-
"""
c6_ipd_4cohorts.py —— #C6 ② 四队列 IPD 拼接补做（Mollin 补齐后重跑）
=========================================================================
派单 ② 原文："Mollin 基线 CST4 尚未统计，② 可先做三队列，Mollin 补齐后再并入重跑"

口径说明：
- 四队列属名交集 = 61 属（Xiao/Plummer/Gustin 三队列 65 属交集 ∩ Mollin 归一后属集）
- richness = 共同属 richness（每样本在 61 交集属中非零属数；Mollin 无 species_observed）
- 基线 CST4 亚组 = base_LB < 0.5
- 统一结局 = 随访 LB >= 50% (康复=1)
- ②a 朴素拼接单模型（logit） + ②b cohort 主效应 + cohort×richness 交互（LR test）

输出：追加到 results_c6_ipd.csv（标注 C6-2a4 / C6-2b4）
"""
import pandas as pd
import numpy as np
import re
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
from scipy.stats import chi2

OUT = '<OUT_DIR>/modeling/results_c6_ipd.csv'
np.random.seed(42)


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


# ================= 1. 加载四队列属级（共同 61 属） =================
def load_log1p(pid, is_mollin=False):
    if is_mollin:
        g = pd.read_csv(f'<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
    else:
        g = pd.read_csv(f'<DATA_DIR>/feat/{pid}_genus_rel.csv', index_col=0)
    g.index = [norm_sid(i) for i in g.index]
    g.columns = [norm_genus(c) for c in g.columns]
    agg = g.T.groupby(level=0).sum().T
    return np.log1p(agg)


# 先算四队列属名交集
G = {'Xiao': load_log1p('PRJNA302078'),
     'Plummer': load_log1p('PRJNA592384'),
     'Gustin': load_log1p('PRJNA691964'),
     'Mollin': load_log1p('', is_mollin=True)}
inter4 = sorted(set.intersection(*[set(g.columns) for g in G.values()]))
print(f'四队列属名交集: {len(inter4)} 属')

# 每样本算共同属 richness（log1p > 0 等价于原始 > 0）
for c in G:
    G[c] = G[c].reindex(columns=inter4).fillna(0)
    G[c]['richness_common'] = (G[c] > 0).sum(axis=1)

# ================= 2. 构建 subject 表 =================
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
G_map = G_map[G_map['control'].isna()]
M_pairs = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')

subj = []

# Xiao
xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
xiao = xiao[xiao['base_cst4'] == 1]
for _, r in xiao.iterrows():
    if r['base'] in G['Xiao'].index:
        subj.append(dict(cohort='Xiao', run=r['base'],
                         richness=G['Xiao'].loc[r['base'], 'richness_common'],
                         richness_sp=None,
                         outcome=int(r['outcome'])))

# Plummer
plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
plum = plum[plum['base_cst4'] == 1]
for _, r in plum.iterrows():
    if r['base'] in G['Plummer'].index:
        subj.append(dict(cohort='Plummer', run=r['base'],
                         richness=G['Plummer'].loc[r['base'], 'richness_common'],
                         richness_sp=None,
                         outcome=int(r['outcome'])))

# Gustin（全队列即 CST4；v1 基线）
pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
for b_run, o_run in pairs:
    # Gustin 结局：修正 LB>=0.5（基于原始相对丰度）
    lb_raw = np.expm1(G['Gustin'].loc[o_run, 'Lactobacillus'])
    subj.append(dict(cohort='Gustin', run=b_run,
                     richness=G['Gustin'].loc[b_run, 'richness_common'],
                     richness_sp=None,
                     outcome=int(lb_raw >= 0.5)))

# Mollin（mollin_pairs.csv 有 base_cst4 + outcome；LB>=0.5 口径）
mol = M_pairs[M_pairs['base_cst4'] == 1].copy()
for _, r in mol.iterrows():
    if r['base_sample'] in G['Mollin'].index:
        subj.append(dict(cohort='Mollin', run=r['base_sample'],
                         richness=G['Mollin'].loc[r['base_sample'], 'richness_common'],
                         richness_sp=None,
                         outcome=int(r['outcome'])))

S = pd.DataFrame(subj)
print('\n=== 四队列 CST4 亚组 subject 表（共同属 richness 口径）===')
print(S.groupby('cohort').agg(n=('outcome', 'size'), rec=('outcome', 'sum'),
                              rich_mean=('richness', 'mean')).round(2).to_string())
print(f'总 n = {len(S)}')

# ================= 3. ②a 朴素拼接 =================
X = S[['richness']].values.astype(float)
yv = S['outcome'].values.astype(int)

lr = LogisticRegression(solver='liblinear', class_weight='balanced', max_iter=5000)
lr.fit(X, yv)
p_all = lr.predict_proba(X)[:, 1]
auc_all = roc_auc_score(yv, p_all)

# bootstrap 95%CI
aucs = []
for _ in range(2000):
    idx = np.random.choice(len(yv), len(yv), replace=True)
    if len(np.unique(yv[idx])) < 2:
        continue
    try:
        lr_b = LogisticRegression(solver='liblinear', class_weight='balanced', max_iter=5000)
        lr_b.fit(X[idx], yv[idx])
        pb = lr_b.predict_proba(X[idx])[:, 1]
        aucs.append(roc_auc_score(yv[idx], pb))
    except Exception:
        pass
aucs = np.array(aucs)
ci_lo, ci_hi = np.percentile(aucs, [2.5, 97.5])

from scipy.stats import spearmanr
rho_p, p_p = spearmanr(S['richness'].values, yv)
print(f'\n②a4 朴素拼接（共同属 richness, 四队列 n={len(S)}）:')
print(f'  AUC {auc_all:.4f} [95%CI {ci_lo:.4f},{ci_hi:.4f}]')
print(f'  spearman rho={rho_p:+.3f} p={p_p:.3f}')

# ================= 4. ②b 交互项 LR test =================
S2 = S.copy()
S2['rich_c'] = StandardScaler().fit_transform(S2[['richness']].values.astype(float))
cohorts_present = sorted(S2['cohort'].unique())
ref = 'Gustin' if 'Gustin' in cohorts_present else cohorts_present[0]
others = [c for c in cohorts_present if c != ref]
dummies = pd.get_dummies(S2['cohort'], prefix='coh').astype(float)
keep = [c for c in dummies.columns if f'coh_{ref}' not in c]
D = dummies[keep]
for c in D.columns:
    S2[f'{c}_rich'] = D[c].values * S2['rich_c'].values
X_design = pd.concat([S2[['rich_c']], D, S2[[f'{c}_rich' for c in D.columns]]], axis=1)
X_design = sm.add_constant(X_design)
y_mod = S2['outcome'].values.astype(float)

interact_cols = [f'{c}_rich' for c in D.columns]
main_cols = ['const', 'rich_c'] + list(D.columns)
X1 = X_design[main_cols + interact_cols]
X0 = X_design[main_cols]
m1 = sm.Logit(y_mod, X1).fit(disp=0, maxiter=500)
m0 = sm.Logit(y_mod, X0).fit(disp=0, maxiter=500)
LR = -2 * (m0.llf - m1.llf)
df_int = len(interact_cols)
p_int = chi2.sf(LR, df_int)
print(f'\n②b4 交互项 LR 检验（四队列, 参考={ref}）:')
print(f'  LR={LR:.3f}, df={df_int}, p={p_int:.4g}')
for c in interact_cols:
    print(f'  {c}: coef={m1.params[c]:+.3f}, p={m1.pvalues[c]:.3f}')

# 各队列 slope
print('\n  队列内单变量 slope:')
for c in cohorts_present:
    sc = S2[S2['cohort'] == c]
    if len(np.unique(sc['outcome'])) < 2:
        continue
    lm = LogisticRegression(solver='liblinear', max_iter=5000).fit(
        sc[['richness']].values, sc['outcome'].values)
    print(f'    {c}: coef={lm.coef_[0][0]:+.4f} (n={len(sc)})')

# ================= 5. 追加写回 results_c6_ipd.csv =================
rows = [
    dict(analysis=f'C6-2a4 naive pooled 4-cohort (common-genus richness)', auc=round(auc_all, 4),
         ci_lo=round(ci_lo, 4), ci_hi=round(ci_hi, 4), n=int(len(yv)),
         note=f'spearman rho={rho_p:+.3f} p={p_p:.3f}'),
    dict(analysis=f'C6-2b4 cohort×richness interaction 4-cohort (LR test)',
         auc=np.nan, ci_lo=np.nan, ci_hi=np.nan, n=int(len(yv)),
         note=f'LR={LR:.3f} df={df_int} p={p_int:.4g} ref={ref}'),
]
# 读旧 CSV，追加新行
try:
    old = pd.read_csv(OUT)
    res = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
except FileNotFoundError:
    res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print(f'\n已追加到 {OUT}')
print(res.tail(4).to_string(index=False))
