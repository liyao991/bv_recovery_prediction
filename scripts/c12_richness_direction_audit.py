# -*- coding: utf-8 -*-
"""
c12_richness_direction_audit.py —— #C12 核实：Xiao richness 方向（四队列版 vs 主分析）不一致的根因拆解
================================================================================
背景：ASUS #C12 请求核实第 8 项遗留——框架 319 行三象限说「Xiao/Plummer 弱正」，
但主分析 Spearman ρ（Xiao 全 -0.054 弱负）方向相反。

三种口径对照，一次性定位差异来源：
  A. 主分析（richness_direction_check.py）：richness = 各队列种级 species_observed，Spearman ρ
  B. 四队列版（c6_ipd_4cohorts.py）：richness = 61 属交集共同属非零数，logistic slope（coef）
  C. 四队列 richness + Spearman ρ（B 的 richness 换成 ρ 统计量，隔离"统计量"因素）
  D. 四队列 richness + 各队列种级 richness 对照（隔离"richness 定义"因素）

只打印，不写回任何 CSV。
"""
import pandas as pd
import numpy as np
import re
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression


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


# ============ 1. 四队列共同属 richness（复现 c6_ipd_4cohorts.py） ============
def load_log1p(pid, is_mollin=False):
    if is_mollin:
        g = pd.read_csv('<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
    else:
        g = pd.read_csv(f'<DATA_DIR>/feat/{pid}_genus_rel.csv', index_col=0)
    g.index = [norm_sid(i) for i in g.index]
    g.columns = [norm_genus(c) for c in g.columns]
    agg = g.T.groupby(level=0).sum().T
    return np.log1p(agg)


G = {'Xiao': load_log1p('PRJNA302078'),
     'Plummer': load_log1p('PRJNA592384'),
     'Gustin': load_log1p('PRJNA691964'),
     'Mollin': load_log1p('', is_mollin=True)}
inter4 = sorted(set.intersection(*[set(g.columns) for g in G.values()]))
print(f'四队列属名交集: {len(inter4)} 属')

for c in G:
    G[c] = G[c].reindex(columns=inter4).fillna(0)
    G[c]['richness_common'] = (G[c] > 0).sum(axis=1)

# ============ 2. 构建 subject 表（与 c6_ipd_4cohorts.py 一致） ============
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
G_map = G_map[G_map['control'].isna()]
M_pairs = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')

subj = []

xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
xiao = xiao[xiao['base_cst4'] == 1]
for _, r in xiao.iterrows():
    if r['base'] in G['Xiao'].index:
        subj.append(dict(cohort='Xiao', run=r['base'],
                         richness=G['Xiao'].loc[r['base'], 'richness_common'],
                         outcome=int(r['outcome'])))

plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
plum = plum[plum['base_cst4'] == 1]
for _, r in plum.iterrows():
    if r['base'] in G['Plummer'].index:
        subj.append(dict(cohort='Plummer', run=r['base'],
                         richness=G['Plummer'].loc[r['base'], 'richness_common'],
                         outcome=int(r['outcome'])))

pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
for b_run, o_run in pairs:
    lb_raw = np.expm1(G['Gustin'].loc[o_run, 'Lactobacillus'])
    subj.append(dict(cohort='Gustin', run=b_run,
                     richness=G['Gustin'].loc[b_run, 'richness_common'],
                     outcome=int(lb_raw >= 0.5)))

mol = M_pairs[M_pairs['base_cst4'] == 1].copy()
for _, r in mol.iterrows():
    if r['base_sample'] in G['Mollin'].index:
        subj.append(dict(cohort='Mollin', run=r['base_sample'],
                         richness=G['Mollin'].loc[r['base_sample'], 'richness_common'],
                         outcome=int(r['outcome'])))

S = pd.DataFrame(subj)
print('\n=== 四队列 CST4 亚组 subject 表 ===')
print(S.groupby('cohort').agg(n=('outcome', 'size'), rec=('outcome', 'sum'),
                              rich_mean=('richness', 'mean')).round(2).to_string())

# ============ 3. 口径 B：四队列 richness + logistic slope ============
print('\n=== 口径 B：共同属 richness + logistic slope（复现三象限）===')
for c in ['Xiao', 'Plummer', 'Gustin', 'Mollin']:
    sc = S[S['cohort'] == c]
    if len(np.unique(sc['outcome'])) < 2:
        print(f'  {c}: 结局无变异性，跳过')
        continue
    lm = LogisticRegression(solver='liblinear', max_iter=5000).fit(
        sc[['richness']].values, sc['outcome'].values)
    print(f'  {c}: logistic coef={lm.coef_[0][0]:+.4f} (n={len(sc)}, rec={int(sc["outcome"].sum())})')

# ============ 4. 口径 C：四队列 richness + Spearman ρ ============
print('\n=== 口径 C：共同属 richness + Spearman ρ（隔离统计量）===')
for c in ['Xiao', 'Plummer', 'Gustin', 'Mollin']:
    sc = S[S['cohort'] == c]
    if len(np.unique(sc['outcome'])) < 2:
        continue
    rho, p = spearmanr(sc['richness'].values, sc['outcome'].values)
    print(f'  {c}: rho={rho:+.4f} p={p:.3f} (n={len(sc)})')

# ============ 5. 口径 A：各队列种级 species_observed + Spearman ρ ============
print('\n=== 口径 A：种级 species_observed + Spearman ρ（主分析口径）===')
alpha = {}
for pid, cname in [('PRJNA302078', 'Xiao'), ('PRJNA592384', 'Plummer'), ('PRJNA691964', 'Gustin')]:
    af = pd.read_csv(f'<DATA_DIR>/feat/{pid}_alpha_family.csv')
    af['run'] = [norm_sid(s) for s in af['sample']]
    af = af.set_index('run')
    alpha[cname] = af

# Xiao 种级 ρ (CST4)
xiao2 = y[y['cohort'] == 'PRJNA302078'].copy()
xiao2['base'] = [norm_sid(s) for s in xiao2['base_sample']]
cst_x = xiao2['base_cst4'].values.astype(bool)
rich_x = alpha['Xiao'].loc[xiao2['base']]['species_observed'].values
yv_x = xiao2['outcome'].values.astype(int)
r_x, p_x = spearmanr(rich_x[cst_x], yv_x[cst_x])
print(f'  Xiao CST4: rho={r_x:+.4f} p={p_x:.3f} (n={int(cst_x.sum())})')
r_xf, p_xf = spearmanr(rich_x, yv_x)
print(f'  Xiao full: rho={r_xf:+.4f} p={p_xf:.3f} (n={len(yv_x)})')

# Plummer 种级 ρ (CST4)
plum2 = y[y['cohort'] == 'PRJNA592384'].copy()
plum2['base'] = [norm_sid(s) for s in plum2['base_sample']]
cst_p = plum2['base_cst4'].values.astype(bool)
rich_p = alpha['Plummer'].loc[plum2['base']]['species_observed'].values
yv_p = plum2['outcome'].values.astype(int)
r_p, p_p = spearmanr(rich_p[cst_p], yv_p[cst_p])
print(f'  Plummer CST4: rho={r_p:+.4f} p={p_p:.3f} (n={int(cst_p.sum())})')
r_pf, p_pf = spearmanr(rich_p, yv_p)
print(f'  Plummer full: rho={r_pf:+.4f} p={p_pf:.3f} (n={len(yv_p)})')

# Gustin 种级 ρ
af3 = alpha['Gustin']
g_pairs2 = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs2.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs2)
G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [norm_sid(i) for i in G_rel.index]
yv_g = (G_rel.loc[list(g_out)]['Lactobacillus'] >= 0.5).values.astype(int)
rich_g = af3.loc[list(g_base)]['species_observed'].values
r_g, p_g = spearmanr(rich_g, yv_g)
print(f'  Gustin all(CST4): rho={r_g:+.4f} p={p_g:.3f} (n={len(yv_g)})')

print('\n=== 结论定位 ===')
print('若 B(logistic coef) 与 C(Spearman, 同 richness) 方向一致 → 差异来自 richness 定义（属级共同 vs 种级）')
print('若 C 与 A 方向一致、B 与 C 方向不同 → 差异来自统计量（logistic vs Spearman）')
