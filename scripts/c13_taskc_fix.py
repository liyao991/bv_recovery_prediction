# -*- coding: utf-8 -*-
"""
c13_taskc_fix.py —— Task C 修正：用主分析同口径 richness（alpha_family 的 species_observed）
复现 -0.454 后再做控制基线 LB 的偏相关。
"""
import pandas as pd
import numpy as np
import re
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

OUT_CSV = '<REPO_DIR>/results/results_c13_taskc_partial.csv'

# Gustin 配对与结局（与 four_cohort_matrix 一致）
G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [re.sub(r'_F\.fastq\.gz$', '', str(i)) for i in G_rel.index]
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)

# 结局：v3 LB>=0.5（33 例无恰好=0.5 者，> 与 >= 等价）
lb_out = G_rel.loc[list(g_out)]['Lactobacillus'].values.astype(float)
baseLB = G_rel.loc[list(g_base)]['Lactobacillus'].values.astype(float)
y_bin = (lb_out >= 0.5).astype(float)

# 主分析同口径 richness：alpha_family 的 species_observed
af = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
print('alpha_family 列:', list(af.columns))
id_col = af.columns[0]
af[id_col] = [re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s)) for s in af[id_col]]
af = af.set_index(id_col)
richness = af.loc[list(g_base)]['species_observed'].values.astype(float)
print(f'n={len(y_bin)} 康复={int(y_bin.sum())} richness 范围 [{richness.min():.0f},{richness.max():.0f}]')

rows = []
def rec(item, value, extra=''):
    rows.append(dict(item=item, value=value, extra=extra))
    print(f'  {item}: {value}  {extra}')

rho_ro, p_ro = stats.spearmanr(richness, y_bin)
rho_lo, p_lo = stats.spearmanr(baseLB, y_bin)
rho_rl, p_rl = stats.spearmanr(richness, baseLB)
rec('rho(richness, outcome)', f'{rho_ro:.3f} (p={p_ro:.4f})', '应复现主分析 -0.454')
rec('rho(baseLB, outcome)', f'{rho_lo:.3f} (p={p_lo:.4f})', '基线 LB 单特征（Gustin CST4 内）')
rec('rho(richness, baseLB)', f'{rho_rl:.3f} (p={p_rl:.4f})', '两者耦合强度')

def rank_resid(v, covar):
    rv = stats.rankdata(v); rc = stats.rankdata(covar)
    X_d = np.column_stack([np.ones(len(rc)), rc])
    beta = np.linalg.lstsq(X_d, rv, rcond=None)[0]
    return rv - X_d @ beta

res_r = rank_resid(richness, baseLB)
res_y = rank_resid(y_bin, baseLB)
r_part = np.corrcoef(res_r, res_y)[0, 1]
n = len(y_bin); k = 1
t_stat = r_part * np.sqrt((n - 2 - k) / (1 - r_part**2))
p_part = 2 * stats.t.sf(abs(t_stat), df=n - 2 - k)
rec('partial rho(richness, outcome | baseLB)', f'{r_part:.3f} (p={p_part:.4f})',
    '控制基线 LB 后')

Xs = StandardScaler().fit_transform(np.column_stack([richness, baseLB]))
lr = LogisticRegression(max_iter=5000).fit(Xs, y_bin)
rec('logistic 双变量 coef', f'richness={lr.coef_[0][0]:.3f}', f'baseLB={lr.coef_[0][1]:.3f}')

# 附加：Mann-Whitney（与主分析一致性交叉验证）
u = stats.mannwhitneyu(richness[y_bin == 1], richness[y_bin == 0])
rec('MWU richness 康复 vs 未康复', f'U={u.statistic:.0f} p={u.pvalue:.4f}', '')

pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f'\n已保存: {OUT_CSV}')
