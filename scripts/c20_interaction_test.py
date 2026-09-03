# -*- coding: utf-8 -*-
"""c20_interaction_test.py —— #C20 任务⑤补充：richness x depth 交互检验"""
import pandas as pd
import numpy as np
import re
import statsmodels.api as sm
from scipy import stats

# 数据（复刻 c20_sensitivity.py 任务⑤加载）
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)

af = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af['run'] = [re.sub(r'_F\.fastq\.gz$', '', str(s)) for s in af['sample']]
af = af.set_index('run')
base_runs = [r for r in g_base if r in af.index]
rich = af.loc[base_runs]['species_observed'].values.astype(float)

counts = pd.read_csv('<OUT_DIR>/PRJNA691964/asv_table.csv', index_col=0)
counts.index = [re.sub(r'_F\.fastq\.gz$', '', str(i)) for i in counts.index]
depth = counts.loc[base_runs].values.astype(float).sum(axis=1)
log_depth = np.log(depth)

G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [re.sub(r'_F\.fastq\.gz$', '', str(i)) for i in G_rel.index]
G_out_lb = G_rel.loc[list(g_out)]['Lactobacillus'].values
y_g = (G_out_lb >= 0.5).astype(int)

gm = {r: i for i, r in enumerate(g_base)}
sel = [gm[r] for r in base_runs]
y_g = y_g[np.array(sel)]
print(f'n={len(y_g)}, rec={y_g.sum()}')

# 标准化避免交互项尺度问题
z_rich = (rich - rich.mean()) / rich.std()
z_depth = (log_depth - log_depth.mean()) / log_depth.std()
inter = z_rich * z_depth

# logistic 含交互项
X = np.column_stack([np.ones(len(y_g)), z_rich, z_depth, inter])
fit = sm.Logit(y_g, X).fit(disp=0)
print('\nlogistic outcome ~ richness + log_depth + richness:log_depth')
print(f'  richness    beta={fit.params[1]:.4f} p={fit.pvalues[1]:.4f}')
print(f'  log_depth   beta={fit.params[2]:.4f} p={fit.pvalues[2]:.4f}')
print(f'  交互项      beta={fit.params[3]:.4f} p={fit.pvalues[3]:.4f}')

# 置换检验：交互项 p 的稳健性（标签置换 5000 次）
rng = np.random.default_rng(42)
obs_b = fit.params[3]
n_perm = 5000
cnt = 0
for _ in range(n_perm):
    yp = rng.permutation(y_g)
    try:
        f2 = sm.Logit(yp, X).fit(disp=0, maxiter=100)
        if abs(f2.params[3]) >= abs(obs_b):
            cnt += 1
    except Exception:
        pass
print(f'\n交互项置换检验（5000 次）: p = {cnt / n_perm:.4f}')
