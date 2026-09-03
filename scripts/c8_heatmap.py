# c8_heatmap.py — #C8 靶标 x richness rho 热图（四队列 CST4）
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

df = pd.read_csv('<REPO_DIR>/results/results_c8_richness_targets.csv')
cst4 = df[df['cohort'].str.contains('CST4')].copy()
order = ['Lactobacillus', 'Gardnerella', 'Atopobium', 'Megasphaera', 'Prevotella_agg',
         'Sneathia', 'Mobiluncus', 'Anaerobes_agg', 'Dysbiosis_index']
cohorts = ['Xiao_CST4', 'Plummer_CST4', 'Gustin_CST4(all)', 'Mollin_CST4']
cname = {'Xiao_CST4': 'Xiao\n(n=54)', 'Plummer_CST4': 'Plummer\n(n=43)',
         'Gustin_CST4(all)': 'Gustin\n(n=33)', 'Mollin_CST4': 'Mollin\n(n=30)'}

mat = np.full((len(order), len(cohorts)), np.nan)
sig = np.full((len(order), len(cohorts)), False, bool)
for i, t in enumerate(order):
    for j, c in enumerate(cohorts):
        row = cst4[(cst4['target'] == t) & (cst4['cohort'] == c)]
        if len(row):
            mat[i, j] = row['rho'].iloc[0]
            sig[i, j] = row['p'].iloc[0] < 0.05

fig, ax = plt.subplots(figsize=(8, 6))
vmax = max(np.nanmax(np.abs(mat)), 0.5)
im = ax.imshow(mat, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
ax.set_xticks(range(len(cohorts)))
ax.set_xticklabels([cname[c] for c in cohorts], fontsize=10)
ax.set_yticks(range(len(order)))
ylabels = {'Lactobacillus': 'Lactobacillus', 'Gardnerella': 'Gardnerella (GV)',
           'Atopobium': 'Atopobium (AV)', 'Megasphaera': 'Megasphaera',
           'Prevotella_agg': 'Prevotella (P.bivia代理)', 'Sneathia': 'Sneathia',
           'Mobiluncus': 'Mobiluncus', 'Anaerobes_agg': '厌氧菌聚合', 'Dysbiosis_index': '失调指数'}
ax.set_yticklabels([ylabels[t] for t in order], fontsize=10)
for i in range(len(order)):
    for j in range(len(cohorts)):
        v = mat[i, j]
        if np.isfinite(v):
            s = '*' if sig[i, j] else ''
            ax.text(j, i, f'{v:.2f}{s}', ha='center', va='center', fontsize=9,
                    color='white' if abs(v) > vmax * 0.55 else 'black')
ax.set_title('richness × 靶标 Spearman ρ（基线 CST4 亚组）\n红=正相关，蓝=负相关，*=p<0.05', fontsize=11)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Spearman ρ')
plt.tight_layout()
plt.savefig('<REPO_DIR>/results/fig_c8_targets_heatmap.png', dpi=150)
print('saved fig_c8_targets_heatmap.png')
