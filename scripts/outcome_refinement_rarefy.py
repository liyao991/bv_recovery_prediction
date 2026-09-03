# outcome_refinement_rarefy.py
# rarefaction 校正版：连续结局（随访 LB 丰度 / 恢复幅度 ΔLB）下的 richness 方向性
# 目的：补 #C4-1 —— 确认"连续结局"口径下 Gustin richness 信号（rho=-0.565/-0.557）经 rarefaction 校正后是否仍稳健，
#       排除"测序深度混杂"对连续结局方向性结果的解释。
# 方法：与 richness_rarefy_check.py 同口径 rarefy（每样本多项抽样 N 次取均值），再与 outcome_refinement.py 的连续结局做 Spearman。
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr

RNG_SEED = 42
N_RArefy = 200
OUT = '<OUT_DIR>/modeling/results_outcome_refinement_rarefy.csv'

def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))

def load_counts(path):
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    df.index = [norm_sid(i) for i in df.index]
    return df

def rarefy_richness(counts, depth, n_repeat=N_RArefy, seed=RNG_SEED):
    total = counts.sum()
    if total < depth:
        return np.nan
    p = counts / total
    rng = np.random.default_rng(seed)
    samp = rng.multinomial(depth, p, size=n_repeat)
    return (samp > 0).sum(axis=1).mean()

def rarefy_series(counts_df, ids, depth):
    base = counts_df.loc[ids].values.astype(int)
    return np.array([rarefy_richness(c, depth) for c in base])

y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
rows = []

# ---- Xiao ----
g = pd.read_csv('<DATA_DIR>/feat/PRJNA302078_genus_rel.csv', index_col=0)
g.index = [norm_sid(i) for i in g.index]
af = pd.read_csv('<DATA_DIR>/feat/PRJNA302078_alpha_family.csv')
af['run'] = [norm_sid(s) for s in af['sample']]; af = af.set_index('run')
xiao = y[y['cohort']=='PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
xiao['out'] = [norm_sid(s) for s in xiao['outcome_sample']]
cx = load_counts('<OUT_DIR>/PRJNA302078/asv_table.csv')
dx = int(cx.loc[xiao['base']].values.sum(axis=1).min())
rx = rarefy_series(cx, xiao['base'], dx)
fu_lb = np.log1p(g.loc[xiao['out']]['Lactobacillus'].values)
base_lb = np.log1p(g.loc[xiao['base']]['Lactobacillus'].values)
r_s, p_s = spearmanr(rx, fu_lb)
r_d, p_d = spearmanr(rx, fu_lb - base_lb)
rows.append(dict(cohort='Xiao', n=len(xiao), rarefy_depth=dx,
    rho_rich_vs_fuLB=r_s, p_rich_vs_fuLB=p_s, rho_rich_vs_deltaLB=r_d, p_rich_vs_deltaLB=p_d))

# ---- Plummer ----
g2 = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_genus_rel.csv', index_col=0)
g2.index = [norm_sid(i) for i in g2.index]
af2 = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_alpha_family.csv')
af2['run'] = [norm_sid(s) for s in af2['sample']]; af2 = af2.set_index('run')
plum = y[y['cohort']=='PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
plum['out'] = [norm_sid(s) for s in plum['outcome_sample']]
cp = load_counts('<OUT_DIR>/PRJNA592384/asv_table.csv')
dp = int(cp.loc[plum['base']].values.sum(axis=1).min())
rp = rarefy_series(cp, plum['base'], dp)
fu_lb2 = np.log1p(g2.loc[plum['out']]['Lactobacillus'].values)
base_lb2 = np.log1p(g2.loc[plum['base']]['Lactobacillus'].values)
r_s2, p_s2 = spearmanr(rp, fu_lb2)
r_d2, p_d2 = spearmanr(rp, fu_lb2 - base_lb2)
rows.append(dict(cohort='Plummer', n=len(plum), rarefy_depth=dp,
    rho_rich_vs_fuLB=r_s2, p_rich_vs_fuLB=p_s2, rho_rich_vs_deltaLB=r_d2, p_rich_vs_deltaLB=p_d2))

# ---- Gustin（三时间点：richness vs D30 LB / D30-D0 ΔLB）----
g3 = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
g3.index = [norm_sid(i) for i in g3.index]
af3 = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af3['run'] = [norm_sid(s) for s in af3['sample']]; af3 = af3.set_index('run')
gm = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
gm = gm[gm['control'].isna()]
pat_v1 = gm[gm['visit'].astype(float)==1].set_index('patient')['run']
pat_v3 = gm[gm['visit'].astype(float)==3].set_index('patient')['run']
patients = sorted(set(pat_v1.index) & set(pat_v3.index))
v1r = pat_v1[patients].values; v3r = pat_v3[patients].values
cg = load_counts('<OUT_DIR>/PRJNA691964/asv_table.csv')
dg = int(cg.loc[v1r].values.sum(axis=1).min())
rg = rarefy_series(cg, v1r, dg)
lb1 = np.log1p(g3.loc[v1r]['Lactobacillus'].values)
lb3 = np.log1p(g3.loc[v3r]['Lactobacillus'].values)
r_s3, p_s3 = spearmanr(rg, lb3)
r_d3, p_d3 = spearmanr(rg, lb3 - lb1)
rows.append(dict(cohort='Gustin', n=len(patients), rarefy_depth=dg,
    rho_rich_vs_fuLB=r_s3, p_rich_vs_fuLB=p_s3, rho_rich_vs_deltaLB=r_d3, p_rich_vs_deltaLB=p_d3))

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print('=== 连续结局 rarefaction 校正版：richness vs 随访 LB / ΔLB ===')
print(f'rarefy 深度: Xiao={dx}, Plummer={dp}, Gustin={dg}; 每样本 {N_RArefy} 次多项抽样')
print(out.round(4).to_string(index=False))
print()
print('【对照】未 rarefy（results_outcome_refinement.csv）：')
print('  Xiao   : rho_vs_fuLB=-0.011(0.932)  rho_vs_ΔLB=+0.221(0.070)')
print('  Plummer: rho_vs_fuLB=+0.102(0.486)  rho_vs_ΔLB=+0.286(0.046)')
print('  Gustin : rho_vs_fuLB=-0.565(0.001)  rho_vs_ΔLB=-0.557(0.001)')
