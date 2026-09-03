# outcome_refinement.py
# 结局精细化：把"二元结局(outcome>=0.5)"换成"连续结局(随访 Lactobacillus 丰度)"
# 目的：检验 richness / 特定菌 与"恢复程度"是否有单调关系（二元分类可能抹掉信号）
# 三队列：Xiao(454) / Plummer(MiSeq) / Gustin(MiSeq, 三时间点)
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr, pearsonr

def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))

y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
rows = []

# ---- 1. Xiao ----
g = pd.read_csv('<DATA_DIR>/feat/PRJNA302078_genus_rel.csv', index_col=0)
g.index = [norm_sid(i) for i in g.index]
af = pd.read_csv('<DATA_DIR>/feat/PRJNA302078_alpha_family.csv')
af['run'] = [norm_sid(s) for s in af['sample']]; af = af.set_index('run')
xiao = y[y['cohort']=='PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
xiao['out'] = [norm_sid(s) for s in xiao['outcome_sample']]
# 连续结局 = 随访 Lactobacillus 丰度（log1p）
fu_lb = np.log1p(g.loc[xiao['out']]['Lactobacillus'].values)
base_lb = np.log1p(g.loc[xiao['base']]['Lactobacillus'].values)
rich = af.loc[xiao['base']]['species_observed'].values
# richness vs 连续结局
r_s, p_s = spearmanr(rich, fu_lb)
r_p, p_p = pearsonr(rich, fu_lb)
# ΔLB = 随访 - 基线
delta = fu_lb - base_lb
r_d, p_d = spearmanr(rich, delta)
# 基线 LB 与随访 LB（自回归基线，看基线菌群能否预测随访）
r_bb, p_bb = spearmanr(base_lb, fu_lb)
rows.append(dict(cohort='Xiao', n=len(xiao),
    rho_rich_vs_fuLB=r_s, p_rich_vs_fuLB=p_s,
    rho_rich_vs_deltaLB=r_d, p_rich_vs_deltaLB=p_d,
    rho_baseLB_vs_fuLB=r_bb, p_baseLB_vs_fuLB=p_bb))

# ---- 2. Plummer ----
g2 = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_genus_rel.csv', index_col=0)
g2.index = [norm_sid(i) for i in g2.index]
af2 = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_alpha_family.csv')
af2['run'] = [norm_sid(s) for s in af2['sample']]; af2 = af2.set_index('run')
plum = y[y['cohort']=='PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
plum['out'] = [norm_sid(s) for s in plum['outcome_sample']]
fu_lb2 = np.log1p(g2.loc[plum['out']]['Lactobacillus'].values)
base_lb2 = np.log1p(g2.loc[plum['base']]['Lactobacillus'].values)
rich2 = af2.loc[plum['base']]['species_observed'].values
r_s2, p_s2 = spearmanr(rich2, fu_lb2)
r_d2, p_d2 = spearmanr(rich2, fu_lb2 - base_lb2)
r_bb2, p_bb2 = spearmanr(base_lb2, fu_lb2)
rows.append(dict(cohort='Plummer', n=len(plum),
    rho_rich_vs_fuLB=r_s2, p_rich_vs_fuLB=p_s2,
    rho_rich_vs_deltaLB=r_d2, p_rich_vs_deltaLB=p_d2,
    rho_baseLB_vs_fuLB=r_bb2, p_baseLB_vs_fuLB=p_bb2))

# ---- 3. Gustin（三时间点，重点）----
g3 = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
g3.index = [norm_sid(i) for i in g3.index]
af3 = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af3['run'] = [norm_sid(s) for s in af3['sample']]; af3 = af3.set_index('run')
gm = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
gm = gm[gm['control'].isna()]
def visit_run(v):
    return gm[gm['visit'].astype(float)==v]['run'].values
v1 = visit_run(1); v2 = visit_run(2); v3 = visit_run(3)
# 按 patient 对齐三时间点
pat_v1 = gm[gm['visit'].astype(float)==1].set_index('patient')['run']
pat_v2 = gm[gm['visit'].astype(float)==2].set_index('patient')['run']
pat_v3 = gm[gm['visit'].astype(float)==3].set_index('patient')['run']
patients = sorted(set(pat_v1.index) & set(pat_v2.index) & set(pat_v3.index))
lb1 = np.log1p(g3.loc[pat_v1[patients].values]['Lactobacillus'].values)
lb2 = np.log1p(g3.loc[pat_v2[patients].values]['Lactobacillus'].values)
lb3 = np.log1p(g3.loc[pat_v3[patients].values]['Lactobacillus'].values)
rich3 = af3.loc[pat_v1[patients].values]['species_observed'].values

# richness vs D30 连续 LB
r_s3, p_s3 = spearmanr(rich3, lb3)
# richness vs D7 LB（短期信号）
r_s3_d7, p_s3_d7 = spearmanr(rich3, lb2)
# richness vs 恢复轨迹（D30 - 基线）
r_d3, p_d3 = spearmanr(rich3, lb3 - lb1)
# D7 到 D30 的回落（复发程度）= lb3 - lb2（负 = 回落）
r_drop, p_drop = spearmanr(rich3, lb3 - lb2)
rows.append(dict(cohort='Gustin', n=len(patients),
    rho_rich_vs_fuLB=r_s3, p_rich_vs_fuLB=p_s3,
    rho_rich_vs_deltaLB=r_d3, p_rich_vs_deltaLB=p_d3,
    rho_baseLB_vs_fuLB=np.nan, p_baseLB_vs_fuLB=np.nan))

out = pd.DataFrame(rows)
out.to_csv('<OUT_DIR>/modeling/results_outcome_refinement.csv', index=False)
print('=== 结局精细化：richness vs 连续随访 LB 丰度 ===')
print(out.round(3).to_string(index=False))
print()
print('【Gustin 三时间点轨迹】')
print(f'  D0 LB median={np.median(lb1):.3f}, D7 LB median={np.median(lb2):.3f}, D30 LB median={np.median(lb3):.3f}')
print(f'  richness vs D7 LB: rho={r_s3_d7:.3f}, p={p_s3_d7:.3f}  (短期信号)')
print(f'  richness vs D30-D7回落: rho={r_drop:.3f}, p={p_drop:.3f}  (复发程度)')
