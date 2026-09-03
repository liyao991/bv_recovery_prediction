# richness_direction_check.py
# 基线种级 richness 与 D30 康复结局的方向性验证（三队列 ρ 是否同号）
# 目标：验证 "richness 高 -> 康复率低"（负相关）是否跨队列稳健，支撑 richness 风险分层作为正向可操作点
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr, mannwhitneyu

def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))

SEED_SAVE = '<OUT_DIR>/modeling/results_richness_direction.csv'
rows = []

# ---- 1. Xiao (PRJNA302078) ----
af = pd.read_csv('<DATA_DIR>/feat/PRJNA302078_alpha_family.csv')
af['run'] = [norm_sid(s) for s in af['sample']]
af = af.set_index('run')
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
a = af.loc[xiao['base']].reset_index(drop=True)
yv = xiao['outcome'].values.astype(int)
cst = xiao['base_cst4'].values.astype(bool)
rich = a['species_observed'].values

for tag, m in [('Xiao_full', ~cst & ~cst), ('Xiao_full_all', np.ones(len(yv), bool)), ('Xiao_CST4', cst)]:
    r, p = spearmanr(rich[m], yv[m])
    med = np.median(rich[m])
    hi = yv[m][rich[m] >= med]; lo = yv[m][rich[m] < med]
    u = mannwhitneyu(rich[m][yv[m] == 1], rich[m][yv[m] == 0]) if (yv[m].sum() > 0 and (1 - yv[m]).sum() > 0) else None
    rows.append(dict(cohort='Xiao', subgroup=tag, n=int(m.sum()),
                     n_rec=float(yv[m].sum()), rho=r, p=p,
                     rec_hi=hi.mean(), rec_lo=lo.mean(),
                     mwu_p=(u.pvalue if u else np.nan)))

# ---- 2. Plummer (PRJNA592384) ----
af2 = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_alpha_family.csv')
af2['run'] = [norm_sid(s) for s in af2['sample']]
af2 = af2.set_index('run')
plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
a2 = af2.loc[plum['base']].reset_index(drop=True)
yv2 = plum['outcome'].values.astype(int)
cst2 = plum['base_cst4'].values.astype(bool)
rich2 = a2['species_observed'].values

for tag, m in [('Plummer_full_all', np.ones(len(yv2), bool)), ('Plummer_CST4', cst2)]:
    if m.sum() < 3:
        continue
    r, p = spearmanr(rich2[m], yv2[m])
    med = np.median(rich2[m])
    hi = yv2[m][rich2[m] >= med]; lo = yv2[m][rich2[m] < med]
    u = mannwhitneyu(rich2[m][yv2[m] == 1], rich2[m][yv2[m] == 0]) if (yv2[m].sum() > 0 and (1 - yv2[m]).sum() > 0) else None
    rows.append(dict(cohort='Plummer', subgroup=tag, n=int(m.sum()),
                     n_rec=float(yv2[m].sum()), rho=r, p=p,
                     rec_hi=hi.mean(), rec_lo=lo.mean(),
                     mwu_p=(u.pvalue if u else np.nan)))

# ---- 3. Gustin (PRJNA691964) ----
af3 = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af3['run'] = [norm_sid(s) for s in af3['sample']]
af3 = af3.set_index('run')
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
yv3 = (G_rel.loc[list(g_out)]['Lactobacillus'] >= 0.5).values.astype(int)  # 正确口径 LB>=0.5（原 >log1p(0.5)=0.4055 为笔误）
rich3 = af3.loc[list(g_base)]['species_observed'].values

r, p = spearmanr(rich3, yv3)
med = np.median(rich3)
hi = yv3[rich3 >= med]; lo = yv3[rich3 < med]
u = mannwhitneyu(rich3[yv3 == 1], rich3[yv3 == 0])
rows.append(dict(cohort='Gustin', subgroup='Gustin_all(CST4)', n=len(yv3),
                 n_rec=float(yv3.sum()), rho=r, p=p,
                 rec_hi=hi.mean(), rec_lo=lo.mean(), mwu_p=u.pvalue))

# ---- 输出 ----
out = pd.DataFrame(rows)
out.to_csv(SEED_SAVE, index=False)
print('=== 基线种级 richness vs D30 康复（方向性验证）===')
print(out[['cohort', 'subgroup', 'n', 'n_rec', 'rho', 'p', 'rec_hi', 'rec_lo', 'mwu_p']].round(3).to_string(index=False))
print()
print('【判读】rho<0 = richness 高->康复率低（与 Gustin 一致）；rho>0 = 方向相反')
print('三队列 rho 是否同号：', '是' if len(set(np.sign([x['rho'] for x in rows if x['rho'] != 0]))) == 1 else '否/有例外')
