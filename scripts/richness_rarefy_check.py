# richness_rarefy_check.py
# rarefaction 校正版：基线种级 richness vs D30 康复结局 的方向性验证
# 目标：排除"测序深度混杂"对 Gustin richness 信号（rho=-0.498 p=0.003）的解释
# 方法：从 ASV 计数表无放回抽样(rarefy)到队列内基线样本最小 reads 数，重算 observed species richness，
#       重复 R 次取均值，再用 rarefied richness 复算 Spearman rho/p，与 results_richness_direction.csv 同格式输出
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr, mannwhitneyu

RNG_SEED = 42
N_RArefy = 200          # 每样本 rarefy 重复次数
OUT = '<OUT_DIR>/modeling/results_richness_rarefy.csv'

def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))

def load_counts(path):
    # ASV 计数表：第 1 行表头(ASV序列)，第 1 列样本名；返回 DataFrame(样本 × ASV) int
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    df.index = [norm_sid(i) for i in df.index]
    return df

def rarefy_richness(counts, depth, n_repeat=N_RArefy, seed=RNG_SEED):
    # counts: 1D int array (某样本 ASV counts)；depth: 抽样深度
    # 多项分布近似无放回抽样，返回该样本 rarefied observed richness（n_repeat 次均值）
    total = counts.sum()
    if total < depth:
        return np.nan
    p = counts / total
    rng = np.random.default_rng(seed)
    samp = rng.multinomial(depth, p, size=n_repeat)  # (n_repeat, n_asv)
    rich = (samp > 0).sum(axis=1)
    return rich.mean()

rows = []

# ================= Xiao (PRJNA302078) =================
counts_x = load_counts('<OUT_DIR>/PRJNA302078/asv_table.csv')
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
yv = xiao['outcome'].values.astype(int)
cst = xiao['base_cst4'].values.astype(bool)

base_counts = counts_x.loc[xiao['base']].values.astype(int)
totals = base_counts.sum(axis=1)
depth = int(totals.min())
rich = np.array([rarefy_richness(c, depth) for c in base_counts])

for tag, m in [('Xiao_full_all', np.ones(len(yv), bool)), ('Xiao_CST4', cst)]:
    r, p = spearmanr(rich[m], yv[m])
    med = np.median(rich[m])
    hi = yv[m][rich[m] >= med]; lo = yv[m][rich[m] < med]
    u = mannwhitneyu(rich[m][yv[m] == 1], rich[m][yv[m] == 0]) if (yv[m].sum() > 0 and (1 - yv[m]).sum() > 0) else None
    rows.append(dict(cohort='Xiao', subgroup=tag, n=int(m.sum()), n_rec=float(yv[m].sum()),
                     rho=r, p=p, rec_hi=hi.mean(), rec_lo=lo.mean(),
                     mwu_p=(u.pvalue if u else np.nan), rarefy_depth=depth))

# ================= Plummer (PRJNA592384) =================
counts_p = load_counts('<OUT_DIR>/PRJNA592384/asv_table.csv')
plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
yv2 = plum['outcome'].values.astype(int)
cst2 = plum['base_cst4'].values.astype(bool)

base_counts2 = counts_p.loc[plum['base']].values.astype(int)
depth2 = int(base_counts2.sum(axis=1).min())
rich2 = np.array([rarefy_richness(c, depth2) for c in base_counts2])

for tag, m in [('Plummer_full_all', np.ones(len(yv2), bool)), ('Plummer_CST4', cst2)]:
    if m.sum() < 3:
        continue
    r, p = spearmanr(rich2[m], yv2[m])
    med = np.median(rich2[m])
    hi = yv2[m][rich2[m] >= med]; lo = yv2[m][rich2[m] < med]
    u = mannwhitneyu(rich2[m][yv2[m] == 1], rich2[m][yv2[m] == 0]) if (yv2[m].sum() > 0 and (1 - yv2[m]).sum() > 0) else None
    rows.append(dict(cohort='Plummer', subgroup=tag, n=int(m.sum()), n_rec=float(yv2[m].sum()),
                     rho=r, p=p, rec_hi=hi.mean(), rec_lo=lo.mean(),
                     mwu_p=(u.pvalue if u else np.nan), rarefy_depth=depth2))

# ================= Gustin (PRJNA691964) =================
counts_g = load_counts('<OUT_DIR>/PRJNA691964/asv_table.csv')
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

g_base_counts = counts_g.loc[list(g_base)].values.astype(int)
depth3 = int(g_base_counts.sum(axis=1).min())
rich3 = np.array([rarefy_richness(c, depth3) for c in g_base_counts])

r, p = spearmanr(rich3, yv3)
med = np.median(rich3)
hi = yv3[rich3 >= med]; lo = yv3[rich3 < med]
u = mannwhitneyu(rich3[yv3 == 1], rich3[yv3 == 0])
rows.append(dict(cohort='Gustin', subgroup='Gustin_all(CST4)', n=len(yv3), n_rec=float(yv3.sum()),
                 rho=r, p=p, rec_hi=hi.mean(), rec_lo=lo.mean(), mwu_p=u.pvalue, rarefy_depth=depth3))

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print('=== rarefaction 校正后：基线种级 richness vs D30 康复（方向性验证）===')
print(f'rarefy 深度: Xiao={depth}, Plummer={depth2}, Gustin={depth3}; 每样本 {N_RArefy} 次多项抽样')
print(out[['cohort', 'subgroup', 'n', 'n_rec', 'rho', 'p', 'rec_hi', 'rec_lo', 'mwu_p', 'rarefy_depth']].round(4).to_string(index=False))
print()
print('【判读】rho<0 = richness 高->康复率低；rho>0 = 方向相反')
signs = set(np.sign([x['rho'] for x in rows if x['rho'] != 0 and not np.isnan(x['rho'])]))
print('三队列 rho 是否同号：', '是' if len(signs) == 1 else '否/有例外')
