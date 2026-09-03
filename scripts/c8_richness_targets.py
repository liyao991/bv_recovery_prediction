# c8_richness_targets.py
# #C8 richness 特质解剖：richness x 核心 8 靶标（四队列，基线 CST4 亚组，相对丰度）
# 目的：验证 "richness 是否 = 乳杆菌塌陷 + BV 厌氧菌扩张" 这一失调状态标量
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr

def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))

# ---------- 靶标定义（属级，可用子集） ----------
SINGLE_GENERA = ['Lactobacillus', 'Gardnerella', 'Atopobium', 'Megasphaera', 'Sneathia', 'Mobiluncus']

def prevotella_cols(cols):
    return [c for c in cols if c.startswith('Prevotella') or c.startswith('Alloprevotella') or c.startswith('Prevotellaceae')]

def build_targets(g):
    """g: 属级相对丰度 DataFrame(index=sample, columns=genus)。返回靶标矩阵 + 厌氧菌聚合 + 失调指数。"""
    cols = list(g.columns)
    m = pd.DataFrame(index=g.index)
    for t in SINGLE_GENERA:
        m[t] = g[t].astype(float) if t in cols else 0.0
    # Lacto 总（含 2020 年 Lactobacillus 拆分出的新属，阴道内通常极低，作敏感性）
    lacto_like = [c for c in cols if c in ('Lactobacillus', 'Limosilactobacillus', 'Ligilactobacillus', 'Liquorilactobacillus', 'Lacticaseibacillus')]
    m['Lacto_total'] = g[lacto_like].sum(axis=1) if lacto_like else m['Lactobacillus']
    pcols = prevotella_cols(cols)
    m['Prevotella_agg'] = g[pcols].sum(axis=1) if pcols else 0.0
    # BV 厌氧菌聚合（核心）
    anaerobes = ['Gardnerella', 'Atopobium', 'Megasphaera', 'Prevotella_agg', 'Sneathia', 'Mobiluncus']
    m['Anaerobes_agg'] = m[anaerobes].sum(axis=1)
    # 失调指数 = 乳杆菌占优度
    denom = m['Lactobacillus'] + m['Anaerobes_agg']
    m['Dysbiosis_index'] = np.where(denom > 0, m['Lactobacillus'] / denom, np.nan)
    return m

def ols_r2(X_raw, rich, log=True):
    """多元回归 richness ~ X，返回 (R2, R2_adj)。log=True 时对丰度 log1p、对 richness log1p。"""
    X = np.column_stack([np.ones(len(rich))] + [X_raw[c].values.astype(float) for c in X_raw.columns])
    y = np.log1p(rich.astype(float)) if log else rich.astype(float)
    Xl = X.copy()
    if log:
        Xl[:, 1:] = np.log1p(Xl[:, 1:])
    coef, _, _, _ = np.linalg.lstsq(Xl, y, rcond=None)
    yhat = Xl @ coef
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    n, k = Xl.shape
    r2_adj = 1 - (1 - r2) * (n - 1) / (n - k - 1) if n > k + 1 else np.nan
    return r2, r2_adj

# ---------- 加载与汇总 ----------
rows = []      # 长表：每队列 x 每靶标一行（Spearman richness x target）
var_rows = []  # 方差分解：每队列一行
consistency = {}  # 靶标 -> 各队列 rho 方向

def analyze(cohort, rich, targets, note=''):
    """rich: Series(sample->richness)；targets: build_targets 输出。"""
    regressors = ['Lactobacillus', 'Gardnerella', 'Atopobium', 'Megasphaera', 'Prevotella_agg', 'Sneathia', 'Mobiluncus']
    # 1. richness x 单靶标 Spearman
    for t in ['Lactobacillus', 'Lacto_total', 'Gardnerella', 'Atopobium', 'Megasphaera', 'Prevotella_agg', 'Sneathia', 'Mobiluncus', 'Anaerobes_agg', 'Dysbiosis_index']:
        vals = targets[t].values.astype(float)
        ok = np.isfinite(vals) & np.isfinite(rich.values)
        if ok.sum() < 3:
            continue
        r, p = spearmanr(rich.values[ok], vals[ok])
        rows.append(dict(cohort=cohort, target=t, n=int(ok.sum()), rho=r, p=p, note=note))
        consistency.setdefault(t, []).append(r)
    # 2. 方差分解：richness ~ 7 靶标
    X = targets[regressors]
    r2, r2_adj = ols_r2(X, rich, log=True)
    r2_raw, r2_raw_adj = ols_r2(X, rich, log=False)
    # 简化版：richness ~ [Lactobacillus, Anaerobes_agg]
    X2 = targets[['Lactobacillus', 'Anaerobes_agg']]
    r2b, r2b_adj = ols_r2(X2, rich, log=True)
    var_rows.append(dict(cohort=cohort, n=int(len(rich)), note=note,
                         R2_log1p=r2, R2_log1p_adj=r2_adj,
                         R2_raw=r2_raw, R2_raw_adj=r2_raw_adj,
                         R2_2var_log1p=r2b, R2_2var_adj=r2b_adj))

# ---- 1. Xiao (PRJNA302078) ----
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
af1 = pd.read_csv('<DATA_DIR>/feat/PRJNA302078_alpha_family.csv')
af1['run'] = [norm_sid(s) for s in af1['sample']]
af1 = af1.set_index('run')
g1 = pd.read_csv('<DATA_DIR>/feat/PRJNA302078_genus_rel.csv', index_col=0)
g1.index = [norm_sid(i) for i in g1.index]
xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
cst = xiao['base_cst4'].values.astype(bool)
for tag, mask in [('CST4', cst), ('full', np.ones(len(xiao), bool))]:
    idx = xiao['base'].values[mask]
    rich = af1.loc[idx]['species_observed'].astype(float)
    targets = build_targets(g1.loc[idx])
    analyze(f'Xiao_{tag}', rich, targets)

# ---- 2. Plummer (PRJNA592384) ----
af2 = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_alpha_family.csv')
af2['run'] = [norm_sid(s) for s in af2['sample']]
af2 = af2.set_index('run')
g2 = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_genus_rel.csv', index_col=0)
g2.index = [norm_sid(i) for i in g2.index]
plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
cst2 = plum['base_cst4'].values.astype(bool)
for tag, mask in [('CST4', cst2), ('full', np.ones(len(plum), bool))]:
    idx = plum['base'].values[mask]
    rich = af2.loc[idx]['species_observed'].astype(float)
    targets = build_targets(g2.loc[idx])
    analyze(f'Plummer_{tag}', rich, targets)

# ---- 3. Gustin (PRJNA691964) —— 基线 = visit 01 ----
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
af3 = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af3['run'] = [norm_sid(s) for s in af3['sample']]
af3 = af3.set_index('run')
g3 = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
g3.index = [norm_sid(i) for i in g3.index]
g_base = G_map[G_map['visit'].astype(float) == 1.0]['run'].tolist()
g_base = [norm_sid(s) for s in g_base]
rich3 = af3.loc[g_base]['species_observed'].astype(float)
targets3 = build_targets(g3.loc[g_base])
analyze('Gustin_CST4(all)', rich3, targets3, note='全队列即 CST4')

# ---- 4. Mollin —— 属级 richness（无种级），基线 CST4 ----
mg = pd.read_csv('<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
mp = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')
m_cst4 = mp[mp['base_cst4'] == 1]
idx_m = m_cst4['base_sample'].tolist()
rich_m = (mg.loc[idx_m] > 0).sum(axis=1).astype(float)  # genus_observed
targets_m = build_targets(mg.loc[idx_m])
analyze('Mollin_CST4', rich_m, targets_m, note='richness=genus_observed(无种级)')

# ---------- 输出 ----------
out = pd.DataFrame(rows)
out.to_csv('<REPO_DIR>/results/results_c8_richness_targets.csv', index=False)

var_df = pd.DataFrame(var_rows)
var_df.to_csv('<REPO_DIR>/results/results_c8_variance_decomp.csv', index=False)

# 跨队列方向一致性（仅 CST4 亚组）
cst4_cohorts = [c for c in out['cohort'].unique() if 'CST4' in c]
cons_rows = []
for t, vals in consistency.items():
    sub = out[(out['target'] == t) & (out['cohort'].isin(cst4_cohorts))]
    if len(sub) == 0:
        continue
    signs = np.sign(sub['rho'].values)
    signs_nz = signs[signs != 0]
    n_neg = int((sub['rho'] < 0).sum())
    n_pos = int((sub['rho'] > 0).sum())
    n_sig = int((sub['p'] < 0.05).sum())
    cons_rows.append(dict(target=t, n_cohorts=len(sub), n_neg=n_neg, n_pos=n_pos, n_sig=n_sig,
                          rho_mean=sub['rho'].mean(), rho_min=sub['rho'].min(), rho_max=sub['rho'].max(),
                          consistent='是' if (n_neg == len(sub) or n_pos == len(sub)) else '否'))
cons_df = pd.DataFrame(cons_rows).sort_values('rho_mean')
cons_df.to_csv('<REPO_DIR>/results/results_c8_direction_consistency.csv', index=False)

# ---------- 打印关键结论 ----------
print('===== #C8 richness x 靶标（基线 CST4 亚组）=====')
print()
piv = out[out['cohort'].isin(cst4_cohorts)].pivot_table(index='target', columns='cohort', values='rho', aggfunc='first')
print('【Spearman rho 矩阵（行=靶标，列=队列 CST4）】')
print(piv.round(3).to_string())
print()
print('【跨队列方向一致性】')
print(cons_df[['target', 'n_cohorts', 'n_neg', 'n_pos', 'n_sig', 'rho_mean', 'consistent']].round(3).to_string(index=False))
print()
print('【方差分解 R2（richness ~ 7 靶标，log1p 口径）】')
print(var_df[['cohort', 'n', 'R2_log1p', 'R2_log1p_adj', 'R2_2var_log1p']].round(3).to_string(index=False))
