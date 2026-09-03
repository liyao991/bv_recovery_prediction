# gustin_fragility.py
# #C2 任务二：Gustin 小样本脆弱性模拟（n=33）
# 目的：证明 rho=-0.498 在 n=33 下是脆弱的，主动坐实"探索性、非稳健"
# 方法：Bootstrap 1000 / Leave-one-out / 置换检验 1000
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr

def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))

rng = np.random.default_rng(42)
SEED_SAVE = '<OUT_DIR>/modeling/results_gustin_fragility.csv'

# ---- 提取 Gustin richness + 结局（复用 richness_direction_check.py 口径）----
af = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af['run'] = [norm_sid(s) for s in af['sample']]
af = af.set_index('run')
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
y = (G_rel.loc[list(g_out)]['Lactobacillus'] > np.log1p(0.5)).values.astype(int)
rich = af.loc[list(g_base)]['species_observed'].values
n = len(y)
rho_obs, p_obs = spearmanr(rich, y)

print(f'Gustin n={n}, 复发={y.sum()}, 康复={(1-y).sum()}')
print(f'观测 rho={rho_obs:.4f}, p={p_obs:.4f}')

# ---- 1. Bootstrap (1000 次) ----
B = 1000
bs_rho = np.empty(B); bs_p = np.empty(B)
for b in range(B):
    idx = rng.integers(0, n, n)
    r, p = spearmanr(rich[idx], y[idx])
    bs_rho[b] = r; bs_p[b] = p
ci_lo, ci_hi = np.percentile(bs_rho, [2.5, 97.5])
frac_sig = (bs_p < 0.05).mean()
print(f'\n[1] Bootstrap {B}次:')
print(f'    rho 95% CI = [{ci_lo:.3f}, {ci_hi:.3f}]  (跨度 {ci_hi-ci_lo:.3f})')
print(f'    p<0.05 比例 = {frac_sig*100:.1f}%')

# ---- 2. Leave-one-out / jackknife ----
loo_rho = np.empty(n); loo_p = np.empty(n)
for i in range(n):
    mask = np.ones(n, bool); mask[i] = False
    r, p = spearmanr(rich[mask], y[mask])
    loo_rho[i] = r; loo_p[i] = p
# 找影响最大的样本
delta = loo_rho - rho_obs
order = np.argsort(-np.abs(delta))
print(f'\n[2] Leave-one-out (jackknife):')
print(f'    LOO rho 范围 = [{loo_rho.min():.3f}, {loo_rho.max():.3f}]')
print(f'    去掉后失显著(p>=0.05)的样本数 = {(loo_p>=0.05).sum()}/{n}')
print(f'    去掉后 rho 翻转符号的样本数 = {((loo_rho*rho_obs)<0).sum()}/{n}')
print(f'    影响最大的 5 个样本（样本索引 -> 去掉后 rho, p）:')
sensitive = []
for i in order[:5]:
    sensitive.append(dict(sample_idx=int(i), loo_rho=round(float(loo_rho[i]),3),
                          loo_p=round(float(loo_p[i]),3),
                          delta_rho=round(float(delta[i]),3),
                          rich=round(float(rich[i]),1), y=int(y[i])))
    print(f'      #{i}: rho={loo_rho[i]:.3f}, p={loo_p[i]:.3f}, Δrho={delta[i]:+.3f}, rich={rich[i]:.0f}, y={y[i]}')

# ---- 3. 置换检验 (1000 次) ----
P = 1000
perm_rho = np.empty(P)
for k in range(P):
    y_perm = rng.permutation(y)
    perm_rho[k], _ = spearmanr(rich, y_perm)
perm_p = (np.abs(perm_rho) >= np.abs(rho_obs)).mean()
print(f'\n[3] 置换检验 {P}次:')
print(f'    置换 p = {perm_p:.4f}  (对照 #B2 的 p={0.037:.3f})')

# ---- 稳定性评级 ----
if ci_lo < 0 < ci_hi:
    rating = '脆弱（95% CI 跨 0，方向不稳定）'
elif frac_sig < 0.5:
    rating = '弱稳健（CI 不跨 0 但 p<0.05 比例 <50%）'
else:
    rating = '稳健（CI 不跨 0 且多数 bootstrap 显著）'

print(f'\n=== 一句话结论 ===')
print(f'rho={rho_obs:.3f} 在 n={n} 下评级：{rating}')
print(f'  - bootstrap 95% CI [{ci_lo:.3f},{ci_hi:.3f}] {"跨0" if ci_lo<0<ci_hi else "不跨0"}')
print(f'  - LOO 失显著 {(loo_p>=0.05).sum()}/{n} 例，符号翻转 {((loo_rho*rho_obs)<0).sum()}/{n} 例')
print(f'  - 置换 p={perm_p:.3f}')

# ---- 输出 ----
out = pd.DataFrame({
    'metric': ['n', 'n_cases', 'n_controls', 'rho_obs', 'p_obs',
               'bs_ci_lo', 'bs_ci_hi', 'bs_frac_p_lt_0.05',
               'loo_n_lose_sig', 'loo_n_flip_sign', 'perm_p', 'rating'],
    'value': [n, int(y.sum()), int((1-y).sum()), rho_obs, p_obs,
              ci_lo, ci_hi, frac_sig,
              int((loo_p>=0.05).sum()), int(((loo_rho*rho_obs)<0).sum()), perm_p, rating]
})
out.to_csv(SEED_SAVE, index=False)
print(f'\n已写入 {SEED_SAVE}')

# 敏感样本清单单独存
pd.DataFrame(sensitive).to_csv('<OUT_DIR>/modeling/results_gustin_loo_sensitive.csv', index=False)
print('已写入 <OUT_DIR>/modeling/results_gustin_loo_sensitive.csv')
