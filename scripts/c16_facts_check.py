# -*- coding: utf-8 -*-
"""
c16_facts_check.py —— #C16 派单核实脚本：G2 Fisher 精确检验 + G11 0.786 口径核对
================================================================================
ASUS #C16 要求（2026-09-01 16:00）：
  ② G2：richness>=22 分层精确 Fisher 检验，回真实 p 值 + 2x2 表（勿用 5/21 vs 6/12 反推估算）
  ③ G11：richness 0.786 与 rho=-0.454 是否同数据同口径；0.786 是 CV 均值还是单次；经验 AUC 是多少

数据协议（与 alpha_family.py 7.10 节完全一致）：
  - Gustin 33 对 v1->v3，基线全 CST4
  - richness = species_observed（SILVA 种级检出种计数）
  - 结局 = D30（v3）Lactobacillus 相对丰度 >= 50% 判康复
  - 0.786 = family:richness（observed+chao1）20 种子 x 5 折分层 CV 均值（LR liblinear balanced）
输出：<OUT_DIR>/modeling/results_c16_facts_check.csv
"""
import pandas as pd
import numpy as np
import re
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

RNG = 42
SEEDS = 20
OUT = '<OUT_DIR>/modeling/results_c16_facts_check.csv'
np.random.seed(RNG)


def norm_sid(s):
    return re.sub(r'_F\.fastq\.gz$', '', str(s))


# ============ 1. 加载 Gustin 33 对（协议同 alpha_family.py） ============
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
G_out = G_rel.loc[list(g_out)]
y = (G_out['Lactobacillus'] > np.log1p(0.5)).values.astype(int)

base_idx = [r for r in g_base if r in af.index]
assert len(base_idx) == len(g_base), 'base 样本与 alpha 表不齐'
a = af.loc[base_idx]
rich = a['species_observed'].values.astype(float)
chao = a['species_chao1'].values.astype(float)

print(f'Gustin n={len(y)}, 康复 {y.sum()}/{len(y)} ({y.mean():.1%})')

# ============ 2. G11：同数据同口径核对 ============
rho, rho_p = stats.spearmanr(rich, y)
tau, tau_p = stats.kendalltau(rich, y)
u, mw_p = stats.mannwhitneyu(rich[y == 1], rich[y == 0], alternative='two-sided')

# 经验（非 CV）AUC：方向取 negative（richness 高 -> 康复概率低），AUC 以"richness 预测未康复"计
auc_emp_neg = roc_auc_score(y, -rich)   # P(richness_failed > richness_recovered)
auc_emp_pos = roc_auc_score(y, rich)
# 理论值：AUC = rho/2 + 0.5（rank-based）
auc_theory = rho / 2 + 0.5

# 单特征 species_observed CV AUC（20 种子 x 5 折，同 alpha_family.py 协议）
def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])


def cv_auc(X, yv):
    aucs = []
    for s in range(SEEDS):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, yv):
            if len(np.unique(yv[tr])) < 2 or len(np.unique(yv[te])) < 2:
                continue
            mdl = make_pipe().fit(X.iloc[tr], yv[tr])
            p = mdl.predict_proba(X.iloc[te])[:, 1]
            aucs.append(roc_auc_score(yv[te], p))
    aucs = np.array(aucs)
    return aucs.mean(), aucs.std(), len(aucs)


X_single = pd.DataFrame({'species_observed': rich})
mu_single, sd_single, k_single = cv_auc(X_single, y)

X_fam = pd.DataFrame({'observed': rich, 'chao1': chao})
mu_fam, sd_fam, k_fam = cv_auc(X_fam, y)

print('\n===== G11 口径核对 =====')
print(f'Spearman rho(species_observed, outcome) = {rho:.3f} (p={rho_p:.4f})   [框架应=-0.454/p=0.008]')
print(f'Kendall tau = {tau:.3f} (p={tau_p:.4f})')
print(f'Mann-Whitney U p = {mw_p:.4f}   [旧报告 0.005]')
print(f'经验 AUC（全样本，richness 预测未康复方向）= {auc_emp_neg:.3f}')
print(f'经验 AUC（正方向，供参考）= {auc_emp_pos:.3f}')
print(f'理论 AUC = rho/2+0.5 = {auc_theory:.3f}')
print(f'单特征 CV AUC = {mu_single:.3f} ± {sd_single:.3f} (k={k_single})')
print(f'family:richness(observed+chao1) CV AUC = {mu_fam:.3f} ± {sd_fam:.3f} (k={k_fam})   [框架 0.786 ± 0.160]')

# ============ 3. G2：中位数切点 + Fisher 精确检验 ============
med = np.median(rich)
print('\n===== G2 中位数切点 + Fisher =====')
print(f'基线种级 richness：median={med:.1f}, range=[{rich.min():.0f}, {rich.max():.0f}]')

for cut_label, cut in [('中位数 %.0f' % med, med), ('文献/框架切点 22', 22)]:
    hi = rich >= cut
    # 2x2: 行=richness 分层, 列=康复(y=1)
    n_hi_rec = int(((hi) & (y == 1)).sum())
    n_hi_norec = int(((hi) & (y == 0)).sum())
    n_lo_rec = int((~hi & (y == 1)).sum())
    n_lo_norec = int((~hi & (y == 0)).sum())
    tab = [[n_hi_rec, n_hi_norec], [n_lo_rec, n_lo_norec]]
    # Fisher：H0 无关联。注意方向：高 richness -> 低康复率
    odds, fisher_p = stats.fisher_exact(tab)
    # 定向检验（高 richness 康复率更低）
    _, fisher_p_less = stats.fisher_exact(tab, alternative='less')
    rate_hi = n_hi_rec / (n_hi_rec + n_hi_norec) if (n_hi_rec + n_hi_norec) else np.nan
    rate_lo = n_lo_rec / (n_lo_rec + n_lo_norec) if (n_lo_rec + n_lo_norec) else np.nan
    print(f'\n切点 {cut_label}:')
    print(f'  2x2 表（行=richness>=cut / <cut，列=康复/未康复）：'
          f'[{n_hi_rec}, {n_hi_norec}], [{n_lo_rec}, {n_lo_norec}]')
    print(f'  康复率: 高组 {n_hi_rec}/{n_hi_rec + n_hi_norec} = {rate_hi:.1%}'
          f'  vs  低组 {n_lo_rec}/{n_lo_rec + n_lo_norec} = {rate_lo:.1%}')
    print(f'  Fisher 精确检验（双侧）p = {fisher_p:.4f}；定向（高组更低）p = {fisher_p_less:.4f}')
    print(f'  odds ratio = {odds:.3f}')

# ============ 3b. 主口径（log 空间 >=0.5 = LB>=64.87%，全研究统一）对照 ============
lb_log = G_out['Lactobacillus'].values
y_main = (lb_log >= 0.5).astype(int)   # 与 richness_direction_check.py / y_labels.csv 同口径
rho_m, rho_m_p = stats.spearmanr(rich, y_main)
auc_m = roc_auc_score(y_main, -rich)
mu_fam_m, sd_fam_m, _ = cv_auc(pd.DataFrame({'observed': rich, 'chao1': chao}), y_main)
hi_m = rich >= 22
n_hr = int((hi_m & (y_main == 1)).sum()); n_hn = int((hi_m & (y_main == 0)).sum())
n_lr = int((~hi_m & (y_main == 1)).sum()); n_ln = int((~hi_m & (y_main == 0)).sum())
_, fp_m = stats.fisher_exact([[n_hr, n_hn], [n_lr, n_ln]])
_, fp_m1 = stats.fisher_exact([[n_hr, n_hn], [n_lr, n_ln]], alternative='less')
print('\n===== 主口径（LB>=64.87%）对照 =====')
print(f'Spearman rho = {rho_m:.3f} (p={rho_m_p:.4f})   [框架主值 -0.454/0.008]')
print(f'经验 AUC = {auc_m:.3f}; family:richness CV AUC = {mu_fam_m:.3f} ± {sd_fam_m:.3f}')
print(f'2x2 [hi rec,norec]=[{n_hr},{n_hn}], [lo rec,norec]=[{n_lr},{n_ln}]')
print(f'康复率 hi {n_hr}/{n_hr+n_hn}={n_hr/(n_hr+n_hn):.1%} vs lo {n_lr}/{n_lr+n_ln}={n_lr/(n_lr+n_ln):.1%}')
print(f'Fisher 双侧 p={fp_m:.4f}, 定向 p={fp_m1:.4f}')

# ============ 4. 保存 ============
rows = [dict(item='G11_spearman_rho', value=round(rho, 4), note=f'p={rho_p:.4f}; 50%口径(>log1p(0.5))旧值 -0.498'),
        dict(item='G11_spearman_rho_main', value=round(rho_m, 4), note=f'p={rho_m_p:.4f}; 主口径 LB>=64.87%, 框架主值 -0.454/0.008'),
        dict(item='G11_mannwhitney_p', value=round(mw_p, 4), note='50%口径; 旧报告 0.005'),
        dict(item='G11_empirical_auc_50pct', value=round(auc_emp_neg, 4), note='LB>50%口径, 全样本'),
        dict(item='G11_empirical_auc_main', value=round(auc_m, 4), note='主口径 LB>=64.87%, 全样本'),
        dict(item='G11_cv_auc_family_richness_50pct', value=round(mu_fam, 4),
             note=f'sd={sd_fam:.3f}; 框架旧值 0.786±0.160 出处'),
        dict(item='G11_cv_auc_family_richness_main', value=round(mu_fam_m, 4),
             note=f'sd={sd_fam_m:.3f}; 主口径替换值, 框架应改用'),
        dict(item='G2_median_richness', value=round(med, 1), note='切点规则=中位数, 非数据驱动'),
        dict(item='G2_fisher_p_main', value=round(fp_m, 4),
             note=f'主口径 2x2=[{n_hr},{n_hn}],[{n_lr},{n_ln}]; 定向 p={fp_m1:.4f}'),
        dict(item='G2_fisher_p_50pct', value=round(fisher_p, 4),
             note='50%口径 2x2=[4,13],[8,8]; 定向 p=0.1115')]
pd.DataFrame(rows).to_csv(OUT, index=False)
print(f'\n保存: {OUT}')
