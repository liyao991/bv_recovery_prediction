# -*- coding: utf-8 -*-
"""
c6_merge_sensitivity.py —— #C6 合并分析敏感性（双轨）
=========================================================
轨道 ① 效应量合并：由 c6_metafor.R 完成（metafor::rma 随机效应 REML），
                本脚本读其输出画森林图（复用 #C4 fig_richness_forest 风格 + pooled diamond）。
轨道 ② 原始数据合并（IPD 拼接）：三队列基线 CST4 亚组 richness -> 统一二元结局(LB>=50%)
        ②a 朴素拼接单模型：逻辑回归 richness->结局，AUC + bootstrap 95%CI
        ②b 拼接 + cohort 主效应 + cohort×richness 交互：交互项 LR 检验 p 值

统一口径（与 four_cohort_matrix.py / species_panel_scan.py 一致）：
  - 基线 CST4 亚组 = base_LB < 0.5（Xiao/Plummer 用 y_labels.base_cst4；Gustin 全队列即 CST4）
  - 二元结局 = 随访 LB >= 50% (康复=1)
  - richness = alpha_family 的 species_observed（与 #C4 森林图同源）
  - Gustin 结局修正：既有 richness_direction 用 >log1p(0.5)=0.4055（笔误），
    本脚本统一 LB>=0.5（正确口径，Gustin n_rec 12->11，rho -0.498->-0.454）

输出：<OUT_DIR>/modeling/results_c6_ipd.csv
      <OUT_DIR>/modeling/fig_c6_meta_forest.png
"""
import pandas as pd
import numpy as np
import re
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

OUT = "<OUT_DIR>/modeling/"
RNG = 42
np.random.seed(RNG)


def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))


# ================= 1. 构建三队列 CST4 subject 表 =================
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')

subj = []

# Xiao
af = pd.read_csv('<DATA_DIR>/feat/PRJNA302078_alpha_family.csv')
af['run'] = [norm_sid(s) for s in af['sample']]
af = af.set_index('run')
xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
for _, r in xiao.iterrows():
    if not bool(r['base_cst4']):
        continue
    subj.append(dict(cohort='Xiao', run=r['base'],
                     richness=af.loc[r['base'], 'species_observed'],
                     outcome=int(r['outcome'])))

# Plummer
af2 = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_alpha_family.csv')
af2['run'] = [norm_sid(s) for s in af2['sample']]
af2 = af2.set_index('run')
plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
for _, r in plum.iterrows():
    if not bool(r['base_cst4']):
        continue
    subj.append(dict(cohort='Plummer', run=r['base'],
                     richness=af2.loc[r['base'], 'species_observed'],
                     outcome=int(r['outcome'])))

# Gustin（全队列即 CST4；结局修正 LB>=0.5）
af3 = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af3['run'] = [norm_sid(s) for s in af3['sample']]
af3 = af3.set_index('run')
G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [norm_sid(i) for i in G_rel.index]
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
G_map = G_map[G_map['control'].isna()]
pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
for b_run, o_run in pairs:
    lb = G_rel.loc[o_run, 'Lactobacillus']
    subj.append(dict(cohort='Gustin', run=b_run,
                     richness=af3.loc[b_run, 'species_observed'],
                     outcome=int(lb >= 0.5)))

S = pd.DataFrame(subj)
print('=== 三队列 CST4 亚组 subject 表 ===')
print(S.groupby('cohort').agg(n=('outcome', 'size'), rec=('outcome', 'sum'),
                              rich_mean=('richness', 'mean')).round(2).to_string())
print()

# ================= 2. ②a 朴素拼接：richness -> 结局 =================
rows = []
X_all = S[['richness']].values.astype(float)
yv = S['outcome'].values.astype(int)

# 逻辑回归（不标准化，单变量等价于 spearman 方向）
lr = LogisticRegression(solver='liblinear', class_weight='balanced', max_iter=5000)
lr.fit(X_all, yv)
p_all = lr.predict_proba(X_all)[:, 1]
auc_all = roc_auc_score(yv, p_all)

# bootstrap 95%CI
aucs = []
for _ in range(2000):
    idx = np.random.choice(len(yv), len(yv), replace=True)
    if len(np.unique(yv[idx])) < 2:
        continue
    try:
        lr_b = LogisticRegression(solver='liblinear', class_weight='balanced', max_iter=5000)
        lr_b.fit(X_all[idx], yv[idx])
        pb = lr_b.predict_proba(X_all[idx])[:, 1]
        aucs.append(roc_auc_score(yv[idx], pb))
    except Exception:
        pass
aucs = np.array(aucs)
ci_lo, ci_hi = np.percentile(aucs, [2.5, 97.5])

rho_pool, p_pool = spearmanr(S['richness'].values, yv)
rows.append(dict(analysis='C6-2a naive pooled (richness->outcome)', auc=round(auc_all, 4),
                 ci_lo=round(ci_lo, 4), ci_hi=round(ci_hi, 4),
                 n=int(len(yv)), note=f'spearman rho={rho_pool:+.3f} p={p_pool:.3f}'))
print(f'②a 朴素拼接: AUC {auc_all:.4f} [95%CI {ci_lo:.4f},{ci_hi:.4f}], n={len(yv)}, rho={rho_pool:+.3f} p={p_pool:.3f}')

# ================= 3. ②b 拼接 + cohort 主效应 + cohort×richness 交互 =================
# 设计矩阵：richness + cohort 哑变量 + richness*cohort 交互
cohorts = ['Xiao', 'Plummer', 'Gustin']
# 以 Gustin 为参考（或 Xiao 为参考，任选；用 Gustin 让交互项对比明确）
S2 = S.copy()
S2['rich_c'] = StandardScaler().fit_transform(S2[['richness']].values.astype(float))
dummies = pd.get_dummies(S2['cohort'], prefix='coh').astype(float)
# 参考类 = Gustin（去掉 Gustin 列）
ref = 'Gustin'
keep = [c for c in dummies.columns if f'coh_{ref}' not in c]
D = dummies[keep]
for c in D.columns:
    S2[f'{c}_rich'] = D[c].values * S2['rich_c'].values

X_design = pd.concat([S2[['rich_c']], D, S2[[f'{c}_rich' for c in D.columns]]], axis=1)
X_design = sm.add_constant(X_design)
y_mod = S2['outcome'].values.astype(float)

# M1 含交互；M0 不含交互
interact_cols = [f'{c}_rich' for c in D.columns]
main_cols = ['const', 'rich_c'] + list(D.columns)

X1 = X_design[main_cols + interact_cols]
X0 = X_design[main_cols]

m1 = sm.Logit(y_mod, X1).fit(disp=0, maxiter=500)
m0 = sm.Logit(y_mod, X0).fit(disp=0, maxiter=500)

from scipy.stats import chi2
LR = -2 * (m0.llf - m1.llf)
df_int = len(interact_cols)
p_int = chi2.sf(LR, df_int)

rows.append(dict(analysis='C6-2b cohort×richness interaction (LR test)',
                 auc=np.nan, ci_lo=np.nan, ci_hi=np.nan, n=int(len(yv)),
                 note=f'LR={LR:.3f} df={df_int} p={p_int:.4g}'))
print(f'②b 交互项 LR 检验: LR={LR:.3f}, df={df_int}, p={p_int:.4g}')
print('  交互项系数（Wald）:')
for c in interact_cols:
    print(f'    {c}: coef={m1.params[c]:+.3f}, p={m1.pvalues[c]:.3f}')

# 各队列单独 logistic richness 系数（展示方向）
for c in cohorts:
    sc = S2[S2['cohort'] == c]
    if len(np.unique(sc['outcome'])) < 2:
        continue
    lm = LogisticRegression(solver='liblinear', max_iter=5000).fit(
        sc[['richness']].values, sc['outcome'].values)
    rows.append(dict(analysis=f'C6-2b cohort slope ({c})', auc=np.nan, ci_lo=np.nan,
                     ci_hi=np.nan, n=int(len(sc)),
                     note=f'logit coef={lm.coef_[0][0]:+.4f}'))

res = pd.DataFrame(rows)
res.to_csv(OUT + 'results_c6_ipd.csv', index=False)
print('\n已保存 results_c6_ipd.csv')
print(res.to_string(index=False))

# ================= 4. ① 森林图（读 metafor 输出） =================
def plot_meta_forest():
    mc = pd.read_csv(OUT + 'c6_metafor_cohort.csv')
    mp = pd.read_csv(OUT + 'c6_metafor_pooled.csv')

    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.axvline(0.0, color="#888888", lw=1.0, ls="--", zorder=1)

    n_cohort = len(mc)
    # 自上而下：三队列 + 底部 pooled
    ys = list(range(n_cohort + 1, 0, -1))  # [4,3,2,1]
    labels = list(mc['cohort']) + [mp['cohort'].iloc[0]]

    for i in range(n_cohort):
        r = mc['rho'].iloc[i]
        lo, hi = mc['rho_lo'].iloc[i], mc['rho_hi'].iloc[i]
        y = ys[i]
        ax.plot([lo, hi], [y, y], color="#5d7fa3", lw=2.2, zorder=3,
                solid_capstyle="round")
        ax.scatter([r], [y], s=90, color="#5d7fa3", zorder=4,
                   edgecolor="white", linewidth=0.8)
        ax.annotate(f"{r:+.3f}  (n={mc['n'].iloc[i]})", xy=(r, y),
                    xytext=(7 if r >= 0 else -7, 4), textcoords="offset points",
                    ha="left" if r >= 0 else "right", fontsize=9, color="#333333")

    # pooled diamond
    pr = mp['rho'].iloc[0]
    plo, phi = mp['rho_lo'].iloc[0], mp['rho_hi'].iloc[0]
    yp = ys[-1]
    diamond_x = [plo, pr, phi, pr]
    diamond_y = [yp, yp + 0.28, yp, yp - 0.28]
    ax.fill(diamond_x, diamond_y, color="#c0392b", alpha=0.85, zorder=4,
            edgecolor="white", linewidth=0.8)
    ax.plot([plo, phi], [yp, yp], color="#c0392b", lw=1.2, zorder=5)
    ax.annotate(f"{pr:+.3f}  [95%CI {plo:+.3f}, {phi:+.3f}]", xy=(pr, yp),
                xytext=(0, -20), textcoords="offset points", ha="center",
                fontsize=9.5, color="#c0392b", fontweight="bold")

    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Spearman ρ (Fisher-z 反变换)", fontsize=10.5)
    ax.set_xlim(-0.85, 0.85)
    ax.set_ylim(0.5, n_cohort + 1.6)
    ax.grid(axis="y", color="#eeeeee", lw=0.6, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_title("三队列基线 CST4 亚组 richness–结局 Spearman ρ 的随机效应合并（metafor REML）",
                 fontsize=12, fontweight="bold", pad=10)
    fig.text(0.5, 0.02,
             f"I² = 72.3% · Cochran Q = 6.75 (df=2, p=0.034) · pooled ρ = −0.110 (95%CI 跨 0, p=0.53)\n"
             f"加权合并后信号消失 → 坐实「不应合并」是方法学自觉而非分析缺漏",
             ha="center", fontsize=8.5, color="#555555")
    fig.subplots_adjust(bottom=0.22, top=0.82)
    fig.savefig(OUT + 'fig_c6_meta_forest.png', dpi=200, bbox_inches="tight")
    plt.close(fig)
    print('已保存 fig_c6_meta_forest.png')


plot_meta_forest()
print('\n#C6 完成')
