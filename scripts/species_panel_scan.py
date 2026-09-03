# -*- coding: utf-8 -*-
"""
species_panel_scan.py —— #C5 任务 2 + 任务 3
=========================================================
任务 2：单物种精细信号扫查（三队列基线 CST4 亚组，逐物种单变量关联二元/连续结局）
任务 3：失调指数比值（乳杆菌占优度 = Lactobacillus / (Lactobacillus + 核心厌氧菌)）

队列：Xiao(454/PRJNA302078) · Plummer(MiSeq/PRJNA592384) · Gustin(MiSeq/PRJNA691964)
统一口径（与 four_cohort_matrix.py / outcome_refinement.py 严格一致）：
  - 二元结局 outcome = 随访 LB >= 50% (康复=1) / <50% (未康复=0)
  - 连续结局 fu_LB = log1p(随访 Lactobacillus 相对丰度)；ΔLB = fu_LB - base_LB
  - 基线 CST4 亚组 = base_LB < 0.5（Gustin 全队列即 CST4）
单变量检验：Spearman（二元/连续/ΔLB）+ ROC-AUC（二元，物种丰度作连续分数）
输出：<OUT_DIR>/modeling/results_species_scan.csv
      <OUT_DIR>/modeling/results_dysbiosis_index.csv
"""
import pandas as pd
import numpy as np
import re
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

OUT1 = '<OUT_DIR>/modeling/results_species_scan.csv'
OUT2 = '<OUT_DIR>/modeling/results_dysbiosis_index.csv'


def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))


# ================= 1. 加载三队列数据 =================
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')

# 各队列种级 + 属级相对丰度矩阵
SP, GE = {}, {}
for pid in ['PRJNA302078', 'PRJNA592384', 'PRJNA691964']:
    sp = pd.read_csv(f'<DATA_DIR>/feat/{pid}_species_rel.csv', index_col=0)
    sp.index = [norm_sid(i) for i in sp.index]
    ge = pd.read_csv(f'<DATA_DIR>/feat/{pid}_genus_rel.csv', index_col=0)
    ge.index = [norm_sid(i) for i in ge.index]
    SP[pid] = sp
    GE[pid] = ge

G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
G_map = G_map[G_map['control'].isna()]  # 去 control


# ================= 2. 构建统一 subject 表 =================
# 每行一个 subject：cohort, base_run, out_run, outcome(二元), cst4, fu_LB, base_LB
subj = []

# Xiao
xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
xiao['out'] = [norm_sid(s) for s in xiao['outcome_sample']]
gx = GE['PRJNA302078']
for _, r in xiao.iterrows():
    fu = gx.loc[r['out'], 'Lactobacillus']
    base = gx.loc[r['base'], 'Lactobacillus']
    subj.append(dict(cohort='Xiao', base_run=r['base'], out_run=r['out'],
                     outcome=int(r['outcome']), cst4=bool(r['base_cst4']),
                     fu_LB=np.log1p(fu), base_LB=np.log1p(base)))

# Plummer
plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
plum['out'] = [norm_sid(s) for s in plum['outcome_sample']]
gp = GE['PRJNA592384']
for _, r in plum.iterrows():
    fu = gp.loc[r['out'], 'Lactobacillus']
    base = gp.loc[r['base'], 'Lactobacillus']
    subj.append(dict(cohort='Plummer', base_run=r['base'], out_run=r['out'],
                     outcome=int(r['outcome']), cst4=bool(r['base_cst4']),
                     fu_LB=np.log1p(fu), base_LB=np.log1p(base)))

# Gustin：v1 -> v3 配对（与 four_cohort_matrix 一致）
gg = GE['PRJNA691964']
pat_v1 = G_map[G_map['visit'].astype(int) == 1].set_index('patient')['run']
pat_v3 = G_map[G_map['visit'].astype(int) == 3].set_index('patient')['run']
patients = sorted(set(pat_v1.index) & set(pat_v3.index))
for p in patients:
    b_run, o_run = pat_v1[p], pat_v3[p]
    fu = gg.loc[o_run, 'Lactobacillus']
    base = gg.loc[b_run, 'Lactobacillus']
    subj.append(dict(cohort='Gustin', base_run=b_run, out_run=o_run,
                     outcome=int(fu > 0.5), cst4=bool(base < 0.5),
                     fu_LB=np.log1p(fu), base_LB=np.log1p(base)))

S = pd.DataFrame(subj)
S['delta_LB'] = S['fu_LB'] - S['base_LB']

print('=== subject 表（全人群）===')
print(S.groupby('cohort').agg(n=('outcome', 'size'), rec=('outcome', 'sum'),
                              cst4=('cst4', 'sum')).to_string())
print()

# ================= 3. 物种/类群定义 =================
# 属级类群（genus_rel 直接取列）
GENUS_TAXA = {
    'Lactobacillus(属)': ['Lactobacillus'],
    'Megasphaera(属)': ['Megasphaera'],
    'Sneathia(属)': ['Sneathia'],
    'Mobiluncus(属)': ['Mobiluncus'],
}
# Prevotella 属聚合：所有含 'Prevotella' 的属级列
# 种级物种（species_rel 精确列名）
SPECIES_TAXA = {
    'G.vaginalis(种)': 'Gardnerella Gardnerella vaginalis',
    'A.vaginae(种)': 'Atopobium Atopobium vaginae DSM 15829',
}

# 缺失物种（SILVA 种级注释无法区分 / 未注释）——单列说明，不参与扫查
MISSING = ['L.iners', 'L.crispatus', 'P.bivia', 'BVAB1-3']


def get_feature(pid, base_runs):
    """返回 DataFrame（index=base_run，列为各 taxa 丰度），全部用相对丰度原始值"""
    ge = GE[pid]
    sp = SP[pid]
    feat = pd.DataFrame(index=base_runs)
    # 属级类群
    for tax, cols in GENUS_TAXA.items():
        use = [c for c in cols if c in ge.columns]
        feat[tax] = ge.loc[base_runs, use].sum(axis=1).values if use else 0.0
    # Prevotella 属聚合
    pv = [c for c in ge.columns if 'Prevotella' in c]
    feat['Prevotella(属聚合)'] = ge.loc[base_runs, pv].sum(axis=1).values if pv else 0.0
    # 种级物种
    for tax, col in SPECIES_TAXA.items():
        if col in sp.columns:
            feat[tax] = sp.loc[base_runs, col].values
        else:
            feat[tax] = 0.0
    return feat


def auc_safe(x, yv):
    if len(np.unique(yv)) < 2 or np.all(x == x[0]):
        return np.nan
    try:
        return roc_auc_score(yv, x)
    except Exception:
        return np.nan


# ================= 4. 任务 2：单物种扫查（CST4 亚组内） =================
rows2 = []
for cohort in ['Xiao', 'Plummer', 'Gustin']:
    sub = S[S['cohort'] == cohort].reset_index(drop=True)
    m = sub['cst4'].values.astype(bool)
    subc = sub[m].reset_index(drop=True)
    if len(subc) < 3:
        continue
    base_runs = subc['base_run'].tolist()
    feat = get_feature({'Xiao': 'PRJNA302078', 'Plummer': 'PRJNA592384',
                        'Gustin': 'PRJNA691964'}[cohort], base_runs)
    yv = subc['outcome'].values.astype(int)
    fu = subc['fu_LB'].values
    dl = subc['delta_LB'].values
    for tax in feat.columns:
        x = feat[tax].values.astype(float)
        r_b, p_b = spearmanr(x, yv)
        r_c, p_c = spearmanr(x, fu)
        r_d, p_d = spearmanr(x, dl)
        a = auc_safe(x, yv)
        rows2.append(dict(cohort=cohort, taxon=tax, n=len(subc),
                          n_rec=int(yv.sum()),
                          rho_binary=round(r_b, 4), p_binary=round(p_b, 4),
                          auc_binary=round(a, 4) if not np.isnan(a) else np.nan,
                          rho_cont=round(r_c, 4), p_cont=round(p_c, 4),
                          rho_delta=round(r_d, 4), p_delta=round(p_d, 4)))

res2 = pd.DataFrame(rows2)

# BH-FDR 校正（对 p_binary / p_cont / p_delta 分别，按全部 taxa×cohort 联合）
from statsmodels.stats.multitest import multipletests
for col, qcol in [('p_binary', 'q_binary'), ('p_cont', 'q_cont'), ('p_delta', 'q_delta')]:
    pvals = res2[col].values
    valid = ~np.isnan(pvals)
    q = np.full(len(pvals), np.nan)
    if valid.sum() > 0:
        q[valid] = multipletests(pvals[valid], method='fdr_bh')[1]
    res2[qcol] = np.round(q, 4)

res2.to_csv(OUT1, index=False)
print('=== 任务2 单物种扫查（三队列基线 CST4 亚组）===')
print(res2.to_string(index=False))
print()

# ================= 5. 任务 3：失调指数（乳杆菌占优度） =================
# 核心厌氧菌属（三队列属级均有）：Gardnerella/Atopobium/Megasphaera/Sneathia/Mobiluncus/
#   Dialister/Porphyromonas/Anaerococcus/Aerococcus + 所有 Prevotella 列
ANAEROBES = ['Gardnerella', 'Atopobium', 'Megasphaera', 'Sneathia', 'Mobiluncus',
             'Dialister', 'Porphyromonas', 'Anaerococcus', 'Aerococcus']

rows3 = []
for cohort in ['Xiao', 'Plummer', 'Gustin']:
    ge = GE[{'Xiao': 'PRJNA302078', 'Plummer': 'PRJNA592384',
             'Gustin': 'PRJNA691964'}[cohort]]
    sub = S[S['cohort'] == cohort].reset_index(drop=True)
    base_runs = sub['base_run'].tolist()
    lacto = ge.loc[base_runs, 'Lactobacillus'].values.astype(float)
    ana = [c for c in ANAEROBES if c in ge.columns]
    pv = [c for c in ge.columns if 'Prevotella' in c]
    ana_sum = ge.loc[base_runs, ana + pv].sum(axis=1).values.astype(float)
    dom = lacto / (lacto + ana_sum + 1e-12)  # 乳杆菌占优度 [0,1]
    yv = sub['outcome'].values.astype(int)
    fu = sub['fu_LB'].values
    dl = sub['delta_LB'].values

    for tag, m in [('full', np.ones(len(yv), bool)), ('CST4', sub['cst4'].values.astype(bool))]:
        if m.sum() < 3:
            continue
        d, yy = dom[m], yv[m]
        r_b, p_b = spearmanr(d, yy)
        r_c, p_c = spearmanr(d, fu[m])
        r_d, p_d = spearmanr(d, dl[m])
        a = auc_safe(d, yy)
        rows3.append(dict(cohort=cohort, subgroup=tag, n=int(m.sum()),
                          n_rec=int(yy.sum()),
                          rho_binary=round(r_b, 4), p_binary=round(p_b, 4),
                          auc_binary=round(a, 4) if not np.isnan(a) else np.nan,
                          rho_cont=round(r_c, 4), p_cont=round(p_c, 4),
                          rho_delta=round(r_d, 4), p_delta=round(p_d, 4)))

res3 = pd.DataFrame(rows3)
res3.to_csv(OUT2, index=False)
print('=== 任务3 乳杆菌占优度（Lactobacillus/(Lactobacillus+核心厌氧菌)）===')
print(res3.to_string(index=False))
print()

# ================= 6. 缺失物种说明 =================
print('=== 种级注释缺失物种（16S 分辨率不足，恰是 PCR+CE 差异化卖点）===')
for sp_ in MISSING:
    print(f'  - {sp_}: 三队列 SILVA v138.1 种级注释均未出现（{sp_}）')
print()
print(f'已保存: {OUT1}')
print(f'已保存: {OUT2}')
