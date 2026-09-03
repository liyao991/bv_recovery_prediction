# -*- coding: utf-8 -*-
"""
four_cohort_matrix.py —— 4 队列跨平台验证矩阵
================================================
队列：Xiao(454/DADA2/SILVA) · Plummer(MiSeq/DADA2/SILVA) · Mollin(MiSeq/Mothur/RDP) · Gustin(MiSeq/DADA2/SILVA)
统一结局：outcome = 随访 LB>=50%(康复,y=1) / LB<50%(未康复,y=0)；基线 CST4 = base_LB<0.5
两组问题：
  A. 原假说（全人群）：治疗前 16S 属级特征能否预测治疗结局（Xiao 训练 -> 各队列外推）
  B. 新叙事（基线 CST4 亚组）：同为失调人群内预测力是否 ≈ 随机（各队列 CST4 亚组 AUC≈0.5 复现）
方法：20 种子 LASSO（内层 5 折 CV 选 C），内部参考用分层 5 折 CV / 重复 holdout
输出：<OUT_DIR>/modeling/results_four_cohort_matrix.csv
"""
import pandas as pd
import numpy as np
import re
from collections import Counter
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV, RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score

RNG = 42
SEEDS = 20
OUT = '<OUT_DIR>/modeling/results_four_cohort_matrix.csv'
np.random.seed(RNG)

# ================= 1. 加载数据 =================
# Xiao（沿用 X_features/y_labels，与先前分析严格一致）
X_all = pd.read_csv('<DATA_DIR>/feat/X_features.csv')
y_all = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
feat_meta = pd.read_csv('<DATA_DIR>/feat/feature_meta.csv')
genus_cols_x = feat_meta[feat_meta['source'] == 'genus']['feature'].tolist()
alpha_cols_x = feat_meta[feat_meta['source'] == 'alpha']['feature'].tolist()
xiao_mask = (y_all['cohort'] == 'PRJNA302078').values
y_xiao = y_all[xiao_mask].reset_index(drop=True)
X_xiao_log1p = X_all[xiao_mask].reset_index(drop=True)

# Mollin
M_rel = pd.read_csv('<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
M_pairs = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')

# Gustin（属矩阵 + run 映射）
G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [re.sub(r'_F\.fastq\.gz$', '', str(i)) for i in G_rel.index]  # SRRxxx_F.fastq.gz -> SRRxxx
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')

# ================= 2. 属名标准化 =================
def norm_genus(name):
    s = str(name).strip()
    if s in ('NA_unclassified', 'Unassigned', 'nan', 'NA', 'NaN', ''):
        return 'NA_unclassified'
    s = s.replace('[', '').replace(']', '')
    if s.lower().startswith('uncultured') or s.lower().startswith('unclassified') or s == 'Bacteria':
        return 'NA_unclassified'
    first = s.split(' ')[0]
    if first in ('Sub', 'sp', 'group', 'incertae', 'Family', 'Order', 'Class', 'Phylum'):
        return 'NA_unclassified'
    return first

def map_cols(cols):
    m = {}
    for c in cols:
        m[c] = norm_genus(c)
    return m

def agg_rel_to_log1p(rel_df, col_map):
    df = rel_df.copy()
    df.columns = [col_map[c] for c in df.columns]
    agg = df.T.groupby(level=0).sum().T
    return np.log1p(agg)

# 各侧相对丰度 -> log1p（属名归一）
xiao_map = map_cols(genus_cols_x)
X_rel = np.expm1(X_xiao_log1p[genus_cols_x]).clip(lower=0)
X_rel = pd.DataFrame(X_rel, columns=genus_cols_x, index=X_xiao_log1p.index)
X_log = agg_rel_to_log1p(X_rel, xiao_map)

mollin_map = map_cols(M_rel.columns)
M_log = agg_rel_to_log1p(M_rel, mollin_map)

gustin_map = map_cols(G_rel.columns)
G_log = agg_rel_to_log1p(G_rel, gustin_map)

inter = sorted(set(X_log.columns) & set(M_log.columns) & set(G_log.columns))
print(f'三队列属名交集: {len(inter)}')
for g in ['Lactobacillus', 'Gardnerella', 'Atopobium', 'Mobiluncus', 'Prevotella',
          'Megasphaera', 'Sneathia', 'Fannyhessea', 'Anaerococcus']:
    print(f'  {g}: {"Y" if g in inter else "N"}')

X_f = X_log[inter].reset_index(drop=True)
M_f = M_log.reindex(columns=inter).fillna(0)
G_f = G_log.reindex(columns=inter).fillna(0)
print(f'特征空间: Xiao {X_f.shape} / Mollin {M_f.shape} / Gustin {G_f.shape}')

# ================= 3. 构建配对与标签 =================
# Xiao（y_labels 已含 outcome/base_cst4）
y_x = y_xiao['outcome'].values.astype(int)
cst4_x = y_xiao['base_cst4'].values.astype(bool)

# Mollin
M_v1 = M_f.loc[M_pairs['base_sample'].tolist()].reset_index(drop=True)
y_m = M_pairs['outcome'].values.astype(int)
cst4_m = M_pairs['base_cst4'].values.astype(bool)

# Gustin：v1 -> v3 配对
g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs) if g_pairs else ([], [])
G_base = G_f.loc[list(g_base)].reset_index(drop=True)
G_out = G_f.loc[list(g_out)].reset_index(drop=True)
# 结局：v3 LB>=0.5 康复
lb_g = 'Lactobacillus'
g_outcome = (G_out[lb_g] > np.log1p(0.5)).values.astype(int)
g_cst4 = (G_base[lb_g] < np.log1p(0.5)).values.astype(bool)
print(f'\nGustin: {len(g_pairs)} 对 v1->v3, 康复率 {g_outcome.mean():.3f}, 基线 CST4 {g_cst4.sum()}/{len(g_cst4)}')
print(f'Xiao: {len(y_x)} 对, 康复率 {y_x.mean():.3f}, 基线 CST4 {cst4_x.sum()}/{len(cst4_x)}')
print(f'Mollin: {len(y_m)} 对, 康复率 {y_m.mean():.3f}, 基线 CST4 {cst4_m.sum()}/{len(cst4_m)}')

# ================= 4. 模型与评估 =================
def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])
GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}

def fit_tuned(Xtr, ytr, seed):
    """内层 5 折 CV 选 C，返回训练好的 pipe"""
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    gs = GridSearchCV(make_pipe(), GRID, cv=skf, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_

def external_auc(Xtr, ytr, Xte, yte, seeds=SEEDS):
    """多种子：train -> predict external，返回 AUC 均值±SD"""
    aucs = []
    for s in range(seeds):
        mdl = fit_tuned(Xtr, ytr, seed=s)
        p = mdl.predict_proba(Xte)[:, 1]
        if len(np.unique(yte)) < 2:
            aucs.append(np.nan); continue
        aucs.append(roc_auc_score(yte, p))
    aucs = np.array(aucs); aucs = aucs[~np.isnan(aucs)]
    return aucs.mean(), aucs.std(), len(aucs)

def internal_cv(X, y, seeds=SEEDS):
    """队列内重复分层 5 折 CV，返回 (AUC 均值, SD, n)"""
    aucs = []
    for s in range(seeds):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, y):
            mdl = fit_tuned(X.iloc[tr], y[tr], seed=s)
            p = mdl.predict_proba(X.iloc[te])[:, 1]
            aucs.append(roc_auc_score(y[te], p))
    return np.mean(aucs), np.std(aucs), len(y)

rows = []
def rec(tag, mu, sd, n, note=''):
    rows.append(dict(analysis=tag, auc_mean=round(float(mu), 4), auc_sd=round(float(sd), 4),
                     n_test=int(n), note=note))
    print(f'  {tag}: AUC {mu:.3f} ± {sd:.3f} (n={n}) {note}')

# ================= 5. 运行矩阵 =================
print('\n===== A. 原假说（全人群）=====')
mu, sd, n = internal_cv(X_f, y_x); rec('A1 Xiao internal CV (full)', mu, sd, n, 'anchor 0.65')
mu, sd, n = external_auc(X_f, y_x, M_v1, y_m); rec('A2 Xiao->Mollin external (full)', mu, sd, n, 'anchor 0.609')
mu, sd, n = external_auc(X_f, y_x, G_base, g_outcome); rec('A3 Xiao->Gustin external (full) [NEW]', mu, sd, n, '原假说第3验证')
mu, sd, n = internal_cv(G_base, g_outcome); rec('A4 Gustin internal CV (full) [NEW]', mu, sd, n)
mu, sd, n = internal_cv(M_v1, y_m); rec('A5 Mollin internal CV (full)', mu, sd, n, 'n=33 小样本参考')

print('\n===== B. 新叙事（基线 CST4 亚组）=====')
mu, sd, n = internal_cv(X_f[cst4_x], y_x[cst4_x]); rec('B1 Xiao CST4 internal CV', mu, sd, n, 'anchor 0.596')
mu, sd, n = external_auc(X_f[cst4_x], y_x[cst4_x], M_v1[cst4_m], y_m[cst4_m]); rec('B2 Xiao->Mollin CST4 external', mu, sd, n, 'anchor 0.527')
mu, sd, n = external_auc(X_f[cst4_x], y_x[cst4_x], G_base[g_cst4], g_outcome[g_cst4]); rec('B3 Xiao->Gustin CST4 external [NEW]', mu, sd, n, '新叙事第3复现')
mu, sd, n = internal_cv(G_base[g_cst4], g_outcome[g_cst4]); rec('B4 Gustin CST4 internal [NEW]', mu, sd, n)
mu, sd, n = internal_cv(M_v1[cst4_m], y_m[cst4_m]); rec('B5 Mollin CST4 internal', mu, sd, n, 'n=30 小样本参考')

# ================= 6. Gustin α 多样性（论文 richness 发现对照） =================
print('\n===== C. Gustin α 多样性（对照论文 richness 发现）=====')
try:
    G_alpha = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha.csv')
    G_alpha['sample'] = [re.sub(r'_F\.fastq\.gz$', '', str(s)) for s in G_alpha['sample']]
    G_alpha_idx = dict(zip(G_alpha['sample'], G_alpha.index))
    ga = G_alpha.set_index('sample').loc[list(g_base)].reset_index(drop=True)
    cols_a = [c for c in ['shannon', 'chao1', 'n_asv'] if c in ga.columns]
    if len(cols_a) >= 2:
        mu, sd, n = internal_cv(ga[cols_a], g_outcome)
        rec('C1 Gustin alpha-only internal CV [NEW]', mu, sd, n, '论文: 基线 richness 高->失败')
    else:
        print('  缺 alpha 列:', cols_a)
except Exception as e:
    print('  alpha 读取失败:', e)

# ================= 7. 保存 =================
res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print(f'\n已保存: {OUT}')
print('\n矩阵总览:')
print(res.to_string(index=False))
