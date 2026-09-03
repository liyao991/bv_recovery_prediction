# -*- coding: utf-8 -*-
"""
c20_sensitivity.py —— #C20 任务二：4 项可选敏感性补充
=================================================================
派单：ASUS #C20（留言板 2026-09-03 18:45）
口径：全部沿用主矩阵协议（83 属交集 / 20 种子 / seed-mean AUC / 嵌套 CV 无泄漏）

任务② 模型族敏感性：elastic-net(l1_ratio=0.5, saga) + ridge(L2)，对照 L1 锚点
任务③ CLR 变换敏感性：伪计数 1e-4 + 中心对数比，替代 log1p，L1 主协议
任务④ LASSO 特征选择频率表：Xiao full/CST4 训练，20 种子非零系数频率
任务⑤ 非抽平 richness + 测序深度协变量（Gustin CST4）

输出：
  <OUT_DIR>/modeling/results_c20_model_family.csv   （任务②）
  <OUT_DIR>/modeling/results_c20_clr.csv            （任务③）
  <OUT_DIR>/modeling/results_c20_feat_freq.csv      （任务④）
  <OUT_DIR>/modeling/results_c20_richness_depth.csv （任务⑤）
"""
import pandas as pd
import numpy as np
import re
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

np.random.seed(42)
SEEDS = 20
PSEUDO = 1e-4   # CLR 伪计数（相对丰度尺度 0.01%）

# ================= 1. 数据加载（复刻 four_cohort_matrix.py） =================
X_all = pd.read_csv('<DATA_DIR>/feat/X_features.csv')
y_all = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
feat_meta = pd.read_csv('<DATA_DIR>/feat/feature_meta.csv')
genus_cols_x = feat_meta[feat_meta['source'] == 'genus']['feature'].tolist()
xiao_mask = (y_all['cohort'] == 'PRJNA302078').values
y_xiao = y_all[xiao_mask].reset_index(drop=True)
X_xiao_log1p = X_all[xiao_mask].reset_index(drop=True)

M_rel = pd.read_csv('<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
M_pairs = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')
G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [re.sub(r'_F\.fastq\.gz$', '', str(i)) for i in G_rel.index]
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')

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

def agg_rel(rel_df):
    df = rel_df.copy()
    df.columns = [norm_genus(c) for c in df.columns]
    return df.T.groupby(level=0).sum().T    # 线性聚合（不 log1p）

# 三队列线性相对丰度（属名归一聚合）
X_rel_lin = pd.DataFrame(np.expm1(X_xiao_log1p[genus_cols_x]).clip(lower=0),
                         columns=genus_cols_x, index=X_xiao_log1p.index)
X_lin = agg_rel(X_rel_lin)
M_lin = agg_rel(M_rel)
G_lin = agg_rel(G_rel)

inter = sorted(set(X_lin.columns) & set(M_lin.columns) & set(G_lin.columns))
print(f'三队列属名交集: {len(inter)}')

# log1p 特征（主协议）
X_f = np.log1p(X_lin[inter]).reset_index(drop=True)
M_f = np.log1p(M_lin.reindex(columns=inter).fillna(0))
G_f = np.log1p(G_lin.reindex(columns=inter).fillna(0))

# CLR 特征（任务③）：行内中心对数比，伪计数 PSEUDO
def clr_transform(lin_df):
    vals = lin_df.reindex(columns=inter).fillna(0).values.astype(float)
    lg = np.log(vals + PSEUDO)
    return pd.DataFrame(lg - lg.mean(axis=1, keepdims=True),
                        columns=inter, index=lin_df.reindex(columns=inter).fillna(0).index)

X_c = clr_transform(X_lin).reset_index(drop=True)
M_c = clr_transform(M_lin)
G_c = clr_transform(G_lin)

# ================= 2. 标签 =================
y_x = y_xiao['outcome'].values.astype(int)
cst4_x = y_xiao['base_cst4'].values.astype(bool)

M_v1 = M_f.loc[M_pairs['base_sample'].tolist()].reset_index(drop=True)
M_v1_c = M_c.loc[M_pairs['base_sample'].tolist()].reset_index(drop=True)
y_m = M_pairs['outcome'].values.astype(int)
cst4_m = M_pairs['base_cst4'].values.astype(bool)

g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)
G_base = G_f.loc[list(g_base)].reset_index(drop=True)
G_base_c = G_c.loc[list(g_base)].reset_index(drop=True)
G_out_df = G_f.loc[list(g_out)].reset_index(drop=True)
g_outcome = (G_out_df['Lactobacillus'] >= np.log1p(0.5)).values.astype(int)
g_cst4 = (G_base['Lactobacillus'] < np.log1p(0.5)).values.astype(bool)

print(f"Xiao n={len(y_x)} rec={y_x.sum()} | Mollin n={len(y_m)} rec={y_m.sum()} | Gustin n={len(g_outcome)} rec={g_outcome.sum()}")

# ================= 3. 模型协议 =================
GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}

def make_pipe(family):
    if family == 'L1':
        clf = LogisticRegression(solver='liblinear', penalty='l1',
                                 class_weight='balanced', max_iter=5000)
    elif family == 'enet':
        clf = LogisticRegression(solver='saga', penalty='elasticnet', l1_ratio=0.5,
                                 class_weight='balanced', max_iter=20000, tol=1e-4)
    elif family == 'ridge':
        clf = LogisticRegression(solver='liblinear', penalty='l2',
                                 class_weight='balanced', max_iter=5000)
    return Pipeline([('sc', StandardScaler()), ('clf', clf)])

def fit_tuned(Xtr, ytr, seed, family):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    gs = GridSearchCV(make_pipe(family), GRID, cv=skf, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_

def internal_cv(X, y, family, seeds=SEEDS):
    aucs = []
    for s in range(seeds):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, y):
            mdl = fit_tuned(X.iloc[tr], y[tr], seed=s, family=family)
            aucs.append(roc_auc_score(y[te], mdl.predict_proba(X.iloc[te])[:, 1]))
    return np.mean(aucs), np.std(aucs)

def external_auc(Xtr, ytr, Xte, yte, family, seeds=SEEDS):
    aucs = []
    for s in range(seeds):
        mdl = fit_tuned(Xtr, ytr, seed=s, family=family)
        aucs.append(roc_auc_score(yte, mdl.predict_proba(Xte)[:, 1]))
    return np.mean(aucs), np.std(aucs)

# 6 个代表 cell 的数据装配器
def get_cells(fs):
    """fs: 特征 DataFrame 集合（Xiao/Mollin/Gustin），返回 cell 列表"""
    Xx = fs['X']; Mm = fs['M']; Gg = fs['G']
    return [
        ('A1 Xiao internal CV (full)', 'internal', Xx, y_x, None, None),
        ('A2 Xiao->Mollin external',   'external', Xx, y_x, Mm, y_m),
        ('A3 Xiao->Gustin external',   'external', Xx, y_x, Gg, g_outcome),
        ('B1 Xiao CST4 internal CV',   'internal', Xx[cst4_x].reset_index(drop=True), y_x[cst4_x], None, None),
        ('B2 Xiao->Mollin CST4',       'external', Xx[cst4_x].reset_index(drop=True), y_x[cst4_x], Mm[cst4_m].reset_index(drop=True), y_m[cst4_m]),
        ('B3 Xiao->Gustin CST4',       'external', Xx[cst4_x].reset_index(drop=True), y_x[cst4_x], Gg[g_cst4].reset_index(drop=True), g_outcome[g_cst4]),
    ]

def run_cells(fs, family, tag):
    rows = []
    for name, kind, Xtr, ytr, Xte, yte in get_cells(fs):
        if kind == 'internal':
            mu, sd = internal_cv(Xtr, ytr, family)
        else:
            mu, sd = external_auc(Xtr, ytr, Xte, yte, family)
        rows.append(dict(cell=name, family=tag, auc=round(mu, 4), sd=round(sd, 4)))
        print(f'  [{tag}] {name:30s} AUC {mu:.3f} ± {sd:.3f}')
    return rows

# ================= 任务②：模型族敏感性 =================
print('\n===== 任务② 模型族敏感性（enet=elastic-net 0.5 / ridge=L2）=====')
feat_log = {'X': X_f, 'M': M_v1, 'G': G_base}
rows2 = run_cells(feat_log, 'enet', 'enet')
rows2 += run_cells(feat_log, 'ridge', 'ridge')
# L1 锚点重跑（同脚本内对照，防 RNG 状态漂移）
rows2 += run_cells(feat_log, 'L1', 'L1')
pd.DataFrame(rows2).to_csv('<OUT_DIR>/modeling/results_c20_model_family.csv', index=False)

# ================= 任务③：CLR 变换敏感性（L1 主协议） =================
print('\n===== 任务③ CLR 变换敏感性（伪计数 1e-4，L1 协议）=====')
feat_clr = {'X': X_c, 'M': M_v1_c, 'G': G_base_c}
rows3 = run_cells(feat_clr, 'L1', 'L1+CLR')
pd.DataFrame(rows3).to_csv('<OUT_DIR>/modeling/results_c20_clr.csv', index=False)

# ================= 任务④：LASSO 特征选择频率表 =================
print('\n===== 任务④ LASSO 特征选择频率（20 种子，L1）=====')
freq_rows = []
for scope, Xtr, ytr in [
    ('Xiao full (A系训练集)', X_f, y_x),
    ('Xiao CST4 (B系训练集)', X_f[cst4_x].reset_index(drop=True), y_x[cst4_x]),
]:
    counts = {}
    nz_list = []
    for s in range(SEEDS):
        mdl = fit_tuned(Xtr, ytr, seed=s, family='L1')
        coef = mdl.named_steps['clf'].coef_.ravel()
        nz = [Xtr.columns[i] for i in np.where(coef != 0)[0]]
        nz_list.append(len(nz))
        for g in nz:
            counts[g] = counts.get(g, 0) + 1
    print(f'  {scope}: 每种子非零属数 min={min(nz_list)} max={max(nz_list)} median={int(np.median(nz_list))}')
    for g, c in sorted(counts.items(), key=lambda kv: -kv[1])[:15]:
        freq_rows.append(dict(train_set=scope, genus=g, freq=f'{c}/{SEEDS}', pct=round(100*c/SEEDS, 1)))
        print(f'    {g:28s} {c}/{SEEDS}')
pd.DataFrame(freq_rows).to_csv('<OUT_DIR>/modeling/results_c20_feat_freq.csv', index=False)

# ================= 任务⑤：非抽平 richness + 测序深度协变量 =================
print('\n===== 任务⑤ Gustin richness + 测序深度（非抽平）=====')
af = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af['run'] = [re.sub(r'_F\.fastq\.gz$', '', str(s)) for s in af['sample']]
af = af.set_index('run')
base_runs = [r for r in g_base if r in af.index]
rich = af.loc[base_runs]['species_observed'].values.astype(float)

# 测序深度：ASV 表行和（基线样本）
counts = pd.read_csv('<OUT_DIR>/PRJNA691964/asv_table.csv', index_col=0)
counts.index = [re.sub(r'_F\.fastq\.gz$', '', str(i)) for i in counts.index]
depth = counts.loc[base_runs].values.astype(float).sum(axis=1)
log_depth = np.log(depth)

# 结局（对应同 33 例基线样本）
gm = {r: i for i, r in enumerate(g_base)}
sel = [gm[r] for r in base_runs]
y_g = g_outcome[np.array(sel)]

r5_rows = []
# ① 非抽平 richness 主关联（复现 ρ=−0.454）
rho, p = stats.spearmanr(rich, y_g)
r5_rows.append(dict(item='spearman_rich_outcome', value=round(rho, 4), note=f'p={p:.4f}; 非抽平 species_observed，主分析同源'))
print(f'  ① 非抽平 richness ρ={rho:.4f} (p={p:.4f})')

# ② depth 与 richness 的共线
rho_dr, p_dr = stats.spearmanr(log_depth, rich)
r5_rows.append(dict(item='spearman_depth_richness', value=round(rho_dr, 4), note=f'p={p_dr:.4f}'))
print(f'  ② Spearman(log depth, richness) ρ={rho_dr:.4f} (p={p_dr:.4f})')

# ③ depth 与结局的边际关联
rho_do, p_do = stats.spearmanr(log_depth, y_g)
r5_rows.append(dict(item='spearman_depth_outcome', value=round(rho_do, 4), note=f'p={p_do:.4f}'))
print(f'  ③ Spearman(log depth, outcome) ρ={rho_do:.4f} (p={p_do:.4f})')

# ④ partial Spearman(richness, outcome | log depth)
def partial_spearman(x, y, z):
    rx = stats.rankdata(x); ry = stats.rankdata(y); rz = stats.rankdata(z)
    # 秩残差
    def resid(a, b):
        A = np.column_stack([np.ones(len(b)), b])
        beta, *_ = np.linalg.lstsq(A, a, rcond=None)
        return a - A @ beta
    ex = resid(rx, rz); ey = resid(ry, rz)
    r, pp = stats.pearsonr(ex, ey)
    return r, pp
pr, ppr = partial_spearman(rich, y_g, log_depth)
r5_rows.append(dict(item='partial_spearman_rich_outcome_ctrl_logdepth', value=round(pr, 4), note=f'p={ppr:.4f}'))
print(f'  ④ partial ρ(rich, outcome | log depth) = {pr:.4f} (p={ppr:.4f})')

# ⑤ logistic: outcome ~ richness + log_depth（richness 系数）
import statsmodels.api as sm
Xlog = np.column_stack([np.ones(len(y_g)), rich, log_depth])
try:
    fit = sm.Logit(y_g, Xlog).fit(disp=0)
    b_rich = fit.params[1]; p_rich = fit.pvalues[1]
    b_dep = fit.params[2]; p_dep = fit.pvalues[2]
    r5_rows.append(dict(item='logit_rich_coef', value=round(b_rich, 4), note=f'p={p_rich:.4f}; 调整 log depth'))
    r5_rows.append(dict(item='logit_logdepth_coef', value=round(b_dep, 4), note=f'p={p_dep:.4f}; 调整 richness'))
    print(f'  ⑤ logistic richness β={b_rich:.4f} (p={p_rich:.4f}) | log depth β={b_dep:.4f} (p={p_dep:.4f})')
except Exception as e:
    print('  ⑤ logistic 拟合失败:', e)

# ⑥ 深度分层双稳态：上半/下半深度各算 ρ
med_d = np.median(depth)
for lab, m in [('depth>=median', depth >= med_d), ('depth<median', depth < med_d)]:
    r, pp = stats.spearmanr(rich[m], y_g[m])
    r5_rows.append(dict(item=f'spearman_rich_outcome_{lab.replace(">=", "_hi_").replace("<", "_lo_")}', value=round(r, 4), note=f'n={m.sum()}, p={pp:.4f}'))
    print(f'  ⑥ {lab}: ρ={r:.4f} (p={pp:.4f}, n={m.sum()})')

pd.DataFrame(r5_rows).to_csv('<OUT_DIR>/modeling/results_c20_richness_depth.csv', index=False)
print('\n全部完成，输出 4 个 CSV。')
