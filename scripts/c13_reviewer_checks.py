# -*- coding: utf-8 -*-
"""
c13_reviewer_checks.py —— #C13 独立统计审核的三项实测验证
==========================================================
站在"最挑剔审稿人 + 独立统计学家"视角，对框架三处软肋做定量验证：

Task A. 外推 AUC 的不确定性被 ±SD 低估
   矩阵协议 B2/B3 报告的 ±SD 只是 20 个训练种子的模型方差，
   不含测试集（n=30/33）的抽样变异。此处：20 种子平均概率 -> 测试集
   bootstrap 5000 次 -> percentile 95% CI。

Task B. Mollin 内部 AUC（B5=0.649 / A5=0.766）的诚实区间
   200 次分层 60/40 holdout（内层 3 折 CV 选 C）-> 均值 + [2.5%, 97.5%]。
   用于判断"三队列复现≈随机"叙事与 Mollin 点估计的张力有多大。

Task C. Gustin richness 信号 vs 基线 LB 耦合（论文自己的批评反身适用）
   论文批评别人"全人群信号=LB 标签耦合"，但自己的唯一阳性信号
   （Gustin richness rho=-0.454）未做控制基线 LB 的偏相关检验。
   此处：Spearman + 秩残差偏 Spearman + 二变量 logistic。

口径与 four_cohort_matrix.py 完全一致（83 属交集、log1p、LASSO-LR、
内层 CV 选 C、class_weight=balanced）。
"""
import pandas as pd
import numpy as np
import re
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV, train_test_split
from sklearn.metrics import roc_auc_score

RNG = 42
np.random.seed(RNG)
OUT_CSV = '<REPO_DIR>/results/results_c13_reviewer_checks.csv'
rows = []

def rec(task, item, value, extra=''):
    rows.append(dict(task=task, item=item, value=value, extra=extra))
    print(f'  [{task}] {item}: {value}  {extra}')

# ================= 数据加载（复刻 four_cohort_matrix 口径） =================
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

def agg_rel_to_log1p(rel_df):
    df = rel_df.copy()
    df.columns = [norm_genus(c) for c in df.columns]
    agg = df.T.groupby(level=0).sum().T
    return np.log1p(agg)

X_rel = np.expm1(X_xiao_log1p[genus_cols_x]).clip(lower=0)
X_rel = pd.DataFrame(X_rel, columns=genus_cols_x, index=X_xiao_log1p.index)
X_log = agg_rel_to_log1p(X_rel)
M_log = agg_rel_to_log1p(M_rel)
G_log = agg_rel_to_log1p(G_rel)

inter = sorted(set(X_log.columns) & set(M_log.columns) & set(G_log.columns))
X_f = X_log[inter].reset_index(drop=True)
M_f = M_log.reindex(columns=inter).fillna(0)
G_f = G_log.reindex(columns=inter).fillna(0)

y_x = y_xiao['outcome'].values.astype(int)
cst4_x = y_xiao['base_cst4'].values.astype(bool)
M_v1 = M_f.loc[M_pairs['base_sample'].tolist()].reset_index(drop=True)
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
G_out = G_f.loc[list(g_out)].reset_index(drop=True)
g_outcome = (G_out['Lactobacillus'] > np.log1p(0.5)).values.astype(int)
g_cst4 = (G_base['Lactobacillus'] < np.log1p(0.5)).values.astype(bool)
print(f'数据: Xiao {len(y_x)}(CST4 {cst4_x.sum()}) Mollin {len(y_m)}(CST4 {cst4_m.sum()}) '
      f'Gustin {len(g_outcome)}(CST4 {g_cst4.sum()}) 交集 {len(inter)} 属')

# ================= 模型 =================
def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])
GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}

def fit_tuned(Xtr, ytr, seed, inner=5):
    gs = GridSearchCV(make_pipe(), GRID,
                      cv=StratifiedKFold(inner, shuffle=True, random_state=seed),
                      scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_

# ================= Task A: 外推 AUC 的测试集 bootstrap CI =================
print('\n===== Task A: 外推 AUC 测试集 bootstrap CI（矩阵协议 ±SD 之外的真实不确定性）=====')
def external_boot(Xtr, ytr, Xte, yte, tag, seeds=20, nboot=5000):
    probs = np.zeros(len(yte))
    for s in range(seeds):
        mdl = fit_tuned(Xtr, ytr, seed=s)
        probs += mdl.predict_proba(Xte)[:, 1]
    probs /= seeds
    auc_pt = roc_auc_score(yte, probs)
    rng = np.random.RandomState(123)
    boots = []
    n = len(yte)
    for _ in range(nboot):
        idx = rng.randint(0, n, n)
        if len(np.unique(yte[idx])) < 2:
            continue
        boots.append(roc_auc_score(yte[idx], probs[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rec('A', tag, f'{auc_pt:.3f}', f'测试集 bootstrap 95%CI [{lo:.3f}, {hi:.3f}] (n_test={n})')
    return auc_pt, lo, hi

# B3: Xiao CST4 -> Gustin CST4（=全部 33 例）
external_boot(X_f[cst4_x], y_x[cst4_x], G_base[g_cst4].reset_index(drop=True),
              g_outcome[g_cst4], 'B3 Xiao->Gustin CST4 (报告值 0.512)')
# B2: Xiao CST4 -> Mollin CST4
external_boot(X_f[cst4_x], y_x[cst4_x], M_v1[cst4_m].reset_index(drop=True),
              y_m[cst4_m], 'B2 Xiao->Mollin CST4 (报告值 0.505±0.020)')
# A3: Xiao full -> Gustin full
external_boot(X_f, y_x, G_base, g_outcome, 'A3 Xiao->Gustin full (报告值 0.556±0.021)')

# ================= Task B: Mollin 内部 holdout 分布 =================
print('\n===== Task B: Mollin 内部 200 次分层 60/40 holdout =====')
def holdout_dist(X, y, tag, n_rep=200):
    aucs = []
    for r in range(n_rep):
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.4, stratify=y,
                                              random_state=1000 + r)
        if len(np.unique(yte)) < 2:
            continue
        mdl = fit_tuned(Xtr, ytr, seed=r, inner=3)
        aucs.append(roc_auc_score(yte, mdl.predict_proba(Xte)[:, 1]))
    aucs = np.array(aucs)
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    rec('B', tag, f'{aucs.mean():.3f}',
        f'holdout 分布 [2.5%,97.5%]=[{lo:.3f},{hi:.3f}] P(>0.5)={np.mean(aucs>0.5):.2f} (n={len(y)})')

holdout_dist(M_v1, y_m, 'Mollin full internal (CV 报告 0.766)')
holdout_dist(M_v1[cst4_m].reset_index(drop=True), y_m[cst4_m],
             'Mollin CST4 internal (CV 报告 0.649)')

# ================= Task C: Gustin richness 偏相关（控制基线 LB） =================
print('\n===== Task C: Gustin richness 信号控制基线 LB（反身应用论文自己的耦合批评）=====')
G_sp = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_species_rel.csv', index_col=0)
G_sp.index = [re.sub(r'_F\.fastq\.gz$', '', str(i)) for i in G_sp.index]
G_sp_base = G_sp.loc[list(g_base)]
richness = (G_sp_base > 0).sum(axis=1).values.astype(float)
baseLB = G_rel.loc[list(g_base)]['Lactobacillus'].values.astype(float)
y_bin = g_outcome.astype(float)

rho_ro, p_ro = stats.spearmanr(richness, y_bin)
rho_lo, p_lo = stats.spearmanr(baseLB, y_bin)
rho_rl, p_rl = stats.spearmanr(richness, baseLB)
rec('C', 'rho(richness, outcome)', f'{rho_ro:.3f} (p={p_ro:.4f})', '应复现 -0.454')
rec('C', 'rho(baseLB, outcome)', f'{rho_lo:.3f} (p={p_lo:.4f})', '基线 LB 单特征')
rec('C', 'rho(richness, baseLB)', f'{rho_rl:.3f} (p={p_rl:.4f})', 'C8 已示 CST4 内两者弱相关')

# 秩残差偏 Spearman：richness~outcome | baseLB
def rank_resid(v, covar):
    rv = stats.rankdata(v); rc = stats.rankdata(covar)
    X_d = np.column_stack([np.ones(len(rc)), rc])
    beta = np.linalg.lstsq(X_d, rv, rcond=None)[0]
    return rv - X_d @ beta

res_r = rank_resid(richness, baseLB)
res_y = rank_resid(y_bin, baseLB)
r_part = np.corrcoef(res_r, res_y)[0, 1]
n = len(y_bin); k = 1
t_stat = r_part * np.sqrt((n - 2 - k) / (1 - r_part**2))
p_part = 2 * stats.t.sf(abs(t_stat), df=n - 2 - k)
rec('C', 'partial rho(richness, outcome | baseLB)', f'{r_part:.3f} (p={p_part:.4f})',
    '控制基线 LB 后是否仍显著')

# 二变量 logistic（标准化系数方向对照）
Xs = StandardScaler().fit_transform(np.column_stack([richness, baseLB]))
lr = LogisticRegression(max_iter=5000).fit(Xs, y_bin)
rec('C', 'logistic coef (richness | baseLB 同模型)', f'{lr.coef_[0][0]:.3f}',
    f'baseLB coef={lr.coef_[0][1]:.3f}')

pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f'\n已保存: {OUT_CSV}')
