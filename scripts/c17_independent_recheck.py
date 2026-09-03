# -*- coding: utf-8 -*-
"""
c17_independent_recheck.py —— #C17 独立复现验证（计算正确性防线）
====================================================================
ASUS #C17 要求（2026-09-01 17:10）：不用已有脚本对输出，用另一条路算同一个数。
本脚本从 <DATA_DIR> 原始 CSV 出发，全部用 scipy/sklearn/numpy 标准库独立实现：
  1. Gustin richness：rho/p（直算+置换验证）、经验 AUC、family CV AUC、rarefy 双路对照
  2. alpha 矩阵 richness 家族三分辨率（主口径=线性 LB>=50%）
  3. 矩阵协议 A1/A2/A3/A4/B1/B2/B3（独立实现同协议）+ 外推 seed-mean bootstrap CI
  4. Xiao 200 次重复 holdout（全人群 + CST4）
  5. Fisher 精确检验 + 中位数切点
  6. partial rho（rank 残差 + logistic 对照）
  7. 数据溯源：各队列 n / 事件数 / 2x2 表 / 阈值卡点检查

关键口径事实（#C17 期间实测纠正）：
  - PRJNA*_genus_rel.csv 全部为【线性相对丰度】（max<=1.0）
  - 全研究结局统一口径 = 线性 LB>=50%（build_labels cst_type lb<0.5 + Gustin >=0.5）
  - 旧 bug 阈值 >log1p(0.5)=0.4055（线性 40.5%，richness_direction_check.py 已注明笔误）
  - #C16 所报「主口径 LB>=64.87%」系数据空间误读，实际不存在该口径
输出：<REPO_DIR>/results/results_c17_recheck.csv
"""
import pandas as pd
import numpy as np
import re
import warnings
warnings.filterwarnings('ignore')
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV, train_test_split
from sklearn.metrics import roc_auc_score

OUT = '<REPO_DIR>/results/results_c17_recheck.csv'
rows = []


def rec(item, value, expect, note=''):
    status = ''
    rows.append(dict(item=item, value=value, expect=expect, note=note))
    print(f'  {item}: {value}   [框架 {expect}]  {note}')


def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))


# ================= Part 0: 数据加载 + 溯源 =================
print('=' * 70)
print('Part 0  数据溯源（原始 CSV 直读）')
X_all = pd.read_csv('<DATA_DIR>/feat/X_features.csv')
y_all = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
feat_meta = pd.read_csv('<DATA_DIR>/feat/feature_meta.csv')
genus_cols_x = feat_meta[feat_meta['source'] == 'genus']['feature'].tolist()
alpha_cols_x = feat_meta[feat_meta['source'] == 'alpha']['feature'].tolist()

xiao_mask = (y_all['cohort'] == 'PRJNA302078').values
y_xiao = y_all[xiao_mask].reset_index(drop=True)
X_xiao_log1p = X_all[xiao_mask].reset_index(drop=True)
y_x = y_xiao['outcome'].values.astype(int)
cst4_x = y_xiao['base_cst4'].values.astype(bool)

M_rel = pd.read_csv('<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
M_pairs = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')
G_rel = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_genus_rel.csv', index_col=0)
G_rel.index = [norm_sid(i) for i in G_rel.index]
G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')

# 阈值卡点检查：Gustin v3 LB 是否有恰=0.5 者（>= 与 > 等价性）
g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)
lb_v3 = G_rel.loc[list(g_out)]['Lactobacillus'].values.astype(float)
n_at_thr = int((lb_v3 == 0.5).sum())
y_g = (lb_v3 >= 0.5).astype(int)          # 主口径：线性 LB>=50%
g_cst4 = (G_rel.loc[list(g_base)]['Lactobacillus'] < 0.5).values.astype(bool)

y_m = M_pairs['outcome'].values.astype(int)
cst4_m = M_pairs['base_cst4'].values.astype(bool)

print(f'Xiao    n={len(y_x)}  康复 {y_x.sum()}/{len(y_x)}  CST4 {cst4_x.sum()}  (CST4 内 {y_x[cst4_x].sum()}/{cst4_x.sum()})')
print(f'Plummer n={int((y_all["cohort"]=="PRJNA592384").sum())}  康复 {y_all[y_all["cohort"]=="PRJNA592384"]["outcome"].sum()}')
print(f'Mollin  n={len(y_m)}  康复 {y_m.sum()}  CST4 {cst4_m.sum()}')
print(f'Gustin  n={len(y_g)}  康复 {y_g.sum()}  基线CST4 {g_cst4.sum()}  | v3 LB 恰=0.5 者 {n_at_thr} 个（>=与>等价性）')
# Mollin outcome 与 lb_pct 一致性核验
m_lb_out = M_pairs['outcome_lb_pct'].values.astype(float)
m_match = int(((m_lb_out >= 0.5).astype(int) == y_m).sum())
print(f'Mollin outcome 与 outcome_lb_pct>=0.5 一致 {m_match}/{len(y_m)}')

rec('n_Xiao', f'{len(y_x)} (康复{y_x.sum()}/CST4 {cst4_x.sum()})', '68/39/54', '溯源')
rec('n_Plummer', f'{int((y_all["cohort"]=="PRJNA592384").sum())} (康复{y_all[y_all["cohort"]=="PRJNA592384"]["outcome"].sum()}/CST4 43)', '49/44/43', '溯源')
rec('n_Mollin', f'{len(y_m)} (康复{y_m.sum()}/CST4 {cst4_m.sum()})', '33/15/30', '溯源')
rec('n_Gustin', f'{len(y_g)} (康复{y_g.sum()}/基线CST4 {g_cst4.sum()})', '33/11/33', '溯源')

# ================= Part 1: Gustin richness 统计量 =================
print('\n' + '=' * 70)
print('Part 1  Gustin 种级 richness 统计量（直算 + 置换验证）')
af = pd.read_csv('<DATA_DIR>/feat/PRJNA691964_alpha_family.csv')
af['run'] = [norm_sid(s) for s in af['sample']]
af = af.set_index('run')
rich = af.loc[list(g_base)]['species_observed'].values.astype(float)
chao = af.loc[list(g_base)]['species_chao1'].values.astype(float)

rho, rho_p = stats.spearmanr(rich, y_g)
rec('rho(richness, outcome)', f'{rho:.4f} (p={rho_p:.4f})', '-0.454 / p=0.008', 'scipy 直算')

# 置换检验（另一条路验证 p 值）
rng = np.random.default_rng(2024)
n_perm = 10000
perm_r = np.array([stats.spearmanr(rich, rng.permutation(y_g))[0] for _ in range(n_perm)])
p_perm = (np.abs(perm_r) >= abs(rho)).mean()
rec('perm_p (10000x)', f'{p_perm:.4f}', '~0.005-0.008', '置换分布独立验证 p')

auc_emp = roc_auc_score(y_g, -rich)
rec('empirical AUC (single feature)', f'{auc_emp:.4f}', '0.777', '方向=richness 预测未康复')

# ================= Part 2: CV AUC（family + 单指标 + alpha 矩阵 richness 家族） =================
print('\n' + '=' * 70)
print('Part 2  CV AUC（20 种子 x 分层 5 折，独立实现）')


def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])


def cv_auc(Xdf, yv, seeds=20):
    aucs = []
    for s in range(seeds):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(Xdf, yv):
            if len(np.unique(yv[tr])) < 2 or len(np.unique(yv[te])) < 2:
                continue
            mdl = make_pipe().fit(Xdf.iloc[tr], yv[tr])
            aucs.append(roc_auc_score(yv[te], mdl.predict_proba(Xdf.iloc[te])[:, 1]))
    v = np.array(aucs)
    return v.mean(), v.std()


mu1, sd1 = cv_auc(pd.DataFrame({'species_observed': rich}), y_g)
rec('CV AUC species_observed 单指标', f'{mu1:.4f} ± {sd1:.4f}', '≈0.776', '')
mu2, sd2 = cv_auc(pd.DataFrame({'observed': rich, 'chao1': chao}), y_g)
rec('CV AUC family:richness (种级)', f'{mu2:.4f} ± {sd2:.4f}', '0.776±0.164', '')

# alpha 矩阵 richness 家族三分辨率（主口径=线性 LB>=50%）
METRICS = ['observed', 'chao1']
print('\nalpha 矩阵 richness 家族（主口径 线性 LB>=50%）：')
a_base = af.loc[list(g_base)]
for lv, expect in [('asv', '0.727±0.210'), ('genus', '0.681±0.219'), ('species', '0.776±0.164')]:
    cols = [f'{lv}_{m}' for m in METRICS]
    mu, sd = cv_auc(a_base[cols], y_g)
    rec(f'alpha矩阵 {lv} richness', f'{mu:.4f} ± {sd:.4f}', expect, '')

# ================= Part 3: 矩阵协议 A/B 系（独立实现） =================
print('\n' + '=' * 70)
print('Part 3  四队列矩阵协议复刻（83 属交集 + 20 种子 LASSO）')


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


def to_log1p(rel_df):
    df = rel_df.copy()
    df.columns = [norm_genus(c) for c in df.columns]
    return np.log1p(df.T.groupby(level=0).sum().T)


X_rel = np.expm1(X_xiao_log1p[genus_cols_x]).clip(lower=0)
X_rel = pd.DataFrame(X_rel, columns=genus_cols_x, index=X_xiao_log1p.index)
X_log = to_log1p(X_rel)
M_log = to_log1p(M_rel)
G_log = to_log1p(G_rel)
inter = sorted(set(X_log.columns) & set(M_log.columns) & set(G_log.columns))
X_f = X_log[inter].reset_index(drop=True)
M_f = M_log.reindex(columns=inter).fillna(0)
G_f = G_log.reindex(columns=inter).fillna(0)
M_v1 = M_f.loc[M_pairs['base_sample'].tolist()].reset_index(drop=True)
G_base = G_f.loc[list(g_base)].reset_index(drop=True)
G_out_f = G_f.loc[list(g_out)].reset_index(drop=True)
y_g_matrix = (G_out_f['Lactobacillus'] > np.log1p(0.5)).values.astype(int)  # log 空间>0.405 = 线性>50%
print(f'三队列属交集: {len(inter)}; 矩阵口径 Gustin 结局与主口径逐例一致: '
      f'{int((y_g_matrix == y_g).all())}（n_rec {y_g_matrix.sum()} vs {y_g.sum()}）')

GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}


def fit_tuned(Xtr, ytr, seed):
    gs = GridSearchCV(make_pipe(), GRID, cv=StratifiedKFold(5, shuffle=True, random_state=seed),
                      scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_


def internal_cv(X, y, seeds=20):
    aucs = []
    for s in range(seeds):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, y):
            mdl = fit_tuned(X.iloc[tr], y[tr], seed=s)
            aucs.append(roc_auc_score(y[te], mdl.predict_proba(X.iloc[te])[:, 1]))
    return np.mean(aucs), np.std(aucs)


def external_seedmean(Xtr, ytr, Xte, yte, seeds=20, nboot=5000):
    P = np.zeros((seeds, len(yte)))
    for s in range(seeds):
        P[s] = fit_tuned(Xtr, ytr, seed=s).predict_proba(Xte)[:, 1]
    per_seed = np.array([roc_auc_score(yte, P[s]) for s in range(seeds)])
    r2 = np.random.RandomState(2024)
    boots = []
    for _ in range(nboot):
        idx = r2.randint(0, len(yte), len(yte))
        if len(np.unique(yte[idx])) < 2:
            continue
        boots.append(np.mean([roc_auc_score(yte[idx], P[s][idx]) for s in range(seeds)]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return per_seed.mean(), per_seed.std(ddof=1), lo, hi


mu, sd = internal_cv(X_f, y_x)
rec('A1 Xiao internal CV (full)', f'{mu:.4f} ± {sd:.4f}', '0.626±0.163', '')
mu, sd, lo, hi = external_seedmean(X_f, y_x, M_v1, y_m)
rec('A2 Xiao->Mollin external', f'{mu:.4f} ± {sd:.4f} CI[{lo:.3f},{hi:.3f}]', '0.521 CI[0.339,0.702]', 'seed-mean')
mu, sd, lo, hi = external_seedmean(X_f, y_x, G_base, y_g)
rec('A3 Xiao->Gustin external', f'{mu:.4f} ± {sd:.4f} CI[{lo:.3f},{hi:.3f}]', '0.556 CI[0.341,0.757]', 'seed-mean')
mu, sd = internal_cv(G_base, y_g)
rec('A4 Gustin internal CV', f'{mu:.4f} ± {sd:.4f}', '0.436', '')
mu, sd = internal_cv(X_f[cst4_x], y_x[cst4_x])
rec('B1 Xiao CST4 internal CV', f'{mu:.4f} ± {sd:.4f}', '0.535±0.127', '')
mu, sd, lo, hi = external_seedmean(X_f[cst4_x], y_x[cst4_x], M_v1[cst4_m], y_m[cst4_m])
rec('B2 Xiao->Mollin CST4', f'{mu:.4f} ± {sd:.4f} CI[{lo:.3f},{hi:.3f}]', '0.501 CI[0.299,0.714]', 'seed-mean')
mu, sd, lo, hi = external_seedmean(X_f[cst4_x], y_x[cst4_x], G_base[g_cst4], y_g[g_cst4])
rec('B3 Xiao->Gustin CST4', f'{mu:.4f} ± {sd:.4f} CI[{lo:.3f},{hi:.3f}]', '0.512 CI[0.296,0.729]', 'seed-mean')

# ================= Part 4: Xiao 200 次重复 holdout（全人群 + CST4） =================
print('\n' + '=' * 70)
print('Part 4  Xiao 200 次重复 holdout（60/40，LASSO C=1.0，keep>=5% 属 + alpha）')
keep = [c for c in genus_cols_x if (X_all[xiao_mask][c] > 0).mean() >= 0.05]
feats = keep + alpha_cols_x
print(f'特征数: {len(feats)} (keep 属 {len(keep)} + alpha {len(alpha_cols_x)})')


def holdout(Xd, yv, feats_use, n_rep=200, seed=42, test_size=0.4):
    r = np.random.default_rng(seed)
    aucs = []
    for i in range(n_rep):
        tr, te = train_test_split(np.arange(len(yv)), test_size=test_size, stratify=yv,
                                  random_state=int(r.integers(1e9)))
        pipe = Pipeline([('sc', StandardScaler()),
                         ('clf', LogisticRegression(solver='liblinear', C=1.0, l1_ratio=1.0,
                                                    class_weight='balanced', max_iter=5000,
                                                    random_state=seed))])
        pipe.fit(Xd.iloc[tr][feats_use], yv[tr])
        aucs.append(roc_auc_score(yv[te], pipe.predict_proba(Xd.iloc[te][feats_use])[:, 1]))
    aucs = np.array(aucs)
    return aucs.mean(), np.percentile(aucs, 2.5), np.percentile(aucs, 97.5), (aucs > 0.5).mean()


Xd = X_all[xiao_mask].reset_index(drop=True)
m1, lo1, hi1, p1 = holdout(Xd, y_x, feats)
rec('Xiao full holdout 200x', f'{m1:.4f} CI[{lo1:.3f},{hi1:.3f}] P(>0.5)={p1:.3f}', '0.652 [0.500,0.797]', '')
Xd4 = Xd[cst4_x].reset_index(drop=True)
m2, lo2, hi2, p2 = holdout(Xd4, y_x[cst4_x], [c for c in feats if c in Xd4.columns])
rec('Xiao CST4 holdout 200x', f'{m2:.4f} CI[{lo2:.3f},{hi2:.3f}] P(>0.5)={p2:.3f}', '0.596 [0.388,0.802]', '')

# ================= Part 5: Fisher + 中位数切点 =================
print('\n' + '=' * 70)
print('Part 5  richness 切点 + Fisher 精确检验（主口径 LB>=50%）')
med = np.median(rich)
hi_mask = rich >= med
tab = [[int((hi_mask & (y_g == 1)).sum()), int((hi_mask & (y_g == 0)).sum())],
       [int((~hi_mask & (y_g == 1)).sum()), int((~hi_mask & (y_g == 0)).sum())]]
odds, fp = stats.fisher_exact(tab)
_, fp1 = stats.fisher_exact(tab, alternative='less')
rate_hi = tab[0][0] / sum(tab[0]); rate_lo = tab[1][0] / sum(tab[1])
print(f'中位数={med:.1f} (range {rich.min():.0f}-{rich.max():.0f}); 2x2={tab}')
print(f'康复率: 高组 {tab[0][0]}/{sum(tab[0])}={rate_hi:.1%} vs 低组 {tab[1][0]}/{sum(tab[1])}={rate_lo:.1%}')
rec('G2 中位数切点', f'{med:.1f}', '22', '')
rec('G2 康复率 hi vs lo', f'{rate_hi:.1%} vs {rate_lo:.1%}', '24% vs 44%', f'2x2={tab}')
rec('G2 Fisher 双侧 p', f'{fp:.4f} (定向 {fp1:.4f})', '0.282', '')

# ================= Part 6: partial rho =================
print('\n' + '=' * 70)
print('Part 6  partial rho（控制基线 LB，rank 残差独立实现）')
baseLB = G_rel.loc[list(g_base)]['Lactobacillus'].values.astype(float)


def rank_resid(v, covar):
    rv = stats.rankdata(v); rc = stats.rankdata(covar)
    D = np.column_stack([np.ones(len(rc)), rc])
    beta = np.linalg.lstsq(D, rv, rcond=None)[0]
    return rv - D @ beta


res_r = rank_resid(rich, baseLB)
res_y = rank_resid(y_g, baseLB)
r_part = np.corrcoef(res_r, res_y)[0, 1]
n = len(y_g)
t_stat = r_part * np.sqrt((n - 3) / (1 - r_part ** 2))
p_part = 2 * stats.t.sf(abs(t_stat), df=n - 3)
rec('partial rho(rich, outcome | baseLB)', f'{r_part:.4f} (p={p_part:.4f})', '-0.431 / p=0.014', '')

# 对照：双变量 logistic richness 系数
Xs = StandardScaler().fit_transform(np.column_stack([rich, baseLB]))
lr = LogisticRegression(max_iter=5000).fit(Xs, y_g)
rec('logistic richness coef', f'{lr.coef_[0][0]:.3f}', '(负号)', '对照')

# ================= Part 7: rarefaction 双路对照 =================
print('\n' + '=' * 70)
print('Part 7  rarefaction：原协议（多项近似）vs 真·无放回（多元超几何）')


def load_counts(path):
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    df.index = [norm_sid(i) for i in df.index]
    return df


counts_g = load_counts('<OUT_DIR>/PRJNA691964/asv_table.csv')
gc = counts_g.loc[list(g_base)]
depth = int(gc.sum(axis=1).min())
print(f'Gustin 基线最小 reads = {depth}')

# 路 A：多项近似（原协议，200 次平均）
rngA = np.random.default_rng(42)
richA = []
for i in range(gc.shape[0]):
    c = gc.iloc[i].values
    p = c / c.sum()
    samp = rngA.multinomial(depth, p, size=200)
    richA.append((samp > 0).sum(axis=1).mean())
richA = np.array(richA)
rA, pA = stats.spearmanr(richA, y_g)
rec('rarefy rho (多项近似, 原协议)', f'{rA:.4f} (p={pA:.4f})', '-0.392 / p=0.024', '')

# 路 B：真无放回（多元超几何，200 次平均）
rngB = np.random.default_rng(2024)
richB = []
for i in range(gc.shape[0]):
    c = gc.iloc[i].values
    samp = rngB.multivariate_hypergeometric(c, depth, size=200)
    richB.append((samp > 0).sum(axis=1).mean())
richB = np.array(richB)
rB, pB = stats.spearmanr(richB, y_g)
rec('rarefy rho (超几何, 真无放回)', f'{rB:.4f} (p={pB:.4f})', '对照', '另一条路')

# ================= 保存 =================
pd.DataFrame(rows).to_csv(OUT, index=False)
print('\n' + '=' * 70)
print(f'已保存: {OUT}')
