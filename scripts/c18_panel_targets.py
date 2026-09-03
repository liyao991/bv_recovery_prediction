# -*- coding: utf-8 -*-
"""
c18_panel_targets.py —— #C18 临床关键菌（检测面板靶标）靶向预测：补全 + 复核
======================================================================
派单 3 项：
  任务1 补全单变量扫查：7 属 + 2 可注释种，三队列基线 CST4 亚组，
        结局=随访线性 LB>=50%（康复=1）；输出 AUC / Spearman rho / p / BH-FDR
  任务2 组合检验：单属层扩到 7、组合层 7 属 / 2 种 / 9 特征 panel，
        对照全集 65 属；LASSO-LR 20 种子 x 5 折内部 CV
  任务3 Plummer AUC=0.5/sd=0 异常诊断：逐折记录预测是否常数、选中的 C、
        补报 PR-AUC / balanced accuracy；附固定 C=1.0 敏感性

口径铁律（与主分析一致）：线性 LB>=50%、CST4=base_LB<0.5、Gustin 全队列即 CST4。
数据空间：genus_rel / species_rel 均为线性相对丰度（#C17 确认）。
Prevotella 聚合修正：先线性求和再 log1p（原 c7 在 log1p 空间求和不严谨）。

输出：<OUT_DIR>/modeling/results_c18_univariate.csv
      <OUT_DIR>/modeling/results_c18_combo.csv
      <REPO_DIR>/results/results_c18_plummer_diag.csv
"""
import pandas as pd
import numpy as np
import re
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from statsmodels.stats.multitest import multipletests

OUT_U = '<OUT_DIR>/modeling/results_c18_univariate.csv'
OUT_C = '<OUT_DIR>/modeling/results_c18_combo.csv'
OUT_D = '<REPO_DIR>/results/results_c18_plummer_diag.csv'
SEEDS = 20
np.random.seed(42)

PROJ = {'Xiao': 'PRJNA302078', 'Plummer': 'PRJNA592384', 'Gustin': 'PRJNA691964'}


def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))


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


# ================= 1. 加载 =================
GE_LIN, SP_LIN = {}, {}
G_LOG = {}
for c, pid in PROJ.items():
    g = pd.read_csv(f'<DATA_DIR>/feat/{pid}_genus_rel.csv', index_col=0)
    g.index = [norm_sid(i) for i in g.index]
    GE_LIN[c] = g
    sp = pd.read_csv(f'<DATA_DIR>/feat/{pid}_species_rel.csv', index_col=0)
    sp.index = [norm_sid(i) for i in sp.index]
    SP_LIN[c] = sp
    # log1p + norm_genus 聚合（供全集 65 属对照，与 c7 协议一致）
    gc = g.copy()
    gc.columns = [norm_genus(col) for col in gc.columns]
    G_LOG[c] = np.log1p(gc.T.groupby(level=0).sum().T)

inter = sorted(set(G_LOG['Xiao'].columns) & set(G_LOG['Plummer'].columns) & set(G_LOG['Gustin'].columns))
print(f'三队列属名交集: {len(inter)} 属')

G_map = pd.read_csv('<DATA_DIR>/meta/gustin_run_map.csv')
G_map = G_map[G_map['control'].isna()]

# ================= 2. subject 表（CST4 亚组） =================
y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
subj = []

xiao = y[y['cohort'] == 'PRJNA302078'].copy()
xiao['base'] = [norm_sid(s) for s in xiao['base_sample']]
for _, r in xiao.iterrows():
    subj.append(dict(cohort='Xiao', base_run=r['base'], cst4=bool(r['base_cst4']),
                     outcome=int(r['outcome'])))

plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
for _, r in plum.iterrows():
    subj.append(dict(cohort='Plummer', base_run=r['base'], cst4=bool(r['base_cst4']),
                     outcome=int(r['outcome'])))

pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
gb, go = zip(*pairs)
gg = GE_LIN['Gustin']
g_outcome = (gg.loc[list(go), 'Lactobacillus'] >= 0.5).values.astype(int)  # 线性 >=50%
g_cst4 = (gg.loc[list(gb), 'Lactobacillus'] < 0.5).values                  # 线性 <50%
for b, o, oc, c4 in zip(gb, go, g_outcome, g_cst4):
    subj.append(dict(cohort='Gustin', base_run=b, cst4=bool(c4), outcome=int(oc)))

S = pd.DataFrame(subj)
print('\n=== subject 表（CST4 亚组）===')
for c in ['Xiao', 'Plummer', 'Gustin']:
    sub = S[(S['cohort'] == c) & S['cst4']]
    print(f'  {c}: CST4 n={len(sub)}, rec={int(sub["outcome"].sum())}, rate={sub["outcome"].mean():.2f}')

# ================= 3. 靶标定义 =================
GENUS7 = ['Lactobacillus', 'Gardnerella', 'Atopobium', 'Megasphaera',
          'Sneathia', 'Mobiluncus']  # + Prevotella 聚合 = 7
SPECIES2 = {'G.vaginalis': 'Gardnerella Gardnerella vaginalis',
            'A.vaginae': 'Atopobium Atopobium vaginae DSM 15829'}
MISSING_SPECIES = ['L.crispatus', 'L.iners', 'P.bivia', 'BVAB2', 'Megasphaera-1']


def genus_features(cohort, base_runs):
    """7 属特征（log1p 空间；Prevotella 先线性求和再 log1p）"""
    ge = GE_LIN[cohort]
    feat = pd.DataFrame(index=base_runs)
    for gname in GENUS7:
        if gname in ge.columns:
            feat[gname] = np.log1p(ge.loc[base_runs, gname].values.astype(float))
        else:
            feat[gname] = 0.0
    pv = [c for c in ge.columns if 'Prevotella' in c]
    feat['Prevotella_agg'] = np.log1p(ge.loc[base_runs, pv].sum(axis=1).values.astype(float)) if pv else 0.0
    return feat


def species_features(cohort, base_runs):
    sp = SP_LIN[cohort]
    feat = pd.DataFrame(index=base_runs)
    for tag, col in SPECIES2.items():
        if col in sp.columns:
            feat[tag] = np.log1p(sp.loc[base_runs, col].values.astype(float))
        else:
            feat[tag] = 0.0
    return feat


# ================= 4. 任务1：单变量扫查（线性丰度，CST4 亚组） =================
def auc_safe(x, yv):
    if len(np.unique(yv)) < 2 or np.all(x == x[0]):
        return np.nan
    try:
        return roc_auc_score(yv, x)
    except Exception:
        return np.nan


rows_u = []
for cohort in ['Xiao', 'Plummer', 'Gustin']:
    sub = S[(S['cohort'] == cohort) & S['cst4']].reset_index(drop=True)
    if len(sub) < 3:
        continue
    base_runs = sub['base_run'].tolist()
    yv = sub['outcome'].values.astype(int)
    ge = GE_LIN[cohort]
    sp = SP_LIN[cohort]
    taxa_x = {}
    for gname in GENUS7:
        taxa_x[f'{gname}(属)'] = ge.loc[base_runs, gname].values.astype(float) if gname in ge.columns else np.zeros(len(sub))
    pv = [c for c in ge.columns if 'Prevotella' in c]
    taxa_x['Prevotella(属聚合)'] = ge.loc[base_runs, pv].sum(axis=1).values.astype(float) if pv else np.zeros(len(sub))
    for tag, col in SPECIES2.items():
        taxa_x[f'{tag}(种)'] = sp.loc[base_runs, col].values.astype(float) if col in sp.columns else np.zeros(len(sub))
    for tax, x in taxa_x.items():
        r_b, p_b = spearmanr(x, yv)
        a = auc_safe(x, yv)
        rows_u.append(dict(cohort=cohort, taxon=tax, n=len(sub), n_rec=int(yv.sum()),
                           auc_binary=round(a, 4) if not np.isnan(a) else np.nan,
                           rho_binary=round(r_b, 4), p_binary=round(p_b, 4)))

res_u = pd.DataFrame(rows_u)
valid = ~np.isnan(res_u['p_binary'].values)
q = np.full(len(res_u), np.nan)
q[valid] = multipletests(res_u['p_binary'].values[valid], method='fdr_bh')[1]
res_u['q_binary'] = np.round(q, 4)
res_u.to_csv(OUT_U, index=False)
print('\n=== 任务1 单变量扫查（7 属 + 2 种，CST4 亚组，线性丰度）===')
print(res_u.to_string(index=False))

# ================= 5. CV 协议（与主分析一致：真 LASSO-L1） =================
# ⚠️ 语义核实（c18_l1ratio_semantics.py, sklearn 1.9.0）：
#   liblinear + l1_ratio=1.0 == 显式 penalty='l1'（系数差 7e-6），主分析一直是真 LASSO。
#   c18 首跑误用默认 L2（Ridge），本版改为显式 L1 主协议 + L2 敏感性对照。
#   Plummer CST4（4 阴性）在 L1 下 C=0.01 收缩为空模型（零系数）-> 常数预测 -> AUC 恰好 0.5。


def make_pipe(penalty='l1'):
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', penalty=penalty,
                                                class_weight='balanced', max_iter=5000))])


GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}


def fit_tuned(Xtr, ytr, seed, penalty='l1'):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    gs = GridSearchCV(make_pipe(penalty), GRID, cv=skf, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_, gs.best_params_['clf__C']


def internal_cv(X, yv, seeds=SEEDS, collect=False, penalty='l1'):
    aucs, details = [], []
    for s in range(seeds):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for k, (tr, te) in enumerate(skf.split(X, yv)):
            if len(np.unique(yv[tr])) < 2 or len(np.unique(yv[te])) < 2:
                continue
            Xtr = X.iloc[tr] if hasattr(X, 'iloc') else X[tr]
            Xte = X.iloc[te] if hasattr(X, 'iloc') else X[te]
            mdl, best_c = fit_tuned(Xtr, yv[tr], seed=s, penalty=penalty)
            p = mdl.predict_proba(Xte)[:, 1]
            aucs.append(roc_auc_score(yv[te], p))
            if collect:
                details.append(dict(seed=s, fold=k, n_test=len(te),
                                    n_pos=int(yv[te].sum()),
                                    auc=round(aucs[-1], 4),
                                    pred_nunique=int(len(np.unique(p))),
                                    best_C=best_c,
                                    pr_auc=round(average_precision_score(yv[te], p), 4),
                                    bal_acc=round(balanced_accuracy_score(yv[te], (p >= 0.5).astype(int)), 4)))
    return np.mean(aucs), np.std(aucs), len(aucs), details


# ================= 6. 任务2：组合检验（L1 主协议 + L2 敏感性） =================
rows_c = []
diag_rows = []
for cohort in ['Xiao', 'Plummer', 'Gustin']:
    sub = S[(S['cohort'] == cohort) & S['cst4']].reset_index(drop=True)
    base_runs = sub['base_run'].tolist()
    yv = sub['outcome'].values.astype(int)
    X7 = genus_features(cohort, base_runs).reset_index(drop=True)
    Xsp = species_features(cohort, base_runs).reset_index(drop=True)
    X9 = pd.concat([X7, Xsp], axis=1)
    Xfull = G_LOG[cohort].reindex(columns=inter).fillna(0).loc[base_runs].reset_index(drop=True)

    layers = []
    for g in list(X7.columns):
        layers.append((f'① 单属:{g}', X7[[g]]))
    for s_ in list(Xsp.columns):
        layers.append((f'① 单种:{s_}', Xsp[[s_]]))
    layers.append(('② 属组合7', X7))
    layers.append(('② 种组合2', Xsp))
    layers.append(('② panel9(属7+种2)', X9))
    layers.append(('③ 全集65属', Xfull))

    for name, X in layers:
        # L1 主协议（与主分析 LASSO-LR 一致）
        collect = (cohort == 'Plummer')  # 任务3：Plummer 逐折诊断
        mu, sd, nf, det = internal_cv(X, yv, collect=collect, penalty='l1')
        rows_c.append(dict(cohort=cohort, layer=name, penalty='L1(LASSO)', n_feat=X.shape[1],
                           auc=round(mu, 4), sd=round(sd, 4), n=len(yv), n_folds=nf))
        if collect:
            for d in det:
                d['layer'] = name
                d['penalty'] = 'L1'
            diag_rows.extend(det)
        # L2 敏感性（Plummer 空模型退化的对照）
        mu2, sd2, nf2, _ = internal_cv(X, yv, collect=False, penalty='l2')
        rows_c.append(dict(cohort=cohort, layer=name, penalty='L2(敏感性)', n_feat=X.shape[1],
                           auc=round(mu2, 4), sd=round(sd2, 4), n=len(yv), n_folds=nf2))

res_c = pd.DataFrame(rows_c)
res_c.to_csv(OUT_C, index=False)
print('\n=== 任务2 组合检验（LASSO-L1 主协议 + L2 敏感性，20 种子 x 5 折 CV）===')
print(res_c.to_string(index=False))

# 三队列均值透视（按正则类型）
print('\n=== 三队列均值透视（L1 主协议）===')
piv = res_c[res_c['penalty'] == 'L1(LASSO)'].pivot_table(index='layer', columns='cohort', values='auc', aggfunc='mean')
print(piv.round(4).to_string())

# ================= 7. 任务3：Plummer 诊断（L1 空模型退化） =================
diag = pd.DataFrame(diag_rows)
diag.to_csv(OUT_D, index=False)
print('\n=== 任务3 Plummer 逐折诊断（CST4 n=43, rec=39, 仅 4 阴性）===')
print('总有效折数:', len(diag))
print('预测为常数(截距-only)的折占比: %.1f%%' % (100 * (diag['pred_nunique'] == 1).mean()))
print('best_C 分布:')
print(diag['best_C'].value_counts().to_string())
nonconst = diag[diag['pred_nunique'] > 1]
print(f'非常数折数: {len(nonconst)}; 其中 AUC 均值={nonconst["auc"].mean():.4f}' if len(nonconst) else '非常数折数: 0')

# 固定 C=1.0 敏感性（L1 无网格，看强正则放宽后是否仍空模型）
print('\n=== 固定 C=1.0 敏感性（Plummer，L1 无内层网格）===')
sub = S[(S['cohort'] == 'Plummer') & S['cst4']].reset_index(drop=True)
yv_p = sub['outcome'].values.astype(int)
X7p = genus_features('Plummer', sub['base_run'].tolist()).reset_index(drop=True)
Xfp = G_LOG['Plummer'].reindex(columns=inter).fillna(0).loc[sub['base_run'].tolist()].reset_index(drop=True)
for nm, X in [('属组合7', X7p), ('全集65属', Xfp)]:
    aucs = []
    for s in range(SEEDS):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, yv_p):
            if len(np.unique(yv_p[tr])) < 2 or len(np.unique(yv_p[te])) < 2:
                continue
            pipe = Pipeline([('sc', StandardScaler()),
                             ('clf', LogisticRegression(solver='liblinear', penalty='l1', C=1.0,
                                                        class_weight='balanced', max_iter=5000))])
            pipe.fit(X.iloc[tr], yv_p[tr])
            aucs.append(roc_auc_score(yv_p[te], pipe.predict_proba(X.iloc[te])[:, 1]))
    print(f'  {nm}: AUC={np.mean(aucs):.4f} ± {np.std(aucs):.4f} ({len(aucs)} 折)')

# Plummer 补报指标（PR-AUC / balanced accuracy，属组合7 与全集）
print('\n=== Plummer 补报指标（逐折均值，PR-AUC 多数类基线 = 39/43 = %.3f）===' % (39 / 43))
for nm in ['② 属组合7', '③ 全集65属']:
    d = diag[diag['layer'] == nm]
    print(f'  {nm}: ROC-AUC={d["auc"].mean():.4f} | PR-AUC={d["pr_auc"].mean():.4f} | bal-acc={d["bal_acc"].mean():.4f}')

print(f'\n已保存: {OUT_U}\n已保存: {OUT_C}\n已保存: {OUT_D}')
print('\n=== 种级 16S 不可得靶标（正式边界声明用）===')
for m in MISSING_SPECIES:
    print(f'  - {m}: SILVA v138.1 种级注释不可得（三队列均未注释到种）')
