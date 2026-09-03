# -*- coding: utf-8 -*-
"""c20_clr_pseudo_robust.py —— #C20 追加：CLR 伪计数选择稳健性快检
伪计数 1e-2 / 1e-3 / 1e-4 / 1e-6 四档，B 系三 cell + A1，L1 主协议。
目的：确认 CLR 结论不依赖伪计数选择。"""
import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

np.random.seed(42)
SEEDS = 20
GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}

# ---- 数据加载（复刻） ----
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
    return df.T.groupby(level=0).sum().T

X_lin = agg_rel(pd.DataFrame(np.expm1(X_xiao_log1p[genus_cols_x]).clip(lower=0),
                             columns=genus_cols_x, index=X_xiao_log1p.index))
M_lin = agg_rel(M_rel)
G_lin = agg_rel(G_rel)
inter = sorted(set(X_lin.columns) & set(M_lin.columns) & set(G_lin.columns))

y_x = y_xiao['outcome'].values.astype(int)
cst4_x = y_xiao['base_cst4'].values.astype(bool)
y_m = M_pairs['outcome'].values.astype(int)
cst4_m = M_pairs['base_cst4'].values.astype(bool)
g_pairs = []
for pat, sub in G_map.groupby('patient'):
    v1 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '01']
    v3 = sub[sub['visit'].astype(int).astype(str).str.zfill(2) == '03']
    if len(v1) == 1 and len(v3) == 1:
        g_pairs.append((v1.iloc[0]['run'], v3.iloc[0]['run']))
g_base, g_out = zip(*g_pairs)
G_out_lb = G_rel.loc[list(g_out)]['Lactobacillus'].values
g_outcome = (G_out_lb >= 0.5).astype(int)
G_base_lb = G_rel.loc[list(g_base)]['Lactobacillus'].values
g_cst4 = (G_base_lb < 0.5).astype(bool)

def clr(lin_df, c):
    v = lin_df.reindex(columns=inter).fillna(0).values.astype(float)
    lg = np.log(v + c)
    return pd.DataFrame(lg - lg.mean(axis=1, keepdims=True), columns=inter,
                        index=lin_df.reindex(columns=inter).fillna(0).index)

def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', penalty='l1',
                                                class_weight='balanced', max_iter=5000))])

def fit_tuned(Xtr, ytr, seed):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    gs = GridSearchCV(make_pipe(), GRID, cv=skf, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_

def internal_cv(X, y):
    aucs = []
    for s in range(SEEDS):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, y):
            mdl = fit_tuned(X.iloc[tr], y[tr], seed=s)
            aucs.append(roc_auc_score(y[te], mdl.predict_proba(X.iloc[te])[:, 1]))
    return np.mean(aucs), np.std(aucs)

def external_auc(Xtr, ytr, Xte, yte):
    aucs = []
    for s in range(SEEDS):
        mdl = fit_tuned(Xtr, ytr, seed=s)
        aucs.append(roc_auc_score(yte, mdl.predict_proba(Xte)[:, 1]))
    return np.mean(aucs), np.std(aucs)

rows = []
for c in [1e-2, 1e-3, 1e-4, 1e-6]:
    Xc = clr(X_lin, c).reset_index(drop=True)
    Mc = clr(M_lin, c)
    Gc = clr(G_lin, c)
    Mv1 = Mc.loc[M_pairs['base_sample'].tolist()].reset_index(drop=True)
    Gb = Gc.loc[list(g_base)].reset_index(drop=True)
    mu_a1, sd_a1 = internal_cv(Xc, y_x)
    mu_b1, sd_b1 = internal_cv(Xc[cst4_x].reset_index(drop=True), y_x[cst4_x])
    mu_b2, _ = external_auc(Xc[cst4_x].reset_index(drop=True), y_x[cst4_x],
                            Mv1[cst4_m].reset_index(drop=True), y_m[cst4_m])
    mu_b3, _ = external_auc(Xc[cst4_x].reset_index(drop=True), y_x[cst4_x],
                            Gb[g_cst4].reset_index(drop=True), g_outcome[g_cst4])
    rows.append(dict(pseudo=c, A1=round(mu_a1, 4), B1=round(mu_b1, 4),
                     B2=round(mu_b2, 4), B3=round(mu_b3, 4)))
    print(f'pseudo={c:g}: A1 {mu_a1:.3f} | B1 {mu_b1:.3f} | B2 {mu_b2:.3f} | B3 {mu_b3:.3f}')

pd.DataFrame(rows).to_csv('<OUT_DIR>/modeling/results_c20_clr_pseudo.csv', index=False)
print('已保存 results_c20_clr_pseudo.csv')
