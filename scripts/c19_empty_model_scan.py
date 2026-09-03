# -*- coding: utf-8 -*-
"""
c19_empty_model_scan.py —— #C19 追加：全部 LOCO cell 空模型种子比例扫描
=================================================================
背景：c19_gustin_ci_probe.py 发现 Gustin->Xiao 20 种子中 13 个空模型
（非零系数=0 -> 常数预测 -> AUC 恒 0.5）。空模型种子的 AUC=0.5 是占位符
而非性能测量（#C18 Plummer 定案同类机制），会人为压窄 seed-mean 的 SD 与
bootstrap CI。本脚本扫描全部 12 个 LOCO cell 的空模型种子数 + 分解
「非空种子子集 AUC」。
"""
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
    return np.log1p(df.T.groupby(level=0).sum().T)

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
G_out_df = G_f.loc[list(g_out)].reset_index(drop=True)
g_outcome = (G_out_df['Lactobacillus'] >= np.log1p(0.5)).values.astype(int)
g_cst4 = (G_base['Lactobacillus'] < np.log1p(0.5)).values.astype(bool)

DATA = {
    'Xiao':   {'full': (X_f, y_x), 'cst4': (X_f[cst4_x].reset_index(drop=True), y_x[cst4_x])},
    'Mollin': {'full': (M_v1, y_m), 'cst4': (M_v1[cst4_m].reset_index(drop=True), y_m[cst4_m])},
    'Gustin': {'full': (G_base, g_outcome), 'cst4': (G_base[g_cst4].reset_index(drop=True), g_outcome[g_cst4])},
}

def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])
GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}

def fit_tuned(Xtr, ytr, seed):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    gs = GridSearchCV(make_pipe(), GRID, cv=skf, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    return gs.best_estimator_

CELLS = [('Xiao', 'Mollin'), ('Xiao', 'Gustin'),
         ('Mollin', 'Xiao'), ('Mollin', 'Gustin'),
         ('Gustin', 'Xiao'), ('Gustin', 'Mollin')]

rows = []
print(f'{"cell":36s} {"空模型种子":>8s} {"全部AUC":>8s} {"非空子集AUC":>10s} {"非空n":>5s} {"平均非零系数":>10s}')
for tr_c, te_c in CELLS:
    for scope in ['full', 'cst4']:
        Xtr, ytr = DATA[tr_c][scope]
        Xte, yte = DATA[te_c][scope]
        aucs, ncoefs, empty = [], [], 0
        for s in range(SEEDS):
            mdl = fit_tuned(Xtr, ytr, seed=s)
            coef = mdl.named_steps['clf'].coef_[0]
            nz = int((coef != 0).sum())
            ncoefs.append(nz)
            if nz == 0:
                empty += 1
                aucs.append(0.5)          # 常数预测 AUC 恒 0.5
            else:
                p = mdl.predict_proba(Xte)[:, 1]
                aucs.append(roc_auc_score(yte, p))
        aucs = np.array(aucs)
        nonempty = [a for a, nz in zip(aucs, ncoefs) if nz > 0]
        tag = f'{tr_c}->{te_c} ({scope})'
        print(f'{tag:36s} {empty:>5d}/20 {aucs.mean():>8.3f} '
              f'{(np.mean(nonempty) if nonempty else float("nan")):>10.3f} {len(nonempty):>5d} {np.mean(ncoefs):>10.1f}')
        rows.append(dict(cell=tag, empty_seeds=empty, auc_all=round(aucs.mean(), 4),
                         auc_nonempty=round(float(np.mean(nonempty)), 4) if nonempty else None,
                         n_nonempty=len(nonempty), mean_nz_coef=round(float(np.mean(ncoefs)), 1)))

pd.DataFrame(rows).to_csv('<OUT_DIR>/modeling/results_c19_empty_scan.csv', index=False)
print('\n已保存: <OUT_DIR>/modeling/results_c19_empty_scan.csv')
