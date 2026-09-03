# -*- coding: utf-8 -*-
# c16_tmp_alpha_matrix_main.py - Gustin alpha family x resolution under MAIN protocol (LB>=64.87%)
import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score


def norm_sid(s):
    return re.sub(r'_F\.fastq\.gz$', '', str(s))


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
y = (G_rel.loc[list(g_out)]['Lactobacillus'] >= 0.5).values.astype(int)  # 主口径
a = af.loc[list(g_base)]

METRICS = ['observed', 'chao1', 'shannon', 'simpson', 'invsimpson', 'pielou']


def make_pipe():
    return Pipeline([('sc', StandardScaler()),
                     ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                class_weight='balanced', max_iter=5000))])


def cv_auc(X, yv):
    aucs = []
    for s in range(20):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X, yv):
            if len(np.unique(yv[tr])) < 2 or len(np.unique(yv[te])) < 2:
                continue
            mdl = make_pipe().fit(X.iloc[tr], yv[tr])
            aucs.append(roc_auc_score(yv[te], mdl.predict_proba(X.iloc[te])[:, 1]))
    v = np.array(aucs)
    return v.mean(), v.std()


print('Gustin 主口径（LB>=64.87%）alpha 家族 x 分辨率（对照 7.10 表）：')
print(f'{"分辨率":<6} {"richness":<18} {"evenness":<12} {"entropy":<12} {"full":<12}')
for lv in ['asv', 'genus', 'species']:
    cols = {m: f'{lv}_{m}' for m in METRICS}
    X_lv = a[[cols[m] for m in METRICS]].copy()
    X_lv['simpson_even'] = X_lv[cols['invsimpson']] / X_lv[cols['observed']].clip(lower=1)
    fams = {
        'richness': [cols['observed'], cols['chao1']],
        'evenness': [cols['pielou']],
        'entropy': [cols['shannon'], cols['simpson'], cols['invsimpson']],
        'full': list(X_lv.columns),
    }
    cells = []
    for fam, ms in fams.items():
        mu, sd = cv_auc(X_lv[ms], y)
        cells.append(f'{mu:.3f}±{sd:.3f}')
    print(f'{lv:<6} {cells[0]:<18} {cells[1]:<12} {cells[2]:<12} {cells[3]:<12}')

# 单指标 species observed/chao1（主口径）
for m in ['observed', 'chao1']:
    mu, sd = cv_auc(a[[f'species_{m}']], y)
    print(f'单指标 species_{m}: CV AUC {mu:.3f} ± {sd:.3f}')
