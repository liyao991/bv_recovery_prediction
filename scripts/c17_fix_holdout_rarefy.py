# -*- coding: utf-8 -*-
# c17_fix_holdout_rarefy.py - 修正两处协议复刻差异后重验：
#   1) holdout keep 属阈值按 diag_cst4_evidence 原样用【全队列 X_features】计算（非 Xiao 子集）
#   2) rarefy 按 richness_rarefy_check 原样【每样本重建 default_rng(42)】（共享流 -> 每样本同流）
import pandas as pd
import numpy as np
import re
import warnings
warnings.filterwarnings('ignore')
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))


X_all = pd.read_csv('<DATA_DIR>/feat/X_features.csv')
y_all = pd.read_csv('<DATA_DIR> placeholder') if False else pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
feat_meta = pd.read_csv('<DATA_DIR>/feat/feature_meta.csv')
genus_cols_x = feat_meta[feat_meta['source'] == 'genus']['feature'].tolist()
alpha_cols_x = feat_meta[feat_meta['source'] == 'alpha']['feature'].tolist()

xiao_mask = (y_all['cohort'] == 'PRJNA302078').values
y_xiao = y_all[xiao_mask].reset_index(drop=True)
Xd = X_all[xiao_mask].reset_index(drop=True)
y_x = y_xiao['outcome'].values.astype(int)
cst4_x = y_xiao['base_cst4'].values.astype(bool)

# --- 修正1：keep 属用全队列 X（diag 原样） ---
keep = [c for c in genus_cols_x if (X_all[c] > 0).mean() >= 0.05]
feats = keep + alpha_cols_x
print(f'特征数: {len(feats)} (keep 属 {len(keep)} + alpha {len(alpha_cols_x)})  [diag 原样全队列口径]')


def holdout(Xd2, yv, feats_use, n_rep=200, seed=42, test_size=0.4):
    r = np.random.default_rng(seed)
    aucs = []
    for i in range(n_rep):
        tr, te = train_test_split(np.arange(len(yv)), test_size=test_size, stratify=yv,
                                  random_state=int(r.integers(1e9)))
        pipe = Pipeline([('sc', StandardScaler()),
                         ('clf', LogisticRegression(solver='liblinear', C=1.0, l1_ratio=1.0,
                                                    class_weight='balanced', max_iter=5000,
                                                    random_state=seed))])
        pipe.fit(Xd2.iloc[tr][feats_use], yv[tr])
        aucs.append(roc_auc_score(yv[te], pipe.predict_proba(Xd2.iloc[te][feats_use])[:, 1]))
    aucs = np.array(aucs)
    return aucs.mean(), np.percentile(aucs, 2.5), np.percentile(aucs, 97.5), (aucs > 0.5).mean()


m1, lo1, hi1, p1 = holdout(Xd, y_x, feats)
print(f'Xiao full  holdout 200x: {m1:.4f} CI[{lo1:.3f},{hi1:.3f}] P(>0.5)={p1:.3f}   [框架 0.652 [0.500,0.797]]')
Xd4 = Xd[cst4_x].reset_index(drop=True)
m2, lo2, hi2, p2 = holdout(Xd4, y_x[cst4_x], [c for c in feats if c in Xd4.columns])
print(f'Xiao CST4  holdout 200x: {m2:.4f} CI[{lo2:.3f},{hi2:.3f}] P(>0.5)={p2:.3f}   [框架 0.596 [0.388,0.802]]')

# --- 修正2：rarefy 每样本重建 rng(42)（原协议） ---
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
y_g = (G_rel.loc[list(g_out)]['Lactobacillus'] >= 0.5).values.astype(int)


def load_counts(path):
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    df.index = [norm_sid(i) for i in df.index]
    return df


counts_g = load_counts('<OUT_DIR>/PRJNA691964/asv_table.csv')
gc = counts_g.loc[list(g_base)]
depth = int(gc.sum(axis=1).min())

richA = []
for i in range(gc.shape[0]):
    c = gc.iloc[i].values
    p = c / c.sum()
    rng = np.random.default_rng(42)          # 原协议：每样本同流
    samp = rng.multinomial(depth, p, size=200)
    richA.append((samp > 0).sum(axis=1).mean())
richA = np.array(richA)
rA, pA = stats.spearmanr(richA, y_g)
print(f'\nrarefy rho (每样本 rng(42), 原协议复刻): {rA:.4f} (p={pA:.4f})   [框架 -0.392 / p=0.024]')
