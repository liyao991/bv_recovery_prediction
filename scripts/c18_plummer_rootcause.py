# -*- coding: utf-8 -*-
"""
c18_plummer_rootcause.py —— 定位 c7_panel_combo 中 Plummer AUC=0.5/sd=0 的确切根因
c18 新协议下 Plummer 不再是 0.5/0.0（单属 Lactobacillus 0.6219），须解释差异来源。
怀疑点：① liblinear 无 random_state 的全局 RNG 蒙特卡洛噪声（#C17 已知）
       ② c7 的 l1_ratio=1.0 参数被 liblinear 忽略（实际 penalty='l2'）
       ③ C=0.01 强收缩下系数数值下溢 -> 预测完全常数 -> AUC 恰好 0.5
"""
import pandas as pd
import numpy as np
import re
import warnings
warnings.filterwarnings('ignore')
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

SEEDS = 20


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


# ==== c7 确切数据路径 ====
def load_genus_log1p(pid):
    g = pd.read_csv(f'<DATA_DIR>/feat/{pid}_genus_rel.csv', index_col=0)
    g.index = [norm_sid(i) for i in g.index]
    g.columns = [norm_genus(c) for c in g.columns]
    agg = g.T.groupby(level=0).sum().T
    return np.log1p(agg)


G = {c: load_genus_log1p(p) for c, p in
     [('Xiao', 'PRJNA302078'), ('Plummer', 'PRJNA592384'), ('Gustin', 'PRJNA691964')]}
inter = sorted(set(G['Xiao'].columns) & set(G['Plummer'].columns) & set(G['Gustin'].columns))

y = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
plum = y[y['cohort'] == 'PRJNA592384'].copy()
plum['base'] = [norm_sid(s) for s in plum['base_sample']]
plum = plum[plum['base_cst4'] == 1]
X_p = G['Plummer'].reindex(columns=inter).fillna(0).loc[plum['base'].tolist()]
yv = plum['outcome'].values.astype(int)
print(f'Plummer CST4: n={len(yv)}, rec={yv.sum()}')

# 验证 c18 的 X7 与 c7 的 G 列数值一致性（单属 Lactobacillus）
XL = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_genus_rel.csv', index_col=0)
XL.index = [norm_sid(i) for i in XL.index]
x18 = np.log1p(XL.loc[plum['base'].tolist(), 'Lactobacillus'].values)
x7c = X_p['Lactobacillus'].values
print('单属 Lactobacillus c7 vs c18 特征最大差:', np.abs(x18 - x7c).max())


def make_pipe(c7_style=True, C=None, rs=None):
    kw = dict(solver='liblinear', class_weight='balanced', max_iter=5000)
    if c7_style:
        kw['l1_ratio'] = 1.0  # c7 原样（liblinear 忽略之）
    if rs is not None:
        kw['random_state'] = rs
    if C is not None:
        kw['C'] = C
    return Pipeline([('sc', StandardScaler()), ('clf', LogisticRegression(**kw))])


GRID = {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}


def run_variant(tag, pipe_factory, use_grid=True):
    aucs, consts, coefs, bestcs = [], 0, [], []
    for s in range(SEEDS):
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        for tr, te in skf.split(X_p, yv):
            if len(np.unique(yv[tr])) < 2 or len(np.unique(yv[te])) < 2:
                continue
            if use_grid:
                gs = GridSearchCV(pipe_factory(), GRID,
                                  cv=StratifiedKFold(5, shuffle=True, random_state=s),
                                  scoring='roc_auc', n_jobs=1)
                gs.fit(X_p.iloc[tr], yv[tr])
                mdl, bc = gs.best_estimator_, gs.best_params_['clf__C']
            else:
                mdl = pipe_factory()
                mdl.fit(X_p.iloc[tr], yv[tr])
                bc = np.nan
            p = mdl.predict_proba(X_p.iloc[te])[:, 1]
            if len(np.unique(p)) == 1:
                consts += 1
            aucs.append(roc_auc_score(yv[te], p))
            coefs.append(mdl.named_steps['clf'].coef_.flatten())
            bestcs.append(bc)
    coefs = np.array(coefs)
    print(f'\n[{tag}] AUC={np.mean(aucs):.4f}±{np.std(aucs):.4f} | 常数预测折 {consts}/{len(aucs)} '
          f'| |coef| max={np.abs(coefs).max():.3e} | best_C mode={pd.Series(bestcs).mode().tolist()}')
    return aucs


# ① c7 原样：无 random_state、l1_ratio=1.0（被忽略）、网格 CV
np.random.seed(42)
run_variant('c7 原样（无 random_state + 网格）', lambda: make_pipe(c7_style=True))

# ② c7 + 固定 random_state=0
run_variant('c7 + random_state=0', lambda: make_pipe(c7_style=True, rs=0))

# ③ 无网格固定 C=0.01（诊断强收缩）
run_variant('固定 C=0.01 无网格', lambda: make_pipe(c7_style=True, C=0.01), use_grid=False)

# ④ 无网格固定 C=1.0
run_variant('固定 C=1.0 无网格', lambda: make_pipe(c7_style=True, C=1.0), use_grid=False)

# ⑤ 单属 Lactobacillus（c7 的 ① 层，无网格 C=0.01，看是否常数预测）
X1 = X_p[['Lactobacillus']]
aucs, consts = [], 0
for s in range(SEEDS):
    skf = StratifiedKFold(5, shuffle=True, random_state=s)
    for tr, te in skf.split(X1, yv):
        if len(np.unique(yv[tr])) < 2 or len(np.unique(yv[te])) < 2:
            continue
        mdl = make_pipe(c7_style=True, C=0.01)
        mdl.fit(X1.iloc[tr], yv[tr])
        p = mdl.predict_proba(X1.iloc[te])[:, 1]
        if len(np.unique(p)) == 1:
            consts += 1
        aucs.append(roc_auc_score(yv[te], p))
print(f'\n[单属 Lactobacillus 固定 C=0.01] AUC={np.mean(aucs):.4f}±{np.std(aucs):.4f} | 常数预测折 {consts}/{len(aucs)}')
print(f'  系数范围示例: {mdl.named_steps["clf"].coef_.flatten()}')

# ⑥ 有效折内阴性数分布（说明 AUC 方差来源）
neg_per_fold = []
for s in range(SEEDS):
    skf = StratifiedKFold(5, shuffle=True, random_state=s)
    for tr, te in skf.split(X_p, yv):
        if len(np.unique(yv[tr])) < 2 or len(np.unique(yv[te])) < 2:
            continue
        neg_per_fold.append(int((yv[te] == 0).sum()))
vc = pd.Series(neg_per_fold).value_counts().sort_index()
print('\n有效测试折内阴性数分布:', dict(vc))
print('=> 每折最多 1 个阴性，AUC 单折取值只能落在 {0, 1/(n_pos+1)...1} 的粗格点上')
