# -*- coding: utf-8 -*-
"""
Xiao(454/DADA2/SILVA) 训练 -> Mollin(MiSeq/Mothur/RDP) 预测
====================================================
三重外推验证：跨测序平台(454->MiSeq) + 跨分类管线(DADA2/SILVA->Mothur/RDP)
            + 跨人群(中国->美国复发性BV)
统一结局：CST4 = LB<50%，outcome = 结局非 CST4(康复)

特征对齐：属级交集（Xiao 属名截断到属级 vs Mollin 属名清洗）
输出：<OUT_DIR>/modeling/results_xiao_mollin_cross.csv
"""
import pandas as pd
import numpy as np
import re
from collections import Counter
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from xgboost import XGBClassifier

RNG = 42
OUT = '<OUT_DIR>/modeling/results_xiao_mollin_cross.csv'
np.random.seed(RNG)

# ---------- 1. 加载 Xiao 特征与标签 ----------
X_all = pd.read_csv('<DATA_DIR>/feat/X_features.csv')
y_all = pd.read_csv('<DATA_DIR>/feat/y_labels.csv')
feat_meta = pd.read_csv('<DATA_DIR>/feat/feature_meta.csv')
genus_cols_x = feat_meta[feat_meta['source'] == 'genus']['feature'].tolist()
alpha_cols_x = feat_meta[feat_meta['source'] == 'alpha']['feature'].tolist()

xiao_mask = (y_all['cohort'] == 'PRJNA302078').values
X_xiao = X_all[xiao_mask].reset_index(drop=True)
y_xiao = y_all[xiao_mask].reset_index(drop=True)
print(f'Xiao: {len(X_xiao)} 例, outcome 分布: {y_xiao["outcome"].value_counts().to_dict()}')
print(f'Xiao 基线 CST4: {y_xiao["base_cst4"].value_counts().to_dict()}')

# ---------- 2. 加载 Mollin 特征与配对 ----------
M_X = pd.read_csv('<DATA_DIR>/meta/mollin_genus_log1p.csv', index_col=0)
M_pairs = pd.read_csv('<DATA_DIR>/meta/mollin_pairs.csv')
print(f'\nMollin 特征矩阵: {M_X.shape} | 配对: {len(M_pairs)}')
print(f'Mollin 配对 outcome 分布: {M_pairs["outcome"].value_counts().to_dict()}')

# ---------- 3. 属名标准化（双侧对齐到属级） ----------
def norm_genus(name):
    """属名标准化：去方括号、取第一个词、去 'sp.'/种级后缀/组群后缀"""
    s = str(name).strip()
    if s in ('NA_unclassified', 'nan', 'NA', 'NaN'):
        return 'NA_unclassified'
    # 去方括号（SILVA 的 [Eubacterium] 形式）
    s = s.replace('[', '').replace(']', '')
    # uncultured/未知 归为 NA
    if s.lower().startswith('uncultured') or s.lower().startswith('unclassified') or s == 'Bacteria':
        return 'NA_unclassified'
    # 取第一个词
    first = s.split(' ')[0]
    # 纯种级残渣（如 pseudomesenteroides、pullorum_sp、Sub）
    if first in ('Sub', 'sp', 'group', 'incertae'):
        return 'NA_unclassified'
    return first

# Xiao 侧属名
xiao_genus_map = {}
for c in genus_cols_x:
    xiao_genus_map[c] = norm_genus(c)
# 检查截断冲突（两个原始列截断到同一属）
mapped = [xiao_genus_map[c] for c in genus_cols_x]
dup = [k for k, v in Counter(mapped).items() if v > 1]
if dup:
    print(f'⚠️ Xiao 截断冲突属: {dup[:10]}')

# Mollin 侧属名
mollin_genus_map = {}
for c in M_X.columns:
    mollin_genus_map[c] = norm_genus(c)

# ---------- 4. 构造交集特征空间（相对丰度域聚合 -> log1p） ----------
def agg_rel_to_log1p(rel_df, col_map):
    """相对丰度矩阵按新属名聚合（相对丰度域求和）-> 转 log1p"""
    df = rel_df.copy()
    df.columns = [col_map[c] for c in df.columns]
    # 相对丰度域聚合
    agg = df.T.groupby(level=0).sum().T
    return np.log1p(agg)

# Xiao 侧：X_features 是 log1p，用 expm1 还原相对丰度
X_rel = np.expm1(X_xiao[genus_cols_x]).clip(lower=0)
X_rel = pd.DataFrame(X_rel, columns=genus_cols_x, index=X_xiao.index)
xiao_agg = agg_rel_to_log1p(X_rel, xiao_genus_map)

# Mollin 侧：直接用相对丰度矩阵
M_rel = pd.read_csv('<DATA_DIR>/meta/mollin_genus_rel.csv', index_col=0)
mollin_agg = agg_rel_to_log1p(M_rel, mollin_genus_map)

# 交集
inter = sorted(set(xiao_agg.columns) & set(mollin_agg.columns))
print(f'\n属名交集: {len(inter)}')
print('交集属示例:', inter[:20])
# 关键 BV 属是否在交集
for g in ['Lactobacillus', 'Gardnerella', 'Atopobium', 'Mobiluncus', 'Prevotella',
          'Megasphaera', 'Sneathia', 'Fannyhessea', 'Anaerococcus', 'Finegoldia',
          'Peptoniphilus', 'Dialister', 'Mycoplasma', 'Ureaplasma', 'Eggerthella',
          'Bifidobacterium', 'Streptococcus', 'Staphylococcus', 'Corynebacterium']:
    print(f'  {g}: {"✅" if g in inter else "❌"}')

Xtr = xiao_agg[inter].reset_index(drop=True)
Mfeat = mollin_agg.reindex(columns=inter).fillna(0)
print(f'\nXiao 训练矩阵: {Xtr.shape} | Mollin 预测矩阵: {Mfeat.shape}')

# ---------- 5. Mollin 配对 V1 特征 ----------
M_Xi = Mfeat.set_index(M_X.index)
mollin_v1 = M_pairs[M_pairs['base_sample'].isin(M_Xi.index)]
print(f'Mollin 配对匹配到特征: {len(mollin_v1)} / {len(M_pairs)}')
M_v1 = M_Xi.loc[mollin_v1['base_sample'].tolist()].reset_index(drop=True)
M_y = mollin_v1['outcome'].reset_index(drop=True).values
print(f'Mollin 测试集: {M_v1.shape}, 阳性率 {M_y.mean():.3f}')

# ---------- 6. 模型工厂 ----------
def make_models():
    return [
        ('LASSO-LR', Pipeline([('sc', StandardScaler()),
                               ('clf', LogisticRegression(solver='liblinear', l1_ratio=1.0,
                                                          class_weight='balanced', max_iter=5000))]),
         {'clf__C': [0.01, 0.05, 0.1, 0.5, 1.0]}),
        ('RF', Pipeline([('sc', StandardScaler()),
                         ('clf', RandomForestClassifier(n_estimators=400, class_weight='balanced',
                                                        random_state=RNG, n_jobs=1))]),
         {'clf__max_depth': [2, 3, 5], 'clf__min_samples_leaf': [2, 5]}),
        ('XGB', Pipeline([('sc', StandardScaler()),
                          ('clf', XGBClassifier(eval_metric='logloss', random_state=RNG, n_jobs=1))]),
         {'clf__max_depth': [2, 3], 'clf__learning_rate': [0.02, 0.05],
          'clf__n_estimators': [150, 300]}),
    ]

# ---------- 7. Xiao 内部 CV（参考） + 外推预测 ----------
from sklearn.pipeline import Pipeline
results = []
ytr = y_xiao['outcome'].values

print('\n=== Xiao 内部 5 折 CV（同特征空间参考）===')
for name, pipe, grid in make_models():
    aucs, prs, briers = [], [], []
    skf = StratifiedKFold(5, shuffle=True, random_state=RNG)
    for tr, te in skf.split(Xtr, ytr):
        gs = GridSearchCV(pipe, grid, cv=3, scoring='roc_auc', n_jobs=1)
        gs.fit(Xtr.iloc[tr], ytr[tr])
        p = gs.predict_proba(Xtr.iloc[te])[:, 1]
        aucs.append(roc_auc_score(ytr[te], p))
        prs.append(average_precision_score(ytr[te], p))
        briers.append(brier_score_loss(ytr[te], p))
    print(f'  {name}: CV AUC {np.mean(aucs):.3f}±{np.std(aucs):.3f} | PR {np.mean(prs):.3f} | Brier {np.mean(briers):.3f}')

print('\n=== Xiao 训练 -> Mollin 外推 ===')
for name, pipe, grid in make_models():
    # 内层 CV 选超参（Xiao 内 5 折）
    skf = StratifiedKFold(5, shuffle=True, random_state=RNG)
    gs = GridSearchCV(pipe, grid, cv=skf, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr)
    p = gs.predict_proba(M_v1)[:, 1]
    auc = roc_auc_score(M_y, p)
    pr = average_precision_score(M_y, p)
    br = brier_score_loss(M_y, p)
    # bootstrap 95% CI
    rng = np.random.RandomState(123)
    aucs_b, prs_b = [], []
    for _ in range(1000):
        idx = rng.randint(0, len(M_y), len(M_y))
        if len(np.unique(M_y[idx])) < 2:
            continue
        aucs_b.append(roc_auc_score(M_y[idx], p[idx]))
        prs_b.append(average_precision_score(M_y[idx], p[idx]))
    lo_a, hi_a = np.percentile(aucs_b, [2.5, 97.5])
    lo_p, hi_p = np.percentile(prs_b, [2.5, 97.5])
    print(f'  {name}: AUC {auc:.3f} [95%CI {lo_a:.3f}-{hi_a:.3f}] | '
          f'PR-AUC {pr:.3f} [{lo_p:.3f}-{hi_p:.3f}] | Brier {br:.3f} | 随机基线 PR={M_y.mean():.3f}')
    results.append(dict(model=name, auc=round(auc, 4), auc_ci_lo=round(lo_a, 4),
                        auc_ci_hi=round(hi_a, 4), pr_auc=round(pr, 4),
                        pr_ci_lo=round(lo_p, 4), pr_ci_hi=round(hi_p, 4),
                        brier=round(br, 4), n_test=len(M_y), pos_rate=round(M_y.mean(), 4)))

# ---------- 8. 反向：Mollin 训练 -> Xiao 预测（对照） ----------
print('\n=== 反向: Mollin 训练 -> Xiao 预测（对照）===')
for name, pipe, grid in make_models():
    if name != 'LASSO-LR':
        continue
    skf = StratifiedKFold(5, shuffle=True, random_state=RNG)
    gs = GridSearchCV(pipe, grid, cv=skf, scoring='roc_auc', n_jobs=1)
    gs.fit(M_v1, M_y)
    p = gs.predict_proba(Xtr)[:, 1]
    auc = roc_auc_score(ytr, p)
    pr = average_precision_score(ytr, p)
    print(f'  {name}: Xiao 测试 AUC {auc:.3f} | PR-AUC {pr:.3f} | Xiao 阳性率 {ytr.mean():.3f}')

# ---------- 9. 保存 + 输出关键属的跨队列丰度对比 ----------
pd.DataFrame(results).to_csv(OUT, index=False)
print(f'\n已保存: {OUT}')

print('\n=== 关键属在两侧队列的平均相对丰度（log1p）===')
for g in ['Lactobacillus', 'Gardnerella', 'Atopobium', 'Mobiluncus', 'Prevotella', 'Megasphaera', 'Sneathia']:
    if g in inter:
        print(f'  {g}: Xiao {Xtr[g].mean():.3f} vs Mollin {M_v1[g].mean():.3f}')
