# dangC_planA_analysis.py
# 档 C 方案 A：基线(M0)特征 vs "1 个月内复发(3) vs 未复发(65)" 结局
# 输出：richness 方向性(ρ/p) + 单特征 AUC(Lactobacillus/richness) + LASSO 全特征 AUC(标注过拟合)
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))

lab = pd.read_csv('<OUT_DIR>/modeling/dangC_labels.csv')  # subject, base_sample, recurrence
lab['base'] = [norm_sid(s) for s in lab['base_sample']]
y = lab['recurrence'].values.astype(int)

af = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_alpha_family.csv')
af['run'] = [norm_sid(s) for s in af['sample']]; af = af.set_index('run')
rich = af.loc[lab['base']]['species_observed'].values

genus = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_genus_rel.csv', index_col=0)
genus.index = [norm_sid(i) for i in genus.index]
G = genus.loc[lab['base']]

species = pd.read_csv('<DATA_DIR>/feat/PRJNA592384_species_rel.csv', index_col=0)
species.index = [norm_sid(i) for i in species.index]
S = species.loc[lab['base']]

print('=== 档 C 方案 A：基线 M0 特征 vs 1 个月内复发（n=68, 复发=3）===')
print()
# 1. richness 方向性
r, p = spearmanr(rich, y)
print(f'richness vs 复发: rho={r:.3f}, p={p:.3f}')
print(f'  复发组 richness median={np.median(rich[y==1]):.1f}, 治愈组 median={np.median(rich[y==0]):.1f}')

# 2. 单特征 AUC
def single_auc(x, y):
    if len(np.unique(y)) < 2: return np.nan
    return roc_auc_score(y, x)
lb_col = 'Lactobacillus' if 'Lactobacillus' in G.columns else [c for c in G.columns if 'Lacto' in c][0]
auc_lb = single_auc(G[lb_col].values, y)
auc_rich = single_auc(rich, y)
print(f'单特征 AUC: Lactobacillus={auc_lb:.3f}, richness={auc_rich:.3f}')

# 3. LASSO 全特征（StratifiedKFold，标注过拟合风险）
def lasso_auc(X, y):
    X = np.asarray(X, float)
    X = np.log1p(X)
    if len(np.unique(y)) < 2: return np.nan, None
    # 3 例事件，3 折 CV 每折只有 1 个正例，用 LOO 风格近似
    clf = LogisticRegression(penalty='l1', solver='liblinear', C=0.1, class_weight='balanced', max_iter=2000)
    try:
        pred = cross_val_predict(clf, X, y, cv=StratifiedKFold(3, shuffle=True, random_state=42), method='decision_function')
        return roc_auc_score(y, pred), pred
    except Exception as e:
        return np.nan, None

auc_genus, _ = lasso_auc(G.values, y)
auc_species, _ = lasso_auc(S.values, y)
print(f'LASSO 全特征 AUC(标注过拟合): 属级={auc_genus:.3f}, 种级={auc_species:.3f}')

# 落盘
out = pd.DataFrame([dict(
    cohort='Plummer(PRJNA592384)', outcome='1mo_recurrence', n=68, n_events=3,
    min_detect_auc=0.921,
    rho_rich=round(r,4), p_rich=round(p,4),
    auc_LB=round(auc_lb,4), auc_richness=round(auc_rich,4),
    auc_genus=round(auc_genus,4) if auc_genus==auc_genus else None,
    auc_species=round(auc_species,4) if auc_species==auc_species else None,
    note='方案A:3例事件,功效崩塌(0.921),AUC不可靠仅参考')])
out.to_csv('<OUT_DIR>/modeling/results_dangC_planA.csv', index=False)
print()
print('已写入 <OUT_DIR>/modeling/results_dangC_planA.csv')
