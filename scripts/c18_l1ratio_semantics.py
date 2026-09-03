# -*- coding: utf-8 -*-
"""c18_l1ratio_semantics.py —— 确定 liblinear + l1_ratio=1.0 的真实正则语义
如果 l1_ratio=1.0 使 liblinear 实际跑 L1（真 LASSO），则全项目协议名副其实，
c18 不带参数的才是偏差（Ridge）；反之亦然。"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression

print('sklearn version:')
import sklearn
print(sklearn.__version__)

rng = np.random.RandomState(0)
n, p = 100, 10
X = rng.randn(n, p)
beta = np.zeros(p)
beta[:3] = [1.5, -1.0, 0.8]
y = (X @ beta + rng.randn(n) * 0.5 > 0).astype(int)

C = 0.05
variants = {
    '默认(penalty=l2, 无l1_ratio)': dict(solver='liblinear', C=C, max_iter=5000),
    'l2 显式': dict(solver='liblinear', C=C, penalty='l2', max_iter=5000),
    'l1 显式': dict(solver='liblinear', C=C, penalty='l1', max_iter=5000),
    'l1_ratio=1.0 (c7风格)': dict(solver='liblinear', C=C, l1_ratio=1.0, max_iter=5000),
    'l1_ratio=0.5': dict(solver='liblinear', C=C, l1_ratio=0.5, max_iter=5000),
    'saga elasticnet l1_ratio=1': dict(solver='saga', C=C, penalty='elasticnet', l1_ratio=1.0, max_iter=20000),
    'saga elasticnet l1_ratio=0.5': dict(solver='saga', C=C, penalty='elasticnet', l1_ratio=0.5, max_iter=20000),
}

coefs = {}
for name, kw in variants.items():
    try:
        m = LogisticRegression(**kw).fit(X, y)
        coefs[name] = m.coef_.flatten()
        z = (np.abs(m.coef_.flatten()) < 1e-10).sum()
        print(f'{name:32s} | 零系数 {z}/{p} | coef[:3]={np.round(m.coef_.flatten()[:3], 4)}')
    except Exception as e:
        print(f'{name:32s} | ERROR: {e}')

print()
ref_l2 = coefs.get('l2 显式')
ref_l1 = coefs.get('l1 显式')
c7 = coefs.get('l1_ratio=1.0 (c7风格)')
if ref_l2 is not None and c7 is not None:
    print('c7风格 vs L2 差异 :', np.abs(ref_l2 - c7).max())
if ref_l1 is not None and c7 is not None:
    print('c7风格 vs L1 差异 :', np.abs(ref_l1 - c7).max())
saga1 = coefs.get('saga elasticnet l1_ratio=1')
if saga1 is not None and ref_l1 is not None:
    print('saga(elasticnet,1) vs liblinear L1 差异 :', np.abs(ref_l1 - saga1).max())
