# -*- coding: utf-8 -*-
"""
Mollin 2022 属级特征矩阵构建（供跨队列外推实验）
====================================================
从 S2 (OTU->Taxonomy) + S4 (865 species 计数表) 构建：
    1. 属级相对丰度矩阵（样本 x 属）
    2. 与 Xiao/SToPBV 相同的 log1p 变换（Xiao X_features.csv 口径）
    3. 属名标准化（Mothur/RDP 属名 vs SILVA 属名对齐用）

输出：
    <DATA_DIR>/meta/mollin_genus_rel.csv — 样本 x 属 相对丰度（0-1）
    <DATA_DIR>/meta/mollin_genus_log1p.csv — log1p 变换后（与 Xiao X_features 同口径）
"""
import pandas as pd
import numpy as np

XLSX = '<DATA_DIR>/meta/mollin2022_S1.xlsx'
OUT_DIR = '<DATA_DIR>/meta/'

# ---------- 1. 解析 S2: Otu -> Taxonomy ----------
s2 = pd.read_excel(XLSX, sheet_name='S2Table')
s2.columns = ['hdr', 'tax', 'seq']
s2 = s2.iloc[1:].reset_index(drop=True)
otu2tax = {}
for _, r in s2.iterrows():
    otu = str(r['hdr']).strip()
    tax = str(r['tax']).strip()
    if not otu.startswith('Otu') or not tax:
        continue
    otu2tax[otu] = tax
print(f'S2 映射解析: {len(otu2tax)} 个 OTU 带分类学')


def tax_to_genus(tax):
    """'>Firmicutes;...;Family;Genus;Species*;OtuXXXX' -> Genus"""
    t = tax.lstrip('>').strip()
    parts = [p for p in t.split(';') if p.strip()]
    if parts and parts[-1].startswith('Otu'):
        parts = parts[:-1]
    if not parts:
        return 'NA'
    if len(parts) >= 2:
        genus = parts[-2].strip()
        if genus.endswith('*') or genus == '':
            genus = parts[-1].strip()
    else:
        genus = parts[0].strip()
    genus = genus.replace('*', '').strip()
    return genus


otu2genus = {o: tax_to_genus(t) for o, t in otu2tax.items()}
from collections import Counter
gc = Counter(otu2genus.values())
print('Mollin 属级分类单元数:', len(gc))
print('top 15:', gc.most_common(15))

# ---------- 2. 读 S4: 865 species 计数表 ----------
s4 = pd.read_excel(XLSX, sheet_name='S4 Table', header=None)
samples = [str(v).strip() for v in s4.iloc[1, 1:].tolist()]
classes = [str(v).strip() for v in s4.iloc[2, 1:].tolist()]
otu_rows = s4.iloc[3:]

data = {}
for _, r in otu_rows.iterrows():
    otu = str(r.iloc[0]).strip()
    if not otu.startswith('Otu'):
        continue
    vals = pd.to_numeric(r.iloc[1:], errors='coerce').fillna(0).values
    data[otu] = vals
counts = pd.DataFrame(data, index=samples).T  # OTU x 样本
counts = counts[counts.index.isin(otu2genus.keys())]
print(f'匹配到分类学的 OTU: {counts.shape[0]} / {len(data)}')

# ---------- 3. 属级聚合 + 相对丰度 ----------
genus_level = counts.groupby(counts.index.map(lambda o: otu2genus.get(o, 'NA'))).sum()
genus_pct = genus_level.div(genus_level.sum(axis=0), axis=1).T  # 样本 x 属
genus_pct = genus_pct.round(6)
print(f'\nMollin 属级相对丰度矩阵: {genus_pct.shape}')
print('列（属）示例:', list(genus_pct.columns[:12]))

# ---------- 4. 保存 ----------
genus_pct.to_csv(OUT_DIR + 'mollin_genus_rel.csv')
log1p = np.log1p(genus_pct)
log1p.to_csv(OUT_DIR + 'mollin_genus_log1p.csv')
print('已保存: mollin_genus_rel.csv (相对丰度), mollin_genus_log1p.csv (log1p)')
print('log1p 矩阵描述:'); print(log1p.describe().loc[['mean', 'max']].round(3).head(8).to_string())
