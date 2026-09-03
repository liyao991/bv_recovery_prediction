# -*- coding: utf-8 -*-
"""
Mollin 2022 (PLoS ONE, Akins-Sobel cohort) — 统一 CST 结局标签构建
====================================================
目标：让 Mollin 套用与 Xiao/SToPBV 完全相同的结局定义
    CST4 = Lactobacillus 相对丰度 < 50%
    outcome = 1 (康复) 当 V2 非 CST4；0 (未康复/CST4) 当 V2 为 CST4

数据源（已下载 mollin2022_S1.xlsx）：
    S1 Table: 元数据（#NAME/Class/Patient/Visit/Nugent/Amsel 等）
    S2 Table: OTU -> Taxonomy 映射（5011 OTU，RDP+NCBI）
    S4 Table: 865 species 计数表（样本 × OTU 计数，含 CLASS 行）

输出：
    <DATA_DIR>/meta/mollin_cst_labels.csv — 样本级 CST 标签（全部样本）
    <DATA_DIR>/meta/mollin_pairs.csv — V1->V2 配对结局
"""
import pandas as pd
import numpy as np
import re

XLSX = '<DATA_DIR>/meta/mollin2022_S1.xlsx'
OUT_DIR = '<DATA_DIR>/meta/'

# ---------- 1. 解析 S2: Otu -> Taxonomy ----------
s2 = pd.read_excel(XLSX, sheet_name='S2Table')
s2.columns = ['hdr', 'tax', 'seq']
s2 = s2.iloc[1:].reset_index(drop=True)  # 去掉表头行
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
    # 去掉末尾 Otu 编号段
    if parts and parts[-1].startswith('Otu'):
        parts = parts[:-1]
    if not parts:
        return 'NA'
    # 属 = 倒数第二段（倒数第一为种）；若只有 5 层则倒数第二为属
    if len(parts) >= 2:
        genus = parts[-2].strip()
        # 若 genus 段其实带 *（种级无属可用），退回倒数第一
        if genus.endswith('*') or genus == '':
            genus = parts[-1].strip()
    else:
        genus = parts[0].strip()
    # 清理
    genus = genus.replace('*', '').strip()
    return genus


otu2genus = {o: tax_to_genus(t) for o, t in otu2tax.items()}
from collections import Counter
gc = Counter(otu2genus.values())
print('映射后属分布 top 15:', gc.most_common(15))
lb_otus = [o for o, g in otu2genus.items() if g == 'Lactobacillus']
print(f'Lactobacillus OTU 数: {len(lb_otus)}')

# ---------- 2. 读 S4: 865 species 计数表 ----------
# row0: 表标题, row1: 样本名, row2: CLASS, row3+: OTU 行
s4 = pd.read_excel(XLSX, sheet_name='S4 Table', header=None)
samples = [str(v).strip() for v in s4.iloc[1, 1:].tolist()]
classes = [str(v).strip() for v in s4.iloc[2, 1:].tolist()]
otu_rows = s4.iloc[3:]
print(f'\nS4 样本列: {len(samples)}, OTU 行: {len(otu_rows)}')
print('CLASS 分布:', pd.Series(classes).value_counts().to_dict())

# 构建 样本 x OTU 计数矩阵
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
genus_pct = genus_level.div(genus_level.sum(axis=0), axis=1)
lb_pct = genus_pct.loc['Lactobacillus'] if 'Lactobacillus' in genus_pct.index else pd.Series(0, index=genus_pct.columns)

# ---------- 4. CST 标签 ----------
cst = pd.DataFrame({
    'sample': samples,
    'class': classes,
    'lactobacillus_pct': lb_pct.values,
})
cst['cst4'] = (cst['lactobacillus_pct'] < 0.50).astype(int)

# 匹配 S1 元数据（Patient/Visit）——列名在第 1 行，数据从第 2 行起
s1 = pd.read_excel(XLSX, sheet_name='S1 Table', header=None)
s1.columns = ['c' + str(i) for i in range(s1.shape[1])]
s1 = s1.iloc[2:].reset_index(drop=True)
s1 = s1.rename(columns={'c0': '#NAME', 'c1': 'Class', 'c2': 'Patient', 'c3': 'Visit'})
s1 = s1.dropna(subset=['Patient']).drop_duplicates(subset=['#NAME'], keep='first').reset_index(drop=True)
s1['#NAME'] = s1['#NAME'].astype(str).str.strip()
s1['Patient'] = s1['Patient'].astype(str).str.strip()
s1['Visit'] = s1['Visit'].astype(str).str.strip()
s1['Class'] = s1['Class'].astype(str).str.strip()
meta_map = s1.set_index('#NAME')[['Patient', 'Visit', 'Class']].to_dict('index')
cst['patient'] = cst['sample'].map(lambda s: meta_map.get(s, {}).get('Patient', np.nan))
cst['visit'] = cst['sample'].map(lambda s: meta_map.get(s, {}).get('Visit', np.nan))
cst['meta_class'] = cst['sample'].map(lambda s: meta_map.get(s, {}).get('Class', np.nan))
cst = cst.dropna(subset=['patient'])

# 去掉重复列名样本（A15V2.1）
cst = cst[~cst['sample'].str.contains(r'\.', na=False)].reset_index(drop=True)

print('\n=== Mollin CST 标签（样本级）===')
print(f'样本数: {len(cst)}, 患者数: {cst["patient"].nunique()}')
print('Visit 分布:', cst['visit'].value_counts().to_dict())
print('\nCST4 (LB<50%) 分布:')
print(cst.groupby('visit')['cst4'].value_counts().unstack(fill_value=0))

# ---------- 5. V1 -> V2 配对结局 ----------
pairs = []
for pid, g in cst.groupby('patient'):
    v1 = g[g['visit'] == 'V1']
    v2 = g[g['visit'] == 'V2']
    if v1.empty or v2.empty:
        continue
    r1, r2 = v1.iloc[0], v2.iloc[0]
    pairs.append({
        'patient': pid,
        'base_sample': r1['sample'], 'outcome_sample': r2['sample'],
        'base_lb_pct': round(r1['lactobacillus_pct'], 4),
        'outcome_lb_pct': round(r2['lactobacillus_pct'], 4),
        'base_cst4': int(r1['cst4']), 'outcome_cst4': int(r2['cst4']),
        'outcome': 1 if r2['cst4'] == 0 else 0,  # 康复 = V2 非 CST4
        'meta_class_base': r1['meta_class'], 'meta_class_outcome': r2['meta_class'],
    })
pairs = pd.DataFrame(pairs)
print(f'\n=== V1->V2 配对: {len(pairs)} 例 ===')
print('outcome (1=康复/非CST4, 0=未康复/CST4) 分布:', pairs['outcome'].value_counts().to_dict())
print('base_cst4 分布:', pairs['base_cst4'].value_counts().to_dict())
print('\n基线状态 × 结局交叉:')
print(pd.crosstab(pairs['base_cst4'], pairs['outcome'], margins=True))

# ---------- 6. 与 S1 硬结局对比 ----------
print('\n=== CST 结局 vs S1 硬结局（meta_class）===')
print('基线 CST4=1 (菌群失调) 的硬结局分布:')
core = pairs[pairs['base_cst4'] == 1]
print(core['meta_class_base'].value_counts().to_dict())
print('\n配对表:')
print(core[['patient', 'base_lb_pct', 'outcome_lb_pct', 'outcome', 'meta_class_base']].head(10).to_string())

# ---------- 7. 保存 ----------
cst.to_csv(OUT_DIR + 'mollin_cst_labels.csv', index=False)
pairs.to_csv(OUT_DIR + 'mollin_pairs.csv', index=False)
print(f'\n已保存: mollin_cst_labels.csv ({cst.shape}), mollin_pairs.csv ({pairs.shape})')
