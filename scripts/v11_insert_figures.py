"""线 B 论文中文草稿 v1 → v1.1 升级
变更：
1. 3.2 节末插入图 2（标签耦合对照）
2. 原图 2（richness 森林图）改编号为图 3
3. 3.4 节末插入图 4（特征家族 α vs 属级对照）
4. 3.6 节末插入图 5（靶标 ρ × 队列热图）
5. 同步改图表清单
"""
from pathlib import Path

p = Path('<REPO_DIR>/results/线B论文中文草稿_v1_2026-09-01.md')
text = p.read_text(encoding='utf-8')

# ========== 1. 原图 2 → 图 3 ==========
old3 = (
    '![图 2 基线 richness 与治疗结局的 Spearman ρ 森林图（三队列五亚组，左=原始 / 右=rarefaction 校正）。Gustin ρ 由 −0.45→−0.39 仍显著（红），Xiao/Plummer 方向相反且不显著（蓝）。Mollin 因 SI 表重建无种级 richness 而未纳入。](fig_richness_forest.png)'
)
new3 = (
    '![图 3 基线 richness 与治疗结局的 Spearman ρ 森林图（三队列五亚组，左=原始 / 右=rarefaction 校正）。Gustin ρ 由 −0.45→−0.39 仍显著（红），Xiao/Plummer 方向相反且不显著（蓝）。Mollin 因 SI 表重建无种级 richness 而未纳入。](fig_richness_forest.png)'
)
assert text.count(old3) == 1, f'原图 2 匹配 {text.count(old3)} 次，应为 1'
text = text.replace(old3, new3)

# 图 2 配套说明 → 图 3 配套说明
old3b = '**图 2 配套说明**'
new3b = '**图 3 配套说明**'
assert text.count(old3b) == 1, f'图2 配套说明匹配 {text.count(old3b)} 次'
text = text.replace(old3b, new3b)

# ========== 2. 3.2 节末插入图 2（标签耦合）==========
# 锚点：3.2 节最后一句"提示基线乳杆菌高者「必康复」的定义层面机械关联贡献了大部分预测力。"
old32_anchor = '全人群去除 Lactobacillus 特征后 AUC 由 0.652 降至 0.567（显著下降），提示基线乳杆菌高者「必康复」的定义层面机械关联贡献了大部分预测力。'
new32 = old32_anchor + '''

![图 2 标签耦合对照：含/去 Lactobacillus 属对预测精度的影响（Xiao 队列）。全人群去 LB 后 AUC 显著下降（0.652→0.567），亚组内几乎不变（0.596→0.584）——表明全人群 0.6 信号很大程度来自 Lactobacillus 属对「康复」标签的强耦合；亚组内已剔除该耦合故不下降。](fig_label_coupling.png)

**图 2 配套说明**：误差棒为 95%CI（去 LB 后区间按 holdout SD×1.96 估）；亚组内"无下降"说明负结果干净、不依赖单一特征；全人群 0.6 信号经此对照被部分归因于定义层面的机械关联。'''
assert text.count(old32_anchor) == 1, f'3.2 锚点匹配 {text.count(old32_anchor)} 次'
text = text.replace(old32_anchor, new32)

# ========== 3. 3.4 节末插入图 4（特征家族对照）==========
old34_anchor = '（「短暂清除陷阱」）。'
new34 = old34_anchor + '''

![图 4 特征家族对照：属级全集 vs α-only richness（多队列对照）。左图显示三队列属级全集 AUC 在 0.44–0.63 之间（≈随机），而 Gustin α-only richness 单独达 0.689 卡边缘——证明「α 强、属级弱」非矛盾而是不同家族的不同判别能力，但前者仅在 Gustin 单队列成立、未跨队列复现。](fig_feature_family.png)

**图 4 配套说明**：(a) 三队列内部 CV 属级全集与 Gustin α-only richness AUC 对比；(b) richness 单变量跨队列 Spearman ρ 方向对比——Gustin 显著（红）、Xiao/Plummer 方向相反不显著（蓝/浅蓝）。两图共同说明"α-only richness 是单队列的局部胜利、非跨队列可泛化信号"。'''
assert text.count(old34_anchor) == 1, f'3.4 锚点匹配 {text.count(old34_anchor)} 次'
text = text.replace(old34_anchor, new34)

# ========== 4. 3.6 节末插入图 5（靶标热图）==========
old36_anchor = '**可预测性的瓶颈不在特征选择（先验靶标 vs 全集），而在结局本身的生物学决定因素**。'
new36 = old36_anchor + '''

![图 5 9 个临床关键菌属（4 队列）Spearman ρ 热图。Xiao/Plummer/Gustin 三队列方向一致性差：仅 Atopobium（蓝）、Anaerobes 聚合、Mobiluncus 方向一致（其中 Atopobium/Mobiluncus 在 ≥2 队列达 p<0.05），其余 6 个靶标方向反复或弱信号；Mollin 因 SI 表重建、种级粒度不可比，整体偏中性。Mollin 列入作"管线异质性"参考。](fig_c8_targets_heatmap.png)

**图 5 配套说明**：红色=正相关（高菌→康复），蓝色=负相关（高菌→未康复），* = p<0.05；色阶对称于 ρ=0；视觉聚类：Gardnerella/Atopobium/Prevotella 在不同队列方向不一致。'''
assert text.count(old36_anchor) == 1, f'3.6 锚点匹配 {text.count(old36_anchor)} 次'
text = text.replace(old36_anchor, new36)

# ========== 5. 改图表清单 ==========
old_list = '''**图 1**：三层信号框架柱状图（横断面诊断 >0.9 / 全人群恢复 ~0.53-0.63 / CST4 亚组 ≈0.5；橙色虚线标随机 0.5）
**图 2**：richness 森林图（三队列五亚组 Spearman ρ；左 raw / 右 rarefaction 校正；Gustin 显著、Xiao/Plummer 方向相反不显著）'''
new_list = '''**图 1**：三层信号框架柱状图（横断面诊断 >0.9 / 全人群恢复 ~0.53-0.63 / CST4 亚组 ≈0.5；橙色虚线标随机 0.5）— 见 Results 3.1
**图 2**：标签耦合对照柱状图（Xiao 队列；全人群 vs 亚组 × 含/去 Lactobacillus 双因素）— 见 Results 3.2
**图 3**：richness 森林图（三队列五亚组 Spearman ρ；左 raw / 右 rarefaction 校正；Gustin 显著、Xiao/Plummer 方向相反不显著）— 见 Results 3.5
**图 4**：特征家族 α vs 属级对照（a：三队列属级全集 vs Gustin α-only AUC；b：richness 跨队列 ρ 方向对比）— 见 Results 3.4
**图 5**：9 临床关键菌属 × 4 队列 Spearman ρ 热图（仅 Atopobium/Anaerobes 聚合/Mobiluncus 方向一致）— 见 Results 3.6'''
assert text.count(old_list) == 1, f'图表清单匹配 {text.count(old_list)} 次'
text = text.replace(old_list, new_list)

# ========== 6. 改标题版本号 ==========
text = text.replace('> **中文稿 v1 · 2026-09-01 · 依据框架 v3.2 撰写**',
                    '> **中文稿 v1.1 · 2026-09-01 · 依据框架 v3.2 撰写；图 2/4/5 新增嵌入**')

# 写回
p.write_text(text, encoding='utf-8')
print('OK: 5 处插入 + 1 处版本号修改完成')
print(f'文件总长 {len(text)} 字符')
