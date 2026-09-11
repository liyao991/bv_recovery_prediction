# 附表 S6 · α 多样性家族 × 分辨率矩阵（Gustin 队列）

> **目的**：系统评估 α 多样性指标（单指标 vs. 家族组合）在不同分类学分辨率（ASV / 属 / 种）下对 BV 治疗结局的预测性能。
> **数据**：PRJNA691964（Gustin 2022, n=33 对）；**结局口径**：随访 Lactobacillus ≥0.5（主口径）。
> **建模**：LASSO-LR（liblinear, l1_ratio=1.0, class_weight='balanced', max_iter=5000）+ 5 折分层 CV × 20 随机种子。

---

## 一、Gustin 主口径矩阵（12 项）

| 分辨率 | 家族 | 构成 | AUC mean | AUC SD | n |
|---|---|---|---|---|---|
| asv | richness | observed + chao1 | 0.727 | 0.210 | 33 |
| asv | evenness | pielou | 0.402 | 0.159 | 33 |
| asv | entropy | shannon + simpson + invsimpson | 0.608 | 0.223 | 33 |
| asv | full | 全部 6 项指标 | 0.691 | 0.211 | 33 |
| genus | richness | observed + chao1 | 0.681 | 0.219 | 33 |
| genus | evenness | pielou | 0.468 | 0.185 | 33 |
| genus | entropy | shannon + simpson + invsimpson | 0.577 | 0.243 | 33 |
| genus | full | 全部 6 项指标 | 0.660 | 0.225 | 33 |
| species | richness | observed + chao1 | **0.776** | 0.164 | 33 |
| species | evenness | pielou | 0.347 | 0.180 | 33 |
| species | entropy | shannon + simpson + invsimpson | 0.504 | 0.231 | 33 |
| species | full | 全部 6 项指标 | 0.747 | 0.184 | 33 |

**要点**：
- **richness 家族**在种级分辨率最高（AUC=0.776±0.164），与主正文中「种级 richness」结论一致；
- 分辨率效应**非单调**：ASV richness 0.727 > 属 richness 0.681，种级最高 0.776；
- evenness/entropy 家族整体较弱，不构成独立预测信号。

---

## 二、「120 项」精确定义

附表 S6 正文仅列出 Gustin 主口径 **12 项**。论文 Methods 中提到的 **120 项** 是以下完整探索矩阵的单元格数：

| 维度 | 水平 | 数量 |
|---|---|---|
| 队列 | PRJNA302078 (Xiao) / PRJNA592384 (Plummer①) / PRJNA691964 (Gustin) / PRJNA398590 (Plummer②) | 4 |
| 分类学分辨率 | ASV / 属 / 种 | 3 |
| 单指标 α 指标 | observed / chao1 / shannon / simpson / invsimpson / pielou | 6 |
| 家族组合 | richness / evenness / entropy / full | 4 |

按每个单元格计 1 项：
- 单指标矩阵 = 4 队列 × 3 分辨率 × 6 指标 = **72 项**
- 家族组合矩阵 = 4 队列 × 3 分辨率 × 4 家族 = **48 项**
- **合计 = 72 + 48 = 120 项**

**可执行状态**：
- 本表给出 Gustin 的 12 项家族组合结果（单指标隐含在家族组合构成中）；
- Xiao / Plummer 队列虽有结局标签，但 α 数据为随访时点（与结局强耦合），不适合作为前瞻性预测证据，未纳入 S6；
- PRJNA398590（Plummer②）已被正文排除（伴侣干预/含男性样本），且无可用结局标签，不具可估算性。

---

## 三、脚本与复现

- 主计算：`line_B/scripts/c16_alpha_matrix_main.py`
- 数据源：`D:/lineB_data/feat/PRJNA691964_alpha_family.csv`
- 配对：`D:/lineB_data/meta/gustin_run_map.csv`
- 输出：`table_S6_alpha_family_matrix.csv`（本文件同路径）
