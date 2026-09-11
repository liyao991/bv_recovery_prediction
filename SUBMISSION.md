# Submission Inventory — Line B BV Recovery Prediction Paper

## 数据可用性 / 代码与结果公开说明

完整代码与复现 CSV / 图均收录于本仓库 `results/` 与 `scripts/`，体积 ~2.5 MB（清理后）。

公开时间：2026-09-11  
清理对象：`refs/` 4 篇 PDF、内部中文草稿 v0–v1.8、框架审阅稿 v2.4–v3.2、审核/核实/统计答辩卡/白话流程图解、大中间文件（PRJNA302078/398590/592384 子目录含 ASV/phyloseq 对象）等共 63 文件 + 4 目录。

---

## 1. 复现主结果（Results §3 与表 / 图）

### 1.1 三层框架（§3.3）

| 论文元素                | 仓库路径                                                       |
| -------------------- | --------------------------------------------------------- |
| 图 1（三层 AUC 主图）       | `results/fig_three_tier_auc.png`                          |
| 表 3 摘要（A1–A4/B1–B4） | `results/results_c4_master_summary.{csv,md}`              |
| 复现 CSV（按 cross-cohort / CST-IV 子表） | `results/results_c4_master_summary.csv`                   |

### 1.2 三档族对比 / 五类模型族（§3.4 与补充 S2）

| 论文元素               | 仓库路径                                  |
| ------------------- | ------------------------------------ |
| 图 2（家族对比）          | `results/fig_feature_family.png`     |
| 图 3（richness 森林图）  | `results/fig_richness_forest.png`    |
| 图 4（标签耦合效应）        | `results/fig_label_coupling.png`     |
| 9 靶标 × 3 队列 × 结局热图 | `results/fig_c18_outcome_heatmap.png` |
| 跨家族/CLR/特征频率敏感性     | `results/results_c20_*.csv`          |
| LOCO 12 cell 矩阵     | `results/results_c19_loco_matrix.csv` |
| 跨队列 ω 一致性          | `results/results_c8_direction_consistency.csv` |

### 1.3 临床界标（§3.5）

| 论文元素                | 仓库路径                                                |
| -------------------- | -------------------------------------------------- |
| 图 5（D7 动态 + 短暂清除陷阱） | `results/fig_S5_d7_dynamic.png`                    |
| 复原 PPV 0.35（7/20）    | `results/results_dysbiosis_index.csv`             |
| 关键菌靶向预测             | `results/附表S7_临床关键菌靶向预测_数字表.md`                      |
| metafor 跨队列 pooled ρ | `results/c6_metafor_pooled.csv` + `fig_c6_meta_forest.png` |

### 1.4 主分析描述性数（正文 Methods/Results 频繁引用）

| 论文元素              | 仓库路径                                                |
| ------------------ | -------------------------------------------------- |
| 三队列基线 CST-IV 数   | `results/table_S2_three_tier_framework.{csv,docx}` |
| 最终 metrics（A1–A4 等） | `results/results_c9_final_metrics.csv`            |
| 41 属 / 65 属 / 83 属矩阵 | `results/results_c8_richness_targets.csv`         |
| 跨家族方差分解           | `results/results_c8_variance_decomp.csv`          |

---

## 2. 补充材料（Supplementary S1–S9）

### 2.1 S1–S5 数字与表

| 元素       | 仓库路径                                            |
| -------- | ---------------------------------------------- |
| Table S1 | `results/table_S1_dada2_params.{csv,docx}`      |
| Table S2 | `results/table_S2_three_tier_framework.{csv,docx}` |
| Table S3 | `results/table_S3_feature_list.{csv,docx}`      |
| Table S4 | `results/table_S4_resolution_annotation.{csv,docx}` |
| Table S5 | `results/table_S5_d7_dynamic.{csv,docx}`       |
| Table S6 | `results/table_S6_alpha_family_matrix.{csv,md}` |
| 附表 S7    | `results/附表S7_临床关键菌靶向预测_数字表.md`                  |
| 附表 S8/S9 | `results/线B论文_补充材料附表S8-S9_2026-09-03.{docx,html,md}` |

### 2.2 S1–S5 图

| 元素    | 仓库路径                          |
| ----- | ---------------------------- |
| Fig S1 | `results/fig_S1_quality_profiles.png` |
| Fig S2 | `results/fig_S2_three_tier.png`       |
| Fig S3 | `results/fig_S3_feature_distribution.png` |
| Fig S4 | `results/fig_S4_resolution_contrast.png` |
| Fig S5 | `results/fig_S5_d7_dynamic.png`         |

### 2.3 Code S1–S2

| 元素            | 仓库路径                                       |
| ------------- | ----------------------------------------- |
| Code S1 DADA2 | `results/code_S1_dada2_pipeline.R`        |
| Code S2 LASSO | `results/code_S2_lasso_pipeline.py`       |

---

## 3. 质量图（README 引用的 SRA fastq 质控核验）

| 队列                 | 仓库路径                                  |
| ------------------ | ------------------------------------ |
| PRJNA302078（Xiao）  | `results/PRJNA302078_quality.pdf` + `..._py.png` |
| PRJNA398590         | `results/PRJNA398590_quality_py.png`  |
| PRJNA592384         | `results/PRJNA592384_quality_py.png`  |
| PRJNA691964         | `results/PRJNA691964_quality_py.png`  |

---

## 4. 复现建模报告（Methods 引用）

- `results/modeling_report_20260829.md` — LASSO-LR / RF / XGB 五类模型族复现报告（含嵌套 CV 协议）

---

## 5. 清理掉的（仅在此文档记录，不上传）

下列文件均只用于作者内部记录，**不应**公开，已物理删除：

- `refs/` 4 篇 PDF（Xiao / Plummer / Gustin / Mollin 原文）
- `results/PRJNA302078`、`results/PRJNA398590`、`results/PRJNA592384` 大型中间对象（ASV/phyloseq RData 等）
- 全部内部中文草稿（v0, v1, v1.1–v1.8 JCM 排版）以及 v1.7 修改项证据链追踪表
- 全部框架审阅稿（v2.4 审阅批注版 → v3.2 审阅稿）
- 统计方法说明、白话答辩卡、白话流程图解（"求背书精简版"含机读检测高风险词）
- 审阅意见待改清单、审阅批注整理、审核意见清单（#C11/#C13/#C15）、核实结果（#C16–#C20）、文献引用内容核对表（#C33）等内部审计产物
- 框架修改待办清单、paper_framework_v1/v2 等迭代草稿
