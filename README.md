# BV 复发预测模型（Line B）

基于公开 16S rRNA 扩增子测序数据，构建细菌性阴道病（BV）治疗后菌群学康复的预测模型。核心结论为**阴性结果**：治疗前 16S 测序对 BV 治疗结局的预测价值有限。

## 数据来源（公开 SRA）

| 队列 | SRA 编号 | 文献 | 测序平台 |
|---|---|---|---|
| Xiao 2016 | PRJNA302078 | Xiao et al., *Sci Rep* 2016; 6:26674 | 454 GS Junior |
| Plummer / SToPBV | PRJNA592384 | Plummer et al., *J Infect Dis* 2023; jiad261 | MiSeq PE |
| StepUp-1 | PRJNA398590 | Plummer et al., *PLoS One* 2018; pone0272012 | MiSeq |
| Gustin 2022 | PRJNA691964 | Gustin et al., *Am J Obstet Gynecol* 2022 | MiSeq |

## 目录结构

- `scripts/` — 分析脚本（Python / R）：SRA 下载、DADA2 去噪、SILVA 分类、特征构建、建模与敏感性分析
- `results/` — 关键中间结果与图表

## 路径占位符说明

脚本中的数据路径统一使用以下占位符，复现时请替换为你本地的实际路径：

- `<DATA_DIR>` — 原始/特征数据根目录（SRA fastq、SILVA 参考库、特征矩阵）
- `<OUT_DIR>` — 计算中间输出目录（DADA2 产物、建模结果 CSV）
- `<REPO_DIR>` — 本仓库根目录
- `<USER_HOME>` — 用户主目录（仅构建工具脚本使用）

## 复现步骤

1. 用 SRA 编号从 ENA / NCBI 下载 fastq（`scripts/download_sra_cohorts.py`）
2. 运行 DADA2 流程生成 ASV 表与分类（`scripts/dada2_pipeline_template.R`）
3. 构建属/种级相对丰度与 α 多样性特征矩阵（`scripts/export_genus_matrix.R`、`scripts/export_alpha_family.R` 等）
4. 运行建模与敏感性分析（`scripts/train_main.py`、`scripts/four_cohort_matrix.py`、`scripts/c20_sensitivity.py` 等）
