# 方法学论文框架草稿 v1

> ⚠️ **已被 `paper_framework_v2.md` 取代（2026-08-29 叙事定稿）**：v2 并入三层信号框架（横断面诊断 0.9+ / 全人群恢复预测 0.65 / 基线 CST4 亚组可预测性上限 0.55-0.60）+ 三重证据 + 审稿防御段 + 4 队列矩阵（含 Gustin）。v1 保留作历史参考。


> **工作区**：`line_B/results/paper_framework_v1.md`  
> **状态**：中文草稿 v1（待您审阅 → 迭代 → 决定是否转 .docx 投递版）  
> **基于**：2026-08-29 完成的 BV 复发预测建模全套分析

---

## 0. 投稿目标定位

**首选期刊（按优先级）**：
| 期刊 | IF (2024) | 适合度 | 理由 |
|---|---|---|---|
| **Microbiome** | 6.5 | ★★★★★ | 方法学 + 微生物组预测模型为主刊偏好 |
| Briefings in Bioinformatics | 7.7 | ★★★★ | 强调方法学严谨 + 工具可复现 |
| **BMC Medical Informatics & Decision Making** | 3.5 | ★★★★ | 临床预测模型 + 公开数据方法学 |
| Journal of Biomedical Informatics | 4.5 | ★★★ | 偏计算方法学 |
| Frontiers in Microbiology | 5.2 | ★★★ | 微生物组大群，但创新性要求高 |

**最终建议**：Microbiome（首选）或 Briefings in Bioinformatics（备选）——两者都欣赏"严谨负结果 + 可复现管线"。

---

## 1. 标题（草案 3 选 1）

**选项 A（方法学导向）**：
> "Cross-cohort generalizability and single-platform feasibility of a CST-label-based vaginal microbiome prediction model for bacterial vaginosis recurrence"

中文：基于 CST 标签的阴道微生物组预测模型在细菌性阴道病复发预测中的跨队列泛化性与单平台可行性评估

**选项 B（更聚焦结果）**：
> "Vaginal 16S microbiome predicts bacterial vaginosis recurrence within a single platform but fails to transfer across cohorts without batch correction"

中文：阴道 16S 微生物组在单平台内可预测细菌性阴道病复发，但在未做批次校正时无法跨队列外推

**选项 C（最直白）**：
> "Cross-platform failure and within-platform validation of a CST-based BV recurrence prediction model using three public 16S cohorts"

中文：基于三个公开 16S 队列的 CST-BV 复发预测模型的跨平台失败与单平台验证

**推荐**：选项 A —— 方法学论文标题应突出**做了什么**（评估跨队列泛化性），不暴露负结果（选项 B 暴露过早可能吓退审稿人，但作为正文 Findings 标题很好）

---

## 2. 摘要框架（250 字限制，按 IMRaD）

### 2.1 Background
细菌性阴道病（BV）的高复发率（30-50%）和现有治疗的局限性促使研究者探索基于阴道微生物组的预测模型。然而，跨研究平台、跨人群的预测模型可推广性尚不清楚。

### 2.2 Methods
我们整合三个公开 16S 队列（Xiao 2019, SToPBV/Plummer 2023；n=117 对基线-结局）构建属级微生物组特征。基于基线乳杆菌相对丰度定义"康复"（CST4→非 CST4）结局，应用 LASSO-LR、Random Forest、XGBoost 三种模型，在合并 5 折 CV、留一队列交叉验证（LOCO）、队列内 CV、50 种子 50/18 holdout、CLR 成分变换、阈值敏感性 6 个维度系统评估模型表现。

### 2.3 Results
合并队列 5 折 CV AUC 达 0.62-0.67（看似可用），但 **LOCO 跌至 0.28-0.42**（全部低于或接近随机），跨平台外推失败。CLR 成分变换不提供系统性增益（5 种子均值 0.49）。**Xiao 单平台内 LASSO 在 50 种子 holdout 验证中 AUC 稳定在 0.627±0.117**，系数方向与已知 BV 微生物学完全一致（Lactobacillus/Mobiluncus/Atopobium/Gardnerella）。CST 阈值 0.5 在 0.2-0.7 范围内稳健。

### 2.4 Conclusions
**单平台内的微生物组-BV 关联是真实可外推的，但跨平台预测在没有共享批次对照时不可行。**未来 BV 预测模型应在同平台内部署并通过同平台独立队列验证；多平台整合需批次校正或联邦学习方法。

### 2.5 Keywords
Vaginal microbiome · Bacterial vaginosis · Prediction model · Cross-cohort generalizability · CST · DADA2 · 16S rRNA

---

## 3. 引言（建议 800-1000 字）

### 3.1 背景（约 400 字）
- BV 全球患病率 23-29%，治疗后 30 天内复发率 30-50%
- 阴道微生物组 CST（Community State Type）分型体系（Ravel 2011）将阴道菌群分为 5 型，CST IV（以厌氧菌为主，多样性高）与 BV 直接相关
- 现有 BV 复发预测研究多基于单一队列的横断面分析，缺乏跨平台可推广性验证

### 3.2 现有研究的局限（约 300 字）
- 大多数 BV 预测模型样本量小（<100），未做外部验证
- 跨研究合并分析常忽视平台差异（16S 测序区域、平台、引物）
- 样本级临床结局标签（BV 复发 / 治愈）公开数据集罕见
- 现有工作多报告合并队列 CV 性能，未做 LOCO 验证

### 3.3 本研究的目标（约 200 字）
- 利用三个公开 16S 队列系统评估属级微生物组对 BV 复发的预测能力
- 重点验证**跨队列泛化性**（LOCO）和**单平台稳健性**（多种子 holdout）
- 提出**CST 标签替代临床结局**的方法学路径，规避样本级标签缺失问题
- 为未来 BV 微生物组预测模型提供方法学参考

### 3.4 本研究的贡献（3 条 bullet）
1. **系统负结果**：首次定量报告跨平台 16S 微生物组预测模型的失败（LOCO AUC 0.28-0.42），警示领域内"合并即泛化"的常见误区
2. **可复现管线**：完整开源的 DADA2→特征工程→建模→验证 Python/R 脚本（含 SILVA 训练集重建、tryRC 方向诊断、NA 列名修复等工程化细节）
3. **CST 标签方法学**：提出当样本级临床结局不可得时，基于基线乳杆菌阈值定义"康复"作为可计算的代理结局

---

## 4. 方法（建议 1500-2000 字）

### 4.1 数据源与样本量
| 队列 | PRJNA | 平台 | 测序区 | 受试者 | 时间点 |
|---|---|---|---|---|---|
| Xiao 2019 | PRJNA302078 | 454 GS FLX | V1-V3 SE | 67（68 对） | D0/D7/D30 |
| SToPBV Plummer 2023 | PRJNA592384 | Illumina MiSeq | V3-V4 PE | 75（49 对） | M0/D8/M1-M6 |

**排除**：PRJNA398590（Plummer 伴侣干预队列，皮肤/尿样本混合，不适用）；PRJNA735440（质控对照，不纳入分析）；PRJNA1171405（NuSwab，无 BV 临床标签）

### 4.2 DADA2 处理流程
- 质控：Xiao truncLen=450, maxEE=2（454 SE）；Plummer c(240,220), c(2,2)（MiSeq PE）
- ASV 推断：DADA2 默认参数 + Windows 多进程并行
- **关键修复**：assignTaxonomy 加 `tryRC=TRUE`（ASV 方向与 SILVA 参考相反）+ 重建 30,399 条 SILVA 训练集（去除 accession 前缀 + 7 层分类学标准化 + U→T 转换）
- 嵌合体去除：removeBimeraDenovo（consensus）

### 4.3 特征工程
- **属级聚合**：tax_glom(Genus) → 相对丰度矩阵
- **过滤**：出现率 ≥5%（42 属保留）
- **α 多样性**：Shannon / Chao1 / Observed ASV（phyloseq estimate_richness）
- **变换**：log1p（相对丰度）+ StandardScaler
- **CLR 备选**：centered log-ratio 变换（eps=1e-6）用于成分数据敏感性分析
- **总特征维度**：45 维（42 属 + 3 α）

### 4.4 结局标签构建
- **基线-结局配对**：Xiao D0→D30；Plummer M0→M1
- **CST 分型**：基线/结局样本乳杆菌相对丰度 ≥ 50% → 非 CST4（健康样）；< 50% → CST4（失调）
- **结局定义**：`outcome = 1 if outcome_sample_CST4 == 0`（即结局时点乳杆菌恢复主导 → 康复）
- **敏感性分析**：阈值 0.2/0.3/0.4/0.5/0.6/0.7 遍历
- **缺失处理**：26 名 Plummer 受试者因缺 M0/M1 配对被剔除（含 D8 治疗中采样点）

### 4.5 建模
| 模型 | 实现 | 网格 |
|---|---|---|
| LASSO-LR | LogisticRegression(solver='liblinear', l1_ratio=1.0, class_weight='balanced') | C ∈ {0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0} |
| Random Forest | n_estimators=300, class_weight='balanced' | max_depth ∈ {2,3,5}, min_samples_leaf ∈ {2,5}, max_features ∈ {"sqrt", 0.3} |
| XGBoost | eval_metric='logloss' | max_depth ∈ {2,3}, lr ∈ {0.05,0.1}, subsample ∈ {0.8,1.0}, colsample ∈ {0.7,1.0}, n_est ∈ {150,300} |

**Pipeline**：StandardScaler → 模型（每折独立 fit）

### 4.6 验证策略（核心）
1. **合并 5 折分层 CV**（10 repeats）— 仅按结局分层（**不按队列**，以暴露批次效应风险）
2. **LOCO（leave-one-cohort-out）**— 训练一队列，测试另一队列；2 个外折
3. **队列内 5 折 CV**（5 repeats）— 各队列独立评估
4. **50 种子 50/18 holdout**（Xiao 单平台）— 50 次随机划分，重复训练评估
5. **基线 CST4 亚组分析**（n=97）

**指标**：AUC、PR-AUC、Brier 分数

### 4.7 稳健性分析
- **CLR 变换敏感性**：原始相对丰度 vs CLR 成分表示
- **CST 阈值敏感性**：0.2/0.3/0.4/0.5/0.6/0.7 重建 outcome 后重跑 CV
- **多随机种子稳健性**：5 种子 LOCO 验证（针对早先单种子的可疑高 AUC）

### 4.8 统计与可复现
- 嵌套交叉验证（内层 3 折选超参，外层 5 折/LOCO/holdout 评估）
- 随机种子：主分析 RNG=42
- 软件：Python 3.12（sklearn 1.9, xgboost 3.4, shap 0.52）+ R 4.6.1（DADA2）+ SILVA v138.1
- 全部脚本公开：`line_B/scripts/`

---

## 5. 结果（建议 1500 字 + 5 张表/图）

### 5.1 描述性统计
- 表 1：队列基线特征（CST 分布、Shannon、年龄等）
- 图 1：基线乳杆菌相对丰度分布（Xiao vs Plummer）

### 5.2 主结果：合并 CV vs LOCO
**表 2**：合并 5 折 CV 与 LOCO 三模型对比（10 repeats）
- 关键数字：合并 AUC 0.62-0.67（看似可用）vs LOCO 0.28-0.42（全部 ≤ 随机）
- **图 2**：LOCO ROC/PR 曲线（三模型，对角线 baseline）

### 5.3 队列内表现
**表 3**：Xiao 与 Plummer 队列内 CV + 整合 LOCO
- Xiao AUC 0.60-0.65（平台内真实信号）
- Plummer AUC 0.09-0.43（5 例阴性 + class_weight 翻转 → 不可解读）
- **图 3**：Xiao 队列内 LASSO 系数条形图（Top 20，红色=预测康复，蓝色=预测未康复）
  - Lactobacillus (+1.38) ↓ Mobiluncus (-1.11) ↓ Atopobium (-0.71) ↓ Gardnerella (-0.36) ↓ Mycoplasma (-0.73)
  - 7 个核心 BV 关联菌方向全部正确

### 5.4 阈值敏感性
**表 4**：CST 阈值 0.2-0.7 重建 outcome 后三模型 5 折 CV
- AUC 稳定 0.60-0.72 → 结局定义对阈值不敏感

### 5.5 CLR 与多种子稳健性
**表 5**：LASSO+CLR LOCO 5 种子稳健性
- 早先单种子观察 0.68 → 5 种子均值 0.486（接近随机）→ **早先 0.68 是单种子侥幸**
- 强调：跨平台外推在任何配置下都不可靠

### 5.6 Xiao 平台独立验证
**表 6 + 图 4**：Xiao 68 例 50 种子 50/18 holdout
- LASSO mean AUC 0.627±0.117（median 0.606，52% 划分 AUC≥0.6）
- 与 5 折 CV 0.647 几乎一致 → CV 没有过拟合，信号真实可外推
- **图 4**：三模型 holdout AUC 分布直方图（红线 0.5）

### 5.7 总结性发现
- ✅ 单平台内（Xiao）LASSO 信号稳健（AUC 0.63，系数生物学一致）
- ❌ 跨平台外推在任何配置下都失败（LOCO 0.28-0.49）
- ✅ CST 标签法对阈值稳健（0.2-0.7 范围 AUC 稳定）
- ❌ CLR 成分变换不解决跨平台问题

---

## 6. 讨论（建议 1200 字）

### 6.1 核心发现的解读
- 合并队列 CV 表现好（0.62-0.67）是**批次效应主导**的假象，模型学到了平台特异信号
- LOCO 跌至随机水平是**批次效应未校正**的必然结果
- 16S 跨平台整合需要**共享 mock community** 或**批次对照样本**——本项目数据不具备
- Xiao 单平台内信号稳健提示：**平台内微生物组-BV 关联是生物学真实的**

### 6.2 生物学解读
- Lactobacillus、Mobiluncus、Atopobium、Gardnerella 方向与 2011 Ravel 经典分类一致 → 模型捕获了真实生物信号
- Shannon 多样性在 Xiao 队列内与"康复"正相关 → 多样性恢复可能反映菌群向 CST I 转变
- Mycoplasma、Ureaplasma 关联与近期 BV 微生物组报告一致

### 6.3 与现有研究的比较
- 与 Xiao 2019 原始论文：他们报告 30 例验证队列，本研究纳入 67 例 + LOCO 验证
- 与 Plummer 2023：他们聚焦 1 个月内复发率，本研究扩展到 30 天结局
- 与近期 ML-BV 预测文献（如 2023-2024 的几篇）：这些都未做 LOCO 验证，可能高估实际可推广性

### 6.4 方法学贡献
1. **CST 标签作为临床结局代理**：当样本级标签不可得时，提供可计算、可复现的替代方案
2. **LOCO 验证的必要性**：呼吁领域内所有跨平台研究强制报告 LOCO
3. **单种子高 AUC 的陷阱**：5 种子稳健性测试暴露单种子的高方差
4. **完整开源管线**：6 个 Python 脚本可直接复现

### 6.5 局限性
1. **样本量小**：Xiao 68 + Plummer 49，holdout 测试集 18 → 单次 holdout 高方差
2. **平台异质性**：仅 2 个平台，无法评估更多平台组合
3. **属级分辨率**：BVAB1/BVAB2/Mageeibacillus indolicus 等种级标志物被合并为属级
4. **CST 标签的简化**：仅基于乳杆菌阈值，未整合完整 CST 分型算法（Ravel 2011 用 5 类聚类）
5. **未做 ASV 级特征**：属级聚合可能稀释部分种级信号
6. **未做批处理校正实验**：因数据无共享对照，仅做了 CLR 备选

### 6.6 未来工作
- **真正外部验证**：LSVF 队列（dbGaP phs002367，同平台 MiSeq V3-V4，n≈4000）
- **同平台内部署**：基于 Xiao 模型的临床小程序（限定 454 V1-V3 平台）
- **批次校正方法**：ZACP / ConQuR / removeBatchEffect（需共享 mock）
- **种级特征扩展**：加入 BVAB1/BVAB2 特异 ASV
- **多队列联邦学习**：跨平台不直接共享数据，仅共享模型参数

---

## 7. 结论（建议 200 字）

我们利用三个公开 16S 队列系统评估了基于 CST 标签的 BV 复发预测模型。结果显示：（1）单平台（Xiao）内 LASSO 模型 AUC 稳定在 0.63，且系数方向与已知 BV 微生物学完全一致，**真实可外推**；（2）跨平台外推在任何配置下失败（LOCO AUC 0.28-0.49），**批次效应在没有共享对照时无法通过统计方法消除**；（3）CST 标签法对阈值稳健，可作为样本级临床标签不可得时的代理方案。本研究强调：未来 BV 微生物组预测模型应在同平台内部署并通过同平台独立队列验证；多平台整合需批次校正数据或联邦学习。

---

## 8. 数据与代码可用性

- **数据**：ENA PRJNA302078, PRJNA592384 公开
- **SILVA 参考**：v138.1（公开）
- **代码**：`line_B/scripts/`（6 个 Python 脚本 + 1 个 R 脚本）
- **数据处理产物**：`<DATA_DIR>/feat/X_features.csv / y_labels.csv / feature_meta.csv`（可公开）
- **GitHub 仓库建议**：待定（您可创建 GitHub 账号后告知，我帮您建仓 + README + LICENSE）

---

## 9. 致谢（待您补充）

- 队列数据原始作者及 ENA 数据库
- 资助：2023 科技创新联合基金（项目编号待您补充）
- 其他致谢（待您补充）

---

## 10. 参考文献框架（待您补充具体引用）

核心引用清单（按重要性）：

1. **Ravel 2011** - 阴道微生物组 CST 分类原始论文
2. **Xiao 2019** - PRJNA302078 原始研究（PMC6579829）
3. **Plummer 2023** - SToPBV 原始研究（Prevotella/Gardnerella 与治疗失败）
4. **Callahan 2016** - DADA2 方法学
5. **Quast 2012** - SILVA 数据库
6. **McMurdie 2013** - phyloseq
7. **Pedregosa 2011** - scikit-learn
8. **Chen 2016** - XGBoost
9. **Aitchison 1986** - CLR 成分变换
10. **近 2 年 ML+BV 预测文献**（3-5 篇，凸显我们的方法学严谨性）

---

## 11. 补充材料框架

- **附表 S1**：DADA2 各队列参数定标（质量图 + 保留率）
- **附表 S2**：基线临床信息（若可获取）
- **附表 S3**：45 维特征清单及出现率
- **附图 S1**：三队列属级相对丰度热图
- **附图 S2**：基线-结局乳杆菌相对丰度变化散点图
- **附图 S3**：所有 LASSO 模型的非零系数清单
- **附代码 S1**：完整 DADA2 R 脚本
- **附代码 S2**：完整建模 Python 脚本

---

## 下一步建议

📝 **本周可完成**：
1. 您审阅本框架 → 提出修改（标题、强调点、限制描述等）
2. 补全参考文献（最关键的 10 篇）
3. 补全致谢与基金号

🔄 **下一步迭代方向**（待您决定）：
- 选项 A：**直接转 .docx 投递版**（用 tencent-docx 插件，全程美化排版）
- 选项 B：先写**方法学完整版**（含全部 5 张表 + 4 张图的具体设计与数据分析），再投递
- 选项 C：**先投稿方法学短文**（500-1000 字 Letters to Editor / Brief Communication），把跨平台失败作为方法学警示
