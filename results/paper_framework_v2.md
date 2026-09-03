# 线B 论文框架（v3.2 定稿）

> 版本 v3.2（2026-09-01）· 中文写作稿 · 待据此撰写论文正文（IMRaD）
> 结局定位：菌群学康复（随访乳杆菌优势度 ≥50%，作为临床治愈的生态代理，非临床硬结局）
> **叙事核心**：治疗前横断面 16S 菌群构成在临床治疗对象（基线 CST IV 亚组）内接近随机，且不被分辨率提升（种级）或时间点后移（D7）改善——指向结局的生物学决定因素；失调总量（richness）与结局关联呈队列异质性（Gustin 内部稳健、跨队列未复现），构成「单队列内部稳健、跨队列消失」的可重复性危机实证。

---

## 0. 投稿目标定位

| 期刊 | IF | 2025中科院 | 新锐分区 | JCR | 适合度 | 理由 |
|---|---|---|---|---|---|---|
| **JCM（J Clin Microbiol）** | 5.9 | **2区（医学）** | **2区** | **Q1** | ★★★★★ | **三轨全达标（2025中科院2区+新锐2区+JCR Q1），收稿范围契合（临床微生物检验+AI/ML 诊断应用），非 OA 无 APC——线 B 首选目标** |
| Microbiome | 12.7 | 1区（生物学） | 1区 | Q1 | ★★★ | 1 区超目标，冲刺选项（难度大） |
| npj Biofilms and Microbiomes | 9.2 | 1区（生物学 Top） | 1区 | Q1 | ★★★ | 1 区超目标 |
| Briefings in Bioinformatics | 7.7 | 1区（生物信息学，待核） | 待核 | Q1 | ★★★ | 方法学向，1 区难度 |
| BMC MIDM | 3.5 | 待核（估计 3-4 区） | 待核 | 待核 | ★★ | 可能不达 2 区 |
| JGAR（**留给线 A**） | 3.2-3.7 | 3区（医学；药学小类2区） | 3区 | Q2 | （线A） | **线 B 不投**（主题+分区都不对口）；耐药主题对口线 A WGS 论文 |

**建议**：线 B 首选 **JCM**（三轨 2 区/Q1 最稳，主题契合临床微生物检验+ML）；Microbiome/npj 为 1 区冲刺备选；**JGAR 留给线 A**（用户原"主攻 JGAR"预期需修正——线 B 投 JGAR 只有 3 区，达不到 2 区目标）。

**AI+医学交叉赛道候选（2026-08 调研补充）**：

| 期刊 | IF | 2025中科院 | JCR | 适配度 | 备注 |
|---|---|---|---|---|---|
| **AIM（Artif Intell Med）** | 7.8 | 2区（医学；计算机AI/医学信息均2区） | Q1 | ★★★★ | 收稿=AI在医学+方法新颖性；需强化"预测精度有限/特征家族对照"方法学叙事 |
| **IEEE JBHI** | 6.8 | 2区（医学） | Q1 | ★★★★ | 国人占比 58.78% 友好，发文量大（659/年）；偏工程/建模，ML+菌群合适 |
| **JBI（J Biomed Inform）** | 4.5-5.9 | 2区（医学） | Q2 | ★★★ | 审稿快（2天初审）、非OA、版面费可减免；稳妥备选 |
| npj Digital Medicine | 15-18 | 1区 | Q1 | ★（够不到） | 要求前瞻验证+算法审计，中国占比仅3%，门槛极高，不投入 |
| Lancet Digital Health | 25.5 | 1区 | Q1 | ★（够不到） | 要求前瞻/RCT/大样本，录用率~12%，不投入 |

> **AI+医学投稿要点（2026 实测趋势）**：① 论文类型四分类（methods/development-validation/prospective/deployment），我们是 **development-validation** 型（TRIPOD+AI 报告标准）；② 审稿人最看重**跨机构外部验证**——我们多队列验证是强项；③ 需补：TRIPOD+AI 清单、校准曲线、Brier、公平性/亚组分析（框架 4.4/5.2 已有规划）；④ 热门导致稿件拥挤+门槛提升，投 AI 刊需报告标准完备。

> **⚠️ 分区体系变更（2026-08-29 核实，重要）**：中科院文献情报中心 2026-03-27 官方声明**自 2026 年起不再更新/发布期刊分区表**（22 年体系停更）；独立第三方"新锐学术"2026-03-24 推出《新锐期刊分区表》（一区=学科前5%，二区=6-20%，含 ESCI/中文期刊，37 种标 under review，每年 3 月发布，有争议）。高校过渡方案：已发表论文按 2025 中科院分区认定（向前就近原则）；新投稿建议综合考虑中科院2025版 + 新锐分区 + JCR + 近5年IF。**选刊策略 = 双轨/三轨保险：目标期刊需 2025中科院 ≥2区 与 新锐 ≥2区 均达标（+JCR Q1 更稳），避免单位口径切换导致降区；候选期刊在新锐分区中的具体位置待查证后填入。**

---

## 1. 标题（定稿 A）

**选项 A（临床问题导向 · ✅ 最终采用）**：
> "Limited predictive value of the pretreatment vaginal microbiome for microbiological recovery (Lactobacillus dominance) after BV treatment: composition, taxonomic resolution, and early (D7) dynamics lack predictive signal; richness–outcome associations are cohort-dependent"

中文：治疗前阴道菌群检测对 BV 治疗后菌群学康复（乳杆菌优势）的预测价值有限：构成、分辨率与早期（D7）动态均缺乏预测信号，richness–结局关联呈队列依赖

**选项 B（负结果导向）**：
> "Pretreatment vaginal microbiome composition near-randomly predicts BV cure within the baseline CST-IV subgroup: a multi-cohort negative-result study"

中文：在基线 CST-IV 亚组内，治疗前阴道菌群构成对治愈的区分接近随机：BV 多队列负结果研究

**选项 C（可重复性警示导向）**：
> "Composition versus magnitude in BV outcome prediction: richness signals are cohort-dependent—a caution against small-cohort overinterpretation"

中文：BV 结局预测中的构成 vs 总量：richness 信号呈队列依赖——对小队列过度解读的警示

**✅ 已定稿**：选项 A —— 同时点出"预测价值有限"（硬负结果）与"richness 队列异质性"（可重复性警示），信息密度最高、最诚实；副标题**显式列出「构成 + 分辨率 + 早期（D7）动态」三个实测维度**（均为负结果），堵住审稿人"只测了构成"的质疑——这三者恰恰是我们系统扫过、逐一证伪的维度，不是"没想到"。（选项 B/C 仅留作历史备选：若日后改投 AIM/JBHI 等方法学刊，可切换到 C 以凸显可重复性/方法学贡献。）

---

## 2. 摘要框架（250 词，IMRaD）

### 2.1 Background
细菌性阴道病（BV）治疗后复发率高（30 天内 30-50%）。既往研究多报告"治疗前菌群可预测结局"，但**预测对象常是全人群（含基线健康者），未聚焦真正的临床治疗对象——基线菌群失调（CST IV）患者**。该人群内的真实可预测性尚不清楚。

### 2.2 Methods
整合四个公开 16S 队列（Xiao 2016 / Plummer 2023 / Mollin 2022 / Gustin 2022；DADA2 或原管线），统一菌群学结局（随访乳杆菌相对丰度跨越 50% 阈值 = 菌群学康复，作为临床治愈的生态代理）。主分析采用 LASSO 逻辑回归（随机森林、XGBoost 作敏感性），以三队列属级交集（83 属）为统一特征空间、20 种子内层交叉验证为主口径；验证：队列内重复分层 holdout、两条跨平台外推（Xiao 训练 → Mollin / Gustin 测试，seed-mean bootstrap 95%CI）、基线 CST IV 亚组分层分析（含去 Lactobacillus 对照与 Spearman 检验）；特征家族对照：α 多样性 vs 属级相对丰度（同一亚组内分别建模）。

### 2.3 Results
全人群层：Xiao 队列内 CV AUC 0.626±0.163（n=68）、跨平台外推 0.521（测试集 bootstrap 95%CI 0.339–0.702）——治疗前菌群在全人群层面有中等预测力，但该信号主要来自结局定义与 Lactobacillus 的耦合（无泄漏 holdout 佐证：去 LB 后全人群 AUC 0.652→0.567）。聚焦临床治疗对象（基线 CST IV 亚组）：**Xiao 内部、Gustin 内部及两条 Xiao 外推路径均 ≈ 随机**（Xiao 内部 0.535±0.127；Xiao→Mollin 0.501[95%CI 0.299–0.714]、Xiao→Gustin 0.512[95%CI 0.296–0.729]；Mollin 内部因 n=30 无结论；Plummer 因结局极度偏斜（事件率 ~90%）不纳入 AUC 建模评估，仅用于 richness 方向与生物学一致性检验）——构成层预测精度不被测量维度改善：**分辨率提升（种级）不优于属级**（CST4 亚组：Xiao 0.511→0.492、Gustin 0.529→0.511）、**时间点后移（D7）未复现**（Gustin D7 单独 0.483、Δ(D7−D0) 0.450，均低于 D0 的 0.530；短暂清除陷阱 D7→D30 持续康复 PPV 0.35）。**失调总量（richness）呈队列异质性**：仅在 Gustin 队列观察到 richness–结局关联（种级 AUC 0.776，ρ=−0.454，p=0.008），该信号在 Gustin 队列内部经多重稳健性检验（rarefaction 校正后 ρ=−0.392/p=0.024、bootstrap 95%CI [−0.713, −0.213] 不跨 0、leave-one-out 33/33 保持显著、置换 p=0.005、混杂排除 partial ρ=−0.431/p=0.014）排除小样本噪声；但 Xiao/Plummer 队列方向相反且不显著，未跨队列复现，且功效分析支持此为真实队列差异（Xiao 功效足以排除 AUC≥0.71 的大效应而信号 ≈0）。**结局精细化（二元 vs 连续）**显示连续型结局下 Gustin richness 信号反而更强（ρ_vs_fuLB=−0.565, p=0.001；ρ_vs_ΔLB=−0.557, p=0.001）——呈剂量-反应增量，结论对结局粒度稳健。

### 2.4 Conclusions
在基线 CST IV 亚组（真正的治疗对象）内，治疗前横断面 16S 检测对菌群学康复（乳杆菌优势）的区分接近随机，且预测精度**不被分辨率提升（种级）或时间点后移（D7 早期反应）改善**——精度受限源于结局的生物学决定因素（复发/再定植主导），超出治疗前构成所携带的信息。**失调总量（richness）呈队列异质性**：Gustin 队列基线 richness 与结局相关且内部稳健（经多重稳健性检验排除小样本噪声与标签耦合），且在连续型结局下信号更强（ρ=−0.565/−0.557，呈剂量-反应增量）——但未跨队列复现，提示关联具队列特异性、须独立验证。本研究发现治疗前横断面 16S 检测对 BV 治疗后菌群学康复（乳杆菌优势）的预测能力有限，并警示以全人群 AUC 代表临床预测力、及对单队列信号过度解读的常见误区。

### 2.5 Keywords
Vaginal microbiome · Bacterial vaginosis · CST IV · Limited predictive value · Treatment outcome · 16S rRNA · Negative result

---

## 3. 引言（800-1000 字）

### 3.1 背景（约 400 字）
- BV 全球患病率高，甲硝唑治疗后 30 天内复发率 30-50%；复发机制与菌群恢复失败相关
- CST 分型（Ravel 2011）：CST IV（厌氧菌主导、多样性高）与 BV 直接对应，是治疗的主要对象
- 近十年大量研究探索"基于菌群的 BV 复发预测"（Xiao 2019 五组分面板 AUC 0.73；Plummer 2023 Prevotella/Gardnerella 关联；近期 ML 诊断 BV 报道横断面准确率高〔Celeste 2023 npj Digit Med〕，但针对「治疗后复发」且限定治疗对象人群的预测研究仍少）——但多未明确限定预测人群

### 3.2 关键空白（约 350 字）
- 多数以预测为目的的研究在全人群评估，**未区分"基线即失调者"（临床治疗对象）与"基线健康者"**——全人群 AUC 可能由标签耦合（基线 LB 高 → 必康复）虚高
- 跨队列验证不足；结局定义异质（临床标签 vs 分子定义）导致可比性差
- 缺乏对"该预测任务本身精度有限"的量化检验——即单次治疗前测序能提供多少预测预后的信息、其"上限"为何，既往从未被定量检验（多数研究默认基线菌群足以决定结局，而未检验这一前提本身）
- **既往文献信号矛盾从未被系统解释**：报告"基线菌群可预测结局"的研究多使用 **α 多样性或 qPCR 浓度**（Gustin 2022 richness/evenness；van de Wijgert 2019 浓度），而使用**属级相对丰度**的研究多报告弱/无信号（Xiao 2016、Mollin 2022）——**特征家族（α 多样性 vs 相对丰度）对可预测性的影响从未在同一队列、同一亚组、同一结局口径下对照检验**，这正是本研究差异化所在

### 3.3 本研究目标（约 200 字）
- 用统一 CST 结局 + 属级特征对齐，在**四个公开队列**上分层评估：全人群 vs 基线 CST IV 亚组
- 量化"治疗前横断面 16S 检测预测菌群学康复（乳杆菌优势）"的预测精度
- 提供跨平台（Xiao 454 → Mollin MiSeq）、跨管线、跨人群的外推证据

### 3.4 贡献（5 条）
1. **首个"预测精度有限"负结果检验**：在临床治疗对象（基线 CST IV 亚组）内，治疗前 16S 属级特征区分结局接近随机（矩阵协议 AUC 0.535±0.127；无泄漏 holdout 独立佐证 0.596 [0.388, 0.802]，CI 含 0.5），多队列（内部 CV 与两条外推路径）一致接近随机（跨平台/管线/人群）
2. **归因证据（全人群 vs 亚组两个层面）**：全人群 0.6 信号来自结局定义与 LB 的耦合（全人群去 LB 0.652→0.567，显著下降）；而亚组内去 LB 几乎不变（0.596→0.584）→ 负结果干净、不依赖 LB 单一特征
3. **预测精度不被测量维度改善**：分辨率提升（种级不优于属级）、时间点后移（D7 未复现）都无法提高"≈随机"的预测精度 → 精度受限源于结局的生物学决定因素（复发/再定植主导），超出治疗前构成所携带的信息
4. **失调总量 vs 构成的精确区分 + 队列异质性警示（差异化核心）**：在同一队列内精确区分"失调总量（richness）与具体构成（属/种谱）"的可预测性差异；并发现 richness–结局关联**呈队列异质性**（Gustin 内部稳健、Xiao/Plummer 方向相反未复现），以"单队列内部稳健、跨队列消失"为微生物组预测领域的可重复性危机提供实证
5. **临床落点（探索性）**：Gustin 队列基线 richness≥22 种者康复率更低（24% vs 44%，Fisher p=0.28 不显著，描述性分层，待独立验证）；"短暂清除陷阱"（D7 转好 PPV 0.35）警示动态监测的价值

---

## 4. 方法（1500-2000 字）

### 4.1 数据源与样本量（4 队列矩阵）

**队列检索与纳入/排除标准**：纳入满足以下条件的公开纵向 16S 队列——① 公开可获取原始 16S 数据且具样本级配对（基线→随访）；② BV 患者接受标准治疗（甲硝唑）并具明确随访结局；③ 基线可判定 CST（可计算乳杆菌相对丰度）；④ 覆盖技术异质性（测序平台 454 vs MiSeq、测序区 V1-V3/V3-V4/V4、生信管线 DADA2 vs 原管线、参考库 SILVA vs RDP），为跨平台/跨管线/跨人群三重外推提供证据基础。**排除**：伴侣干预或含男性样本（StepUp-1 PRJNA398590 / StepUp-2 PRJNA735440）、无结局标签（NuSwab PRJNA1171405）、水对照。最终纳入 4 个队列（见下矩阵）。

| 队列 | PRJNA | 平台 | 测序区 | 配对样本 | 基线 CST4 | 结局定义 |
|---|---|---|---|---|---|---|
| Xiao 2016 (Sci Rep) | PRJNA302078 | 454 GS Junior | V1-V3 SE | 68 对（D0→D30） | **54/68** | CST（LB 阈值 50%） |
| Plummer 2023 (SToPBV) | PRJNA592384 | Illumina MiSeq | V3-V4 PE | 49 对（M0→M1） | **43/49** | CST（LB 阈值 50%） |
| Mollin 2022 | SI 表重建 | Illumina MiSeq | V4 | 33 对（V1→V2） | 30/33 | CST（LB 阈值 50%） |
| Gustin 2022 (= CONRAD BV 试验) | PRJNA691964 | Illumina MiSeq | V4 2×150 | 33 例（V1→V2→V3；101 runs = 33×3 + 2 水对照） | **33/33 全 CST4**（专门治疗队列） | CST（LB 主导；4 分类结局可映射） |

> ⚠️ **关键修正（2026-08-29 已核对）**：基线 CST4 = Xiao **54/68**、Plummer **43/49**（旧记录 47/39 已作废）。

### 4.2 DADA2 与特征
- Xiao truncLen=450/maxEE=2（454 SE）；Plummer① c(240,220)/c(2,2)；Gustin 2×150 试跑 merge 校准
- assignTaxonomy tryRC=TRUE；SILVA v138.1（重建 30,399 条训练集）
- **richness 定义**：主分析 richness = species_observed（SILVA 种级注释的检出种计数）。定义敏感：naive 全种计数（丰度>0 即计）得 ρ=−0.342(p=0.051)，与主分析 ρ=−0.454(p=0.008) 方向一致但显著性跨 0.05 边界——正文须报定义口径并说明该敏感性。
- 属级聚合（tax_glom）→ 出现率 ≥5% 过滤（42 属）+ α 多样性 4 项（Shannon / Chao1 / ASV richness / Pielou evenness）；log1p + StandardScaler；**α 多样性家族单独建模以做特征家族对照（见 5.6，特征已现成）**
- 跨队列对齐：SILVA vs RDP 属名标准化 → 131 属交集（Xiao→Mollin 实验）
- **全集口径**：Xiao/Plummer/Gustin 三队列属名交集 = **65 属**（靶标组合全集口径）；与 Xiao 单队列 5% 过滤子集（42 属）结论一致——纳入 Plummer/Gustin 后交集收窄但不改变核心负结果（详见 5.6 ② 分辨率对照）。
- **特征空间口径（审稿人必问）**：主分析（5.1 三队列复现）特征空间 = **Xiao/Mollin/Gustin 三队列交集 83 属**，非各队列全集。特征空间对 AUC 影响显著——Gustin 内部 CV 用 83 属交集得 **0.436**、用 179 属全集得 **0.520**，故 Methods 必须写明「特征空间 = 三队列交集 83 属」，否则审稿人会问「为何不用全集」。（注意：本条 83 属 = Xiao/Mollin/Gustin 交集，与上条 65 属 = Xiao/Plummer/Gustin 交集属不同队列组合，两者都保留并标注队列组合即可）

### 4.3 结局（统一 CST 口径）
- `outcome = 1 if 结局时点 LB 相对丰度 ≥ 阈值（默认 0.5）`（= **菌群学康复**；本研究结局为菌群学代理结局，作为临床治愈的生态学代理，非临床硬结局）；敏感性分析阈值 0.2-0.7 遍历
- **核心分层**：基线 CST IV 亚组（基线 LB<50%）＝ 临床治疗对象，为主分析人群
- **结局精细化（见 5.8）**：主分析用**菌群学二元**（LB≥50%）；敏感性用**菌群学连续**（随访 LB 丰度 + 恢复幅度 ΔLB）→ 检验结论对结局定义的**粒度**（二元/连续）是否稳健。**临床对照**：二元结局仅回答"是否达到康复"，连续型结局则量化"恢复的程度"——后者保留恢复幅度的完整信息，用以确认"预测精度有限"并非人为二分阈值所致。临床硬结局（Amsel/Nugent 复发）因样本级标签不可得（Plummer 原文无逐例清单、SRA/BioSample 无临床字段、run_alias 的 R/RV 为母研究「任意随访复发」而非「1 个月复发」），仅作 Methods 局限带过（6.5 第 6 条）

### 4.4 模型与验证
- 模型：主分析 LASSO-LR（L1 正则，liblinear）；RF / XGB 作敏感性（结论一致，见补充材料）
- 验证：① 队列内重复分层 holdout（200 次 60/40，无泄漏）；② Xiao→Mollin 跨平台外推（20 种子）；③ 基线 CST IV 亚组分层分析（含去 LB 对照 + Spearman）；④ 合并 5 折 CV（仅作诊断，不进正文）
- 指标：AUC（内部 CV 报 ±SD；外推报 seed-mean bootstrap 95%CI）

### 4.5 统计与可复现
- 随机种子主分析 RNG=42；稳健性多种子
- 外推 AUC 的 95%CI 反映测试集抽样变异、条件于训练队列（Xiao），训练集抽样变异未纳入（双层 bootstrap 留作未来）
- Python 3.12（sklearn 1.9 / xgboost 3.4）+ R 4.6.1（DADA2）+ SILVA v138.1
- 脚本：`line_B/scripts/`（train_main / diag_cst4_evidence / xiao_mollin_cross / xiao_holdout 等）
- **检验分类（confirmatory vs exploratory）**：confirmatory（先验假设检验）= Gustin richness 方向性复现——对 Gustin 2022（AJOG 226:225.e1-15，PMID 34560047）标题级已发表方向性先验（Fig 4C：持续清除组基线 richness/evenness 显著低于失败组）的独立再分析确认（directional replication）；exploratory（探索性）= 其余全部（α×分辨率 120 项矩阵、构成层等），其显著性经多重比较校正评估（Bonferroni 0.05/120、BH-FDR），探索性结果仅作描述不作推断

### 4.6 分析范式（对照传统 meta 分析）
> **本研究不是传统 meta 分析，而是个体水平多队列整合（individual-participant-data, IPD）**：对 Xiao/Plummer/Gustin 三个队列的 fastq **统一重跑 DADA2、统一 CST 结局口径、统一模型协议**后再分别建模（Mollin 因原始数据限制采用原文 SI 表重建，作为管线异质性的天然纳入）。与传统"提取已发表效应量→加权合并"的 aggregate meta 分析相比，IPD 有三点不同：① 数据来自原始测序而非二手效应量，可自行控制混杂（结局口径、注释管线一致）；② **首要终点是跨队列复现性（方向一致性），而非合并效应量**；③ 当队列异质性真实存在时，加权合并会用一个 pooled 估计**掩盖**"单队列稳健、跨队列消失"的异质性结构——本研究刻意不做合并，是方法学自觉而非分析缺漏（若需回应对"为何不合并"的质疑，可补 `metafor` 随机效应合并作敏感性分析，预期 I² 高、pooled 不显著，反证"不应合并"）。

---

**实测验证（支撑"刻意不合并"的方法学立场）**：
- **`metafor` REML 随机效应合并**（示意性敏感性分析；注：k=3 下 I²/τ² 估计功效不足，仅展示 pooled 跨 0 的事实，不作决定性论据）：pooled ρ=**−0.110**（95%CI [−0.425, +0.229] 跨 0、p=0.53）；I²=72.3%、Cochran Q=6.75 (df=2, p=0.034) → **异质性显著、合并估计无意义**（见 `fig_c6_meta_forest.png`；DL 版 −0.109 结论稳健）。
- **IPD 朴素拼接**（三队列 n=130）：AUC=**0.620**（[0.518, 0.725]）、ρ=−0.204（p=0.020）——**此 0.62 非真实信号，而是 Gustin 队列结构主导的混杂假象**（Gustin 33/33 全 CST4，richness–LB 强负相关把拼接结果拉高）。与下方交互项**绑定解读**，避免审稿人误读为"合并有信号"。
- **异质性检验（cohort×richness 交互项）**：LR=9.387（df=2），**p=0.0092** → 队列间 richness 效应方向相反，**异质性真实**，不能用 pooled 掩盖。
- **四队列版**（n=160，统一共同属 richness 61 属，纳入 Mollin）：②a AUC=**0.591**（[0.511, 0.683]）、ρ=−0.157（p=0.048）；②b 交互 LR=9.480（df=3，**p=0.0235**）——**4 队列仍稳健：朴素拼接 ≈0.59、交互项仍显著**，进一步坐实"队列异质性真实"。

→ **结论**：合并估计无意义（I² 高、pooled 跨 0）；朴素拼接的"0.62 信号"是队列结构混杂假象，**必须与交互项 p=0.009 绑定解读**——避免审稿人误读；异质性真实，本研究"刻意不合并"的方法学立场得到实证支持。

## 5. 结果（分层报告，勿合并报告）

### 5.0 三层信号框架（论文核心图表）

| 层次 | 任务 | 代表 AUC | 说明 |
|---|---|---|---|
| ① 横断面诊断 | 你现在是不是 BV | **0.9+** | 菌群现状高度区分（文献一致） |
| ② 全人群恢复预测 | 会不会康复 | **~0.53-0.63**（Xiao 内 0.626 / Xiao→Mollin 0.521 / Xiao→Gustin 0.556） | 3 条外推路径全 >0.5 但信号弱，主要来自 LB 标签耦合 |
| ③ **基线 CST4 亚组（真任务）** | **同为失调，谁更难恢复** | **≈0.5**（Xiao 0.535 / Xiao→Mollin 0.501 / **Xiao→Gustin 0.512**） | **多队列稳健复现 ≈ 随机 = 预测精度有限** |

> **图 1 · 三层信号框架柱状图**（`line_B/results/fig_three_tier_auc.png`）：纵轴 AUC，横轴三层；内部 CV 报 ±SD、外推报 seed-mean bootstrap 95%CI（图需同步重绘）；橙色虚线标 0.5（随机）。**图 2 · richness 森林图**（`line_B/results/fig_richness_forest.png`）：左=raw、右=rarefaction 校正；Gustin −0.454→−0.392 显著稳健（红），Xiao/Plummer 方向相反不显著（蓝）。

### 5.1 三队列复现证据（基线 CST4 亚组；矩阵协议 = Xiao/Mollin/Gustin 交集 83 属 + 20 种子内层 CV 选 C，为最终采用口径；早期单队列锚点降级为 CI 佐证；外推 AUC 点估计与 95%CI 同源自 seed-mean 复现管线单次 bootstrap）

> **队列纳入说明**：本复现覆盖 Xiao / Mollin / Gustin 三队列。Plummer 队列因基线 CST IV 亚组结局极度偏斜（43 例中 39 例康复、仅 4 例未康复），LASSO 模型退化为空模型（AUC 恒 0.5，系占位符而非性能测量），无法可靠评估，故不纳入本复现，仅用于 richness 方向与生物学一致性检验（见 5.2 报告矩阵）。

**证据 1 · 队列内（Xiao，无泄漏 200 次分层 holdout）**
- n=54（27+/27-）：矩阵协议口径 AUC **0.535±0.127**；无泄漏 200 次分层 holdout 独立佐证 **0.596 [0.388, 0.802]**，P(AUC>0.5)=0.81 → CI 含 0.5
- 去 LB 对照（亚组层面）：0.596→0.584（几乎不变，负结果干净）；对照全人群层面去 LB：0.652→0.567（显著下降 = 标签耦合，见摘要 2.3）
- 单特征 LB：Spearman 0.306（p=0.025）——"模型预测力弱"而非"生物学无关"

**证据 2 · 跨平台外推（Xiao 454 → Mollin MiSeq，20 种子）**
- 基线 CST4 亚组外推 AUC 矩阵协议口径 **0.501（95%CI 0.299–0.714）**（≈随机；早期锚点 holdout 佐证 0.527±0.027）；Mollin 内部 0.649（n=30，holdout [0.36,0.92] 覆盖 0.5，统计无结论）

**证据 3 · 第 3 独立队列（Xiao → Gustin，矩阵协议，2026-08-29 新增）**
- 外推 AUC **0.512（95%CI 0.296–0.729）**（≈随机）；Gustin 内部 0.436
- **关键结构事实**：Gustin（CONRAD BV 试验）33/33 基线全为 CST4（专门治疗队列，所有人基线 LB<50%）→ 其"全人群"= "CST4 亚组"（A4=B4=0.436）；A3(0.556)>B3(0.512) 之差 = CST 分层/标签耦合贡献的量化

**证据 4 · 文献佐证（Xiao 2016 原文 + Mollin 2022）**
- Xiao 原文：D0/D7 不按治愈/失败聚类，仅 D30 与结局同期相关；Mollin 原文："治疗前组成在结局组间无总体区分" → **3 个独立团队一致观察**

→ **叙事升级：多队列稳健复现（跨 454/MiSeq、DADA2/Mothur、SILVA/RDP），审稿人"偶然性"质疑基本打消**

### 5.1b 临床关键菌靶向预测（先验靶标 vs 全集对照）

> **动机**：Plummer 2023 报告治疗前 Prevotella（普雷沃菌）、治疗后 Gardnerella（加德纳菌）与复发相关，提示**关键菌**可能携带预后信息，而临床检测面板正是围绕这一类关键菌设计的。这一文献线索引出一个必答问题——"构成层 ≈ 随机，是否因用了全菌群特征、而非这些关键菌？"本节以先验关键菌靶向特征正面检验。

**靶标**：7 个临床关键菌属（Lactobacillus、Gardnerella、Atopobium、Megasphaera、Sneathia、Mobiluncus、Prevotella 属聚合）+ 2 个可注释种（G. vaginalis、A. vaginae）。**边界声明**：8 靶标中 5 个种级靶标（L. crispatus / L. iners / P. bivia / BVAB2 / Megasphaera-1）在 SILVA v138.1 种级注释不可得（三队列均未注释到种）——属级代理即 16S 技术上限，为 16S 与 qPCR/CE 分子检测面板的差异化点（属级代理 vs 种/型级靶标，两法互补）。

**结果（基线 CST4 亚组，LASSO-L1 主协议，与主分析同协议）**：
- **单变量（9 靶标 × 3 队列，27 项）**：Xiao 队列 Prevotella（ρ=+0.346，p=0.010）与 Lactobacillus（ρ=+0.306，p=0.025）为 nominal 信号，但 27 项 BH-FDR 联合校正后最小 q=0.28，**全部不显著**（完整表见附表 S7）；Gardnerella 属级在 Xiao 方向为负（高 GV→不康复，ρ=−0.259）但 p=0.059 边际、校正后消失。
- **组合**：Xiao 队列 7 属靶标组合 CV AUC **0.7197±0.156**，高于 65 属全集（0.5283）与主矩阵 B1（0.535）——靶标组合在单队列内部确有中等信号；**但跨队列完全不复现**：Gustin 0.4612（<0.5 方向反转）、Plummer 因仅 4 例未康复退化为空模型（AUC 不可估计，勿把 0.5 当数字）。

**核心句**：**临床关键菌靶向特征与全特征集共享同一模式——单队列内可有中等信号、跨队列均不可迁移**。提示可预测性的瓶颈不在特征选择（先验靶标 vs 全集），而在结局本身的生物学决定因素。

### 5.2 报告矩阵（写作铁律）

- ✅ **可报告**：Xiao 队列内 CV 0.626±0.163（n=68）/ 无泄漏 holdout 0.652 [0.500, 0.797] / 跨平台外推 0.521（→Mollin）、0.556（→Gustin） —— 原假说"治疗前 16S 预测结局"在全人群成立
- ⚠️ **必须声明**：0.6 信号主要来自 CST 标签与 LB 的定义耦合；基线 CST4 亚组内预测力 ≈ 随机（任务本身预测精度有限）
- ❌ **不进正文**：合并 CV 0.62-0.67（批次效应假象，LOCO 0.28-0.42 实证）；Plummer 队列内 CV 0.43（90% 阳性率 + 5 例阴性的统计精度极限，**非 class_weight 反转**——已实测排除）

### 5.3 生物学一致性（Xiao LASSO 系数）
Lactobacillus (+1.38，预测康复)、Mobiluncus (-1.11)、Atopobium (-0.71)、Gardnerella (-0.36)、Mycoplasma (-0.73)、Ureaplasma、Anaerococcus —— **7 个核心 BV 关联菌方向全部正确** → 模型捕获真实生物学信号（但不足以支撑个体预测）

### 5.4 稳健性
- CST 阈值 0.2-0.7 遍历：全人群 AUC 稳定（0.60-0.72，早期单队列锚点口径，仅作稳健性说明）
- CLR 变换无系统性增益（5 种子 0.486）
- 合并 CV 0.62-0.67 = 批次假象（不进正文，Methods 中说明）

### 5.5 D7/D8 动态预测（已填结果：D7 早期反应不能预测 30 天结局）
> 数据基础：Gustin（PRJNA691964）= CONRAD BV 试验，33 例 v1(D0)→v2(D7)→v3(D30) 三访视，结局定义与 CST 同构（LB 主导）；Xiao 亦有 D7（67 例）、Plummer① 有 D8 可佐证。

**结果（Gustin，33 例）**：

| 特征组 | 全人群 AUC | CST4 亚组 AUC |
|---|---|---|
| D0 only | 0.530 | 0.528 |
| D7 only | **0.483（最低）** | 0.481 |
| D0+D7 | **0.586（最佳）** | 0.576 |
| Δ(D7-D0) | 0.450 | 0.443 |

**关键发现 · 短暂清除陷阱**：D7 转好 ≠ D30 康复——D7 时 LB≥50% 的 20 例中仅 **7 例** D30 持续康复（**PPV 0.35**，Fisher OR 1.21，p=1.0），与 Gustin 论文"短暂清除（transient clearance）35.7%"类别吻合。

**结论**：
1. **与 Xiao 原文一致（非矛盾）**：Xiao 2016 原文即观察到"D0/D7 菌群不按治愈/失败聚类、仅 D30 与结局同期相关"（即 D7 早期菌群无法区分结局）；本研究在 Gustin 队列独立复现——D7 单独预测力反而最低（0.483），两团队结论一致；
2. D0+D7 组合仅微弱增益（0.586 vs 0.530，不显著）；
3. **短期清除后多复发，D7 早期反应不能预测 30 天结局** → "预测精度有限"叙事延伸至时间维度（结局受复发/再定植主导）；
4. n=33 局限，所有差异不显著，措辞用"未观察到 / 微弱增益"。

> **⚠️ 与原作者结论的关系（防止"相同数据、不同结论"的误读）**：本研究的 D7 结论与 Gustin / Xiao 原文**不矛盾**，三者提问层面不同——① Gustin 报告的"短暂清除 35.7%"是**现象描述**（D7 转好者中部分 D30 复发），本研究"短暂清除陷阱 PPV 0.35"是对同一现象的**定量复现**（互相印证）；② Gustin 的"richness 关联复发"用的是 **α 多样性标量**（见 5.6），本研究 5.5 的 D7 预测用的是**属级相对丰度**——两者是不同"尺子"，正是本研究核心论点"信号在总量、不在具体构成"；③ Xiao 原文结论本就是"D0/D7 不按结局聚类、仅 D30 同期相关"，与本研究一致。故"相同数据得出不同结论"实为"不同问题、不同尺子"，而非结论冲突。

→ **叙事影响**：D7 这条"正向路径"未兑现，但贡献了一个新的临床洞察——**"短暂清除陷阱"**（治疗后 1 周菌群转好不代表 1 月后康复），呼应论文 4 分类结局（sustained/transient/delayed/no clearance），成为 Discussion 的独立亮点。

### 5.6 特征家族对照：α 多样性 vs 属级相对丰度（预留分析——差异化核心，用户指定优先级）

> **定位**：本研究与 Gustin 2022 的关系 = **再分析确认 + 跨队列独立检验**，而非对立。
> - **再分析确认（一致点，confirmatory 先验假设检验）**：Gustin 报告"基线 richness/evenness 高者更易治疗失败"（α 多样性有信号）——这是本研究的**唯一 confirmatory 假设**（对 Gustin 2022 AJOG 标题级已发表方向性先验 Fig 4C 的独立再分析确认，directional replication），我们**同意**这一判断，并预期在 Xiao 队列复现；若复现，两个独立队列互证，"α 多样性有价值"更硬。
> - **补全（独有点）**：Gustin **只测了 α 多样性这一把尺子**，从未在同一人群里测过属级相对丰度。本研究补上两块 Gustin 提供不了的信息：① 属级相对丰度在治疗对象内 ≈ 随机（负结果）；② 同一队列内两家族对照，证明"有信号/没信号的分野在尺子"。
> - 三者合起来 = "**可预测性取决于特征家族**"——把过去"有人说能预测（α/qPCR/种级）、有人说不能（属级丰度）"的文献矛盾统一解释，**容纳而非推翻 Gustin 的发现**。

**设计**（与 5.3 同口径）：
1. 人群：基线 CST4 亚组（Xiao n=54；Gustin n=33）+ 全人群对照
2. 家族 A（α 多样性）：Shannon / Chao1 / ASV richness / Pielou evenness —— 已含在特征矩阵，**无需重算特征，仅按家族单独建模**
3. 家族 B（属级相对丰度）：现有 42 属 log1p（对照基准）
4. 模型：LASSO-LR（同主分析）；报告两家族 AUC ± 95%CI、去 LB 对照、Spearman 方向（与 Gustin "richness 高→失败"方向对比）
5. 复现检验：Xiao（发现）→ Gustin（独立复现，且是 Gustin 原文"正向"队列）——α 在 Gustin 复现其 richness 信号 → 特征家族对照成立；不复现 → Gustin 信号不稳健（n=28 小样本），本研究负结果更全面
6. **【新扩展·分辨率维度】属级 vs 种级**（特征现成：SILVA v138.1 已注释到种，仅需种级聚合，不必 tax_glom 到属）——与家族对照同框架，形成**"家族 × 分辨率"矩阵**（α / 属级 / 种级）；文献佐证：Mollin 种级组合（M. lornae + Gsp07 + F. magna）、FLAST 全长 16S（L. iners→愈、P. bivia→未愈）均显示种级信号更强；若复现"种级 > 属级"→"分辨率"成为第二条可预测性分界，论文再升一级

**结果（α 家族 × 分辨率矩阵，120 项 × 20 种子 CV）**：

**① α 家族（Gustin，n=33，治疗队列；注意：跨队列仅此队列显著，见下方"队列异质性"）**：
- **richness 家族独强**：ASV 0.727 / 属 0.681 / **种 0.776**（主口径 LB≥50%，与全研究统一结局口径一致）；evenness 0.35-0.47 倒挂无信号；entropy 0.50-0.61
- **方向显著（仅 Gustin）**：基线种级 richness 与 D30 康复 **ρ=−0.454 (p=0.008)**，richness≥22 种组康复率 24% vs 低组 44%（切点=基线 richness 中位数 22 种；Fisher 双侧 p=0.28，组间差异不显著）→ **基线越杂越难康复**（描述性分层）（与 Gustin 论文方向一致；但 Xiao/Plummer 队列方向相反且不显著，见下方"队列异质性"）
- Xiao 全人群全家族无信号；Xiao CST4 亚组 evenness 弱信号（pielou 0.624，rho=+0.222 p=0.106 不显著）；Plummer 结局偏斜不可评估

**② 分辨率对照（任务 1，属级 vs 种级）**：

| 队列 | 属级 full | 种级 full | 属级 CST4 | 种级 CST4 |
|---|---|---|---|---|
| Xiao | 0.621 | 0.575 | 0.511 | 0.492 |
| Gustin | 0.527 | 0.499 | 0.529 | 0.511 |
| Xiao→Gustin 种级外推 | **0.582**（93 种交集） | — | — | — |

- **"种级 > 属级"假设不成立**：种级不优于属级（Xiao/Gustin 略降），CST4 亚组内种级同样 ≈0.5 → **分辨率提升无法提高预测精度，排除"属级粒度不足"解释**，反哺叙事更硬（精度受限源于结局的生物学决定因素（复发/再定植主导），超出治疗前构成所携带的信息）
- **意外发现：种级跨平台外推 0.582 > 属级 0.556**——种级特征跨平台迁移更稳健（SILVA 种名跨队列一致），作为独立方法学报告点

**叙事精确化（关键升级）**：
- 从"可预测性取决于特征家族"进一步精确为——**"失调的总量（richness）与具体构成（属/种谱）的可预测性不同"**：α 标量在治疗对象内（Gustin 全 CST4）有信号 0.73-0.78，而高维构成特征同亚组 ≈0.5 → 预测精度受限针对"构成信息"
- **分辨率梯度非单调（主口径修正）**：richness 信号种级最高（0.776）、属级最低（0.681）、ASV 居中（0.727）——主口径（LB≥50%）下"随分辨率单调增强"不成立，须如实报告为"种级最高、非单调"；α 为标量聚合、种级滤噪不丢信息的机制解释不再适用

**队列异质性警示**：
- richness–结局关联**未跨队列复现**：Spearman ρ Xiao 全/CST4 = −0.054(p=0.662)/+0.045(p=0.745)、Plummer = +0.100(p=0.492)/+0.071(p=0.650)，均方向相反且不显著；**仅 Gustin −0.454(p=0.008) 显著**
- **rarefaction 校正后仍稳健**：Gustin ρ 由 −0.454(p=0.008) → −0.392(p=0.024)；控制测序深度的 partial ρ=−0.389(p=0.025)、logistic 多变量调整后 richness p=0.040（#C20 补充）。**但 richness×深度交互显著（置换 p=0.004，低深度组 ρ=−0.77 强、高深度组 −0.02 消失）——测序深度混杂不能完全排除，以局限呈现，勿写「排除」**（rarefy 仅控制深度对 richness 估计的影响）。Xiao/Plummer 校正后仍不显著（−0.033/+0.109、+0.095/+0.090）
- **连续型结局下信号更强（剂量-反应增量，见 5.8）**：Gustin richness 换连续型结局后 ρ_vs_fuLB −0.565(p=0.001)、ρ_vs_ΔLB −0.557(p=0.001)，强于二元 −0.454；rarefy 后仍稳健（−0.378/−0.383，p=0.030/0.028）
- **脆弱性模拟——结果指向稳健而非脆弱**：bootstrap 1000 次 ρ 95%CI = [−0.713, −0.213] 不跨 0（p<0.05 比例 86.3%）；leave-one-out ρ 范围 [−0.551, −0.461]，33/33 保持显著、0/33 翻转符号（无单样本驱动）；置换检验 1000 次 p=0.005 → **Gustin 内部稳健，排除小样本噪声/过拟合**
- **功效分析——异质性为真实而非功效差异**：Xiao（n=54）最小可检测 AUC 0.711，其 richness 观测 ≈0（CST4 ρ=+0.045；full ρ=−0.054）→ **功效足以排除 AUC≥0.71 的大效应；中等效应不能排除（见 6.5 局限）**；Plummer（n=43，事件率 0.91 失衡）最小可检测 AUC 0.869，设计上几乎检不出中等效应；Gustin（n=33）最小可检测 AUC 0.778，richness 信号（0.776）恰好卡在检出线 → 解释仅 Gustin 报出信号
- **定性（升级）**：richness = "Gustin 队列观察到、内部稳健（多重稳健性检验：rarefaction 校正 + bootstrap/LOO/置换 + 混杂排除）、但未跨队列复现"的**队列特异信号** → 作为"单队列内部稳健、跨队列消失"的可重复性危机实证（比"小队列脆弱"更精确：脆弱点在跨队列，不在队列内部）

**临床可操作点（探索性）**：Gustin 队列基线 richness ≥22 种者 D30 康复率更低（24% vs 44%，Fisher 双侧 p=0.28，描述性分层）——潜在风险分层指标（**待独立验证**）

> **核心图 2（位置）· richness 森林图**（`line_B/results/fig_richness_forest.png`）：三队列五亚组 Spearman ρ，**左=raw / 右=rarefaction 校正**；Gustin ρ=−0.454→−0.392 仍显著（红），Xiao/Plummer 方向相反不显著（蓝）。图下注释一行："Gustin 信号 rarefaction 校正后仍稳健；Xiao/Plummer 均不显著且方向相反 → 队列异质性；Mollin 因 SI 表重建无种级 richness 而未纳入。"

**剂量-反应增量（结局精细化的正向证据，详见 5.8）**：将 Gustin richness 关联从二元结局扩展为连续型结局（随访 LB / ΔLB），信号反而增强——ρ_vs_fuLB = −0.565（p=0.001）、ρ_vs_ΔLB = −0.557（p=0.001），均强于二元的 −0.454；rarefy 后仍稳健（−0.378/−0.383，p=0.030/0.028）。基线 richness 越高，随访 LB 越低、ΔLB 越负，呈**单调剂量-反应**——是"结局精细化"维度最有力的正向量化证据，但**仍仅在 Gustin 成立**（Xiao/Plummer ΔLB 方向相反）。

**写作注意**：① "构成 vs 失调总量"对比是论文新论点（5.7 叙事的精确版）；② 跨队列 α 家族信号不一致（Gustin=richness / XiaoCST4=evenness），仅 Gustin 显著 → 措辞用"Gustin 队列观察到、未跨队列复现、待独立验证"；③ C1(0.689) 实为被 shannon 稀释的 richness 信号，论文建议报告纯 richness；④ 已补 rarefaction 校正（深度混杂排除）+ 脆弱性模拟/功效分析（小样本噪声排除）；n=33 的局限已从"小样本脆弱"转为"跨队列未复现"，措辞升级。

### 5.7 正向路径盘点（论文叙事 = 负结果 + 正向出路，防止被读成"菌群检测无用论"）

| # | 正向路径 | 时点 | 特征类型 | 状态（2026-08-29 定稿） |
|---|---|---|---|---|
| 1 | **失调总量（richness）预测菌群学康复** | 治疗前（基线） | α 多样性·richness（标量） | ⚠️ **队列异质性（内部稳健、跨队列未复现）**：Gustin 种级 richness 0.776，ρ=−0.454 (p=0.008，rarefaction 后 p=0.024；**多重稳健性检验**=bootstrap/LOO/置换 + rarefaction 校正 + 混杂排除)，Xiao/Plummer 未复现；**连续型结局下信号更强（ρ=−0.565/−0.557）呈剂量-反应增量** → **队列特异探索性线索，待独立验证** |
| 2 | 治疗后早期 16S 动态预测 | D7/D8（治后） | 属级丰度 + D0→D7 Δ | ❌ **本研究未复现**：Gustin D7 单独最低(0.483)，D0+D7 仅微弱增益(0.586)；**短暂清除陷阱**(PPV 0.35) |
| 3 | 种级标志物（分辨率提升） | 治疗前 | 种级 | ❌ **排除**：种级不优于属级(Xiao/Gustin 略降)；但**种级跨平台外推更稳健**(0.582>0.556) |
| 4 | 横断面诊断（三层框架①） | 现状 | 属级丰度 | 文献共识 AUC 0.9+（非本研究新结果，证明菌群-失调关联极强） |
| 5 | 生物学一致性佐证 | — | LASSO 系数方向 | ✅ **本研究已实证**：7 个核心 BV 关联属方向全对（LB+、Gardnerella- 等），支撑"预测力弱 ≠ 生物学无关" |

> 说明（最终版）：三条"正向路径"验证结果各不相同——① richness（失调总量）**仅在 Gustin 队列兑现**（内部稳健但跨队列未复现），降为队列特异探索性线索；② D7 动态**未兑现**，反而贡献"短暂清除陷阱"洞察；③ 分辨率提升**被排除**（种级不优于属级）。**论文叙事收敛为**：治疗前横断面菌群构成（属级/种级）与治疗后早期反应（D7）都无法预测 30 天结局；失调总量（richness）仅在 Gustin 队列观察到关联（内部稳健、跨队列未复现）——负结果三件套（构成/分辨率/时间点）+ 一个"队列异质性"警示，临床落点清晰（探索性）。

### 5.8 结局精细化：二元 vs 连续

> **目的**：检验核心结论「预测精度有限」是否对结局定义的**粒度**（二元 vs 连续）稳健。临床硬结局（Amsel/Nugent 复发）因样本级标签不可得而无法纳入，仅作 Methods 局限（见 6.5 第 6 条）。具体结果数字已写入 5.6「队列异质性」段与图 2（forest plot），本段仅作方法学收口。

| 档 | 结局定义 | 队列 | 结论（详见 5.6 段） |
|---|---|---|---|
| 主分析·二元 | 随访 LB≥50% = 康复 | 三队列 | 构成层 ≈0.5（5.1/5.6）；Gustin richness ρ=−0.454(p=0.008) |
| 敏感性·连续 | 随访 LB 丰度（log1p）+ 恢复幅度 ΔLB | 三队列 | 构成层同样 ≈0.5；Gustin richness ρ=−0.565/−0.557 反而更强 |

**方法学收口**：

1. **结论对结局粒度稳健**：横断面构成谱 ≈ 随机在二元与连续两档下都成立（Xiao/Plummer/Gustin 全 0.44–0.60），不因结局换粒度而改变 → 「预测精度有限」是结局的客观特性，非特定阈值假象。
2. **Gustin richness 在连续型结局下呈剂量-反应增量**（详见 5.6 段）：二元 ρ=−0.454 → 连续 ρ=−0.565/−0.557 反而增强，说明基线 richness 与随访 LB 丰度、恢复幅度存在**单调剂量-反应**——比二元截断更能捕捉关联强度，但**仍仅在 Gustin 成立**（Xiao/Plummer ΔLB 方向相反），未跨队列复现。
3. **rarefaction 校正后信号仍稳健**（详见 5.6 段）：Gustin ρ=−0.565→−0.378(p=0.030)、−0.557→−0.383(p=0.028)，排除测序深度混杂。
4. **「连续型结局救回信号」假设不成立**：三队列 ΔLB 方向不一致（Xiao 边缘正、Plummer 弱正、Gustin 强负），连续型结局没把 Xiao/Plummer 的 richness 信号"救"回来，反而进一步坐实队列异质性。

**叙事落点**：结局定义的粒度（二元/连续）**不改变核心结论** → 支撑 Discussion「对结局定义稳健」段；连续型结局的剂量-反应增量是「结局精细化」维度最有力的正向量化证据。

---

## 6. 讨论（1200 字）

### 6.1 核心发现
- **预测精度有限——横断面构成谱 ≈ 随机在临床治疗对象内多队列稳健复现**：基线 CST IV 亚组（真正的治疗对象）内，治疗前横断面 16S 属级检测区分结局的 AUC 跨队列稳定在 0.44–0.60（Xiao 内 0.535、Gustin 内 0.436、两条 Xiao 外推 0.501/0.512，测试集 bootstrap 95%CI 均含 0.5）。Mollin 队列内点估计 0.60–0.65、但 n=30 小样本、holdout 区间宽（[0.36,0.92]）覆盖 0.5——**统计上无结论**，既不支持也不反驳随机，故不纳入「≈随机复现」主证据。三层信号框架（图 1）表明：菌群-BV 关联在横断面（诊断）层面极强（AUC 0.9+），但"从失调基线预测恢复"这一任务本身信息量有限——**关键意义是：预测精度受限针对的是构成信息（高维特征），不是整个预测问题**。
- **预测精度不被测量维度改善**：分辨率提升（种级不优于属级，CST4 亚组：Xiao 0.511→0.492、Gustin 0.529→0.511）、时间点后移（D7 早期反应未复现，Gustin D7 单独 AUC 0.483、Δ(D7−D0) 0.450，均低于 D0 的 0.530；D0+D7 仅微弱增益 0.586）都无法提高"≈随机"的预测精度——指向**精度受限来自结局的生物学决定因素（复发/再定植主导），而非测量方式**。这一发现将"预测精度有限"叙事从分辨率与时间点两个维度同时锁定。
- **失调总量（richness）——单队列内部稳健、跨队列消失的队列异质性**：Gustin 队列种级 richness 与二元结局 ρ=−0.454（p=0.008），基线 richness≥22 种者 D30 康复率更低（24% vs 44%，Fisher p=0.28）。该信号在 Gustin 内部经**多重稳健性检验**：① rarefaction 校正后 ρ=−0.392（p=0.024）排除测序深度混杂；② bootstrap 1000 次 95%CI [−0.713, −0.213] 不跨 0；③ leave-one-out 33/33 保持显著、0/33 翻转符号（无单样本驱动）；④ 置换检验 p=0.005；⑤ 混杂排除：控制基线乳杆菌丰度后 richness–结局关联几乎不衰减（partial ρ=−0.431, p=0.014），排除该信号为 LB 标签耦合假象。但 Xiao/Plummer 队列方向相反（Xiao full −0.054 / CST4 +0.045、Plummer full +0.100 / CST4 +0.071）且不显著，**未跨队列复现**——**精确读为"单队列内部稳健、跨队列消失"，是微生物组预测领域可重复性危机的典型缩影**。功效分析进一步支持此为真实队列差异：Xiao 最小可检测 AUC 0.711 而观测 ≈0，功效足以排除 AUC≥0.71 的大效应，中等效应不能排除（见 6.5 局限）。
- **连续型结局下的剂量-反应增量（结局精细化的正向证据）**：将 Gustin richness 关联从二元结局（康复/未康复）扩展为连续型结局（随访 LB 丰度 / 恢复幅度 ΔLB），关联强度反而增强——ρ_vs_fuLB = −0.565（p=0.001）、ρ_vs_ΔLB = −0.557（p=0.001），均强于二元的 −0.454；rarefy 后仍稳健（−0.378/−0.383，p=0.030/0.028）。基线 richness 越高，随访 LB 越低、恢复幅度越小，呈**单调剂量-反应**——比二元"康复/未康复"更能捕捉关联强度，是"结局精细化"维度最有力的正向量化证据。但该剂量-反应同样**未跨队列复现**（Xiao/Plummer ΔLB 方向相反），故仍是队列特异。
- **0.6 的全人群信号 = 标签耦合，不等于临床预测力**：全人群 AUC ~0.63 主要源于基线 LB 高→必康复的定义耦合（全人群去 LB 对照 0.652→0.567，显著下降 = 标签耦合）；而亚组内去 LB 几乎不变（0.596→0.584），说明亚组负结果干净、不依赖 LB 单一特征——两层对照共同支撑"全人群信号来自耦合、治疗对象内无信号"的结论。

- **richness 队列效应三象限（四队列版，方向性观察）**：Mollin（n=30）richness 效应 slope ≈ 0（−0.004），与 Gustin 负（−0.092）、Xiao/Plummer 弱正形成"三象限"分布（注：此系四队列版 61 属共同 richness 的回归 slope，与正文种级 Spearman ρ 量纲/口径不同，仅看方向；Xiao 弱正对应 CST4 亚组口径 +0.032，与主分析 CST4 ρ=+0.045 同向，主分析另报 full −0.054 均 ≈0 不显著，无方向矛盾）——**方向性观察**：richness 与结局的关联高度依赖队列结构（治疗类型/随访时长/测序平台），**不能作为通用预测因子**。n=30 小样本，Mollin 中性 slope 可能只是噪声，故措辞用"方向性观察"而非"三象限实锤"；但即便如此，"跨队列效应方向不一致"本身已坐实队列异质性。

### 6.2 与既往研究呼应
- **Xiao 2016 原文"动态观测"结论**（D0/D7 不聚类、仅 D30 与结局同期相关）与本研究"预测精度有限"独立互证；其"D7 引入补强"的设计线索在本研究 Gustin 队列**未复现**（D7 单独最低 0.483）——两个独立团队在不同人群、平台、时点下得到一致的"治疗前构成无法预测 30 天结局"结论，**是负结果稳健性的最强外部证据**。
- **Plummer 2023 的 Prevotella 治前/Gardnerella 治后关联**：属关联性证据（方向），未做亚组内预测力检验。本研究的基线 CST4 亚组内属级 AUC ≈0.5 与该报告的"属级信号弱"一致，但**精确化了边界条件**——亚组内（排除标签耦合）后，关联进一步塌缩至随机。
- **近期报道高 AUC 的预测模型（0.7–0.93）**：多在全人群/横断面评估，未限定治疗对象人群。本研究去 LB 对照证明：全人群 0.6 信号几乎完全来自标签耦合——**若这类模型在全人群评估，其 AUC 可能包含类似的标签耦合成分**，建议在治疗对象亚组内复评。（近期 ML+BV 预测具体文献在投稿前按 PMID 核准后补入，见第 10 节 F 组）
- **与 Gustin 2022 的特征家族对话（关键定位）**：Gustin 报告基线 richness/evenness 高者更易失败——本研究在其队列内**再分析确认 richness 信号且方向一致**（ρ=−0.454，p=0.008），并通过**多重稳健性检验**（rarefaction + bootstrap/LOO/置换 + 混杂排除）证明该信号在 Gustin 内部不来自小样本噪声/过拟合。但该信号在 Xiao/Plummer 队列**未复现**，且功效分析显示 Xiao 功效足够（最小可检测 AUC 0.711）而信号为零——指向**真实的队列特异性**（可能源于治疗/人群/混杂差异），而非统计噪声。这一"单队列内部稳健、跨队列消失"的模式，是微生物组预测领域**可重复性危机**的典型缩影。同时本研究补上 Gustin 未测的属级/种级构成（≈0.5），与 Gustin 的 α 标量信号形成**"具体构成不可预测、失调总量可观察但队列特异"**的完整图景——**容纳而非推翻** Gustin 的发现。
- **"短暂清除陷阱"与 Gustin 论文吻合**：本研究 D7 转好 PPV 0.35（即 D7 时 LB≥50% 的 20 例中仅 7 例 D30 持续康复）与 Gustin 论文"短暂清除（transient clearance）35.7%"类别吻合——**两个独立团队独立观察到 D7 早期反应不能预测 D30 持续康复**，为这一现象提供外部佐证，提示治疗后动态监测需重新评估其临床价值。
- **种级跨平台迁移更稳健（0.582 > 属级 0.556）**：SILVA 种名跨队列一致性优于属级归并（不同 SILVA/RDP 命名方案的属级映射差异）——**作为独立方法学报告点**，提示在跨平台外推场景下种级特征反而更稳健。

### 6.3 方法学贡献
1. **预测精度负结果检验范式**：分层（全人群 vs 基线失调亚组）+ 去标签耦合对照（去 LB）应成为 BV 预测模型的标准检查——本研究发现即便相对低 AUC（如 0.6）也可能是标签耦合假象，需在临床治疗对象（基线 CST IV）亚组内复评。
2. **统一 CST 结局 + 结局精细化（粒度×层面）范式**：解决样本级临床标签缺失 + 结局异质，使跨队列可比；进一步以**二元+连续**两档结局检验结论对结局粒度的稳健性（Gustin richness 在连续型结局下信号反而更强，呈剂量-反应增量，验证结论不依赖二元截断）——为结局定义敏感性问题提供可复用方法学模板。
3. **归因修正**：16S 跨平台外推失败（LOCO 0.28-0.42）主因是特征未对齐 + 结局异质，统一后外推可行（0.521，Xiao→Mollin）——平台本身并非障碍，关键在特征与结局的**对齐质量**。
4. **特征家族对照范式**：在同一队列/亚组/结局口径下对照 α 多样性与属级相对丰度两特征家族的可预测性，为"菌群预测文献中正负结果并存"提供统一解释框架——具体构成≈0.5 而 α 标量在治疗队列内（Gustin）有信号 0.73-0.78，预测力差异来源于特征家族而非生物学差异（结果见 5.6）。
5. **可重复性警示范式**：以 richness 的队列异质性为实证，示范"单队列显著信号须经**跨队列方向性检验 + 深度校正 + 内部稳健性模拟（bootstrap/LOO/置换）**三层验证后方可用于临床解释"，呼应微生物组预测领域的可重复性危机；并借功效分析区分"真无信号（Xiao）"与"设计检不出（Plummer）"——为负结果研究提供方法学范本：**报告负结果时应同时报告功效**，让审稿人能够区分"未检出"与"不存在"。
6. **种级跨平台外推更稳健的方法学点**：属级粒度的命名学异质（不同 SILVA/RDP 命名方案映射差异）是跨平台外推的隐藏障碍，种级反而跨队列一致——提示跨平台外推场景下种级特征或为更优选择。
7. **"IPD 多队列整合 + 复现性优先"范式（对照传统 meta）**：区别于"提取已发表效应量→加权合并"的 aggregate meta 分析，本研究基于原始测序数据统一再分析（IPD，见 4.6），并以**跨队列复现性（方向一致性）而非合并效应量**为首要终点。该选择的方法学价值在于——当队列异质性真实存在时，随机效应合并会用一个 pooled 估计掩盖"单队列稳健、跨队列消失"的异质性结构；本研究以**森林图 + 异质性定性 + 功效归因**替代合并，为微生物组预测的"多队列证据综合"提供一种比 meta 合并更尊重异质性的范式（注：此定位的价值取决于审稿人是否认可，故措辞不作过度拔高，仅作为"多队列证据综合的替代范式"陈述）。

### 6.4 临床含义与未来方向
- **单靠治疗前一次横断面 16S 检测（属级/种级构成）不足以指导个体化治疗决策**——本研究在临床治疗对象（基线 CST IV）亚组内接近随机（AUC 0.44–0.60），排除了"属级粒度不足"和"时间点过晚"两个潜在解释，结论直接指向结局的生物学决定因素（复发/再定植主导），超出治疗前构成所携带的信息。
- **治疗后 D7 早期反应同样不可靠（短暂清除陷阱 PPV 0.35）**——D7 时 LB≥50% 的 20 例中仅 7 例 D30 持续康复，与 Gustin 论文"短暂清除 35.7%"独立吻合——提示**单时间点动态监测的临床价值需重新评估**，建议转向多时间点连续监测或更晚时点评估（如 D14/D21）。
- **探索性临床信号 = 失调总量（richness）**：Gustin 队列基线 richness≥22 种者 D30 康复率更低（24% vs 44%，描述性分层）——潜在风险分层指标，但**未跨队列复现**且属队列特异，**须独立验证后方可临床应用**；可作为"谁需更积极随访/延长治疗"的假设。
- **基线失调程度与恢复程度的连续关系（有限讨论，回应临床"严重程度"关切）**：主分析以二元结局（康复/未康复）呈现，但结局本身可还原为更细的连续度量——随访乳杆菌丰度与恢复幅度（ΔLB）。以连续结局分析时，基线 richness 与恢复程度呈**剂量-反应**（ρ=−0.565/−0.557）：基线菌群越失调（richness 越高），随访乳杆菌越低、恢复幅度越小。这提示基线失调程度与预后之间是"连续"而非"二元"的关系——**基线越失调，恢复越差**。临床含义：基线失调最重者或需更积极的随访与延长治疗；但此梯度仅在 Gustin 单队列观察、且基于菌群学结局（非 Amsel/Nugent 临床评分），须独立验证后方可外推。
- **特征选择含义——临床检测应优先纳入失调度标量**：本研究证明**失调总量（α 标量）有可观察信号而具体构成（属/种）无**，提示临床分子检测在物种分型之外，应优先纳入**失调度标量**（如乳杆菌相对占优度）作为辅助读数。Gustin 在 α 标量上的可观察信号（种 richness 0.776，ρ=−0.454）提示：即便队列特异，**失调度方向**仍是有意义的探索性指标。
- **临床关键菌靶向预测同样不优于随机（划清"诊断关联"与"预后预测"边界）**：7 个临床关键菌属的靶标组合在 Xiao 队列内 CV AUC 0.72、但跨队列不复现（Gustin 0.46）、单变量 27 项 FDR 校正后全不显著——提示**临床检测面板的关键菌靶标，同样无法在基线预测菌群学康复**。这一结论并不否定关键菌在 BV 诊断中的作用（横断面诊断 AUC 0.9+），而是划清"诊断关联"与"预后预测"的边界：关键菌与 BV 状态强关联（诊断层面），但其基线丰度无法预测谁康复、谁复发（预后层面）。
- **未来方向**：构成/分辨率/时间点三维都触顶 → 需转向**功能基因/菌株级**维度（金标准验证队列 LSVF 获批后补做）；近期可考虑**菌株水平（strain-level）的复发 vs 再感染判别**（Xiao 2016 原文提到 9 疑似复发 + 8 再感染的判别需求），是构成分类的下一粒度。

### 6.5 局限
1. **亚组样本量小**（Xiao CST4 n=54；Gustin n=33）——两条外推路径 CI 上界分别为 0.71（→Mollin）与 0.73（→Gustin）、Xiao CST4 holdout CI 上界 0.80，**中等效应（0.6–0.7）仍不能排除，必须写局限**。但功效分析显示 Xiao 功效足以排除 AUC≥0.71 的大效应，其构成层负结果为"真无信号"而非样本不足；Gustin n=33 最小可检测 AUC 0.778，richness 信号（0.776）恰卡检出线——**该阳性结果须经独立队列复现后方可外推**。
2. **richness 信号仅在 Gustin 队列统计显著；该信号在 Gustin 内部稳健（多重稳健性检验：bootstrap/LOO/置换 + rarefaction 校正 p=0.024 + 混杂排除）但未跨队列复现**（Xiao/Plummer 方向相反不显著）→ 队列特异性的生物学来源（治疗方案/人群特征/微生物生态差异）未明，须独立队列验证。另：Gustin 观测效应（0.776）恰卡最小可检测线（0.778），真实效应或受赢家诅咒（winner's curse）高估、可能更小，异质性结论系针对观测效应量而言、对真实效应量留有余地。
3. **Plummer 结局偏斜（事件率 0.91）无法可靠评估**；Mollin 来自 SI 表重建（小样本 n=30）；种级注释率参差（Xiao 69.9% / Plummer 62.7% / Gustin 49.1%）影响分辨率对照的可比性——也是种级跨平台外推 0.582 > 0.556 的反向解释（注释率高的队列贡献更多）。
4. **结局为 CST 代理（LB 阈值），非临床硬结局（Amsel/Nugent）**——需声明。已用**结局精细化（5.8）**检验结论对结局**粒度**（二元 vs 连续）与**层面**（菌群学 vs 临床）的稳健性，**二元+连续一致**（构成≈随机、仅 Gustin richness 有信号且连续下更强），结论对粒度稳健；**临床复发层面因样本级标签不可得而无法纳入**（见第 6 条）。
5. **Gustin 队列 n 差异（本研究 33 vs 原文 28）**：本研究 n=33（基线全 CST4，含全部 SRA 样本去水对照 101 runs = 33×3 + 2），原文报告 n=28（口服甲硝唑治疗试验完成者；治疗失败 25% / 短暂清除 35.7% / 延迟清除 7.1% / 持续清除 32.1%）——**差异源于访视完整性与样本-受试者映射**（33 例均有完整 D0/D7/D30 访视，原文 28 例为治疗完成分析集），不改变"全 CST4 治疗队列"判定。
6. **临床硬结局（Amsel/Nugent 复发）无法纳入分析**：Plummer 原文"入组后 1 个月内复发"的样本级标签**不可得**（ENA/BioSample 无临床字段、补充材料仅汇总统计、正文无逐例清单；SRA run_alias 的 R/RV 后缀为母研究"任意随访复发 M1–M6"标记而非"1 个月复发"——16 个 R/RV 中仅 3 例 ≤M1，事件数极少）——故本研究结局为菌群学 CST 代理（二元 + 连续敏感性），临床复发结局的检验留作未来工作（如 LSVF 队列获批后）。
7. **richness 的跨队列比较受 16S 测序区混杂**：Xiao 用 V1-V3、Gustin/Mollin 用 V4，捕获的 richness 谱本就不同；rarefaction 校正测序深度、不校正测序区，故「Xiao 功效足以排除大效应」的推断假设 richness 跨队列同度量，测序区差异可能造成测量衰减、部分解释 Xiao ≈0。

---

## 7. 结论（200 字）

本研究通过四个公开 16S 队列、统一菌群学结局（乳杆菌优势 ≥50%，作为临床治愈的生态代理）与分层分析，量化了"治疗前阴道菌群预测 BV 治疗后菌群学康复"的预测精度：全人群层面中等信号（AUC ~0.63，主要来自标签耦合）；**聚焦临床治疗对象（基线 CST IV 亚组）后接近随机（AUC 0.44–0.60；外推路径可排除 AUC≥0.71 的大效应、中等效应不能排除）**，且预测精度不被分辨率提升（种级）或时间点后移（D7）改善。**失调总量（richness）与结局的关联呈队列异质性**：仅 Gustin 队列显著且内部稳健（rarefaction 校正 + bootstrap/LOO/置换 + 混杂排除多重稳健性检验排除小样本噪声），**连续型结局下信号更强（ρ=−0.565/−0.557）呈剂量-反应增量**——结论对结局粒度稳健；但 Xiao/Plummer 未跨队列复现，提示关联具队列特异性、须独立验证。结论：菌群学康复的预测精度受限源于结局的生物学决定因素（复发/再定植主导），超出治疗前构成所携带的信息；临床不应过度依赖治疗前横断面菌群构成，失调总量（richness）作为探索性风险指标须独立验证。

---

## 8. 数据与代码可用性
- 数据：ENA PRJNA302078 / PRJNA592384 / PRJNA691964（Gustin）+ Mollin SI 表
- SILVA v138.1；代码：`line_B/scripts/`
- 数据处理产物：`<DATA_DIR>/feat/`（X_features / y_labels / feature_meta）
- GitHub 仓库：待创建（账号 229691084@qq.com 已就绪，建 repo 后回填链接；邮箱信息不入正文）

## 9. 致谢（待补）
- 队列原始作者与 ENA；**2023 年福建省科技创新联合资金项目（2023Y9353）**

## 10. 参考文献（正文用，约 42 篇；按主题分组，编号为正文引用序）

**A. BV 流行病学、诊断与复发（6）**
1. Amsel R, Totten PA, Spiegel CA, Chen KC, Eschenbach D, Holmes KK. Nonspecific vaginitis: diagnostic criteria and microbial and epidemiologic associations. Am J Med. 1983;74(1):14-22.
2. Nugent RP, Krohn MA, Hillier SL. Reliability of diagnosing bacterial vaginosis is improved by a standardized method of Gram stain interpretation. J Clin Microbiol. 1991;29(2):297-301.
3. Allsworth JE, Peipert JF. Prevalence of bacterial vaginosis: 2001-2004 National Health and Nutrition Examination Survey data. Obstet Gynecol. 2007;109(1):114-120.
4. Koumans EH, Sternberg M, Bruce C, et al. The prevalence of bacterial vaginosis in the United States, 2001-2004; associations with symptoms, sexual behaviors, and reproductive health. Sex Transm Dis. 2007;34(11):864-869.
5. Peebles K, Velloza J, Balkus JE, McClelland RS, Barnabas RV. High global burden and costs of bacterial vaginosis: a systematic review and meta-analysis. Sex Transm Dis. 2019;46(5):304-311.
6. Bradshaw CS, Morton AN, Hocking J, et al. High recurrence rates of bacterial vaginosis over the course of 12 months after oral metronidazole therapy and factors associated with recurrence. J Infect Dis. 2006;193(11):1478-1486.

**B. 阴道微生态与 CST 分型（6）**
7. Ravel J, Gajer P, Abdo Z, et al. Vaginal microbiome of reproductive-age women. Proc Natl Acad Sci U S A. 2011;108(Suppl 1):4680-4687.
8. Fredricks DN, Fiedler TL, Marrazzo JM. Molecular identification of bacteria associated with bacterial vaginosis. N Engl J Med. 2005;353(18):1899-1911.
9. Gajer P, Brotman RM, Bai G, et al. Temporal dynamics of the human vaginal microbiota. Sci Transl Med. 2012;4(132):132ra52.
10. Srinivasan S, Hoffman NG, Morgan MT, et al. Bacterial communities in women with bacterial vaginosis: high resolution phylogenetic analyses reveal relationships of microbiota to clinical criteria. PLoS One. 2012;7(6):e37818.
11. Onderdonk AB, Delaney ML, Fichorova RN. The human microbiome during bacterial vaginosis. Clin Microbiol Rev. 2016;29(2):223-238.
12. France MT, Ma B, Gajer P, et al. VALENCIA: a nearest centroid classification method for vaginal microbial communities based on composition. Microbiome. 2020;8(1):166.

**C. 本研究四个原始队列 + 直接相关文献（6）**
13. Xiao B, Niu X, Han N, et al. Predictive value of the composition of the vaginal microbiota in bacterial vaginosis, a dynamic study to identify recurrence-related flora. Sci Rep. 2016;6:26674.（PRJNA302078）
14. Xiao B, Wu C, Song W, et al. Association analysis on recurrence of bacterial vaginosis revealed microbes and clinical variables important for prediction outcome. Front Cell Infect Microbiol. 2019;9:189. DOI: 10.3389/fcimb.2019.00189.（五组分面板，基于 D7 数据）
15. Plummer EL, Vodstrcil LA, Danielewski JA, et al. Combined oral and topical antimicrobial therapy for male partners of women with bacterial vaginosis: acceptability, tolerability and impact on the genital microbiota of couples—a pilot study. PLoS One. 2018;13(2):e0192260.（StepUp-1，PRJNA398590）
16. Plummer EL, Vodstrcil LA, Murray GL, et al. Gardnerella vaginalis clade distribution is associated with bacterial vaginosis recurrence in a cohort of women who have sex with women. mBio. 2021.（StepUp-2，PRJNA735440；卷期待核）
17. Plummer EL, Vodstrcil LA, Fairley CK, et al. Sexual partners, bacterial vaginosis, and the vaginal microbiota: a multicohort pooled analysis（SToPBV/StepUp 合并） J Infect Dis. 2023.（jiad261；卷期待核）
18. Vodstrcil LA, et al. SToPBV「oral contraceptive pill」母研究（PRJNA592384；2019；完整题录待核）
19. Mollin DL, et al. PLoS One. 2022.（复发性 BV 队列，SI 表重建数据源；完整题录待核）
20. Gustin AT, et al. Recurrent bacterial vaginosis following metronidazole treatment is associated with microbiota richness at diagnosis. Am J Obstet Gynecol. 2022;226(2):225.e1-225.e15. DOI: 10.1016/j.ajog.2021.09.018; PMID: 34560047.（PRJNA691964）

**D. 16S 生信与微生物组方法（7）**
21. Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods. 2016;13(7):581-583.
22. Quast C, Pruesse E, Yilmaz P, et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013;41(D1):D590-D596.
23. McMurdie PJ, Holmes S. phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PLoS One. 2013;8(4):e61217.
24. Bolyen E, Rideout JR, Dillon MR, et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol. 2019;37(8):852-857.
25. Bokulich NA, Kaehler BD, Rideout JR, et al. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin. Microbiome. 2018;6(1):90.
26. Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011;17(1):10-12.
27. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.

**E. 机器学习与微生物组预测方法（8）**
28. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. J Mach Learn Res. 2011;12:2825-2830.
29. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. Proc 22nd ACM SIGKDD. 2016:785-794.
30. Topçuoğlu BD, Lesniak NA, Ruffin MT 4th, Wiens J, Schloss PD. A framework for effective application of machine learning to microbiome-based classification problems. mBio. 2020;11(3):e00434-20. DOI: 10.1128/mBio.00434-20.
31. Knights D, Costello EK, Knight R. Supervised classification of human microbiota. FEMS Microbiol Rev. 2011;35(2):343-359.
32. Pasolli E, Truong DT, Malik F, Waldron L, Segata N. Machine learning meta-analysis of large metagenomic datasets: tools and biological insights. PLoS Comput Biol. 2016;12(7):e1004977.
33. Wirbel J, Zych K, Essex M, et al. Microbiome meta-analysis and cross-disease comparison enabled by the SIAMCAT machine learning toolbox. Genome Biol. 2021;22(1):93.
34. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. Adv Neural Inf Process Syst. 2017;30:4765-4774.
35. Chekroud AM, Hawrilenko M, Loho H, et al. Illusory generalizability of clinical prediction models. Science. 2024;383(6679):164-167. DOI: 10.1126/science.adg8538.

**F. 近期 ML+BV 预测与公平性（2）**
36. Celeste C, Ming D, Broce J, et al. Ethnic disparity in diagnosing asymptomatic bacterial vaginosis using machine learning. NPJ Digit Med. 2023;6(1):211. DOI: 10.1038/s41746-023-00953-1; PMID: 37978250.
37. （近期 ML+BV 复发预测高 AUC 报道 1-2 篇，待补具体文献）

**G. 统计与软件（3）**
38. R Core Team. R: A language and environment for statistical computing. Vienna: R Foundation for Statistical Computing; 2024.
39. Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48.
40. van der Maaten L, Hinton G. Visualizing data using t-SNE. J Mach Learn Res. 2008;9:2579-2605.

**H. BV 病理与治疗机制（2，备用）**
41. Muzny CA, Schwebke JR. Pathogenesis of bacterial vaginosis: discussion of current hypotheses. J Infect Dis. 2016;214(Suppl 1):S1-S5.
42. van de Wijgert JHHM, Verwijs MC, et al.（阴道菌群 qPCR 浓度与结局关联；完整题录待核）

> 说明：编号 13-19 为本研究四个队列的直接来源文献，其中 16/17/18/19/42 的完整卷期页码在正式投稿前需按 PMID 逐一核准；正文引用编号与本节一一对应。

## 11. 补充材料
- 附表 S1：各队列 DADA2 参数与质量概况
- 附表 S2：三层信号框架完整数字（4 队列矩阵）
- 附表 S3：特征清单（42 属 + α 多样性 4 项，含 α 家族/属级家族分组）
- **附表 S4（新增）**：分辨率对照（属级 vs 种级，4 队列 full/CST4）+ 种级注释率（Xiao 69.9% / Plummer 62.7% / Gustin 49.1%）
- **附表 S5（新增）**：D7 动态四组特征 AUC（D0/D7/D0+D7/Δ，全人群+CST4）+ D7 CST 转换表（短暂清除陷阱）
- **附表 S6（新增）**：α 家族 × 分辨率矩阵（120 项，ASV/属/种 × 6 指标）
- **附表 S7（新增）**：临床关键菌靶向预测——单变量扫查（9 靶标 × 3 队列，27 项 AUC/ρ/p/q）+ 组合检验（13 层 × 3 队列 × L1/L2 双协议）。**实际数字表已单独成文**：`line_B/results/附表S7_临床关键菌靶向预测_数字表.md`（S7a 单变量 27 项 / S7b 组合检验 L1 主协议含 Plummer 空模型脚注 / S7c 5 个种级靶标 16S 不可得边界声明）
- 附图 S1：Xiao→Mollin 外推 ROC；附图 S2：holdout AUC 分布直方图；附图 S3：特征家族对照（α vs 属级）AUC 对比
- **附图 S4（新增）**：richness 剂量反应（richness≥22 种 vs 康复率）；**附图 S5（新增）**：D7→D30 状态转换图（短暂清除陷阱）
- 附代码 S1/S2：DADA2 R + 建模 Python
