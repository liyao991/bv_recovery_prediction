# -*- coding: utf-8 -*-
"""框架 v3.2：G8/G10/G13/G15 + 图表编号统一（精确替换，带断言）"""
import io, sys

PATH = r"<REPO_DIR>\results\paper_framework_v2.md"

with io.open(PATH, "r", encoding="utf-8") as f:
    text = f.read()

edits = [
    # G8：摘要 Methods 与最终协议脱节 -> 重写为矩阵协议主口径
    ("LASSO-LR / RF / XGB 三模型；验证：队列内重复分层 holdout、Xiao→Mollin 跨平台外推、基线 CST IV 亚组分层分析（含去 Lactobacillus 对照与 Spearman 检验）；**特征家族对照：α 多样性 vs 属级相对丰度（同一亚组内分别建模）**。",
     "主分析采用 LASSO 逻辑回归（随机森林、XGBoost 作敏感性），以三队列属级交集（83 属）为统一特征空间、20 种子内层交叉验证为主口径；验证：队列内重复分层 holdout、两条跨平台外推（Xiao 训练 → Mollin / Gustin 测试，seed-mean bootstrap 95%CI）、基线 CST IV 亚组分层分析（含去 Lactobacillus 对照与 Spearman 检验）；特征家族对照：α 多样性 vs 属级相对丰度（同一亚组内分别建模）。", 1),
    # G8：模型明确主分析 LASSO-LR
    ("- 模型：LASSO-LR / RF / XGB（同 v1 网格）",
     "- 模型：主分析 LASSO-LR（L1 正则，liblinear）；RF / XGB 作敏感性（结论一致，见补充材料）", 1),
    # G8：指标删 PR-AUC/Brier/P(AUC>0.5)（正文零兑现）
    ("- 指标：AUC（+95% CI）、P(AUC>0.5)、PR-AUC、Brier",
     "- 指标：AUC（内部 CV 报 ±SD；外推报 seed-mean bootstrap 95%CI）", 1),
    # G10：bootstrap CI 条件于训练集声明
    ("- 随机种子主分析 RNG=42；稳健性多种子",
     "- 随机种子主分析 RNG=42；稳健性多种子\n- 外推 AUC 的 95%CI 反映测试集抽样变异、条件于训练队列（Xiao），训练集抽样变异未纳入（双层 bootstrap 留作未来）", 1),
    # G13：45 维 vs 46 维对齐（α 多样性 4 项）
    ("出现率 ≥5% 过滤（42 属）+ α 多样性（Shannon/Chao1/ASV，可加 Pielou evenness）= 45 维",
     "出现率 ≥5% 过滤（42 属）+ α 多样性 4 项（Shannon / Chao1 / ASV richness / Pielou evenness）", 1),
    # G13：附表 S3
    ("- 附表 S3：45 维特征清单（含 α 家族/属级家族分组）",
     "- 附表 S3：特征清单（42 属 + α 多样性 4 项，含 α 家族/属级家族分组）", 1),
    # G13：图编号统一（核心图1/2 -> 图1/2）+ 误差棒 SD->CI
    ("**核心图 1 · 三层信号框架柱状图**（`line_B/results/fig_three_tier_auc.png`）：纵轴 AUC，横轴三层；误差棒=队列内 CV/外推 20 种子 SD；橙色虚线标 0.5（随机）。**核心图 2 · richness 森林图**（`line_B/results/fig_richness_forest.png`）：左=raw、右=rarefaction 校正；Gustin −0.454→−0.392 显著稳健（红），Xiao/Plummer 方向相反不显著（蓝）。",
     "**图 1 · 三层信号框架柱状图**（`line_B/results/fig_three_tier_auc.png`）：纵轴 AUC，横轴三层；内部 CV 报 ±SD、外推报 seed-mean bootstrap 95%CI（图需同步重绘）；橙色虚线标 0.5（随机）。**图 2 · richness 森林图**（`line_B/results/fig_richness_forest.png`）：左=raw、右=rarefaction 校正；Gustin −0.454→−0.392 显著稳健（红），Xiao/Plummer 方向相反不显著（蓝）。", 1),
    # G13：森林图缺 Mollin 图注说明
    ("图下注释一行：\"Gustin 信号 rarefaction 校正后仍稳健；Xiao/Plummer 均不显著且方向相反 → 队列异质性。\"",
     "图下注释一行：\"Gustin 信号 rarefaction 校正后仍稳健；Xiao/Plummer 均不显著且方向相反 → 队列异质性；Mollin 因 SI 表重建无种级 richness 而未纳入。\"", 1),
    # G13：6.1 图 S1 -> 图 1
    ("三层信号框架（图 S1）表明",
     "三层信号框架（图 1）表明", 1),
    # G13：5.4 旧锚点数字补口径
    ("- CST 阈值 0.2-0.7：AUC 稳定 0.60-0.72（结局定义稳健）",
     "- CST 阈值 0.2-0.7 遍历：全人群 AUC 稳定（0.60-0.72，早期单队列锚点口径，仅作稳健性说明）", 1),
    # G15：高 AUC 声称标文献待补
    ("（Xiao 2019 五组分面板 AUC 0.73；Plummer 2023 Prevotella/Gardnerella 关联；2026 年联合模型 AUC 0.928）",
     "（Xiao 2019 五组分面板 AUC 0.73；Plummer 2023 Prevotella/Gardnerella 关联；近期机器学习联合模型报道 AUC 达 0.93〔待补具体文献〕）", 1),
    # G15：6.2 批评软化
    ("**近期高 AUC 预测模型（0.7–0.93）**：多在全人群/横断面评估，未限定治疗对象人群。本研究去 LB 对照证明：全人群 0.6 信号几乎完全来自标签耦合——**这类报告的高 AUC 可能存在同样的耦合假象**，建议在亚组内复评。",
     "**近期报道高 AUC 的预测模型（0.7–0.93〔待补具体文献〕）**：多在全人群/横断面评估，未限定治疗对象人群。本研究去 LB 对照证明：全人群 0.6 信号几乎完全来自标签耦合——**若这类模型在全人群评估，其 AUC 可能包含类似的标签耦合成分**，建议在治疗对象亚组内复评。", 1),
]

failed = []
for old, new, expect in edits:
    n = text.count(old)
    if n < expect:
        failed.append((old, "匹配 %d 次（预期 >=%d）" % (n, expect)))
    else:
        text = text.replace(old, new)
        print("OK  (%d 处) %s" % (n, old[:36]))

if failed:
    print("\n[FAILED]")
    for old, msg in failed:
        print("  -", msg, "->", old[:60])
    sys.exit(1)

with io.open(PATH, "w", encoding="utf-8") as f:
    f.write(text)

print("\n全部替换完成，已写回文件。")
