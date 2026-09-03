# -*- coding: utf-8 -*-
"""
#C4-4 收口：全结果总表 + richness ρ 森林图 + 三档 AUC 柱状图
数据源：<OUT_DIR>/modeling/ 下各 results_*.csv（口径已在 paper_framework_v2.md §5 定稿）
输出：<OUT_DIR>/modeling/
  - fig_richness_forest.png          richness–结局 Spearman ρ 森林图（原始 + rarefy 校正，5 亚组）
  - fig_three_tier_auc.png           三层信号框架三档 AUC 柱状图
  - results_c4_master_summary.md     全结果总表（Markdown）
  - results_c4_master_summary.csv    全结果总表（结构化，供后续引用）
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

OUT = "<OUT_DIR>/modeling/"

# ---------------------------------------------------------------------------
# 数据（硬编码自 results_*.csv，逐条注释来源）
# ---------------------------------------------------------------------------

# richness–结局 Spearman ρ（results_richness_direction.csv / results_richness_rarefy.csv）
# (label, n, rho_raw, p_raw, rho_rarefy, p_rarefy, rarefy_depth)
RICH = [
    ("Xiao · 全人群",       68, -0.054, 0.662, -0.033, 0.787, 4254),
    ("Xiao · 基线CST4",     54,  0.045, 0.745,  0.109, 0.431, 4254),
    ("Plummer · 全人群",    49,  0.100, 0.492,  0.095, 0.515, 1341),
    ("Plummer · 基线CST4",  43,  0.071, 0.650,  0.090, 0.565, 1341),
    ("Gustin · 全CST4",     33, -0.454, 0.008, -0.3916, 0.024, 1014),
]

def fisher_ci(r, n):
    """Spearman ρ 的 Fisher z 近似 95%CI（SE = 1/sqrt(n-3)）"""
    r = np.clip(r, -0.999, 0.999)
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    lo = np.tanh(z - 1.96 * se)
    hi = np.tanh(z + 1.96 * se)
    return lo, hi

# 三层信号框架（results_four_cohort_matrix.csv，矩阵协议口径）
# tier, [(sub, auc, sd), ...]
THREE_TIER = {
    "① 横断面诊断\n（你现在是不是 BV）":  [("文献共识", 0.90, None)],
    "② 全人群恢复预测\n（会不会康复）":   [("Xiao 队列内", 0.6256, 0.1633),
                                            ("Xiao→Mollin", 0.5278, 0.0101),
                                            ("Xiao→Gustin", 0.5564, 0.0206)],
    "③ 基线CST4亚组\n（同为失调谁难恢复）": [("Xiao 队列内", 0.5353, 0.1267),
                                            ("Xiao→Mollin", 0.5049, 0.0204),
                                            ("Xiao→Gustin", 0.5124, 0.0000)],
}

# ---------------------------------------------------------------------------
# 图 1：richness ρ 森林图（原始 + rarefy 两面板）
# ---------------------------------------------------------------------------
def plot_forest():
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True,
                             gridspec_kw={"wspace": 0.06})
    panels = [
        (axes[0], "raw", "原始（未校正深度）", "rho_raw"),
        (axes[1], "rarefy", "rarefaction 校正后", "rho_rarefy"),
    ]
    yticks = list(range(len(RICH)))[::-1]  # 自上而下
    labels = [r[0] for r in RICH]

    for ax, key, title, _ in panels:
        ax.axvline(0.0, color="#888888", lw=1.0, ls="--", zorder=1)
        for i, (label, n, r_raw, p_raw, r_rarefy, p_rarefy, depth) in enumerate(RICH):
            y = yticks[i]
            if key == "raw":
                r, p = r_raw, p_raw
            else:
                r, p = r_rarefy, p_rarefy
            lo, hi = fisher_ci(r, n)
            sig = p < 0.05
            color = "#c0392b" if sig else "#5d7fa3"
            ax.plot([lo, hi], [y, y], color=color, lw=2.2, zorder=3,
                    solid_capstyle="round")
            ax.scatter([r], [y], s=72, color=color, zorder=4,
                       edgecolor="white", linewidth=0.8)
            # 标注数值 + p
            ax.annotate(f"{r:+.2f}  (p={p:.3f})", xy=(r, y),
                        xytext=(6 if r >= 0 else -6, 4),
                        textcoords="offset points",
                        ha="left" if r >= 0 else "right",
                        fontsize=8.5, color=color, fontweight="bold" if sig else "normal")
        ax.set_yticks(yticks)
        ax.set_yticklabels(labels if key == "raw" else [])
        ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
        ax.set_xlabel("Spearman ρ", fontsize=10)
        ax.set_xlim(-1.05, 1.05)
        ax.grid(axis="y", color="#eeeeee", lw=0.6, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # 图例
    legend_handles = [
        Patch(facecolor="#c0392b", label="p < 0.05（显著）"),
        Patch(facecolor="#5d7fa3", label="p ≥ 0.05（不显著）"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("基线 richness 与治疗结局的关联（Spearman ρ）· 三队列五亚组",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.text(0.5, -0.08,
             "虚线 ρ=0。误差棒 = Fisher z 近似 95%CI。Gustin 信号 rarefaction 校正后仍稳健（−0.454→−0.392，p=0.008→0.024）；\n"
             "Xiao/Plummer 均不显著且方向相反 → 队列异质性。",
             ha="center", fontsize=8.5, color="#555555")
    fig.tight_layout()
    fig.savefig(OUT + "fig_richness_forest.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------------------
# 图 2：三档 AUC 柱状图
# 设计：3 个 tier 组，组内 1-3 个子柱；子标签写在柱顶上方、数值写在柱顶；
#       整组 tier 名（粗体）作为该组中心 x 轴刻度；图例放在图顶（fig-level）避免遮挡
# ---------------------------------------------------------------------------
def plot_three_tier():
    tier_names = list(THREE_TIER.keys())
    # 三个 tier 组的中心 x 坐标
    group_centers = [0, 1.6, 3.2]
    sub_offset   = [-0.22, 0.0, 0.22]  # 组内三个子柱的相对偏移
    group_colors = ["#5b9bd5", "#e8a33d", "#c0392b"]
    group_labels = ["① 横断面诊断\n（你现在是不是 BV）",
                    "② 全人群恢复预测\n（会不会康复）",
                    "③ 基线CST4亚组\n（同为失调谁难恢复）"]

    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    ax.axhline(0.5, color="#888888", lw=1.3, ls="--", zorder=1)
    ax.text(3.55, 0.503, "随机 (AUC=0.5)", fontsize=9, color="#666666", va="bottom")

    for gi, (tier, items) in enumerate(THREE_TIER.items()):
        m = len(items)
        base = group_centers[gi]
        if m == 1:
            xs = [base]
        else:
            xs = [base + off for off in sub_offset[:m]]
        for (sub, auc, sd), x in zip(items, xs):
            hatch = "///" if gi == 0 else None
            ax.bar(x, auc, width=0.20,
                   color=group_colors[gi], edgecolor="white",
                   hatch=hatch, alpha=0.93, zorder=3)
            if sd is not None and sd > 0:
                ax.errorbar([x], [auc], yerr=[sd], fmt="none",
                            ecolor="#222222", lw=1.1, capsize=3, zorder=4)
            # 数值（柱顶）
            ax.text(x, auc + 0.022, f"{auc:.2f}", ha="center",
                    fontsize=9.5, fontweight="bold", color="#222222", zorder=5)
            # 子标签（再上一行）
            ax.text(x, auc + 0.065, sub, ha="center", fontsize=8.2,
                    color="#444444", zorder=5)

    # x 轴：组中心作为刻度，整组 tier 名作为刻度文本（两行，加粗）
    ax.set_xticks(group_centers)
    ax.set_xticklabels(group_labels, fontsize=10.5, fontweight="bold")
    ax.set_xlim(-0.7, 3.7)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("AUC", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#eeeeee", lw=0.6, zorder=0)

    # 图例：放在图顶部（fig-level 水平），不遮挡任何柱
    legend_handles = [
        Patch(facecolor="#5b9bd5", hatch="///", edgecolor="white", label="① 横断面诊断（文献共识）"),
        Patch(facecolor="#e8a33d", label="② 全人群恢复预测（矩阵协议）"),
        Patch(facecolor="#c0392b", label="③ 基线CST4亚组（真任务 ≈ 随机）"),
    ]
    ax.legend(handles=legend_handles, loc="upper center",
              bbox_to_anchor=(0.5, 1.10), ncol=3,
              frameon=False, fontsize=9.5)

    ax.set_title("三层信号框架：菌群–结局关联的可预测性上限",
                 fontsize=13, fontweight="bold", pad=36)
    fig.text(0.5, -0.02,
             "误差棒 = 队列内 CV / 外推 20 种子的 SD。Gustin 全人群 = 全 CST4（专门治疗队列，33/33）。",
             ha="center", fontsize=8.5, color="#555555")
    fig.subplots_adjust(top=0.80, bottom=0.14, left=0.07, right=0.97)
    fig.savefig(OUT + "fig_three_tier_auc.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------------------
# 主汇总表（Markdown + CSV）
# ---------------------------------------------------------------------------
def build_master_summary():
    lines = []
    a = lines.append
    a("# #C4 全结果总表（收口版）")
    a("")
    a("> 口径出处：`<OUT_DIR>/modeling/results_*.csv` + `paper_framework_v2.md §5`。")
    a("> 生成时间：2026-08-30。矩阵协议口径 = 三队列交集 83 属 + 20 种子内层 CV 选 C。")
    a("")

    a("## 一、三层信号框架（论文核心）")
    a("")
    a("| 层次 | 任务 | 代表 AUC | 说明 |")
    a("|---|---|---|---|")
    a("| ① 横断面诊断 | 你现在是不是 BV | **0.9+** | 菌群现状高度区分（文献一致，非本研究新结果） |")
    a("| ② 全人群恢复预测 | 会不会康复 | **~0.53–0.63**（Xiao 内 0.626 / →Mollin 0.528 / →Gustin 0.556） | 3 条外推全 >0.5 但信号弱，主要来自 LB 标签耦合 |")
    a("| ③ 基线 CST4 亚组（真任务） | 同为失调，谁更难恢复 | **≈0.5**（Xiao 0.535 / →Mollin 0.505 / →Gustin 0.512） | 三队列稳健复现 ≈ 随机 = 可预测性上限 |")
    a("")

    a("## 二、四队列矩阵（`results_four_cohort_matrix.csv`）")
    a("")
    a("| 分析 | AUC | SD | n | 注 |")
    a("|---|---|---|---|---|")
    a("| A1 Xiao 队列内 CV（全人群） | 0.626 | 0.163 | 68 | 锚点 0.65 |")
    a("| A2 Xiao→Mollin 外推（全人群） | 0.528 | 0.010 | 20 | 锚点 0.609 |")
    a("| A3 Xiao→Gustin 外推（全人群） | 0.556 | 0.021 | 20 | 原假说第 3 验证 |")
    a("| A4 Gustin 队列内 CV（全人群） | 0.436 | 0.172 | 33 | Gustin 全人群=全 CST4 |")
    a("| B1 Xiao CST4 队列内 CV | 0.535 | 0.127 | 54 | 锚点 0.596 |")
    a("| B2 Xiao→Mollin CST4 外推 | 0.505 | 0.020 | 20 | 锚点 0.527 |")
    a("| B3 Xiao→Gustin CST4 外推 | 0.512 | 0.000 | 20 | 新叙事第 3 复现 |")
    a("| B4 Gustin CST4 队列内 | 0.437 | 0.172 | 33 | =A4 |")
    a("| C1 Gustin α-only 队列内 | 0.689 | 0.201 | 33 | richness 信号（被 shannon 稀释） |")
    a("")

    a("## 三、richness–结局 Spearman ρ（三队列五亚组，含 rarefy 校正）")
    a("")
    a("| 亚组 | n | ρ_raw | p_raw | ρ_rarefy | p_rarefy | 深度 | 结论 |")
    a("|---|---|---|---|---|---|---|---|")
    for label, n, r_raw, p_raw, r_rarefy, p_rarefy, depth in RICH:
        verdict = "**显著·稳健**" if (p_raw < 0.05 and p_rarefy < 0.05) else "不显著"
        a(f"| {label} | {n} | {r_raw:+.3f} | {p_raw:.3f} | {r_rarefy:+.3f} | {p_rarefy:.3f} | {depth} | {verdict} |")
    a("")
    a("- **Gustin**：−0.454(p=0.008)→rarefy −0.392(p=0.024)，方向不变、显著性保持 → 排除深度混杂。bootstrap 95%CI [−0.713,−0.213] 不跨 0、LOO 0/33 失显著、置换 p=0.005 → 内部稳健。")
    a("- **Xiao/Plummer**：方向相反（Xiao 全 −0.054、Plummer 全 +0.100）且均不显著，rarefy 后仍不显著 → 未跨队列复现。")
    a("")

    a("## 四、结局精细化（连续结局 fuLB / ΔLB）")
    a("")
    a("| 队列 | n | ρ_rich·fuLB | p | ρ_rich·ΔLB | p | ρ_baseLB·fuLB | p |")
    a("|---|---|---|---|---|---|---|---|")
    a("| Xiao | 68 | −0.011 | 0.932 | +0.221 | 0.070 | +0.374 | 0.002 |")
    a("| Plummer | 49 | +0.102 | 0.486 | +0.286 | 0.046 | +0.046 | 0.752 |")
    a("| Gustin | 33 | −0.565 | 0.001 | −0.557 | 0.001 | — | — |")
    a("")
    a("- Gustin 剂量–反应增量：基线 richness 越高 → 随访 LB 越低 / ΔLB 越负（恢复越差），−0.565/−0.557 强于二元结局的 −0.454。")
    a("- rarefy 版（`results_outcome_refinement_rarefy.csv`）：Gustin ρ_vs_fuLB −0.378(p=0.030)、ρ_vs_ΔLB −0.383(p=0.028) 仍稳健。")
    a("")

    a("## 五、档 C 二元复发（方案 A，3 例；方案 C 已终止）")
    a("")
    a("| 队列 | 结局 | n | 事件 | min AUC | ρ_rich | p | AUC_LB | AUC_rich | AUC_属 | AUC_种 |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    a("| Plummer(592384) | 1mo 复发 | 68 | 3 | 0.921 | +0.119 | 0.334 | 0.503 | 0.667 | 0.500 | 0.500 |")
    a("")
    a("- 方案 C 终止原因：① 样本级「1 个月复发」标签不可得（ENA/BioSample 空、补充材料仅汇总）；② 即便取到 StepUp 也仅 +4 例（论文 Table 2：StepUp 4 例 + SToPBV 12 例 = 16），功效仍崩塌。")
    a("")

    a("## 六、功效分析（`results_power_analysis.csv`）")
    a("")
    a("| 队列 | n | 事件率 | min detectable AUC | target 0.65 需 n | target 0.70 需 n |")
    a("|---|---|---|---|---|---|")
    a("| Xiao CST4 | 54 | 0.50 | 0.711 | 110 | 62 |")
    a("| Plummer CST4 | 43 | 0.91 | 0.869 | 318 | 178 |")
    a("| Gustin | 33 | 0.36 | 0.778 | 120 | 67 |")
    a("")
    a("- **异质性为真实而非功效差异**：Xiao 功效足够（检出 0.711）而信号为零 → 真无信号；Plummer 设计上几乎检不出中等效应；Gustin 信号 0.786 恰卡检出线 → 解释仅 Gustin 报出信号。")
    a("")

    a("## 七、分辨率对照（属级 vs 种级）与 D7 动态（要点）")
    a("")
    a("- **种级不优于属级**：Xiao 属 0.621/种 0.575、CST4 属 0.511/种 0.492；Gustin 属 0.527/种 0.499。种级跨平台外推 0.582 > 属级 0.556（方法学点）。")
    a("- **D7 动态未复现**：Gustin D0 0.530 / D7 0.483（最低）/ D0+D7 0.586 / Δ 0.450；短暂清除陷阱 PPV 0.35。")
    a("")

    a("## 八、一句话结论")
    a("")
    a("> 治疗前 16S **构成**快照在临床治疗对象（基线 CST IV 亚组）内接近随机（三队列 0.51–0.60，CI 含 0.5），且不被分辨率（种级）或时间点（D7）突破 = 结局定义层面的可预测性上限；**失调总量（richness）** 仅在 Gustin 队列观察到关联且内部稳健（rarefaction + bootstrap/LOO/置换三重验证），但未跨队列复现 → 队列特异信号，须独立验证。")
    a("")
    return "\n".join(lines)

def build_master_csv():
    rows = []
    rows.append(["section", "analysis", "auc_mean", "auc_sd", "n", "note"])
    four = [
        ("四队列矩阵", "A1 Xiao 队列内CV(全)", 0.6256, 0.1633, 68, "锚点0.65"),
        ("四队列矩阵", "A2 Xiao→Mollin(全)", 0.5278, 0.0101, 20, "锚点0.609"),
        ("四队列矩阵", "A3 Xiao→Gustin(全)", 0.5564, 0.0206, 20, "第3验证"),
        ("四队列矩阵", "A4 Gustin 队列内CV(全)", 0.4355, 0.1716, 33, "全CST4"),
        ("四队列矩阵", "B1 Xiao CST4 队列内", 0.5353, 0.1267, 54, "锚点0.596"),
        ("四队列矩阵", "B2 Xiao→Mollin CST4", 0.5049, 0.0204, 20, "锚点0.527"),
        ("四队列矩阵", "B3 Xiao→Gustin CST4", 0.5124, 0.0, 20, "第3复现"),
        ("四队列矩阵", "B4 Gustin CST4", 0.437, 0.1718, 33, "=A4"),
        ("四队列矩阵", "C1 Gustin α-only", 0.6885, 0.201, 33, "richness信号"),
    ]
    for r in four:
        rows.append(list(r))
    for label, n, r_raw, p_raw, r_rarefy, p_rarefy, depth in RICH:
        rows.append(["richness", f"{label} ρ_raw", r_raw, None, n, f"p={p_raw}"])
        rows.append(["richness", f"{label} ρ_rarefy", r_rarefy, None, n, f"p={p_rarefy} depth={depth}"])
    import csv
    with open(OUT + "results_c4_master_summary.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerows(rows)

if __name__ == "__main__":
    plot_forest()
    plot_three_tier()
    md = build_master_summary()
    with open(OUT + "results_c4_master_summary.md", "w", encoding="utf-8") as f:
        f.write(md)
    build_master_csv()
    print("DONE")
    print(" - fig_richness_forest.png")
    print(" - fig_three_tier_auc.png")
    print(" - results_c4_master_summary.md")
    print(" - results_c4_master_summary.csv")
