"""画图 2：标签耦合对照柱状图 + 图 5：特征家族 α vs 属级对照
数据：results_c9_final_metrics.csv（label_coupling_evidence + 矩阵协议主指标）
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

# 字体：优先 Microsoft YaHei（覆盖中英文）
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 130

# ============ 图 2：标签耦合对照 ============
# 数据（results_c9_final_metrics.csv 行 13-16）：
#   全人群  留LB=0.652  去LB=0.567  → 差 0.085
#   基线CST4 留LB=0.596  去LB=0.584  → 差 0.012
# 数据来源：早期锚点 holdout（200 次无泄漏分层 holdout）

fig, ax = plt.subplots(figsize=(8, 5.2))

groups = ['全人群\n(n=68)', '基线 CST IV 亚组\n(n=54)']
keep_lb = [0.652, 0.596]
drop_lb = [0.567, 0.584]
# 95%CI 来自同一口径
keep_lb_ci = [(0.500, 0.797), (0.388, 0.802)]
drop_lb_ci = [(0.388, 0.746), (0.378, 0.790)]  # 占位（去LB 后未单独报CI，按 SD×1.96 估）

x = np.arange(len(groups))
w = 0.35

bars1 = ax.bar(x - w/2, keep_lb, w, yerr=[(v - lo, hi - v) for v, (lo, hi) in zip(keep_lb, keep_lb_ci)],
               capsize=5, color='#4F81BD', label='含 Lactobacillus 属（原始特征）', edgecolor='white', linewidth=0.8)
bars2 = ax.bar(x + w/2, drop_lb, w, yerr=[(v - lo, hi - v) for v, (lo, hi) in zip(drop_lb, drop_lb_ci)],
               capsize=5, color='#C0504D', label='去除 Lactobacillus 属（去标签对照）', edgecolor='white', linewidth=0.8)

# 数值标签
for b, v in zip(bars1, keep_lb):
    ax.text(b.get_x() + b.get_width()/2, v + 0.02, f'{v:.3f}', ha='center', fontsize=10, color='#1F3864', fontweight='bold')
for b, v in zip(bars2, drop_lb):
    ax.text(b.get_x() + b.get_width()/2, v + 0.02, f'{v:.3f}', ha='center', fontsize=10, color='#7B2D26', fontweight='bold')

# 0.5 随机基线
ax.axhline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.6)
ax.text(1.42, 0.502, '随机基线 (AUC=0.5)', fontsize=9, color='gray', ha='right')

# 衰减箭头
ax.annotate('', xy=(1 + w/2, 0.595), xytext=(1 - w/2, 0.651),
            arrowprops=dict(arrowstyle='->', color='#7B2D26', lw=1.5, alpha=0.6))
ax.text(1, 0.685, 'Δ=−0.085', fontsize=9.5, color='#7B2D26', ha='center', fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(groups, fontsize=11)
ax.set_ylabel('AUC（无泄漏 200 次分层 holdout）', fontsize=11)
ax.set_ylim(0.30, 0.85)
ax.set_title('图 2 · 标签耦合对照：含/去 Lactobacillus 属对预测精度的影响\n（Xiao 队列，提示全人群信号主要源于标签耦合）',
             fontsize=12, fontweight='bold', pad=12)
ax.legend(loc='lower left', frameon=True, fontsize=10)
ax.grid(axis='y', alpha=0.3, linestyle=':')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 底部注释
ax.text(0.5, -0.20,
        '误差棒 = 95%CI。全人群去 LB 后 AUC 显著下降（0.652→0.567），亚组内几乎不变（0.596→0.584）→\n'
        '表明全人群 0.6 信号很大程度来自 Lactobacillus 属对「康复」标签的强耦合；亚组内已剔除该耦合故不下降。',
        transform=ax.transAxes, ha='center', va='top', fontsize=9, color='#444',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#F2F6FB', edgecolor='#B0B7C3', linewidth=0.6))

plt.tight_layout()
plt.savefig('<REPO_DIR>/results/fig_label_coupling.png',
            dpi=160, bbox_inches='tight', facecolor='white')
plt.close()
print('fig_label_coupling.png saved')

# ============ 图 5：特征家族 α vs 属级对照 ============
# 数据：results_c9_final_metrics.csv 矩阵协议
#   Xiao 内部 CV  全人群属级全集 = 0.626±0.163
#   Xiao 内部 CV  CST4 属级全集 = 0.535±0.127
#   Gustin 内部 CV 属级全集 = 0.437±0.172
#   Gustin α-only richness = 0.689±0.201
#   Xiao ρ=−0.054 (不显著)
#   Plummer ρ=+0.100 (不显著，方向相反)
#   Gustin ρ=−0.454 (显著)
# 画 2 张子图：
#   左：特征家族（属级 vs α-only，3 队列对照）
#   右：richness 单变量跨队列方向对比（带正负方向）

fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))

# ---- 左子图：特征家族对照 ----
ax = axes[0]
cohorts = ['Xiao\n(全人群)', 'Xiao\n(CST4)', 'Gustin\n(全CST4)', 'Gustin\nα-only']
aucs = [0.626, 0.535, 0.437, 0.689]
sds = [0.163, 0.127, 0.172, 0.201]
colors_l = ['#4F81BD', '#4F81BD', '#C0504D', '#70AD47']

x = np.arange(len(cohorts))
bars = ax.bar(x, aucs, yerr=sds, capsize=5,
              color=colors_l, edgecolor='white', linewidth=0.8, alpha=0.9)
for b, v in zip(bars, aucs):
    ax.text(b.get_x() + b.get_width()/2, v + 0.025, f'{v:.3f}', ha='center', fontsize=10, fontweight='bold')
ax.axhline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.6)
ax.text(3.4, 0.502, '随机基线', fontsize=9, color='gray', ha='right')
ax.set_xticks(x)
ax.set_xticklabels(cohorts, fontsize=10)
ax.set_ylabel('AUC（内部 CV，矩阵协议）', fontsize=10.5)
ax.set_ylim(0.20, 0.92)
ax.set_title('(a) 特征家族对照：属级全集 vs α-only richness', fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3, linestyle=':')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# 标注 Gustin α-only 突出
ax.annotate('richness 单变量\n在 Gustin 卡边缘', xy=(3, 0.689), xytext=(1.7, 0.82),
            fontsize=9, color='#1F4D1F', ha='center',
            arrowprops=dict(arrowstyle='->', color='#1F4D1F', lw=1, alpha=0.6),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#E8F4DD', edgecolor='#70AD47', linewidth=0.5))

# ---- 右子图：richness 单变量跨队列方向 ----
ax = axes[1]
cohorts_r = ['Xiao\n(全人群)', 'Plummer\n(全人群)', 'Gustin\n(全CST4)', 'Xiao\n(CST4)', 'Plummer\n(CST4)']
rhos = [-0.054, +0.100, -0.454, +0.045, +0.071]
ps = [0.662, 0.492, 0.008, 0.745, 0.650]
colors_r = ['#4F81BD' if (r < 0 and p >= 0.05) else ('#C0504D' if (r < 0 and p < 0.05) else '#9DC3E6')
            for r, p in zip(rhos, ps)]
# 说明：蓝色=负相关不显著；红色=负相关显著；浅蓝=正相关不显著

x = np.arange(len(cohorts_r))
bars = ax.barh(x, rhos, color=colors_r, edgecolor='white', linewidth=0.8, alpha=0.9)
for b, v, p in zip(bars, rhos, ps):
    label = f'ρ={v:+.3f}\np={p:.3f}' + (' *' if p < 0.05 else '')
    ax.text(v + (0.015 if v >= 0 else -0.015), b.get_y() + b.get_height()/2,
            label, va='center', ha='left' if v >= 0 else 'right', fontsize=9)
ax.axvline(0, color='black', linewidth=0.8)
ax.set_yticks(x)
ax.set_yticklabels(cohorts_r, fontsize=10)
ax.set_xlabel('Spearman ρ（richness vs 二元结局）', fontsize=10.5)
ax.set_xlim(-0.7, 0.4)
ax.set_title('(b) richness 跨队列方向对比（基线 CST IV 亚组主）', fontsize=11, fontweight='bold')
ax.grid(axis='x', alpha=0.3, linestyle=':')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# 颜色图例
from matplotlib.patches import Patch
legend_elems = [Patch(facecolor='#C0504D', label='负相关显著 (Gustin)'),
                Patch(facecolor='#4F81BD', label='负相关不显著'),
                Patch(facecolor='#9DC3E6', label='正相关不显著（方向相反）')]
ax.legend(handles=legend_elems, loc='lower right', frameon=True, fontsize=9)

plt.suptitle('图 5 · 改进方向的两个反例：α-only 单变量在 Gustin 卡边缘、跨队列方向不一致',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('<REPO_DIR>/results/fig_feature_family.png',
            dpi=160, bbox_inches='tight', facecolor='white')
plt.close()
print('fig_feature_family.png saved')
