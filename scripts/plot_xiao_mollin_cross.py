# -*- coding: utf-8 -*-
"""Xiao 内部 CV vs Xiao→Mollin 外推 AUC 对比图"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# 数据（20 种子均值±SD）
groups = ['全量人群', '基线 CST4 亚组']
internal = [0.650, 0.546]
internal_sd = [0.055, 0.044]
external = [0.609, 0.527]
external_sd = [0.043, 0.027]

x = np.arange(len(groups))
w = 0.32

fig, ax = plt.subplots(figsize=(8, 5.5))
b1 = ax.bar(x - w/2, internal, w, yerr=internal_sd, capsize=4,
            color='#4C72B0', alpha=0.85, label='Xiao 内部 5 折 CV（20 种子）')
b2 = ax.bar(x + w/2, external, w, yerr=external_sd, capsize=4,
            color='#DD8452', alpha=0.85, label='Xiao 训练 → Mollin 外推（20 种子）')

ax.axhline(0.5, color='red', ls='--', lw=1.2, label='随机（AUC=0.5）')
ax.set_xticks(x)
ax.set_xticklabels(groups, fontsize=12)
ax.set_ylabel('AUC', fontsize=13)
ax.set_ylim(0.35, 0.80)
ax.set_title('跨平台外推几乎不损失性能：\n真正瓶颈是基线 CST4 亚组任务本身（内部仅 0.546）', fontsize=13)
ax.legend(fontsize=10, loc='upper right')

# 标注差值
for i in range(2):
    d = internal[i] - external[i]
    ax.annotate(f'Δ={d:.3f}', xy=(x[i], max(internal[i], external[i]) + 0.035),
                ha='center', fontsize=10, color='#333333')

# 标注 CST4 亚组
ax.annotate('基线 CST4 内（LB 已压缩）\n内部与外推均 ≈ 随机',
            xy=(1, 0.545), xytext=(0.55, 0.42),
            arrowprops=dict(arrowstyle='->', color='#888888'),
            fontsize=10, color='#555555')

plt.tight_layout()
plt.savefig('<OUT_DIR>/modeling/fig_xiao_mollin_cross.png', dpi=150)
print('已保存 fig_xiao_mollin_cross.png')
