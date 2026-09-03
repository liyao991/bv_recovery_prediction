# -*- coding: utf-8 -*-
import io

path = r'<REPO_DIR>\results\paper_framework_v2.md'
with io.open(path, 'r', encoding='utf-8') as f:
    text = f.read()

edits = [
    ("Gustin 队列基线 richness≥22 种者康复率减半（待独立验证）",
     "Gustin 队列基线 richness≥22 种者康复率更低（24% vs 44%，Fisher p=0.28 不显著，描述性分层，待独立验证）"),
    ("richness≥22 种组康复率 24% vs 低组 50% → **基线越杂越难康复**",
     "richness≥22 种组康复率 24% vs 低组 44%（切点=基线 richness 中位数 22 种；Fisher 双侧 p=0.28，组间差异不显著）→ **基线越杂越难康复**（描述性分层）"),
    ("→ **功效足够而真无信号**；Plummer",
     "→ **功效足以排除 AUC≥0.71 的大效应；中等效应不能排除（见 6.5 局限）**；Plummer"),
    ("D30 康复率减半（24% vs 50%）——潜在风险分层指标（**待独立验证**）",
     "D30 康复率更低（24% vs 44%，Fisher 双侧 p=0.28，描述性分层）——潜在风险分层指标（**待独立验证**）"),
    ("Gustin 种级 richness 0.786，ρ=−0.454 (p=0.008",
     "Gustin 种级 richness 0.776，ρ=−0.454 (p=0.008"),
    ("基线 richness≥22 种者 D30 康复率减半（24% vs 50%）",
     "基线 richness≥22 种者 D30 康复率更低（24% vs 44%，Fisher p=0.28）"),
    ("Xiao 最小可检测 AUC 0.711 而观测 ≈0，功效足够而真无信号。",
     "Xiao 最小可检测 AUC 0.711 而观测 ≈0，功效足以排除 AUC≥0.71 的大效应，中等效应不能排除（见 6.5 局限）。"),
    ("有信号 0.72-0.79，预测力差异",
     "有信号 0.73-0.78，预测力差异"),
    ("D30 康复率减半（24% vs 50%）——潜在风险分层指标，但**未跨队列复现**",
     "D30 康复率更低（24% vs 44%，描述性分层）——潜在风险分层指标，但**未跨队列复现**"),
    ("（种 richness 0.786，ρ=−0.454）",
     "（种 richness 0.776，ρ=−0.454）"),
    ("1. **亚组样本量小**（Xiao CST4 n=54；Gustin n=33）——CI 上界 0.80，**不能排除中等效应，必须写局限**。但功效分析显示 Xiao 功效足够检出 AUC≥0.711，其构成层负结果为\"真无信号\"而非样本不足；Gustin n=33 最小可检测 AUC 0.778，richness 信号（0.786）恰卡检出线——**该阳性结果须经独立队列复现后方可外推**。",
     "1. **亚组样本量小**（Xiao CST4 n=54；Gustin n=33）——两条外推路径 CI 上界分别为 0.71（→Mollin）与 0.73（→Gustin）、Xiao CST4 holdout CI 上界 0.80，**中等效应（0.6–0.7）仍不能排除，必须写局限**。但功效分析显示 Xiao 功效足以排除 AUC≥0.71 的大效应，其构成层负结果为\"真无信号\"而非样本不足；Gustin n=33 最小可检测 AUC 0.778，richness 信号（0.776）恰卡检出线——**该阳性结果须经独立队列复现后方可外推**。"),
    ("另：Gustin 观测效应（0.786）恰卡最小可检测线（0.778）",
     "另：Gustin 观测效应（0.776）恰卡最小可检测线（0.778）"),
    ("故「Xiao 功效足够而真无信号」的推断",
     "故「Xiao 功效足以排除大效应」的推断"),
    ("聚焦临床治疗对象（基线 CST IV 亚组）后接近随机（AUC 0.51–0.60），且预测精度不被分辨率提升",
     "聚焦临床治疗对象（基线 CST IV 亚组）后接近随机（AUC 0.51–0.60；外推路径可排除 AUC≥0.71 的大效应、中等效应不能排除），且预测精度不被分辨率提升"),
]

for old, new in edits:
    cnt = text.count(old)
    if cnt != 1:
        print("[SKIP count=%d] %s" % (cnt, old[:50]))
    else:
        text = text.replace(old, new)
        print("[OK] %s" % old[:50])

with io.open(path, 'w', encoding='utf-8') as f:
    f.write(text)

print("DONE")
