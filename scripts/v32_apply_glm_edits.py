# -*- coding: utf-8 -*-
"""框架 v3.1 -> v3.2：落 GLM 剩余措辞项 G3/G5/G6/G7/G9（精确字符串替换，带断言）"""
import io, sys

PATH = r"<REPO_DIR>\results\paper_framework_v2.md"

with io.open(PATH, "r", encoding="utf-8") as f:
    text = f.read()

# (old, new, 预期至少匹配次数)
edits = [
    # G3：LOO "0/33 失显著" 反向误读 -> 全部保持显著
    ("leave-one-out 0/33 失显著、置换 p=0.005",
     "leave-one-out 33/33 保持显著、置换 p=0.005", 1),
    ("leave-one-out ρ 范围 [−0.551, −0.461]，0/33 失显著、0/33 翻转符号（无单样本驱动）",
     "leave-one-out ρ 范围 [−0.551, −0.461]，33/33 保持显著、0/33 翻转符号（无单样本驱动）", 1),
    ("③ leave-one-out 0/33 失显著、0/33 翻转符号（无单样本驱动）",
     "③ leave-one-out 33/33 保持显著、0/33 翻转符号（无单样本驱动）", 1),
    # G6：D7 0.483 不是"全队列最低"（Δ 组 0.450 更低）
    ("Gustin D7 单独 0.483 全队列最低；短暂清除陷阱",
     "Gustin D7 单独 0.483、Δ(D7−D0) 0.450，均低于 D0 的 0.530；短暂清除陷阱", 1),
    ("Gustin D7 单独 AUC 0.483，全队列最低；D0+D7 仅微弱增益 0.586",
     "Gustin D7 单独 AUC 0.483、Δ(D7−D0) 0.450，均低于 D0 的 0.530；D0+D7 仅微弱增益 0.586", 1),
    # G7：三队列"独立"复现（循环性）-> 多队列
    ("三队列独立复现（跨平台/管线/人群）",
     "多队列（内部 CV 与两条外推路径）一致接近随机（跨平台/管线/人群）", 1),
    # G9：Plummer"仅作全人群评估"失实 -> 不纳入建模、仅用于 richness 方向
    ("Plummer 因结局极度偏斜仅作全人群评估）",
     "Plummer 因结局极度偏斜（事件率 ~90%）不纳入 AUC 建模评估，仅用于 richness 方向与生物学一致性检验）", 1),
    ("故不纳入本复现，仅在全人群评估中报告（见 5.2 报告矩阵）",
     "故不纳入本复现，仅用于 richness 方向与生物学一致性检验（见 5.2 报告矩阵）", 1),
    # G5：区间 0.51–0.60 与 0.436 冲突 -> 0.44–0.60
    ("（Xiao/Plummer/Gustin 全 0.51–0.60）",
     "（Xiao/Plummer/Gustin 全 0.44–0.60）", 1),
    ("AUC 跨队列稳定在 0.51–0.60（Xiao 内 0.535",
     "AUC 跨队列稳定在 0.44–0.60（Xiao 内 0.535", 1),
    ("亚组内接近随机（AUC 0.51–0.60），排除了",
     "亚组内接近随机（AUC 0.44–0.60），排除了", 1),
    ("接近随机（AUC 0.51–0.60；外推路径",
     "接近随机（AUC 0.44–0.60；外推路径", 1),
]

failed = []
for old, new, expect in edits:
    n = text.count(old)
    if n < expect:
        failed.append((old, "匹配 %d 次（预期 >=%d）" % (n, expect)))
    else:
        text = text.replace(old, new)
        print("OK  (%d 处) %s" % (n, old[:50]))

if failed:
    print("\n[FAILED]")
    for old, msg in failed:
        print("  -", msg, "->", old[:60])
    sys.exit(1)

with io.open(PATH, "w", encoding="utf-8") as f:
    f.write(text)

print("\n全部替换完成，已写回文件。")
