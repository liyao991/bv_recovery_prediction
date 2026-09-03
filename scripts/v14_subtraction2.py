# -*- coding: utf-8 -*-
"""
v1.4 减法（第二轮）：把「作者讲不清」的 3 个统计点再收窄
- #9  双重校正（Bonferroni 0.05/120 + BH-FDR）→ 单一 BH-FDR
- #19 交互项似然比检验 → 三队列方向对比（白话）
- #25 「赢家诅咒」术语 → 白话解释
带断言：每处 old 必须恰好出现 1 次，否则退出不写。
"""
import io, sys

path = r"<REPO_DIR>/results/线B论文中文草稿_v1_2026-09-01.md"

with io.open(path, encoding="utf-8") as f:
    s = f.read()

edits = [
    # 1. #9 双重校正 → 单一 FDR（Methods 2.5）
    (
        "多重比较校正：α×分辨率 120 项矩阵采用 Bonferroni（0.05/120）与 BH-FDR 双重校正评估。",
        "多重比较校正：α×分辨率 120 项矩阵采用 BH-FDR 校正评估。",
    ),
    # 2. #19 交互项（Methods 2.5）→ 删除该句
    (
        "跨队列异质性检验：cohort×richness 交互项似然比检验。软件：",
        "软件：",
    ),
    # 3. #19 交互项（Results 3.5）→ 三队列方向对比白话
    (
        "跨队列异质性检验显著（cohort×richness 交互项似然比 9.387，p=0.0092）。",
        "三队列方向不一致（Gustin 显著负相关，Xiao/Plummer 不显著且方向相反），提示关联存在队列特异性。",
    ),
    # 4. #25 赢家诅咒 → 白话
    (
        "真实效应或受赢家诅咒高估。",
        "真实效应可能被高估（多队列中唯一显著的信号，天然存在「挑出的最显著结果偏大」的偏倚）。",
    ),
]

ok = True
for i, (old, new) in enumerate(edits, 1):
    n = s.count(old)
    if n != 1:
        print(f"[FAIL] edit {i}: expected 1 match, got {n}: {old[:50]}...")
        ok = False
    else:
        s = s.replace(old, new)
        print(f"[OK]   edit {i}: replaced")

if not ok:
    print("ABORT: some edits did not match exactly; nothing written.")
    sys.exit(1)

with io.open(path, "w", encoding="utf-8") as f:
    f.write(s)

print(f"\nDone. {len(edits)} edits applied.")
