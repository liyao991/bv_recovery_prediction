# -*- coding: utf-8 -*-
"""v1.3 减法改稿：压缩 Methods 2.4/2.5 与 Results 3.5 的稳健性方法清单，删 metafor。
带 count 断言，任何一处未精确匹配即失败退出。"""
import io

path = r"<REPO_DIR>/results/线B论文中文草稿_v1_2026-09-01.md"
with io.open(path, "r", encoding="utf-8") as f:
    text = f.read()

MINUS = "\u2212"   # 负号 −
SUP2  = "\u00b2"   # 上标²
TIMES = "\u00d7"   # 乘号×
TAU   = "\u03c4"   # τ

edits = [
    # 1. Methods 2.4：删 "200 次"
    ("队列内重复分层 holdout（200 次 60/40，无泄漏）",
     "队列内重复分层 holdout（60/40，无泄漏）", 1),

    # 2. Methods 2.5：richness 稳健性 5 项 -> 3 项（删 bootstrap 1000、偏相关）
    ("richness 内部稳健性检验：rarefaction 校正（排除测序深度混杂）、bootstrap 1000 次、leave-one-out、置换检验 1000 次、混杂排除（控制基线乳杆菌丰度后计算偏相关）。",
     "richness 内部稳健性检验：rarefaction 校正（排除测序深度混杂）、leave-one-out、置换检验 1000 次。", 1),

    # 3. Methods 2.5：删 metafor 句
    ("跨队列异质性检验：cohort" + TIMES + "richness 交互项似然比检验，并以 metafor 随机效应合并作示意性敏感性分析（k=3 下 I" + SUP2 + "/" + TAU + SUP2 + " 估计功效不足，仅展示 pooled 跨 0 的事实，不作决定性论据）。",
     "跨队列异质性检验：cohort" + TIMES + "richness 交互项似然比检验。", 1),

    # 4. Methods 2.5 软件列表：删 metafor
    ("R 4.6.1（DADA2、phyloseq、metafor）",
     "R 4.6.1（DADA2、phyloseq）", 1),

    # 5. Results 3.5 第 121 行：稳健性数字压缩（删 bootstrap CI、偏相关）
    ("该信号在 Gustin 队列内部经多重稳健性检验排除小样本噪声与深度混杂：rarefaction 校正后 \u03c1=" + MINUS + "0.392（p=0.024）；bootstrap 1000 次 95%CI [" + MINUS + "0.713, " + MINUS + "0.213] 不跨 0；leave-one-out 33/33 保持显著、0/33 翻转符号（无单样本驱动）；置换检验 p=0.005；控制基线乳杆菌丰度后偏相关 \u03c1=" + MINUS + "0.431（p=0.014）几乎不衰减。",
     "该信号在 Gustin 队列内部经稳健性检验排除小样本噪声与深度混杂：rarefaction 校正后 \u03c1=" + MINUS + "0.392（p=0.024）与置换检验（p=0.005）仍稳健，33/33 leave-one-out 无单样本驱动。", 1),

    # 6. Results 3.5 第 123 行：删 metafor 句
    ("跨队列异质性检验显著（cohort" + TIMES + "richness 交互项似然比 9.387，p=0.0092）；metafor 随机效应合并 pooled \u03c1=" + MINUS + "0.110（95%CI [" + MINUS + "0.425, +0.229]，p=0.53），I" + SUP2 + "=72.3%，提示合并估计无意义。",
     "跨队列异质性检验显著（cohort" + TIMES + "richness 交互项似然比 9.387，p=0.0092）。", 1),

    # 7. 图引用笔误修正：图 2 -> 图 3（该处指 richness 森林图）
    ("**但该关联未跨队列复现**（图 2）",
     "**但该关联未跨队列复现**（图 3）", 1),

    # 8. 讨论 4.1：删"随机效应合并无意义"
    ("跨队列异质性检验显著、随机效应合并无意义。",
     "跨队列异质性检验显著。", 1),
]

fails = 0
for old, new, expect in edits:
    c = text.count(old)
    if c != expect:
        print("[FAIL] count=%d (expect %d): %s" % (c, expect, old[:45]))
        fails += 1
        continue
    text = text.replace(old, new)
    print("[OK] (%d) %s" % (c, old[:45]))

if fails:
    print("\n>>> 有 %d 处未匹配，未写入文件。" % fails)
    raise SystemExit(1)

with io.open(path, "w", encoding="utf-8") as f:
    f.write(text)
print("\n>>> 全部 %d 处替换完成并写入。" % len(edits))
