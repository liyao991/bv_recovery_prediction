# -*- coding: utf-8 -*-
"""框架 v3.2：剥离工作日志痕迹（状态行重写 + 删"下一步"待办 + 清理内部标签）"""
import io, sys

PATH = r"<REPO_DIR>\results\paper_framework_v2.md"

with io.open(PATH, "r", encoding="utf-8") as f:
    text = f.read()

# ---- A. 行级结构处理 ----
lines = text.split("\n")

# 1. 标题
assert lines[0].startswith("# 线B 论文框架"), "标题行异常: %r" % lines[0]
lines[0] = "# 线B 论文框架（v3.2 定稿）"

# 2. 状态行（index 2=工作区, 3=状态超长, 4=叙事核心）-> 简洁 header
new_header = [
    "> 版本 v3.2（2026-09-01）· 中文写作稿 · 待据此撰写论文正文（IMRaD）",
    "> 结局定位：菌群学康复（随访乳杆菌优势度 ≥50%，作为临床治愈的生态代理，非临床硬结局）",
    "> **叙事核心**：治疗前横断面 16S 菌群构成在临床治疗对象（基线 CST IV 亚组）内接近随机，且不被分辨率提升（种级）或时间点后移（D7）改善——指向结局的生物学决定因素；失调总量（richness）与结局关联呈队列异质性（Gustin 内部稳健、跨队列未复现），构成「单队列内部稳健、跨队列消失」的可重复性危机实证。",
]
assert lines[2].startswith("> **工作区**"), "工作区行异常: %r" % lines[2][:30]
assert lines[3].startswith("> **状态**"), "状态行异常: %r" % lines[3][:30]
assert lines[4].startswith("> **叙事核心**"), "叙事核心行异常: %r" % lines[4][:30]
lines[2:5] = new_header

# 3. 删"下一步"待办章节（连前面的 --- 分隔线）
cut = None
for i, l in enumerate(lines):
    if l.strip() == "## 下一步":
        j = i - 1
        while j >= 0 and lines[j].strip() == "":
            j -= 1
        cut = j if (j >= 0 and lines[j].strip() == "---") else i
        break
if cut is not None:
    lines = lines[:cut]

text = "\n".join(lines)

# ---- B. 内部标签清理（保留科学内容，去掉 #C 系列/回传/拍板标签）----
edits = [
    ("## 1. 标题（✅ 定稿 A · 2026-08-30 用户拍板）", "## 1. 标题（定稿 A）"),
    ("**✅ 已定稿（2026-08-30 用户拍板）**：选项 A", "**✅ 已定稿**：选项 A"),
    ("**richness 定义（#C13 补充）**", "**richness 定义**"),
    ("（#C7 靶标组合全集口径）", "（靶标组合全集口径）"),
    ("**特征空间口径（#C9 补充·审稿人必问）**", "**特征空间口径（审稿人必问）**"),
    ("**实测验证（2026-08-30 #C6 回传，支撑“刻意不合并”的方法学立场）**", "**实测验证（支撑“刻意不合并”的方法学立场）**"),
    ("**结果（Gustin，33 例，2026-08-29 回传）**", "**结果（Gustin，33 例）**"),
    ("**结果（2026-08-29 回传 · α 家族 × 分辨率矩阵，120 项 × 20 种子 CV）**", "**结果（α 家族 × 分辨率矩阵，120 项 × 20 种子 CV）**"),
    ("；#C17 独立复现逐位一致", ""),
    ("**队列异质性警示（#B2 + #C2 回传，2026-08-30 定稿）**", "**队列异质性警示**"),
    ("**脆弱性模拟（#C2）——结果指向稳健而非脆弱**", "**脆弱性模拟——结果指向稳健而非脆弱**"),
    ("**功效分析（#C2）——异质性为真实而非功效差异**", "**功效分析——异质性为真实而非功效差异**"),
    ("### 5.8 结局精细化：二元 vs 连续（#C3/#C4，2026-08-30 定稿）", "### 5.8 结局精细化：二元 vs 连续"),
    ("【已回传，见下方②分辨率对照】", ""),
    ("（#B2，深度混杂排除）", "（深度混杂排除）"),
    ("（#C2，小样本噪声排除）", "（小样本噪声排除）"),
    ("**richness 队列效应三象限（#C6 ② 四队列版，方向性观察，措辞谨慎）**", "**richness 队列效应三象限（四队列版，方向性观察）**"),
]

failed = []
for old, new in edits:
    n = text.count(old)
    if n < 1:
        failed.append(old)
    else:
        text = text.replace(old, new)
        print("OK  (%d 处) %s" % (n, old[:40]))

if failed:
    print("\n[WARN 未匹配]")
    for old in failed:
        print("  -", old[:60])

with io.open(PATH, "w", encoding="utf-8") as f:
    f.write(text)

print("\n工作日志痕迹剥离完成，已写回文件。")
