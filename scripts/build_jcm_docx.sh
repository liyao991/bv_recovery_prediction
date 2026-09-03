#!/usr/bin/env bash
# 一键：md → HTML → docx，按 JCM（ASM）投稿排版：Letter 页面 + 1 英寸边距 + 12pt + 双倍行距 + 行号
# 用法：bash build_jcm_docx.sh <md路径> <输出docx路径>
set -e

# ---- 本项目固定路径 ----
VENV="<USER_HOME>/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
PLUGIN="<USER_HOME>/.workbuddy/plugins/cache/workbuddy-builtin/tencent-docx/0.4.1"
SKILL_SCRIPTS="$PLUGIN/skills/html-to-docx/scripts"
BUILD_PY="<REPO_DIR>/scripts/build_jcm_docx.py"
LINE_NUM_PY="<REPO_DIR>/scripts/add_line_numbers.py"

# ---- 参数 ----
MD="$1"
OUT_DOCX="$2"
if [ -z "$MD" ] || [ -z "$OUT_DOCX" ]; then
  echo "用法: bash build_jcm_docx.sh <md路径> <输出docx路径>"
  exit 2
fi

# ---- POSIX 路径转 Windows 路径 ----
if command -v cygpath >/dev/null 2>&1; then
  MD=$(cygpath -m "$MD")
  OUT_DOCX=$(cygpath -m "$OUT_DOCX")
fi

# ---- 1. 依赖自检 ----
if ! "$VENV" -c "import markdown, bs4, lxml, docx, click, html4docx" >/dev/null 2>&1; then
  echo "[自愈] 检测到依赖缺失，正在补装 ..."
  "$VENV" -m pip install markdown beautifulsoup4 lxml python-docx click html-for-docx -i https://pypi.tuna.tsinghua.edu.cn/simple 2>&1 | tail -5
fi

# ---- 2. md → html（JCM CSS）----
echo "[1/3] md → html: $MD"
"$VENV" "$BUILD_PY" "$MD" "$OUT_DOCX"

# ---- 3. html → docx（Letter + 1 英寸边距）----
echo "[2/3] html → docx (Letter / 1-inch margins)"
HTML="${OUT_DOCX%.docx}.html"
PYTHONPATH="$SKILL_SCRIPTS" "$VENV" -m html_to_docx convert "$HTML" -o "$OUT_DOCX" \
  --page-size Letter \
  --margin-top 2.54 --margin-bottom 2.54 --margin-left 2.54 --margin-right 2.54

# ---- 4. 行号注入（JCM 硬性要求）----
echo "[3/3] 注入行号"
"$VENV" "$LINE_NUM_PY" "$OUT_DOCX"

echo "完成: $OUT_DOCX"
