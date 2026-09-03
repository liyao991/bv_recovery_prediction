#!/usr/bin/env bash
# 一键：md → HTML → docx（自愈式：自动检查并补装依赖）
# 用法：bash build_docx.sh <md路径> <输出docx路径>
# 解决痛点：venv 依赖反复被清 + PYTHONPATH 手工设置 + 两步转换繁琐
set -e

# ---- 本项目固定路径 ----
VENV="<USER_HOME>/.venv-html-to-docx/Scripts/python.exe"
PLUGIN="<USER_HOME>/.workbuddy/plugins/cache/workbuddy-builtin/tencent-docx/0.4.1"
SKILL_SCRIPTS="$PLUGIN/skills/html-to-docx/scripts"
BUILD_PY="<REPO_DIR>/scripts/build_framework_docx.py"
export PATH="<USER_HOME>/.local/bin:$PATH"

# ---- 参数 ----
MD="$1"
OUT_DOCX="$2"
if [ -z "$MD" ] || [ -z "$OUT_DOCX" ]; then
  echo "用法: bash build_docx.sh <md路径> <输出docx路径>"
  exit 2
fi

# ---- POSIX 路径(/d/...) 转 Windows 路径(D:/...)，供 Python 使用 ----
if command -v cygpath >/dev/null 2>&1; then
  MD=$(cygpath -m "$MD")
  OUT_DOCX=$(cygpath -m "$OUT_DOCX")
fi

# ---- 1. 依赖自检，缺了就自动补装（含 html4docx=html-for-docx，勿漏；固定清华镜像加速）----
if ! "$VENV" -c "import markdown, bs4, lxml, docx, click, html4docx" >/dev/null 2>&1; then
  echo "[自愈] 检测到依赖缺失，正在补装 markdown/bs4/lxml/python-docx/click/html-for-docx ..."
  uv pip install --python "$VENV" markdown beautifulsoup4 lxml python-docx click html-for-docx -i https://pypi.tuna.tsinghua.edu.cn/simple 2>&1 | tail -5
fi

# ---- 2. md → html ----
echo "[1/2] md → html: $MD"
"$VENV" "$BUILD_PY" "$MD" "$OUT_DOCX"

# ---- 3. html → docx ----
echo "[2/2] html → docx"
HTML="${OUT_DOCX%.docx}.html"
PYTHONPATH="$SKILL_SCRIPTS" "$VENV" -m html_to_docx convert "$HTML" -o "$OUT_DOCX"

echo "完成: $OUT_DOCX"
