# -*- coding: utf-8 -*-
"""把 paper_framework_v2.md 转成审阅稿 docx。
用法：<venv>/Scripts/python.exe build_framework_docx.py <md路径> <输出docx路径>
依赖：markdown, beautifulsoup4（venv 内）
"""
import sys
import markdown

CSS = """
<style>
body { font-family: 'Microsoft YaHei', 'SimSun', sans-serif; font-size: 11pt; color: #1a1a1a; line-height: 1.6; }
h1 { font-size: 18pt; color: #1f3864; border-bottom: 2px solid #1f3864; padding-bottom: 6px; }
h2 { font-size: 15pt; color: #1f3864; border-bottom: 1px solid #b0b7c3; padding-bottom: 4px; margin-top: 18px; }
h3 { font-size: 13pt; color: #2e5395; margin-top: 14px; }
h4 { font-size: 12pt; color: #333; }
table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 10pt; }
th, td { border: 1px solid #999; padding: 5px 8px; text-align: left; vertical-align: top; }
th { background-color: #dce6f1; font-weight: bold; }
blockquote { border-left: 4px solid #2e5395; padding: 6px 12px; margin: 10px 0; background-color: #f2f6fb; color: #333; }
blockquote blockquote { background-color: #fbf3e0; border-left-color: #c08a00; }
code { background-color: #f0f0f0; padding: 1px 4px; font-family: Consolas, monospace; font-size: 9.5pt; }
pre { background-color: #f5f5f5; padding: 10px; border: 1px solid #ddd; overflow-x: auto; }
hr { border: none; border-top: 1px solid #ccc; margin: 16px 0; }
strong { color: #000; }
img { max-width: 100%; display: block; margin: 12px auto; }
p img { margin: 12px auto; }
</style>
"""

TEMPLATE = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
{css}
</head>
<body>
{body}
</body>
</html>
"""


def main():
    if len(sys.argv) < 3:
        print("用法: build_framework_docx.py <md> <out.docx>")
        sys.exit(2)
    md_path, out_docx = sys.argv[1], sys.argv[2]
    with open(md_path, encoding="utf-8") as f:
        md_text = f.read()

    body = markdown.markdown(
        md_text,
        extensions=["tables", "fenced_code", "sane_lists"],
        output_format="html5",
    )
    html = TEMPLATE.format(css=CSS, body=body)

    html_path = out_docx.rsplit(".", 1)[0] + ".html"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)
    print("HTML written:", html_path)


if __name__ == "__main__":
    main()
