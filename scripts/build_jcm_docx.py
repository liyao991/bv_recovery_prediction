# -*- coding: utf-8 -*-
"""把论文 md 按 JCM（Journal of Clinical Microbiology / ASM）投稿排版转成 HTML。
规范依据：JCM Instructions to Authors（ASM 系列稳定版）
  - 字号 ≥12pt（正文 12pt）
  - 全篇双倍行距（double-spaced，行距 ≥6mm），含图注、表脚注、参考文献
  - 页面 Letter（8.5×11 in），边距 1 英寸（由 build_jcm_docx.sh 传参）
  - 连续页码、行号（行号由 add_line_numbers.py 后处理注入）
  - 斜体词斜体、段落 lead-in 粗体
用法：<venv>/Scripts/python.exe build_jcm_docx.py <md路径> <输出docx路径>
"""
import sys
import markdown

# JCM 投稿排版 CSS：黑字、双倍行距、12pt、Times New Roman + 宋体、无彩色装饰
CSS = """
<style>
@page { @bottom-center { content: counter(page); } }
body {
  font-family: 'Times New Roman', 'SimSun', serif;
  font-size: 12pt;
  line-height: 2;
  color: #000000;
  text-align: justify;
}
h1 {
  font-family: 'Times New Roman', 'SimSun', serif;
  font-size: 12pt;
  font-weight: bold;
  line-height: 2;
  text-align: center;
  margin: 12pt 0;
}
h2 {
  font-family: 'Times New Roman', 'SimSun', serif;
  font-size: 12pt;
  font-weight: bold;
  line-height: 2;
  margin: 12pt 0 6pt 0;
}
h3 {
  font-family: 'Times New Roman', 'SimSun', serif;
  font-size: 12pt;
  font-weight: bold;
  font-style: italic;
  line-height: 2;
  margin: 10pt 0 4pt 0;
}
h4 {
  font-family: 'Times New Roman', 'SimSun', serif;
  font-size: 12pt;
  font-weight: bold;
  line-height: 2;
  margin: 8pt 0 2pt 0;
}
p {
  line-height: 2;
  margin: 0 0 6pt 0;
  text-indent: 0;
}
table {
  border-collapse: collapse;
  width: 100%;
  margin: 12pt 0;
  font-size: 10pt;
  line-height: 1.3;
}
th, td {
  border: 1px solid #000000;
  padding: 4px 6px;
  text-align: left;
  vertical-align: top;
}
th { font-weight: bold; }
blockquote {
  font-size: 12pt;
  line-height: 2;
  margin: 12pt 0;
  padding: 0 0 0 18pt;
  border-left: none;
  color: #000000;
}
code {
  font-family: 'Courier New', monospace;
  font-size: 11pt;
}
pre {
  font-family: 'Courier New', monospace;
  font-size: 11pt;
  line-height: 1.4;
  margin: 12pt 0;
  padding: 6pt 10pt;
  border: 1px solid #000000;
}
hr { border: none; border-top: 1px solid #000000; margin: 12pt 0; }
strong { font-weight: bold; }
em { font-style: italic; }
img { max-width: 100%; display: block; margin: 12pt auto; }
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
        print("用法: build_jcm_docx.py <md> <out.docx>")
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
