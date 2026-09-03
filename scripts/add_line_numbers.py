# -*- coding: utf-8 -*-
"""给 docx 的每个 section 注入行号（JCM/ASM 硬性要求：无行号会被编辑退回）。
OOXML 中 <w:lnNumType> 须按顺序放在 <w:pgMar> 之后、<w:pgNumType> 之前。
用法：<venv>/Scripts/python.exe add_line_numbers.py <docx路径>
"""
import sys
from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

# sectPr 子元素规范顺序（lnNumType 前的关键节点）
_ORDER_BEFORE = [
    "w:headerReference", "w:footerReference", "w:footnotePr", "w:endnotePr",
    "w:type", "w:pgSz", "w:pgMar", "w:paperSrc", "w:pgBorders",
]
# lnNumType 之后应接这些（用于确定插入锚点）
_ORDER_AFTER = ["w:pgNumType", "w:cols", "w:formProt", "w:vAlign", "w:noEndnote",
                "w:titlePg", "w:textDirection", "w:bidi", "w:rtlGutter",
                "w:docGrid", "w:printerSettings"]


def add_line_numbers(docx_path: str) -> int:
    doc = Document(docx_path)
    inserted = 0
    for section in doc.sections:
        sectPr = section._sectPr
        # 已存在则跳过
        if sectPr.find(qn("w:lnNumType")) is not None:
            continue
        lnNumType = OxmlElement("w:lnNumType")
        lnNumType.set(qn("w:countBy"), "1")
        lnNumType.set(qn("w:restart"), "continuous")

        # 找锚点：lnNumType 应插在第一个 "之后" 节点前；找不到就 append
        anchor = None
        for tag in _ORDER_AFTER:
            node = sectPr.find(qn(tag))
            if node is not None:
                anchor = node
                break
        if anchor is not None:
            anchor.addprevious(lnNumType)
        else:
            sectPr.append(lnNumType)
        inserted += 1

    doc.save(docx_path)
    return inserted


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: add_line_numbers.py <docx>")
        sys.exit(2)
    n = add_line_numbers(sys.argv[1])
    print(f"行号已注入 {n} 个 section")
