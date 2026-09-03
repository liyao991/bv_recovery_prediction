# -*- coding: utf-8 -*-
"""v3.1 第一批落改：M3「五重」→「多重」措辞降级（机械替换）"""
import io

path = r'<REPO_DIR>\results\paper_framework_v2.md'
with io.open(path, encoding='utf-8') as f:
    txt = f.read()

repls = [
    ('五重稳健性检验', '多重稳健性检验'),
    ('五重稳健性验证', '多重稳健性检验'),
    ('五重验证', '多重稳健性检验'),
]

for old, new in repls:
    n = txt.count(old)
    txt = txt.replace(old, new)
    print(f'{old!r} -> {new!r}: {n} 处')

# 残留检查
left = txt.count('五重')
print('残留「五重」:', left)

with io.open(path, 'w', encoding='utf-8') as f:
    f.write(txt)
print('done')
