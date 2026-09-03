# build_dangC_labels.py
# 档 C 结局精细化：构建"入组后 1 个月内复发 vs 未复发"的样本级标签（Plummer PRJNA592384）
# 按 ASUS 07:58/08:19 定夺的原文口径：
#   - recurrence = 1 个月内复发（M1 复发，3 例：1038/1051/1057）
#   - cure = 达到 1 个月未复发（其余有 M0 基线者；含 M2-M5 复发者，因 1 个月时未复发）
# 输出：<OUT_DIR>/modeling/dangC_labels.csv（subject, base_sample, recurrence_label, first_recur_tp, note）
import pandas as pd, re

def norm_sid(s):
    return re.sub(r'(_F|_R1|_R2)?\.fastq\.gz$', '', str(s))

al = pd.read_csv('<DATA_DIR>/meta/PRJNA592384_run_alias.csv')
# 排除对照
CTRL = ['RNAlater', 'Taq', 'NTC']
al = al[~al['alias'].str.contains('|'.join(CTRL), na=False)]

def parse(a):
    a = str(a)
    m = re.search(r'^(\d+)_(M\d+)(.*)$', a)
    if not m:
        return None, None, ''
    suffix = m.group(3)
    recur = 'R' if ('R' in suffix) else ''
    return m.group(1), m.group(2), recur

parsed = al['alias'].apply(parse)
al['subject'] = [p[0] for p in parsed]
al['tp'] = [p[1] for p in parsed]
al['recur'] = [p[2] for p in parsed]
al = al[al['subject'].notna()].copy()

# 规范化时间点：M0A/M0B/M0_D -> M0；M1R -> M1 + R 标记；M5RV -> M5 + RV
def norm_tp(tp):
    t = str(tp)
    m = re.match(r'^(M\d+)', t)
    return m.group(1) if m else t

al['tp_norm'] = al['tp'].apply(norm_tp)

# 每个 subject 的复发时点（最早）
recur_subj = al[al['recur'] != ''].copy()
recur_subj['tp_norm'] = recur_subj['tp'].apply(norm_tp)
first_recur = recur_subj.groupby('subject')['tp_norm'].apply(lambda x: sorted(x, key=lambda t: int(t[1:]))[0])

# 每个 subject 的 M0 基线样本
m0 = al[al['tp_norm'] == 'M0'].copy()
# 标准 M0 优先排序，再 drop_duplicates 取每 subject 第一行
m0['_is_std'] = (m0['tp'] == 'M0').astype(int)
m0 = m0.sort_values('_is_std', ascending=False)
base = m0.drop_duplicates(subset='subject', keep='first').reset_index(drop=True)

rows = []
for _, r in base.iterrows():
    s = r['subject']
    fr = first_recur.get(s, '')
    # 方案 A：复发 = 最早复发时点 == M1（1 个月内）
    label = 1 if fr == 'M1' else 0
    note = f'recur@{fr}' if fr else 'never'
    rows.append(dict(subject=s, base_sample=r['run'], recurrence=label, first_recur_tp=fr, note=note))

out = pd.DataFrame(rows)
out.to_csv('<OUT_DIR>/modeling/dangC_labels.csv', index=False)

print('=== 档 C 方案 A 标签（1 个月内复发 vs 未复发）===')
print(f'总 subject（有 M0 基线）: {len(out)}')
print(f'复发组(≤M1): {out["recurrence"].sum()} 例')
print(f'治愈组(其余): {(out["recurrence"]==0).sum()} 例')
print()
print('复发组明细:')
print(out[out['recurrence']==1][['subject','base_sample','first_recur_tp']].to_string(index=False))
print()
print('治愈组中"M2-M5 复发但 1 个月未复发"的 subject（按原文归 cure）:')
m2plus = out[(out['recurrence']==0) & (out['first_recur_tp']!='')]
print(m2plus[['subject','first_recur_tp']].to_string(index=False))
print(f'共 {len(m2plus)} 例（M2-M5 复发者归入 cure）')
