# -*- coding: utf-8 -*-
"""diag_kmer.py - simulate DADA2 8-mer bootstrap matching for top ASVs."""
import csv, gzip, collections, sys

asv_file = r"<OUT_DIR>/PRJNA592384/asv_table.csv"
ref_file = r"<DATA_DIR>/silva/silva_138.1_train_set_bac_arch.fa.gz"

with open(asv_file) as f:
    rd = csv.reader(f)
    header = next(rd)
    asv_cols = header[1:]
    rows = list(rd)
tot = collections.Counter()
for r in rows:
    for i, v in enumerate(r[1:]):
        if v:
            tot[asv_cols[i]] += int(v)
seq, cnt = tot.most_common(1)[0]
print("top ASV len:", len(seq), "abundance:", cnt, flush=True)

refs = []
f = gzip.open(ref_file, "rt")
hdr = None; cur = []
for line in f:
    line = line.strip()
    if not line: continue
    if line.startswith(">"):
        if hdr is not None: refs.append((hdr, "".join(cur)))
        hdr = line[1:]; cur = []
    else:
        cur.append(line)
if hdr is not None: refs.append((hdr, "".join(cur)))
f.close()
print("refs loaded:", len(refs), flush=True)

# build 8-mer -> set of ref indices (use frozenset-free approach: dict of sets)
kmer_map = collections.defaultdict(set)
n_done = 0
for ri, (h, s) in enumerate(refs):
    for i in range(len(s) - 7):
        k = s[i:i+8]
        if "N" in k or any(c not in "ACGT" for c in k):
            continue
        kmer_map[k].add(ri)
    n_done += 1
    if n_done % 5000 == 0:
        print("  indexed", n_done, "refs, unique kmers:", len(kmer_map), flush=True)
print("total unique ACGT 8-mers in refs:", len(kmer_map), flush=True)

def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

for label, q in [("fwd", seq), ("rc", revcomp(seq))]:
    kmers = [q[i:i+8] for i in range(len(q)-7)]
    hits = [kmer_map.get(k, set()) for k in kmers]
    n_hit = sum(1 for h in hits if h)
    genus_ct = collections.Counter()
    for h in hits:
        for ri in h:
            tax = refs[ri][0].split(" ", 1)[1]
            genus = tax.split(";")[5]
            genus_ct[genus] += 1
    top_g = genus_ct.most_common(5)
    print(f"{label}: {n_hit}/{len(kmers)} kmers hit refs; top genera: {top_g}", flush=True)
print("== DIAG_DONE ==", flush=True)
