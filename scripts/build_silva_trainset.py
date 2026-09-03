# -*- coding: utf-8 -*-
"""
build_silva_trainset.py - Build DADA2-compatible SILVA v138.1 training set
Problem: raw SILVA NR99 trunc file has headers with 2-20 semicolon levels
(Eukaryota deep taxonomy), which breaks DADA2 assignTaxonomy
("Invalid map from references to genus").
Fix: keep only Bacteria/Archaea, standardize to 7 levels
(Domain;Phylum;Class;Order;Family;Genus;Species), pad missing with NA.
Output: gzipped FASTA, same format as DADA2 official silva_nr99_v138.1_train_set
"""
import gzip
import re

SRC = "<DATA_DIR>/silva/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
DST = "<DATA_DIR>/silva/silva_138.1_train_set_bac_arch.fa.gz"
def main():
    n_total = 0
    n_keep = 0
    n_std7 = 0
    n_short = 0
    n_long = 0
    bad = 0
    # per-genus representative selection: keep at most MAX_PER_GENUS longest seqs
    MAX_PER_GENUS = 10
    records = {}  # genus -> list of (len, hdr, seq)

    with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as fi:
        acc = None
        tax = None
        seq = []
        def flush():
            nonlocal n_total, n_keep
            if acc is None:
                return
            n_total += 1
            kd = tax.split(";")[0] if tax else ""
            if kd not in ("Bacteria", "Archaea"):
                return
            n_keep += 1
            parts = tax.split(";")
            if len(parts) == 7:
                pass
            elif len(parts) < 7:
                parts = parts + ["NA"] * (7 - len(parts))
            else:
                parts = parts[:7]
            genus = parts[5] if len(parts) >= 6 else "NA"
            seq_dna = "".join(seq).strip().replace("U", "T")
            lst = records.setdefault(genus, [])
            lst.append((len(seq_dna), acc, ";".join(parts), seq_dna))
        for line in fi:
            if line.startswith(">"):
                flush()
                hdr = line[1:].strip()
                m = re.match(r"^(\S+)\s+(.*)$", hdr)
                if not m:
                    bad += 1
                    acc = None
                    continue
                acc, tax = m.group(1), m.group(2)
                seq = []
            else:
                seq.append(line.strip())
        flush()

    # write representatives
    out_total = 0
    seen_seqs = set()
    dup_skipped = 0
    with gzip.open(DST, "wt", encoding="ascii") as fo:
        for genus, lst in records.items():
            lst.sort(key=lambda x: x[0], reverse=True)
            for item in lst[:MAX_PER_GENUS]:
                _, acc, tax_str, seq_dna = item
                # dedup identical (taxonomy, seq) pairs
                key = (tax_str, seq_dna)
                if key in seen_seqs:
                    dup_skipped += 1
                    continue
                seen_seqs.add(key)
                # header is PURE taxonomy (no accession prefix): on this R box
                # ShortRead::readFasta id() returns the WHOLE header line, so an
                # accession prefix would leak into the Kingdom column.
                fo.write(">" + tax_str + "\n")
                fo.write(seq_dna + "\n")
                out_total += 1
    print(f"total headers: {n_total}")
    print(f"kept (Bacteria/Archaea): {n_keep}")
    print(f"unique genera: {len(records)}  (NA genus: {records.get('NA', 0)} seqs kept)")
    print(f"output refs (<= {MAX_PER_GENUS}/genus, dedup): {out_total} (dup skipped: {dup_skipped})")
    print(f"output: {DST}")

if __name__ == "__main__":
    main()
