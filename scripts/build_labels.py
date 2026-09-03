# -*- coding: utf-8 -*-
"""
build_labels.py - Build CST-based outcome labels for BV recurrence prediction.
Pipeline:
  1. Load genus-level relative abundance + alpha diversity (from phyloseq)
  2. Map run -> subject & timepoint via ENA sample_alias
  3. CST typing: CST4 (dysbiosis) = Lactobacillus < LB_THRESHOLD & BV enrichment
     else CST non-4 (Lactobacillus-dominated, healthy-like)
  4. Pair baseline (D0/M0) -> follow-up (D30/M1) outcome per subject
  5. Save sample_meta.csv for modelling
"""
import os
import re
import csv
import pandas as pd
import numpy as np

META = "<DATA_DIR>/meta"
FEAT = "<DATA_DIR>/feat"
OUT = os.path.join(META, "sample_meta.csv")

BV_GENERA = ["Gardnerella", "Atopobium", "Prevotella", "Sneathia",
             "Megasphaera", "Mobiluncus", "Fannyhessea", "Dialister",
             "Prevotella_7", "Prevotella_9", "Aerococcus", "Fastidiosipila"]

def load_alias(proj):
    d = {}
    with open(os.path.join(META, f"{proj}_run_alias.csv")) as f:
        for r in csv.DictReader(f):
            d[r["run"]] = r["alias"]
    return d

def parse_xiao(alias):
    # P30D7 -> (30, D7); P118D0 -> (118, D0)
    m = re.match(r"P(\d+)D(\d+)", alias)
    if m:
        return int(m.group(1)), f"D{m.group(2)}"
    return None, None

def parse_stopbv(alias):
    # 1087_M0 -> (1087, M0); 1087_D8 -> (1087, D8)
    m = re.match(r"(\d+)_(M\d+|D\d+)$", alias)
    if m:
        return int(m.group(1)), m.group(2)
    return None, None

def cst_type(lb_pct, gard_pct, bv_pct, shannon, lb_threshold=0.5):
    """CST typing: return 1 if CST4-like (dysbiosis), 0 if Lactobacillus-dominated."""
    if lb_pct < lb_threshold:
        return 1
    return 0

def main():
    all_rows = []

    for proj, parse_fn in [("PRJNA302078", parse_xiao), ("PRJNA592384", parse_stopbv)]:
        alias_map = load_alias(proj)
        gen = pd.read_csv(os.path.join(FEAT, f"{proj}_genus_rel.csv"), index_col=0)
        alpha = pd.read_csv(os.path.join(FEAT, f"{proj}_alpha.csv"))
        alpha = alpha.set_index("sample")

        # sample names in genus table look like "SRRxxxx_F.fastq.gz"; extract run id
        def run_of(sname):
            m = re.match(r"(SRR\d+)", sname)
            return m.group(1) if m else None

        rows = []
        for sname in gen.index:
            run = run_of(sname)
            alias = alias_map.get(run, "")
            subj, tp = parse_fn(alias) if parse_fn else (None, None)
            if subj is None:
                print(f"  [warn] {proj} {sname} alias={alias!r} unparsed")
                continue
            lb = gen.loc[sname].get("Lactobacillus", 0.0)
            gard = gen.loc[sname].get("Gardnerella", 0.0)
            bv = sum(gen.loc[sname].get(g, 0.0) for g in BV_GENERA if g in gen.columns)
            sh = alpha.loc[sname, "shannon"] if sname in alpha.index else np.nan
            rows.append({
                "cohort": proj, "sample": sname, "run": run, "alias": alias,
                "subject": subj, "timepoint": tp,
                "Lactobacillus_pct": lb, "Gardnerella_pct": gard, "BV_genera_pct": bv,
                "shannon": sh,
            })

        df = pd.DataFrame(rows)
        # CST typing with threshold diagnostics
        lb_vals = df["Lactobacillus_pct"]
        print(f"=== {proj}: {len(df)} samples ===")
        print(f"  subjects: {df['subject'].nunique()}")
        print(f"  timepoints: {df['timepoint'].value_counts().to_dict()}")
        print(f"  Lactobacillus%: median {lb_vals.median():.3f}, "
              f"q25 {lb_vals.quantile(0.25):.3f}, q75 {lb_vals.quantile(0.75):.3f}")
        print(f"  samples with LB<30%: {(lb_vals<0.30).sum()}, <50%: {(lb_vals<0.50).sum()}, <70%: {(lb_vals<0.70).sum()}")

        df["cst4"] = df.apply(lambda r: cst_type(r["Lactobacillus_pct"], r["Gardnerella_pct"],
                                                 r["BV_genera_pct"], r["shannon"]), axis=1)
        all_rows.append(df)

    meta = pd.concat(all_rows, ignore_index=True)

    # ---- build baseline->outcome pairs ----
    pairs = []
    for cohort, base_tp, out_tp in [("PRJNA302078", "D0", "D30"), ("PRJNA592384", "M0", "M1")]:
        sub = meta[meta["cohort"] == cohort]
        base = sub[sub["timepoint"] == base_tp].set_index("subject")
        outc = sub[sub["timepoint"] == out_tp].set_index("subject")
        common = base.index.intersection(outc.index)
        print(f"=== {cohort}: {base_tp}->{out_tp} pairs: {len(common)} subjects ===")
        for sid in common:
            b, o = base.loc[sid], outc.loc[sid]
            pairs.append({
                "cohort": cohort, "subject": sid,
                "base_sample": b["sample"], "outcome_sample": o["sample"],
                "base_cst4": int(b["cst4"]), "outcome_cst4": int(o["cst4"]),
                "outcome": 1 if int(o["cst4"]) == 0 else 0,  # 1 = recovered (Lacto-dominated)
                "base_LB": float(b["Lactobacillus_pct"]),
                "base_Gardnerella": float(b["Gardnerella_pct"]),
                "base_BV": float(b["BV_genera_pct"]),
                "base_shannon": float(b["shannon"]) if not pd.isna(b["shannon"]) else np.nan,
            })

    pairs_df = pd.DataFrame(pairs)
    print(f"\nTOTAL pairs: {len(pairs_df)}")
    print(f"outcome distribution:\n{pairs_df['outcome'].value_counts()}")
    print(f"by cohort:\n{pairs_df.groupby(['cohort','outcome']).size()}")

    pairs_df.to_csv(OUT, index=False)
    print(f"\nsaved -> {OUT}")

    # also save full per-sample CST table for diagnostics
    meta.to_csv(os.path.join(META, "sample_cst_all.csv"), index=False)
    print("saved ->", os.path.join(META, "sample_cst_all.csv"))

if __name__ == "__main__":
    main()
