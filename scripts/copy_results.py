# -*- coding: utf-8 -*-
"""
copy_results.py - copy DADA2 deliverables from ASCII compute dir back to workspace
Compute runs on <OUT_DIR>/<cohort>/ (R cannot open Chinese paths);
workspace results live at <REPO_DIR>/results/<cohort>/
Usage: python copy_results.py [cohort...]   (default: all finished cohorts)
"""
import os
import shutil
import sys

SRC_ROOT = "<OUT_DIR>"
DST_ROOT = "<REPO_DIR>/results"
COHORTS = ["PRJNA302078", "PRJNA592384", "PRJNA398590"]

def main():
    targets = sys.argv[1:] if len(sys.argv) > 1 else COHORTS
    for c in targets:
        src = os.path.join(SRC_ROOT, c)
        if not os.path.isdir(src):
            print(f"[skip] {c}: no compute dir")
            continue
        dst = os.path.join(DST_ROOT, c)
        os.makedirs(dst, exist_ok=True)
        n = 0
        for fn in os.listdir(src):
            fp = os.path.join(src, fn)
            if os.path.isfile(fp):
                shutil.copy2(fp, os.path.join(dst, fn))
                n += 1
        print(f"[ok] {c}: copied {n} files -> {dst}")
    print("done")

if __name__ == "__main__":
    main()
