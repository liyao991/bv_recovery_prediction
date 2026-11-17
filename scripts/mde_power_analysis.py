#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Minimum detectable AUC (MDE) for the baseline CST IV subgroups.

Solves, for each cohort,

    AUC - 0.5 = (z_{1-alpha/2} + z_{1-beta}) * SE_HM(AUC)

for AUC, where SE_HM is the Hanley-McNeil standard error evaluated at the AUC
being tested, alpha = 0.05 (two-sided) and power 1-beta = 0.80. This is the
calculation described in the Methods of the manuscript.

Output: results/results_mde_power.csv

Usage: python scripts/mde_power_analysis.py
"""

import csv
import os
from statistics import NormalDist

Z = NormalDist().inv_cdf(0.975) + NormalDist().inv_cdf(0.80)  # 1.95996 + 0.84162


def se_hanley_mcneil(auc, n_pos, n_neg):
    """Hanley-McNeil standard error of an AUC estimate.

    n_pos = number of samples with the outcome (microbiological recovery),
    n_neg = number without it.
    """
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    var = (auc * (1.0 - auc)
           + (n_pos - 1) * (q1 - auc * auc)
           + (n_neg - 1) * (q2 - auc * auc)) / (n_pos * n_neg)
    return var ** 0.5


def mde(n_pos, n_neg, lo=0.5001, hi=0.999, tol=1e-10):
    """Smallest AUC that is detectable at the configured alpha and power."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (mid - 0.5) - Z * se_hanley_mcneil(mid, n_pos, n_neg) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


# cohort, total n in the CST IV subgroup, recovered, not recovered
COHORTS = [
    ("Xiao",    54, 27, 27),
    ("Plummer", 43, 39, 4),
    ("Gustin",  33, 11, 22),
]

rows = []
for name, n, n_pos, n_neg in COHORTS:
    rows.append({
        "cohort": name,
        "n_cst4": n,
        "n_recovered": n_pos,
        "n_not_recovered": n_neg,
        "recovery_rate": "{:.3f}".format(n_pos / n),
        "mde_auc": "{:.3f}".format(mde(n_pos, n_neg)),
        "alpha": 0.05,
        "power": 0.80,
    })

out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "results_mde_power.csv")

with open(out_path, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

for r in rows:
    print("{cohort:<8} n={n_cst4:>3}  {n_recovered:>2}/{n_not_recovered:<2}  "
          "MDE={mde_auc}".format(**r))
print("\nwritten:", out_path)
