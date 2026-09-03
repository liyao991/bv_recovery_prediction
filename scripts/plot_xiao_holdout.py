# -*- coding: utf-8 -*-
"""Plot Xiao holdout AUC distribution."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "<OUT_DIR>/modeling"
res = pd.read_csv(os.path.join(OUT, "results_xiao_holdout_raw.csv"))

fig, ax = plt.subplots(figsize=(8, 5))
colors = {"LASSO-LR": "#2c3e50", "RF": "#e67e22", "XGB": "#27ae60"}
for model in ["LASSO-LR", "RF", "XGB"]:
    sub = res[res["model"] == model]["AUC"].dropna()
    ax.hist(sub, bins=12, alpha=0.5, label=f"{model} (mean={sub.mean():.3f}, sd={sub.std():.3f})",
            color=colors[model], edgecolor="black", linewidth=0.4)
ax.axvline(0.5, color="red", linestyle="--", lw=1, label="chance (0.5)")
ax.set_xlabel("Holdout test AUC (18 samples)", fontsize=11)
ax.set_ylabel("Frequency (of 50 random splits)", fontsize=11)
ax.set_title("Xiao single-cohort holdout AUC distribution (n=68, 50 train / 18 test, 50 seeds)")
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_xiao_holdout_auc_dist.png"), dpi=150)
print("saved -> fig_xiao_holdout_auc_dist.png")
