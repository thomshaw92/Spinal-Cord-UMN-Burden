"""Figure S3: per-level pallor against the upper and lower limb Penn UMN Scale scores."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import LIMB, OUT, WIDTH, save, style

style()
per = pd.read_csv(f"{OUT}/pallor_per_level.csv")
x = np.arange(len(per))
fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.8))
panels = (("le_marginal", "ue_marginal", "PUMNS-LE", "PUMNS-UE", "Spearman rho",
           "a. Per-level pallor vs regional PUMNS\n(marginal)", "*p<0.05, **p<0.01 versus zero", 0.7),
          ("le", "ue", "LE | UE", "UE | LE", "Partial Spearman rho",
           "b. Per-level pallor vs PUMNS region\n(partial, adjusted for the other region)", "*p<0.05, **p<0.01", 0.5))
for ax, (le, ue, lle, lue, ylabel, title, note, top) in zip(axes, panels):
    for k, (col, lab, c) in enumerate(((le, lle, LIMB["LE"]), (ue, lue, LIMB["UE"]))):
        bars = ax.bar(x + (k - 0.5) * 0.36, per[f"{col}_rho"], 0.36, color=c, edgecolor="black", linewidth=0.5, label=lab)
        for b, p in zip(bars, per[f"{col}_p"]):
            star = "**" if p < 0.01 else "*" if p < 0.05 else ""
            if star:
                ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.02, star, ha="center", va="bottom",
                        fontweight="bold", color=c)
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_xticks(x, per.level)
    ax.set_xlabel("Cervical level")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="normal")
    ax.set_ylim(-0.15, max(top, per[[f"{le}_rho", f"{ue}_rho"]].max().max() + 0.20))
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    ax.text(0.02, 0.98, note, transform=ax.transAxes, ha="left", va="top", style="italic", color="gray")
plt.tight_layout()
save(fig, "FigS3")
