"""Figure S1: longitudinal lateral CST MTR and FA by phenotype (needs bayes.R)."""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from common import COLOUR, OUT, WIDTH, save, style

style()
post = pd.read_csv(f"{OUT}/longitudinal_posteriors.csv")
me = pd.read_csv(f"{OUT}/longitudinal_marginal_effects.csv")
GROUPS = ["classic-ALS", "slow-ALS", "UMN-dominant", "PLS"]
CONTRASTS = [("PLS", "PLS"), ("UMN-dominant", "UMN-dominant"), ("slow-ALS", "slow-ALS")]
TITLE = {"mtr": "Lateral CST MTR", "fa": "Lateral CST FA"}

fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 5.6), gridspec_kw={"width_ratios": [1.0, 0.92]})
for row, m in enumerate(["mtr", "fa"]):
    ax, axf = axes[row]
    t = me[(me.model == "group") & (me.metric == m)]
    for g in GROUPS:
        x = t[t.level == g].sort_values("time_mo")
        ax.plot(x.time_mo, x.fit, color=COLOUR[g], lw=2.2, zorder=3)
        ax.fill_between(x.time_mo, x.lo, x.hi, color=COLOUR[g], alpha=0.11, lw=0)
    ax.set_xlabel("Months from first scan")
    ax.set_ylabel(f"{TITLE[m]} (z)")
    ax.set_title(f"{'ac'[row]}. {TITLE[m]}: trajectories", loc="left", pad=6)
    p = post[(post.model == "group") & (post.metric == m)]
    for j, (g, _) in enumerate(CONTRASTS):
        r = p[p.variable == f"b_time_mo:group{g.replace('-', 'M')}"].iloc[0]
        y = len(CONTRASTS) - 1 - j
        excl = r["q2.5"] > 0 or r["q97.5"] < 0
        axf.plot([r["q2.5"], r["q97.5"]], [y, y], color=COLOUR[g], lw=2.8 if excl else 1.8, solid_capstyle="round", zorder=3)
        axf.plot(r["mean"], y, "o", color=COLOUR[g], ms=8 if excl else 6.5, zorder=4,
                 markeredgecolor="#3f3f3f", markeredgewidth=0.7)
        axf.annotate(f"pd {r.pd:.3f}", (r["q97.5"], y), xytext=(5, 0), textcoords="offset points",
                     va="center", color="0.25" if excl else "0.45")
    axf.axvline(0, ls="--", color="grey", lw=1)
    axf.set_yticks(range(len(CONTRASTS)), [c for _, c in CONTRASTS][::-1])
    axf.set_ylim(-0.6, len(CONTRASTS) - 0.4)
    axf.set_xlabel("Slope difference vs classical ALS\n(z-score per month)")
    axf.set_title(f"{'bd'[row]}. {TITLE[m]}: fixed effects", loc="left", pad=6)
    axf.margins(x=0.30)
leg = fig.legend(handles=[Line2D([0], [0], color=COLOUR[g], lw=2.6) for g in GROUPS], labels=GROUPS,
                 loc="lower center", ncol=4, frameon=False, title="Phenotype", columnspacing=1.8,
                 handletextpad=0.6, bbox_to_anchor=(0.5, -0.035))
leg.get_title().set_fontweight("bold")
fig.tight_layout()
save(fig, "FigS1")
