"""Figure S2: brain-cord coupling by phenotype (needs bayes.R)."""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from common import COLOUR, OUT, WIDTH, save, style

style()
post = pd.read_csv(f"{OUT}/coupling_posteriors.csv")
post = post[(post.brain == "fdc") & (post.cord == "fa")]
me = pd.read_csv(f"{OUT}/coupling_marginal_effects.csv")
me = me[(me.brain == "fdc") & (me.cord == "fa")]
GROUPS = ["classic-ALS", "slow-ALS", "UMN-dominant", "PLS"]
ORDER = ["PLS", "UMN-dominant", "slow-ALS"]

fig, (a, b) = plt.subplots(1, 2, figsize=(WIDTH, 3.6), gridspec_kw={"width_ratios": [1.0, 1.05]})
for g in GROUPS:
    t = me[me.group == g].sort_values("brain_z")
    a.plot(t.brain_z, t.fit, color=COLOUR[g], lw=2.4, zorder=3)
    a.fill_between(t.brain_z, t.lo, t.hi, color=COLOUR[g], alpha=0.13, lw=0)
a.set_xlabel("Brain CST FDC (z)")
a.set_ylabel("Cord lateral CST FA (z)")
a.set_title("a. Marginal effects", loc="left", pad=6)
for i, g in enumerate(ORDER):
    r = post[post.variable == f"b_brain_z:group{g.replace('-', 'M')}"].iloc[0]
    y = len(ORDER) - 1 - i
    b.plot([r["q2.5"], r["q97.5"]], [y, y], color=COLOUR[g], lw=2.2)
    b.plot(r.q50, y, "o", color=COLOUR[g], ms=9, markeredgecolor="black", markeredgewidth=0.7)
    b.annotate(f"pd = {r.pd:.2f}", (r["q97.5"], y), xytext=(6, 0), textcoords="offset points",
               va="center", fontweight="bold", color=COLOUR[g])
b.axvline(0, ls="--", color="grey", lw=1)
b.set_yticks(range(len(ORDER)), ORDER[::-1])
b.set_xlabel("Coupling slope difference\nvs classical ALS")
b.set_title("b. Slope differences", loc="left", pad=6)
b.set_xlim(-1.6, 2.1)
leg = fig.legend(handles=[Line2D([0], [0], color=COLOUR[g], lw=2.6) for g in GROUPS], labels=GROUPS,
                 loc="lower center", ncol=4, frameon=False, title="Phenotype", columnspacing=1.8,
                 handletextpad=0.6, bbox_to_anchor=(0.5, -0.06))
leg.get_title().set_fontweight("bold")
plt.tight_layout()
save(fig, "FigS2")
