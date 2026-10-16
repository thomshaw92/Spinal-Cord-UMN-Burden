"""Figure 3: cord and brain imaging against the Penn UMN Scale."""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.stats import rankdata

from common import COLOUR, GROUPS, WIDTH, partial_spearman, phenotype_legend, save, sessions, style

style()
d = sessions()


def ptext(p):
    return "p < 0.0001" if p < 1e-4 else f"p = {p:.4f}"


def scatter(ax, x, y, xlabel, ylabel, title, label=None):
    for g in GROUPS:
        s = label == g
        ax.scatter(x[s], y[s], c=COLOUR[g], s=55, edgecolor="black", linewidth=0.4, alpha=0.85, zorder=3)
    slope, icept, *_ = stats.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, icept + slope * xs, color="black", linewidth=1.4, linestyle="--", alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=16)


def annotate(ax, text):
    ax.text(0.5, 1.015, text, transform=ax.transAxes, ha="center", va="bottom")


fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 5.8))
for ax, (m, ylabel, title) in zip(axes[0], (("fa", "Cervical lateral CST FA", "a. Lateral CST FA vs PUMNS"),
                                            ("fdc", "Brain CST FDC", "b. Brain CST FDC vs PUMNS"))):
    s = d[["phenotype", "pumns_total", m]].dropna()
    scatter(ax, s.pumns_total.values, s[m].values, "Penn UMN Scale (total)", ylabel, title, s.phenotype.values)
    r, p = stats.spearmanr(s.pumns_total, s[m])
    annotate(ax, f"rho = {r:+.2f}, {ptext(p)}, n = {len(s)}")

s = d[["phenotype", "fa", "pumns_total", "duration_rank"]].dropna()
rz = rankdata(s.duration_rank)


def resid(v):
    v = rankdata(v)
    slope, icept, *_ = stats.linregress(rz, v)
    return v - (slope * rz + icept)


scatter(axes[1, 0], resid(s.pumns_total), resid(s.fa), "PUMNS (rank residual after duration adjustment)",
        "Lateral CST FA (rank residual)", "c. FA-PUMNS, duration-adjusted", s.phenotype.values)
r, p, n = partial_spearman(s, "fa", "pumns_total", "duration_rank")
annotate(axes[1, 0], f"partial rho = {r:+.2f}, {ptext(p)}, n = {n}")

s = d[["phenotype", "fdc", "fa"]].dropna()
scatter(axes[1, 1], s.fdc.values, s.fa.values, "Brain CST FDC", "Cervical lateral CST FA",
        "d. Brain-cord coupling", s.phenotype.values)
r, p = stats.spearmanr(s.fdc, s.fa)
annotate(axes[1, 1], f"rho = {r:+.2f}, {ptext(p)}, n = {len(s)}")

plt.tight_layout()
leg = fig.add_axes([0, -0.035, 1, 0.04])
leg.axis("off")
phenotype_legend(leg, None)
save(fig, "Fig3")
