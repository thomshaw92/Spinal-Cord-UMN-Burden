"""Figure 4: visual lateral-column pallor."""
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from common import (COLOUR, DATA, GROUPS, LIMB, SHORT, WIDTH, pallor, phenotype_legend, save, sessions,
                    style, with_pallor)

style()
cells = pallor()
m = with_pallor(sessions(), cells)
base = m[(m.session == 1) & m.pallor_mean.notna()]
LEVELS = [f"C{i}" for i in range(1, 8)]

fig = plt.figure(figsize=(WIDTH * 0.955, 8.55))
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.42], hspace=0.58,
                      left=0.075, right=0.985, top=0.905, bottom=0.075)
row1 = gs[0].subgridspec(1, 2, wspace=0.26)
row2 = gs[1].subgridspec(1, 2, wspace=0.34)
row3 = gs[2].subgridspec(1, 2, width_ratios=[0.62, 1.0], wspace=0.14)


def scatter(ax, x, xlabel, title):
    s = m[[x, "pallor_mean", "phenotype"]].dropna()
    for g in GROUPS:
        q = s[s.phenotype == g]
        ax.scatter(q[x], q.pallor_mean, s=42, color=COLOUR[g], edgecolor="black", linewidth=0.4, alpha=0.9, zorder=3)
    r, p = stats.spearmanr(s[x], s.pallor_mean)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Mean pallor grade")
    ax.set_title(title, pad=16)
    ptxt = "p < 0.0001" if p < 1e-4 else f"p = {p:.4f}"
    ax.text(0.5, 1.015, f"rho = {r:+.2f}, {ptxt}, n = {len(s)}", transform=ax.transAxes, ha="center", va="bottom")


scatter(fig.add_subplot(row1[0, 0]), "pumns_total", "Penn UMN Scale (total)", "a. Pallor vs PUMNS")

ax = fig.add_subplot(row1[0, 1])
rng = np.random.default_rng(7)
for i, g in enumerate(GROUPS):
    v = base.loc[base.phenotype == g, "pallor_mean"].values
    ax.scatter(rng.normal(i, 0.07, len(v)), v, s=42, color=COLOUR[g], edgecolor="black", linewidth=0.4, alpha=0.9, zorder=3)
    ax.hlines(np.median(v), i - 0.28, i + 0.28, color="black", lw=1.8, zorder=4)
p = stats.mannwhitneyu(base.loc[base.spectrum == "UMN", "pallor_mean"],
                       base.loc[base.spectrum == "ALS", "pallor_mean"], alternative="two-sided").pvalue
ax.set_xticks(range(4), [SHORT[g] for g in GROUPS])
ax.set_ylabel("Mean pallor grade")
ax.set_title("b. Mean pallor by phenotype", pad=16)
ax.text(0.5, 1.015, f"UMN- vs ALS-spectrum p = {p:.3f}", transform=ax.transAxes, ha="center", va="bottom")

ax = fig.add_subplot(row2[0, 0])
for k, (col, lab, c) in enumerate((("pumns_le", "LE", LIMB["LE"]), ("pumns_ue", "UE", LIMB["UE"]),
                                   ("pumns_bulbar", "Bulbar", "#8795A6"))):
    vals, stars = [], []
    for pal in ["pallor_max", "pallor_mean"]:
        s = m[[pal, col]].dropna()
        r, p = stats.spearmanr(s[pal], s[col])
        vals.append(r)
        stars.append("**" if p < 0.01 else "*" if p < 0.05 else "")
    x = np.arange(2) + (k - 1) * 0.26
    ax.bar(x, vals, 0.26, color=c, label=lab)
    for xi, vi, si in zip(x, vals, stars):
        ax.text(xi, vi + 0.012, si, ha="center", va="bottom", fontweight="bold")
ax.set_ylim(0, 0.54)
ax.set_xlim(-0.48, 2.25)
ax.set_xticks(range(2), ["Max", "Mean"])
ax.set_xlabel("Pallor metric")
ax.set_ylabel("Spearman rho")
ax.set_title("c. Pallor somatotopy: LE > UE", pad=16)
lg = ax.legend(loc="center left", bbox_to_anchor=(0.70, 0.58), frameon=False, title="PUMNS",
               handlelength=1.1, handletextpad=0.5, labelspacing=0.6)
lg.get_title().set_fontweight("bold")

ax = fig.add_subplot(row2[0, 1])
scatter(ax, "fa", "Cervical lateral CST FA", "d. Pallor vs CST FA")
ax.set_yticks([2, 3, 4])

ax = fig.add_subplot(row3[0, 0])
c = cells.merge(sessions()[["participant", "session", "phenotype"]], on=["participant", "session"])
hm = (c[c.session == 1].dropna(subset=["grade"]).groupby(["phenotype", "level"]).grade.mean()
      .unstack().reindex(GROUPS)[LEVELS])
lo, hi = np.floor(hm.values.min() * 10) / 10, np.ceil(hm.values.max() * 10) / 10
im = ax.imshow(hm.values, cmap="Reds", vmin=lo, vmax=hi, aspect="auto")
ax.set_xticks(range(7), LEVELS)
ax.set_yticks(range(4), [SHORT[g] for g in GROUPS])
ax.set_xlabel("Vertebral level")
ax.set_title("e. Per-level mean pallor", pad=8)
fig.colorbar(im, ax=ax, fraction=0.022, pad=0.012).set_ticks(np.round(np.linspace(lo, hi, 4), 1))

ax = fig.add_subplot(row3[0, 1])
ax.axis("off")
ax.set_title("f. Lateral-column pallor examples, C2\u2013C5", loc="left", pad=8)
images = f"{DATA}/images.npz"
if os.path.exists(images):
    z = np.load(images)
    bb = ax.get_position()
    w_in, h_in = bb.width * fig.get_figwidth(), bb.height * fig.get_figheight()
    row_h = (0.825 * w_in / (40 * 4 / 26)) / h_in
    gap = max(0.02, (0.97 - 0.11 - 4 * row_h) / 3)
    for i, g in enumerate(GROUPS):
        y0 = 0.97 - row_h - i * (row_h + gap)
        ax.text(0.0, y0 + row_h / 2, f"{i + 1}. {SHORT[g]}", transform=ax.transAxes, fontweight="bold", va="center")
        a = ax.inset_axes([0.175, y0, 0.825, row_h])
        img = z[f"f4_strip{i}"]
        a.imshow(img, cmap="gray", vmin=0, vmax=1, origin="lower", interpolation="bilinear", aspect="auto")
        a.set_xticks([])
        a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
        if i == 3:
            H, W = img.shape
            for k, L in enumerate(["C2", "C3", "C4", "C5"]):
                a.text((k + 0.5) * W / 4, -0.14 * H, L, ha="center", va="top")

leg = fig.add_axes([0.02, 0.004, 0.96, 0.032])
leg.axis("off")
phenotype_legend(leg)
save(fig, "Fig4")
