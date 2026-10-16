"""Figure 2: UMN-spectrum versus ALS-spectrum imaging at enrolment."""
import os

import matplotlib.pyplot as plt
import numpy as np

from common import COLOUR, DATA, WIDTH, mann_whitney, phenotype_legend, save, sessions, style

style()
d = sessions()
base = d[d.session == 1]
PANELS = {"gma": ("Cervical grey-matter area", "Cervical grey-matter area\nC2-C5 (mm²)"),
          "fa": ("Lateral CST FA", "Lateral CST FA"),
          "fdc": ("Brain CST FDC", "Brain CST FDC"),
          "csa": ("Cord CSA C2-C3", "Cord CSA C2-C3 (mm²)"),
          "mtr": ("Lateral CST MTR", "Lateral CST MTR (a.u.)")}

fig = plt.figure(figsize=(WIDTH, 7.0))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.10], hspace=0.62, wspace=0.38)
for i, (m, (title, ylabel)) in enumerate(PANELS.items()):
    ax = fig.add_subplot(gs[i // 3, i % 3])
    b = base.dropna(subset=[m])
    _, p, r = mann_whitney(b.loc[b.spectrum == "UMN", m], b.loc[b.spectrum == "ALS", m])
    rng = np.random.default_rng(7)
    for j, s in enumerate(["UMN", "ALS"]):
        x = b[b.spectrum == s]
        ax.scatter(j + rng.uniform(-0.18, 0.18, len(x)), x[m], c=x.phenotype.map(COLOUR), s=48,
                   edgecolor="black", linewidth=0.4, alpha=0.85, zorder=3)
        ax.hlines(x[m].median(), j - 0.3, j + 0.3, color="black", linewidth=1.8, zorder=5)
    ax.set_xticks([0, 1], ["UMN-\nspectrum", "ALS-\nspectrum"])
    ax.set_ylabel(ylabel)
    ax.set_title(f"{chr(97 + i)}. {title}", pad=16)
    ax.text(0.5, 1.015, f"p = {p:.3f}, r = {r:+.2f}{' *' if p < 0.05 else ''}",
            transform=ax.transAxes, ha="center", va="bottom")

ax = fig.add_subplot(gs[1, 2])
ax.axis("off")
ax.set_title("f. Metric legend", pad=16)
images = f"{DATA}/images.npz"
if os.path.exists(images):
    z = np.load(images)
    img, (lo, hi) = z["f2_img"], z["f2_window"]
    H, W = img.shape
    colours = {"f2_gm": (1, .83, 0), "f2_rcst": (1, .5, .06), "f2_lcst": (.88, .13, .63)}
    for y0, overlay in ((0.600, False), (0.185, True)):
        a = ax.inset_axes([0, y0, 1, 0.40])
        a.imshow(img, cmap="gray", origin="lower", vmin=lo, vmax=hi, interpolation="bilinear")
        if overlay:
            for k, c in colours.items():
                rgba = np.zeros((H, W, 4))
                rgba[z[k]] = (*c, .38 if k == "f2_gm" else .75)
                a.imshow(rgba, origin="lower", interpolation="nearest")
            a.contour(z["f2_cord"].astype(float), levels=[0.5], colors=[(0, .9, 1)], linewidths=1.6)
        x0, x1 = -.5 - W * 0.34, W - .5 + W * 0.34
        a.set_xlim(x0, x1)
        a.set_ylim(-.5, H - .5)
        a.axis("off")
        a.add_patch(plt.Rectangle((x0, -.5), x1 - x0, H, facecolor="black", zorder=-5))
        for s, x, y, ha, va in (("A", W / 2, H - 1, "center", "top"), ("P", W / 2, 1, "center", "bottom"),
                                ("R", x0 + W * 0.03, H / 2, "left", "center"),
                                ("L", x1 - W * 0.03, H / 2, "right", "center")):
            a.text(x, y, s, color="w", fontweight="bold", ha=ha, va=va)
y = 0.02
for c, txt in (((1, .83, 0), "Grey matter: GMA"), ((1, .5, .06), "Right lateral CST"),
               ((.88, .13, .63), "Left lateral CST"), ((0, .9, 1), "Whole cord: CSA")):
    ax.add_patch(plt.Rectangle((0.02, y - 0.026), 0.070, 0.052, color=c, transform=ax.transAxes, clip_on=False))
    ax.text(0.12, y, txt, va="center", transform=ax.transAxes)
    y -= 0.092

leg = fig.add_subplot(gs[2, :])
leg.axis("off")
phenotype_legend(leg, "Phenotype (point colour)")
save(fig, "Fig2")
