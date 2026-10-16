"""Shared paths, colours, data loading and statistics."""
import os

import matplotlib as mpl
import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("DATA", os.path.join(ROOT, "data"))
OUT = os.path.join(ROOT, "output")
os.makedirs(OUT, exist_ok=True)

GROUPS = ["PLS", "UMN-dominant", "slow-ALS", "classic-ALS"]
LABEL = {"PLS": "PLS", "UMN-dominant": "UMN-dominant", "slow-ALS": "Slow-ALS", "classic-ALS": "Classic-ALS"}
SHORT = {"PLS": "PLS", "UMN-dominant": "UMN-d", "slow-ALS": "S-ALS", "classic-ALS": "C-ALS"}
COLOUR = {"PLS": "#a50026", "UMN-dominant": "#f46d43", "slow-ALS": "#313695", "classic-ALS": "#74add1"}
LIMB = {"LE": "#8c510a", "UE": "#35978f"}
WIDTH = 174 / 25.4


def sessions():
    d = pd.read_csv(f"{DATA}/sessions.csv")
    d["spectrum"] = np.where(d.phenotype.isin(GROUPS[:2]), "UMN", "ALS")
    d["pumns_ue"] = d.pumns_ue_left + d.pumns_ue_right
    d["pumns_le"] = d.pumns_le_left + d.pumns_le_right
    return d


def pallor():
    return pd.read_csv(f"{DATA}/pallor.csv")


def pallor_by_session(cells):
    return (cells.dropna(subset=["grade"]).groupby(["participant", "session"]).grade
            .agg(pallor_mean="mean", pallor_max="max").reset_index())


def with_pallor(d, cells):
    return d.merge(pallor_by_session(cells), on=["participant", "session"], how="left")


def spearman(d, x, y):
    d = d[[x, y]].dropna()
    r, p = stats.spearmanr(d[x], d[y])
    return r, p, len(d)


def partial_spearman(d, x, y, z):
    d = d[[x, y, z]].dropna()
    r = pg.partial_corr(data=d, x=x, y=y, covar=z, method="spearman")
    return r["r"].iloc[0], r["p_val" if "p_val" in r else "p-val"].iloc[0], len(d)


def mann_whitney(a, b):
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return u, p, 2 * u / (len(a) * len(b)) - 1


def participant_means(d):
    return d.groupby("participant").mean(numeric_only=True).reset_index()


def style():
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.labelsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
        "legend.title_fontsize": 10, "axes.linewidth": 0.8, "axes.spines.top": False,
        "axes.spines.right": False, "savefig.dpi": 600})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=600, bbox_inches="tight", facecolor="white")


def phenotype_legend(ax, title="Phenotype"):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", linestyle="", markersize=10, markerfacecolor=COLOUR[g],
                markeredgecolor="black", markeredgewidth=0.5, label=LABEL[g]) for g in GROUPS]
    lg = ax.legend(handles=h, loc="center", frameon=False, ncol=4, title=title,
                   handletextpad=0.5, columnspacing=2.0)
    lg.get_title().set_fontweight("bold")
