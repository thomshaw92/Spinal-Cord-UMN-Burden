"""Statistics: results.csv, Supplementary Tables S1, S3, S4 and S6-S8, and per-level pallor."""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control as bh

from common import (GROUPS, OUT, mann_whitney, pallor, partial_spearman, participant_means,
                    sessions, spearman, with_pallor)

d = sessions()
cells = pallor()
d = with_pallor(d, cells)
base = d[d.session == 1]
means = participant_means(d)
rows = []


def add(section, analysis, n=np.nan, estimate=np.nan, p=np.nan, q=np.nan):
    rows.append(dict(section=section, analysis=analysis, n=n, estimate=estimate, p=p, q=q))


METRICS = {"gma": "Cervical grey-matter area, C2-C5 (mm2)", "fa": "Lateral CST FA",
           "fdc": "Brain CST FDC", "csa": "Cord CSA, C2-C3 (mm2)", "mtr": "Lateral CST MTR (a.u.)"}


def enrolment(m):
    return d[d.session == 1].dropna(subset=[m])


def first_usable(m):
    return d.dropna(subset=[m]).sort_values(["participant", "session"]).groupby("participant").head(1)


def iqr(v):
    return f"{np.median(v):.2f} [{np.percentile(v, 25):.2f}-{np.percentile(v, 75):.2f}] ({len(v)})"


s3, s7 = [], []
for rule, pick in (("enrolment", enrolment), ("first usable", first_usable)):
    block = []
    for m, label in METRICS.items():
        x = pick(m)
        a, b = x.loc[x.spectrum == "UMN", m], x.loc[x.spectrum == "ALS", m]
        u, p, r = mann_whitney(a, b)
        h, pk = stats.kruskal(*[g[m] for _, g in x.groupby("phenotype")])
        n = {g: int((x.phenotype == g).sum()) for g in GROUPS}
        block.append(dict(rule=rule, metric=label, n_umn=len(a), n_als=len(b), U=u, mw_p=p, r=r,
                          **{f"n_{g}": n[g] for g in GROUPS}, kw_H=h, kw_p=pk,
                          umn=iqr(a), als=iqr(b)))
    for row, q1, q2 in zip(block, bh([x["mw_p"] for x in block]), bh([x["kw_p"] for x in block])):
        row["mw_q"], row["kw_q"] = q1, q2
    s7 += block
s7 = pd.DataFrame(s7)
s7.drop(columns=["umn", "als"]).to_csv(f"{OUT}/table_S7.csv", index=False)
s7[s7.rule == "enrolment"][["metric", "umn", "als", "U", "mw_p", "r", "mw_q"]].to_csv(
    f"{OUT}/table_S3.csv", index=False)
for _, r in s7[s7.rule == "enrolment"].iterrows():
    add("3.2", f"Kruskal-Wallis, {r.metric}", r[[f"n_{g}" for g in GROUPS]].sum(), r.kw_H, r.kw_p, r.kw_q)
    add("3.2", f"UMN- vs ALS-spectrum, {r.metric}", r.n_umn + r.n_als, r.r, r.mw_p, r.mw_q)

x = enrolment("gma")
a, b = x.loc[x.phenotype == "UMN-dominant", "gma"], x.loc[x.phenotype != "UMN-dominant", "gma"]
u, p, r = mann_whitney(a, b)
add("3.2", f"UMN-dominant vs rest, grey-matter area (medians {a.median():.2f} vs {b.median():.2f})",
    len(x), r, p)

s1 = []
for m in ["csa", "gma", "mtr", "fa", "fdc", "fd", "logfc", "alsfrsr_total", "pumns_total"]:
    row = {"metric": m}
    for g in GROUPS + ["all"]:
        y = d[d[m].notna() & ((d.phenotype == g) | (g == "all"))]
        row[g] = f"{len(y)} ({y.participant.nunique()})"
    s1.append(row)
pd.DataFrame(s1).to_csv(f"{OUT}/table_S1.csv", index=False)

PAIRS = [("FA vs PUMNS", "fa", "pumns_total"), ("FDC vs PUMNS", "fdc", "pumns_total"),
         ("MTR vs PUMNS", "mtr", "pumns_total"), ("GMA vs PUMNS", "gma", "pumns_total"),
         ("CSA vs PUMNS", "csa", "pumns_total"), ("FA vs PUMNS-LE", "fa", "pumns_le"),
         ("FA vs PUMNS-UE", "fa", "pumns_ue"), ("FDC vs FA (brain-cord)", "fdc", "fa"),
         ("pallor_max vs FA", "pallor_max", "fa"), ("pallor_max vs FDC", "pallor_max", "fdc"),
         ("pallor_mean vs PUMNS", "pallor_mean", "pumns_total"),
         ("pallor_mean vs PUMNS-LE", "pallor_mean", "pumns_le"),
         ("pallor_mean vs PUMNS-UE", "pallor_mean", "pumns_ue")]
SCHEMES = {"pooled": (d, "duration_rank"), "participant mean": (means, "duration_mean_rank"),
           "baseline": (base, "duration_rank")}
s6 = []
for label, x, y in PAIRS:
    for scheme, (df, cov) in SCHEMES.items():
        r, p, n = spearman(df, x, y)
        row = dict(association=label, scheme=scheme, n=n, rho=r, p=p)
        if label.startswith(("FA vs P", "FDC vs P", "MTR", "GMA", "CSA", "pallor_mean")):
            row["partial_rho"], row["partial_p"], _ = partial_spearman(df, x, y, cov)
        s6.append(row)
s6 = pd.DataFrame(s6)
s6.to_csv(f"{OUT}/table_S6.csv", index=False)
s6[(s6.scheme == "pooled") & s6.association.str.endswith("vs PUMNS")
   & ~s6.association.str.startswith("pallor")].to_csv(f"{OUT}/table_S4.csv", index=False)

for m in ["fa", "fdc", "mtr"]:
    for sub in ["pumns_le", "pumns_ue"]:
        r, p, n = spearman(d, m, sub)
        pr, pp, _ = partial_spearman(d, m, sub, "duration_rank")
        add("3.4", f"{m} vs {sub} (partial {pr:+.2f}, p {pp:.3g})", n, r, p)
for m in ["csa", "fdc", "gma", "fa", "mtr"]:
    r, p, n = spearman(d, m, "alsfrsr_total")
    add("3.4", f"{m} vs ALSFRS-R total", n, r, p)
for m in ["fa", "mtr"]:
    r, p, n = spearman(d, m, "alsfrsr_bulbar")
    add("3.4", f"{m} vs ALSFRS-R bulbar", n, r, p)

for label, x, y, df in (("FDC vs FA", "fdc", "fa", d), ("FD vs FA, baseline", "fd", "fa", base),
                        ("log FC vs CSA, baseline", "logfc", "csa", base),
                        ("FDC vs MTR, baseline", "fdc", "mtr", base),
                        ("FD vs FA", "fd", "fa", d), ("log FC vs CSA", "logfc", "csa", d)):
    r, p, n = spearman(df, x, y)
    add("3.5", label, n, r, p)

q = cells.quality
add("3.6.1", "rated level-side cells", int(cells.grade.notna().sum()))
add("3.6.1", "non-diagnostic cells", int((q == "fail").sum()))
add("3.6.1", "cells outside the acquisition", int(q.isna().sum()))
pb = base.dropna(subset=["pallor_max"])
add("3.6.1", "baseline sessions with maximum grade below 2", int((pb.pallor_max < 2).sum()))
for g in GROUPS:
    add("3.6.1", f"median maximum grade at baseline, {g}", int((pb.phenotype == g).sum()),
        pb.loc[pb.phenotype == g, "pallor_max"].median())
for m in ["pallor_mean", "pallor_max"]:
    a, b = pb.loc[pb.spectrum == "UMN", m], pb.loc[pb.spectrum == "ALS", m]
    u, p, r = mann_whitney(a, b)
    add("3.6.1", f"UMN- vs ALS-spectrum, {m} (medians {a.median():.2f} vs {b.median():.2f})", len(pb), r, p)
    h, p = stats.kruskal(*[g[m] for _, g in pb.groupby("phenotype")])
    add("3.6.1", f"Kruskal-Wallis, {m}", len(pb), h, p)

grid = [(x, y) for x in ["mtr", "fa", "fdc", "pumns_total", "duration_rank"] for y in ["pallor_max", "pallor_mean"]]
res = [spearman(d, x, y) for x, y in grid]
for (x, y), (r, p, n), qv in zip(grid, res, bh([p for _, p, _ in res])):
    add("3.6.2", f"{y} vs {x}", n, r, p, qv)
for m in ["pallor_mean", "pallor_max"]:
    r, p, n = partial_spearman(d, m, "pumns_total", "duration_rank")
    add("3.6.2", f"{m} vs PUMNS, adjusted for duration", n, r, p)

s8 = []
for tag, keep in (("all rated", cells.grade.notna()),
                  ("good or acceptable", cells.grade.notna() & cells.quality.isin(["good", "acceptable"]))):
    c = cells[keep]
    x = with_pallor(sessions(), c)
    for label, a, b in (("mean pallor vs PUMNS total", "pallor_mean", "pumns_total"),
                        ("maximum pallor vs lateral CST FA", "pallor_max", "fa"),
                        ("mean pallor vs lateral CST FA", "pallor_mean", "fa")):
        r, p, n = spearman(x, a, b)
        s8.append(dict(analysis=label, subset=tag, n=n, estimate=r, p=p))
    xb = x[(x.session == 1) & x.pallor_mean.notna()]
    _, p, r = mann_whitney(xb.loc[xb.spectrum == "UMN", "pallor_mean"], xb.loc[xb.spectrum == "ALS", "pallor_mean"])
    s8.append(dict(analysis="UMN- vs ALS-spectrum, mean pallor at baseline", subset=tag, n=len(xb), estimate=r, p=p))
    c3 = (c[c.level == "C3"].groupby(["participant", "session"]).grade.max().rename("P").reset_index()
          .merge(sessions(), on=["participant", "session"]))
    r, p, n = partial_spearman(c3, "P", "pumns_le", "pumns_ue")
    s8.append(dict(analysis="C3 pallor vs PUMNS-LE, partialling PUMNS-UE", subset=tag, n=n, estimate=r, p=p))
pd.DataFrame(s8).to_csv(f"{OUT}/table_S8.csv", index=False)

lvl = (cells.dropna(subset=["grade"]).groupby(["participant", "session", "level"]).grade.max()
       .rename("P").reset_index().merge(sessions(), on=["participant", "session"]))
per = []
for L in [f"C{i}" for i in range(1, 8)]:
    x = lvl[lvl.level == L]
    le = partial_spearman(x, "P", "pumns_le", "pumns_ue")
    ue = partial_spearman(x, "P", "pumns_ue", "pumns_le")
    mle, mue = spearman(x, "P", "pumns_le"), spearman(x, "P", "pumns_ue")
    per.append(dict(level=L, n=le[2], le_rho=le[0], le_p=le[1], ue_rho=ue[0], ue_p=ue[1],
                    le_marginal_rho=mle[0], le_marginal_p=mle[1], ue_marginal_rho=mue[0], ue_marginal_p=mue[1]))
per = pd.DataFrame(per)
per["le_q"] = bh(per.le_p)
per.to_csv(f"{OUT}/pallor_per_level.csv", index=False)
for _, r in per.iterrows():
    add("3.6.3", f"{r.level} pallor vs PUMNS-LE | UE (UE | LE {r.ue_rho:+.2f}, p {r.ue_p:.2g})",
        r.n, r.le_rho, r.le_p, r.le_q)

ps = d.dropna(subset=["pallor_max"])
both = ps.groupby("participant").filter(lambda g: len(g) == 2).sort_values(["participant", "session"])
for m in ["pallor_max", "pallor_mean"]:
    w = both.pivot(index="participant", columns="session", values=m)
    add("3.6.3", f"within-participant change, {m} (Wilcoxon)", len(w), np.nan,
        stats.wilcoxon(w[1], w[2]).pvalue)

pd.DataFrame(rows).to_csv(f"{OUT}/results.csv", index=False)
print(f"wrote {len(rows)} results and tables S1, S3, S4, S6, S7, S8 to {OUT}")
