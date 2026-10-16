"""Table 1: participant characteristics at enrolment."""
import pandas as pd

from common import DATA, GROUPS, OUT, sessions

d = sessions().merge(pd.read_csv(f"{DATA}/demographics.csv"), on=["participant", "session"])
base = d[d.session == 1]


def median_iqr(s, dec=1):
    s = s.dropna()
    return f"{s.median():.{dec}f} [{s.quantile(.25):.{dec}f}, {s.quantile(.75):.{dec}f}]"


def count(s, values):
    s = s.dropna()
    n = s.astype(str).str.strip().str.lower().isin([v.lower() for v in values]).sum()
    return f"{n} ({100 * n / len(s):.0f}%)"


ROWS = {"n": lambda b: str(len(b)),
        "Sex, female": lambda b: count(b.sex, ["Female"]),
        "Age at baseline scan (years)": lambda b: median_iqr(b.age_at_scan),
        "Disease duration at scan (months)": lambda b: median_iqr(b.months_since_onset),
        "Diagnostic delay (months)": lambda b: median_iqr(b.diagnostic_delay_months),
        "King's stage at baseline": lambda b: median_iqr(pd.to_numeric(
            b.kings_stage.astype(str).str.extract(r"(\d)")[0], errors="coerce"), 0),
        "Onset region, limb": lambda b: count(b.onset_location, ["LL", "UL"]),
        "Onset region, bulbar": lambda b: count(b.onset_location, ["bulbar"]),
        "ALSFRS-R total (/48)": lambda b: median_iqr(b.alsfrsr_total),
        "ALSFRS-R bulbar (/12)": lambda b: median_iqr(b.alsfrsr_bulbar),
        "Penn UMN Scale total (/32)": lambda b: median_iqr(b.pumns_total),
        "ECAS total (/136)": lambda b: median_iqr(b.ecas_total),
        "Deceased at analysis": lambda b: str(int(b.deceased.sum())),
        "Months from onset to death or last follow-up": lambda b: median_iqr(b.months_onset_to_death_or_followup),
        "Sessions contributed": lambda b: str(int(d.participant.isin(b.participant).sum()))}

cols = {g: base[base.phenotype == g] for g in GROUPS} | {"All": base}
pd.DataFrame({k: [f(b) for f in ROWS.values()] for k, b in cols.items()}, index=list(ROWS)).to_csv(f"{OUT}/table1.csv")
