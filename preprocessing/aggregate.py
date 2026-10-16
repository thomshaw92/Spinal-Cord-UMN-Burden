"""Collect the per-session SCT tables into one long CSV per measure."""
import glob
import os
import sys

import pandas as pd

deriv = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ["BIDS_DIR"], "derivatives", "sct")
out = sys.argv[2] if len(sys.argv) > 2 else "sct_tables"
os.makedirs(out, exist_ok=True)
for name in ["csa_T2w", "gma_T2starw", "cord_csa_T2starw", "wm_csa_T2starw", "mtr_tracts",
             "dti_FA_tracts", "dti_MD_tracts", "dti_AD_tracts", "dti_RD_tracts"]:
    parts = []
    for f in sorted(glob.glob(f"{deriv}/sub-*/ses-*/*/{name}_perlevel.csv")):
        ses_dir = os.path.dirname(os.path.dirname(f))
        parts.append(pd.read_csv(f).assign(sub=os.path.basename(os.path.dirname(ses_dir)),
                                           ses=os.path.basename(ses_dir)))
    if parts:
        pd.concat(parts).to_csv(f"{out}/{name}.csv", index=False)
        print(f"{name}: {len(parts)} sessions")
