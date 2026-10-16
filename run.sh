#!/usr/bin/env bash
# ./run.sh            statistics, Figures 2-4 and S3, Supplementary Tables S1, S3, S4, S6-S8
# ./run.sh --bayes    also Table 1, the Bayesian models and Figures S1-S2
set -euo pipefail
cd "$(dirname "$0")/src"
python3 stats.py
for f in fig2 fig3 fig4 figS3; do python3 "$f.py"; done
if [ "${1:-}" = "--bayes" ]; then
  python3 table1.py
  Rscript bayes.R
  Rscript als_sensitivity.R
  python3 figS1.py
  python3 figS2.py
fi
