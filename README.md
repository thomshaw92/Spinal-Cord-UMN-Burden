# Spinal-Cord-UMN-Burden

Analysis code for cervical spinal cord and brain MRI, and visual lateral-column pallor, as measures
of upper motor neuron burden in motor neuron disease.

## Layout

```
preprocessing/   Spinal Cord Toolbox 6.5: disc labelling, T2w, T2*w, diffusion and MT processing
src/             statistics, figures and tables
data/            input files (see data/README.md)
output/          written by run.sh
```

## Running

Python 3.10 or later with `pip install -r requirements.txt`. The Bayesian models also need R with
`tidyverse`, `brms` and `cmdstanr`.

```
./run.sh            # statistics, Figures 2-4 and S3, Supplementary Tables S1, S3, S4, S6-S8
./run.sh --bayes    # also Table 1, Table S2 posteriors and Figures S1-S2
```

| Output | Article |
|---|---|
| `results.csv` | every statistic quoted in Results 3.2 and 3.4 to 3.6 |
| `table_S1.csv` ... `table_S8.csv` | Supplementary Tables S1, S3, S4, S6, S7, S8 |
| `pallor_per_level.csv` | Results 3.6.3 and Fig. S3 |
| `Fig2`, `Fig3`, `Fig4`, `FigS1`-`FigS3` | the figures |
| `table1.csv` | Table 1 |
| `als_sensitivity_posteriors.csv` | Results 3.4 sensitivity analysis and Supplementary Table S2D |
| `longitudinal_posteriors.csv`, `coupling_posteriors.csv` | Supplementary Table S2 |

Brain fixel-based analysis followed Chang et al., BMC Medicine 24:295 (2026); the tract-selection
code is at https://doi.org/10.18112/openneuro.ds005874.v1.1.0.

## Licence

MIT
