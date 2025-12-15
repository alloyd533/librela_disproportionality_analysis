\# Bedinvetmab (Librela) Pharmacovigilance Analysis

Disproportionality signal detection and ADR rate estimation for bedinvetmab using EudraVigilance veterinary reports.

\## Overview

This repository contains three analysis scripts:

\| Script \| Purpose \|

\|--------\|---------\|

\| \`1_cleaning.R\` \| Standardises EudraVigilance exports, detects polypharmacy, maps reactions to VeDDRA hierarchy \|

\| \`2_dp_analysis.R\` \| Calculates PRR, EBGM, and IC for each preferred term; flags consensus signals (≥2 methods) \|

\| \`3_rate_sim_hypers.R\` \| Simulates individual dog treatment histories to estimate exposure denominators and ADR rates with uncertainty \|

\## Methods

\### Signal Detection

Three standard pharmacovigilance methods, following FDA study:

\- \*\*PRR\*\* — Proportional Reporting Ratio with χ² test (signal: PRR ≥ 2, χ² ≥ 4)

\- \*\*MGPS/EBGM\*\* — Empirical Bayes Geometric Mean via \`openEBGM\` (signal: EB05 ≥ 2)

\- \*\*BCPNN/IC\*\* — Information Component with Dirichlet-based posterior simulation (signal: IC₀.₀₂₅ \> 0)

Comparator drugs: meloxicam, carprofen, firocoxib, enflicoxib, grapiprant, robenacoxib.

\### Rate Simulation

Simulates individual dog treatment trajectories (entry time via logistic uptake, survival from breed-size life tables, exponential dropout) in batches until cumulative doses match assumed global sales, then divides observed AE counts by simulated exposed population. Dog weights and ages were manually mapped from the values provided in Monteiro. Age \~ Beta(4.45, 2.69), Weight \~ LogNormal(3.26, 0.57)

Hyperparameters sampled via Latin Hypercube:

\| Parameter \| Range \| Distribution \|

\|-----------\|-------\|--------------\|

\| Age–weight correlation \| −0.3 to 0 \| Uniform \|

\| Monthly dropout rate \| 0.5–3% \| Uniform \|

\| Dose wastage buffer \| 5–30% \| Uniform \|

\| Underreporting rate \| 0–95% \| Beta(2.5, 7) scaled \|

\| Total doses sold \| 30–31M \| Uniform \|

\| Reporting rate drift \| 0.7–1.3× \| Uniform \|

\## Repository Structure

\`\`\`

├── scripts/

│ ├── 1_cleaning.R

│ ├── 2_dp_analysis.R

│ └── 3_rate_simulation.R

├── data/

│ ├── raw/ \# EudraVigilance exports (not tracked)

│ ├── reference/ \# VeDDRA terminology, life tables

│ └── processed/ \# Derived datasets - created by running scripts

└── output/

├── tables/ \# HTML tables (gt)

├── figures/ \# PNG plots

└── simulation/ \# Simulation results (RDS/CSV)

\`\`\`

\## Data Requirements

\| File \| Location \| Source \|

\|------\|----------\|--------\|

\| Annual EudraVigilance exports \| \`data/raw/YYYY_eudra_DRUG.csv\` \| EudraVigilance veterinary download \|

\| VeDDRA terminology \| \`data/reference/\*.xlsx\` \| EMA VeDDRA clinical terms \|

\## Usage

\`\`\`r

\# Install dependencies

install.packages("renv")

renv::restore()

\# Run sequentially

source("scripts/1_cleaning.R")

source("scripts/2_dp_analysis.R")

source("scripts/3_rate_simulation.R") \# \~150 min with 8 cores

\`\`\`

\### Configuration

Key parameters in \`3_rate_simulation.R\`:

\`\`\`r

N_SIMS \<- 1000

N_WORKERS \<- 8

CUTOFF_DATE \<- ymd("2025-06-30")

\`\`\`

\## Dependencies

tidyverse, readxl, here, janitor, scales, gt, gtExtras, jsonlite, openEBGM, survival, truncnorm, viridis, tictoc, rogressr, htmltools, future, furrr, lhs, MASS, Hmisc, ggtext, glue

\## References

\- Montoya et al. (2023). Life expectancy tables for dogs and cats. \*Frontiers in Veterinary Science\*, 10, 1082102.

\- Monteiro et al. (2025). Global pharmacovigilance reporting of bedinvetmab. \*Frontiers in Veterinary Science\*, 12, 1558222.

\- Zoetis global sales data: [zoetispetcare.com/products/librela](%5Bhttps://www.zoetispetcare.com/products/librela)](<https://www.zoetispetcare.com/products/librela>))
