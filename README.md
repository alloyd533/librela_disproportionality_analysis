# Bedinvetmab (Librela) Pharmacovigilance Analysis

Disproportionality signal detection and ADR rate estimation for bedinvetmab using EudraVigilance veterinary reports.

## Scripts

| Script | Purpose |
|--------|---------|
| 1_cleaning.R | Standardises EudraVigilance exports, detects polypharmacy, maps reactions to VeDDRA hierarchy |
| 2_dp_analysis.R | Calculates PRR, EBGM, and IC for each preferred term; flags consensus signals (≥2 methods) |
| 3_rate_sim_hypers.R | Simulates individual dog treatment histories to estimate exposure denominators and ADR rates with uncertainty |

## Methods

### Signal Detection

Three standard pharmacovigilance methods following EMA/FDA conventions:

- **PRR** — Proportional Reporting Ratio with chi-squared test (signal: PRR >= 2, chi-squared >= 4)
- **MGPS/EBGM** — Empirical Bayes Geometric Mean via openEBGM (signal: EB05 >= 2)
- **BCPNN/IC** — Information Component with Dirichlet-based posterior simulation (signal: IC025 > 0)

Comparators: meloxicam, carprofen, firocoxib, enflicoxib, grapiprant, robenacoxib.

### Rate Simulation

Simulates individual dog treatment trajectories (entry time via logistic uptake, survival from breed-size life tables, exponential dropout) in batches until cumulative doses match assumed global sales, then divides observed AE counts by simulated exposed population.

**Hyperparameters** sampled via Latin Hypercube (n = 1000):

| Parameter | Range | Distribution |
|-----------|-------|--------------|
| Age-weight correlation | -0.3 to 0 | Uniform |
| Monthly dropout rate | 0.5-3% | Uniform |
| Dose wastage buffer | 5-30% | Uniform |
| Underreporting multiplier | 1-20x | Beta(2.5, 7) scaled |
| Total doses sold | 30-31M | Uniform |
| Reporting rate drift | 0.7-1.3x | Uniform |

**Fixed assumptions** from literature: age at treatment as Beta(4.45, 2.69) scaled to 1-18 years, weight as LogNormal(3.26, 0.57) truncated at 100kg (both Monteiro 2025), survival from size-stratified life tables (Montoya 2023).

## Data

Place EudraVigilance exports in data/raw/ as YYYY_eudra_DRUG.csv. Place VeDDRA terminology file in data/reference/.

## Usage
```r
install.packages("renv")
renv::restore()

source("scripts/1_cleaning.R")
source("scripts/2_dp_analysis.R")
source("scripts/3_rate_simulation.R")  # ~150 min with 8 cores
```

## Dependencies

tidyverse, readxl, here, janitor, scales, gt, gtExtras, jsonlite, openEBGM, survival, truncnorm, viridis, tictoc, progressr, htmltools, future, furrr, lhs, MASS, Hmisc, ggtext, glue

## References

Montoya et al. (2023). Life expectancy tables for dogs and cats. Frontiers in Veterinary Science 10:1082102.

Monteiro et al. (2025). Global pharmacovigilance reporting of bedinvetmab. Frontiers in Veterinary Science 12:1558222.