# Bayesian DTA-NMA Code for NCCT-Based AI LVO Triage

This folder contains a cleaned RStudio-ready code package for the Bayesian diagnostic test accuracy network meta-analysis.

## Recommended Use

Open `Bayesian_DTA_NMA_RStudio.R` in RStudio and click `Source`.

The script is self-contained and expects the `data/` folder to be in the same directory. Users can control which analyses run by editing the switches at the top of the script:

- `RUN_PRIMARY_ANALYSIS`
- `RUN_META_REGRESSION`
- `RUN_SENSITIVITY_ANALYSES`
- `RUN_PRIOR_SENSITIVITY`

The MCMC settings in the script are the full analysis settings used for reproducible results with the supplied data.

## Contents

- `Bayesian_DTA_NMA_RStudio.R`: standalone RStudio script for the full analysis workflow.
- `data/`: original analysis input CSV files.

The script regenerates the primary Bayesian DTA-NMA, covariate meta-regression, SA-1 to SA-5 scenario sensitivity analyses, and the half-normal prior sensitivity analysis. SA-5 excludes the Sunwoo et al. (2026) US cohort. Half-normal prior sensitivity analysis, corresponding to SA-6.

## Requirements

Install JAGS and the R package `R2jags` before running the script.
