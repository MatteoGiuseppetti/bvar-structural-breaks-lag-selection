# Bayesian Vector Autoregressions with Structural Breaks and Regime-Specific Lag Selection: A Marginal Likelihood Approach to Model Uncertainty

R code accompanying the master's dissertation *Bayesian Vector Autoregressions with Structural Breaks and Regime-Specific Lag Selection: A Marginal Likelihood Approach to Model Uncertainty* (Matteo Giuseppetti, Università Cattolica del Sacro Cuore, Data Analytics for Business, A.Y. 2025/2026, supervisor: Prof. Peluso).

This repository contains the code used for the simulation study and the empirical application developed in the dissertation.

## Overview

The dissertation develops a Bayesian vector autoregressive framework for identifying structural breaks while accounting for regime-specific uncertainty over the lag structure.

In particular, the framework determines:

- **where** structural breaks are located, through model-marginalized evidence over admissible break configurations;
- **which lags** are relevant in each selected regime, through posterior model probabilities over a binary lag-inclusion space.

Regime-specific marginal likelihoods are available in closed form under a Matrix-Normal-Inverse-Wishart (MNIW) prior. Uncertainty over lag configurations is handled by marginalizing over the candidate lag models using a uniform prior.

The methodology extends the structural-break framework with regime-specific model uncertainty to Bayesian vector autoregressions in which the lag configuration can differ across regimes.

For a small number of breaks, admissible break configurations can be evaluated exactly. The code also implements an exact dynamic-programming procedure and a collapsed Metropolis-within-Gibbs MCMC sampler for higher-dimensional break configurations. Whenever exact results are computationally available, they are used to validate the stochastic algorithm.

## Repository structure

```text
BVAR_thesis_full_code.R   # complete R code used in the dissertation
PCECTPI.csv               # FRED: Personal Consumption Expenditures Price Index
UNRATE.csv                # FRED: Civilian Unemployment Rate
FEDFUNDS.csv              # FRED: Effective Federal Funds Rate
README.md
```

`BVAR_thesis_full_code.R` contains the full implementation in a single script. It includes:

1. **Simulation study** (dissertation Chapters 4-5): exact break-search procedures, dynamic programming, the general MCMC sampler, and its validation.
2. **Empirical application** (dissertation Chapter 6): application to quarterly U.S. PCE inflation, unemployment, and the federal funds rate, including structural-break selection, regime-specific lag selection, diagnostics, robustness analysis, and figure generation.

## Data

The empirical application uses publicly available U.S. macroeconomic series from the Federal Reserve Economic Data (FRED) database maintained by the Federal Reserve Bank of St. Louis:

- `PCECTPI`: Personal Consumption Expenditures Price Index;
- `UNRATE`: Civilian Unemployment Rate;
- `FEDFUNDS`: Effective Federal Funds Rate.

The CSV files included in the repository contain the data used by the empirical code.

## Citation

If you use this code, please cite the dissertation:

> Giuseppetti, M. (2026). *Bayesian Vector Autoregressions with Structural Breaks and Regime-Specific Lag Selection: A Marginal Likelihood Approach to Model Uncertainty*. Master's dissertation, Università Cattolica del Sacro Cuore.

## Author

Matteo Giuseppetti  
Master's Degree in Data Analytics for Business  
Università Cattolica del Sacro Cuore
