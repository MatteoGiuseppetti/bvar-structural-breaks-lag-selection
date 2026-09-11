# Regime-Specific Bayesian VAR with Structural Breaks and Lag-Model Uncertainty

R code accompanying the master's dissertation *"Regime-Specific Bayesian VAR
with Structural Breaks and Lag-Model Uncertainty"* (Matteo Giuseppetti,
Università Cattolica del Sacro Cuore, Data Analytics for Business, A.Y.
2025/2026, supervisor: Prof. Peluso).

This repository contains the implementation used for the simulation study and
the empirical application developed in the dissertation.

## What this is

This repository implements a Bayesian VAR framework that identifies structural
breaks while accounting for regime-specific uncertainty over the lag structure.

In particular, the framework determines:

- **where** the structural breaks are, through model-marginalized evidence over
  admissible break configurations;
- **which lags** matter in each selected regime, through posterior model
  probabilities over a binary lag-inclusion space.

The regime-specific marginal likelihoods are available in closed form under a
Matrix-Normal-Inverse-Wishart (MNIW) prior. Model uncertainty over lag
configurations is handled by marginalizing over the candidate lag models with a
uniform prior.

The methodology builds on Chib and Smith's structural-break framework with
regime-specific model uncertainty and extends it to Bayesian vector
autoregressions with regime-specific lag-model uncertainty.

For a small number of breaks, exact evaluation over admissible break
configurations is feasible. The code also implements an exact dynamic-programming
maximization procedure and a general collapsed Metropolis-within-Gibbs MCMC
sampler for higher-dimensional break configurations. The stochastic algorithm
is validated against exact results whenever the exact posterior is computationally
available.

## Repository structure

```text
BVAR_thesis_full_code.R   # single R script, run top to bottom
PCECTPI.csv               # FRED: Personal Consumption Expenditures Price Index
UNRATE.csv                # FRED: Civilian Unemployment Rate
FEDFUNDS.csv              # FRED: Effective Federal Funds Rate
README.md
```

`BVAR_thesis_full_code.R` is organized as a single linear script and contains
two main parts:

1. **Simulation study** (dissertation Chapters 4-5) -- one-break and two-break
   exact evaluations, dynamic-programming search, the general MCMC sampler, and
   simulation-based validation.
2. **Empirical application** (dissertation Chapter 6) -- application to
   quarterly U.S. PCE inflation, unemployment, and the federal funds rate,
   including residual diagnostics, structural-break selection,
   regime-specific lag selection, robustness analysis, and figure generation.

## Requirements

- R (>= 4.0)
- [`posterior`](https://mc-stan.org/posterior/) package, used for R-hat and ESS
  in the MCMC convergence checks:

```r
install.packages("posterior")
```

All remaining computations use base R.

## Running the code

Run the script from a working directory containing the three CSV files:

```r
source("BVAR_thesis_full_code.R")
```

The script reruns the simulation study and the empirical application. Some
simulation sections include long MCMC chains and may therefore require
substantial computation time.

The empirical section:

- constructs the quarterly macroeconomic data set;
- estimates the baseline specification and the robustness specification;
- generates the figures used in the dissertation;
- saves the final empirical results in an R data file.

The code is intended to be read and executed as a single script rather than as
an R package.

## Data

The empirical application uses publicly available U.S. macroeconomic series
from the Federal Reserve Economic Data (FRED) database maintained by the
Federal Reserve Bank of St. Louis:

- `PCECTPI`: Personal Consumption Expenditures Price Index;
- `UNRATE`: Civilian Unemployment Rate;
- `FEDFUNDS`: Effective Federal Funds Rate.

The CSV files included in the repository contain the data used by the
replication code.

## Treatment of the COVID-19 period

The empirical application treats the extreme observations associated with the
COVID-19 period separately in the likelihood construction. The quarterly time
index is preserved, so the surrounding observations are not compressed or
re-dated. The exact treatment implemented in the code follows the specification
described in Section 6.1.3 of the dissertation.

## Citation

If you use this code, please cite the dissertation:

> Giuseppetti, M. (2026). *Regime-Specific Bayesian VAR with Structural
> Breaks and Lag-Model Uncertainty*. Master's dissertation, Università
> Cattolica del Sacro Cuore.

## Author

Matteo Giuseppetti  
Master's Degree in Data Analytics for Business  
Università Cattolica del Sacro Cuore
