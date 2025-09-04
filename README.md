# Survival Simulations for Penalized Cox Models

This repository contains simulation studies, real-data applications, 
and benchmarking experiments for **Exclusive Lasso**, **Elastic Net**, 
**Group Lasso**, and **IPF-Lasso**.

## Structure

- `R/` – reusable R functions (data generation, preprocessing, fitting, evaluation)
- `configs/` – configuration files for different simulation scenarios
- `scripts/` – drivers (`run_all.R`, `time_compare.R`)
- `results/` – output CSVs and RDS files (ignored in git)

## Usage

Run all simulations:

```r
source("scripts/run_all.R")

Run the fair wall-clock budget comparison:

```r
source("scripts/time_compare.R")
