# Bayesian Monotonic Regression using JAGS

This repository contains the code and models for performing Bayesian monotonic regression using JAGS to analyze experimental data on the effects of vinclozolin on androgen receptor (AR) activity. The data were collected from lab-grown Chinese Hamster ovary cells exposed to different levels of the fungicide vinclozolin. The project aims to model the relationship between vinclozolin concentration and AR activity, assuming a monotonic decreasing effect.

## Project Description

The dataset `vin.txt` contains data from five experimental runs, where:

- `conc` represents the concentration of vinclozolin in nano-moles.
- `effect` is the measured androgen receptor (AR) activity.

## Repository Structure

The repository contains the following files:

- **M1s2596860.jags**: JAGS model file (Model 1) for the baseline monotonic regression model with a single latent variable sequence for all experimental runs.
- **M2s2596860.jags**: JAGS model file (Model 2) for a modified monotonic regression model where each experimental run has its own decreasing sequence of latent variables.
- **README.md**: This file, containing an overview of the repository and its contents.
- **s2596860.R**: R code file containing the analysis, model fitting, and diagnostics for the monotonic regression models. It includes the setup for the JAGS models, data loading, and model diagnostics.

  
### Analysis

The R code file (`s2596860.R`) performs the following tasks:
1. Loads the data from `vin.txt`.
2. Defines and runs the JAGS models (Model 1 and Model 2).
3. Conducts convergence diagnostics (including trace plots and effective sample size checks).
4. Generates the required plots, including the effect vs. dose index plot with model predictions and credible intervals.
5. Compares the models using the Deviance Information Criterion (DIC) to assess which model fits the data better.

## How to Run the Analysis

1. Ensure you have R and JAGS installed on your system.
2. Download the data file `vin.txt` and place it in the working directory.
3. Run the R script (`s2596860.R`) to execute the analysis. The script will automatically load the data, run the models in JAGS, and produce the necessary plots.

## Required R Packages

- `rjags` for interfacing with JAGS.
- `coda` for MCMC diagnostics.
- `ggplot2` for generating the plots.

Make sure these packages are installed before running the analysis.

## Plots

The analysis will produce the following plots:

1. **Effect vs. Dose Plot**: A plot showing the observed AR activity (`effect`) against the vinclozolin concentration (`conc`), color-coded by experimental run. The plot will also show the model's predicted mean effect with 95% credible intervals.
   
2. **Convergence Diagnostics**: Plots showing the trace of the parameter with the highest and lowest effective sample size.

## Model Comparison

The models will be compared using the Deviance Information Criterion (DIC). A comment on which model is preferable based on DIC will be included in the final report.


