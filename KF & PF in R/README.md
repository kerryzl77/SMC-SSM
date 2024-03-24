# Project README

This repository contains a collection of R scripts for modeling population dynamics through state-space models and Sequential Monte Carlo (SMC) methods. Tailored to specific model formulations detailed the paper, these scripts utilize various probabilistic models to estimate and analyze population states influenced by birth, death, and migration data.

## Directory Structure

The main directory, `KF & PF in R`, is divided into two subdirectories for organizing scripts related to Kalman Filters (`KF`) and Particle Filters (`PF`), each containing further subdivisions as follows:

### Kalman Filter Scripts
Located within the `kalman_filter` folder, the following scripts apply Kalman filter techniques:

- **Kalman Filter - Model 1.R**: Fits a linear Gaussian state-space model to population data using the Kalman filter, starting with a 'diffuse' covariance matrix and optimizing noise parameters.

- **optimizeKalmanParams.R**: Optimizes the Kalman Filter's noise parameters using the BFGS algorithm.

### Particle Filter Scripts
The `particle_filter` subfolder includes scripts for various particle filter models:

- **SMC - PoisLogN.R** (Model 2 - Poisson-Lognormal case): Implements a particle filter with a Poisson process for births and a log-normal observation model.

- **BF-Coded.R** (Model 3): Implements a hard-coded Bootstrap Filter for a non-linear state-space model.

#### Within NimbleSMC Folder
Nested inside `particle_filter`, the `NimbleSMC` folder contains scripts that utilize the NimbleSMC package for advanced SMC models:

- **SMC - Bin.R**: Simulates survival and birth events using binomial processes.

- **SMC bin dynamic rate.R**: Showcases a dynamic birth and death rate model, highlighting issues with compilation due to stochasticity.

- **SMC-Normal-phi/beta.R**: Establishes a vague prior and a normal model.

- **SMC-Normal - sdo/sdp.R**: Employs exact birth and death data to estimate dynamic process and observation errors.

- **SMC-Normal-Exact.R**: Utilizes exact birth and death data with a process equation identical to the Kalman filter model.

## Data Sources and Methods

### Population Data
The dataset includes national and subnational population estimates for the UK, incorporating birth, death, and migration statistics, provided by ONS, NRS, and NISRA.

### Migration Data
Migration data are derived from various sources to represent long-term migration figures, compensating for the UK's lack of a population registration system.
