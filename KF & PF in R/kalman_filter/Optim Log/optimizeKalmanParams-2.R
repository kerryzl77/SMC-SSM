library(ggplot2)
library(data.table)
library(dplyr)
library(FKF.SP)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data_path <- "../data/combined_Data.csv"
Combined_Data <- read.csv(data_path)
Combined_Data <- Combined_Data[,-1] # Remove Index

#####################################################################################################
# Optimizing Kalman Filter Parameters  process error (HHt) and observation noise (GGt) using L-BFGS-B 
#####################################################################################################
# Initial state
N_initial <- Combined_Data$Population[!is.na(Combined_Data$Population)][1]

# Initialize vectors
years <- Combined_Data$Year

# Constant Factor 1x1 matrix in measurement equation
Zt <- matrix(1, nrow = 1, ncol = 1)

# Initial state and covariance matrix
a0 <- as.double(N_initial) # Ensure initial state is of type 'double'
P0 <- matrix(1) # Initial variance

# Define the state transition matrix (Tt)
Tt <- matrix(1, nrow = 1, ncol = 1) # 1x1 Transition matrix

# Process error
HHt_Default<- matrix(0.03 * Combined_Data$Population[length(Combined_Data$Population) - 1], nrow = 1, ncol = 1)

# Observation noise
# Approximate with the variance of the yearly differences in population
yearly_differences <- diff(na.omit(Combined_Data$Population))
process_noise_estimate <- sqrt(var(yearly_differences))
GGt_Default <- matrix(process_noise_estimate)

########################
# Observed population
y <- Combined_Data$Population
yt <- rbind(y)
yt[is.na(yt)] <- 0 # Replace NA with 0 for Kalman filtering

dt <- ct <- matrix(0) # intercepts of the transition and measurement equations

#################################################################################################################

# Process error defaults
deviation <- (0.03 * Combined_Data$Population[length(Combined_Data$Population) - 1])
GGt_Default <- log(deviation)  # Log-transformed default
HHt_Default <- log(deviation)  # Log-transformed default

# Define the objective function for optimization
objective_fn <- function(params, yt, a0, P0, dt, ct, Zt, Tt) {
  HHt <- matrix(exp(params[1])) # Measurement noise variance, exponentiated
  GGt <- matrix(exp(params[2])) # Process noise variance, exponentiated
  
  # Run the Kalman filter with the current parameters
  kf <- fkf.SP(HHt = HHt, GGt = GGt, Tt = Tt, Zt = Zt, a0 = a0, P0 = P0, yt = yt, dt = dt, ct = ct, verbose = TRUE)
  
  # Negative log-likelihood
  return(-kf$logLik)
}

# Run the optimization to estimate HHt and GGt
optim_results <- optim(
  par = c(HHt = HHt_Default, GGt = GGt_Default), # Starting values (log-transformed)
  fn = objective_fn, # Objective function
  yt = rbind(yt), # Observed data
  a0 = a0, # Initial state
  P0 = P0, # Initial variance
  dt = dt, # Transition intercept
  ct = ct, # Measurement intercept
  Zt = Zt, # Measurement matrix
  Tt = Tt, # Transition matrix
  method = "L-BFGS-B", # Optimization method
  lower = c(-Inf, -Inf), # Lower bounds for log-transformed HHt and GGt
  upper = c(Inf, Inf) # Upper bounds for log-transformed HHt and GGt
)

# Print the optimization results
print(optim_results)

# Extract optimized values in real scale
optimized_HHt <- exp(optim_results$par[1])
optimized_GGt <- exp(optim_results$par[2])
optimized_GGt;optimized_HHt
