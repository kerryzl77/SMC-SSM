library(ggplot2)
library(data.table)
library(readxl)
library(dplyr)
library(nimble)
library(nimbleSMC)
library(coda)

Combined_Data <- read.csv('/Users/liuzikai/Desktop/Population Model/combined_Data.csv')
Combined_Data <- Combined_Data[,-1]
deviation <- 0.015 * Combined_Data$Population[length(Combined_Data$Population) - 1]

################################################################################################################

stateSpaceModelCode <- nimbleCode({
  # Priors
  phi ~ dunif(0, 0.1)   # phi (death rate)
  beta ~ dunif(0, 0.1)  # beta (birth rate)
  
  # Initial state
  S[1] <- N_initial # Initial surviving population
  Y[1] ~ dnorm((S[1]), sd = sdo) # Observation for initial state
  
  
  # Process model for subsequent years
  for (t in 2:T) {
    b[t] ~ dbin(beta, S[t-1]) # Binomial birth process
    d[t] ~ dbin(phi, S[t-1])  # Binomial death process
    S[t] <- S[t-1] + b[t] - d[t]
    
    # Observation process
    Y[t] ~ dnorm((S[t]), sd = sdo)
    
  }
  
})

data <- list(
  Y = Combined_Data$Population  # Observed population data
)

constants <- list(T = nrow(Combined_Data), N_initial = Combined_Data$Population[1], 
                  sdo = deviation)

# Initial values for the MCMC
phi_initial <- mean((Combined_Data$Death_Count / Combined_Data$Population))
beta_initial <- mean(Combined_Data$Birth_Count / Combined_Data$Population)
initial_values <- list(phi = phi_initial, beta = beta_initial)
inits <- list(
  phi = phi_initial, 
  beta = beta_initial,
  N_initial = Combined_Data$Population[1] # First year's population
)

# Build the model
stateSpaceModel <- nimbleModel(stateSpaceModelCode, data = data, constants = constants, inits = inits)

# Compile the model
compiledModel <- compileNimble(stateSpaceModel)

# Define MCMC configuration with specific monitors
mcmcConf <- configureMCMC(stateSpaceModel, monitors = c("phi", "beta", paste("S[", 1:constants$T, "]", sep = "")))

# Build the MCMC object
mcmcObj <- buildMCMC(mcmcConf)

# Compile the MCMC object
compiledMCMC <- compileNimble(mcmcObj, project = stateSpaceModel)

# Run MCMC
niter <- 1000  # Total number of iterations
nburnin <- 200  # Number of burn-in iterations
mcmcRun <- runMCMC(compiledMCMC, niter = niter, nburnin = nburnin)

# Summary of the MCMC output
summary(mcmcRun)

timeSeriesEstimates <- apply(mcmcRun[, 1:(ncol(mcmcRun) - 2)], 2, mean)
years = Combined_Data$Year  
timeSeriesEstimates_df = data.frame(Year = years, Estimated_Population = timeSeriesEstimates)

# Estimation as line and Obs as points
plot(timeSeriesEstimates_df$Year, timeSeriesEstimates_df$Estimated_Population, type = "l", col = "red",
     xlab = "Year", ylab = "Population", main = "MCMC - Population Estimates vs Actual")
points(Combined_Data$Year, Combined_Data$Population, col = "blue", pch = 16)



