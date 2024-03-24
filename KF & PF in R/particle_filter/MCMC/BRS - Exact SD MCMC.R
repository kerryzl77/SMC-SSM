library(ggplot2)
library(data.table)
library(readxl)
library(dplyr)
library(nimble)
library(nimbleSMC)
library(coda)

Combined_Data <- read.csv('/Users/liuzikai/Desktop/Population Model/combined_Data.csv')
Combined_Data <- Combined_Data[,-1]

################################################################################################################

stateSpaceModelCode <- nimbleCode({
  # Priors 
  sdp ~ dunif(1e-04, 1e06)
  sdo ~ dunif(1e-04, 1e06)
  
  # Initial state
  S[1] ~ dnorm(N_initial, sd = sdp) # Initial surviving population
  Y[1] ~ dnorm((S[1]), sd = sdo) # Observation for initial state
  
  # Process model for subsequent years
  for (t in 2:T) {
    phi <- phi_list[t]
    beta <- beta_list[t]
    S[t] ~ dnorm(S[t-1] * (phi) * (1 + beta), sd = sdp) # Net surviving population
    
    # Normal observation process
    Y[t] ~ dnorm((S[t]), sd = sdo)
  }
})


data <- list(
  Y = Combined_Data$Population,  # Observed population data
  phi_list = 1 - Combined_Data$Death_rate,
  beta_list = Combined_Data$Birth_rate
)

constants <- list(T = nrow(Combined_Data),N_initial = Combined_Data$Population[1])

sdo_start <- 0.015 * Combined_Data$Population[length(Combined_Data$Population) - 1]
sdp_start <- sdo_start

inits <- list(sdp = sdp_start, sdo = sdo_start)

# Build the model
stateSpaceModel <- nimbleModel(stateSpaceModelCode, data = data, constants = constants, inits = inits)

# Compile the model
compiledModel <- compileNimble(stateSpaceModel)

# Define MCMC configuration with specific monitors
mcmcConf <- configureMCMC(stateSpaceModel, monitors = c("sdp", "sdo"))

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

# > summary(mcmcRun)
# sdo              sdp        
# Min.   :891838   Min.   :894311  
# 1st Qu.:895747   1st Qu.:895837  
# Median :895849   Median :895963  
# Mean   :895673   Mean   :896530  
# 3rd Qu.:895943   3rd Qu.:897186  
# Max.   :898396   Max.   :900274 


# Convert MCMC samples to data frame
sdp_samples_df <- as.data.frame(mcmcRun[, "sdp"])
sdo_samples_df <- as.data.frame(mcmcRun[, "sdo"])

# Plot for sdp
hist(sdp_samples_df$`mcmcRun[, "sdp"]`, main = "Posterior Distribution of sdp", xlab = "sdp", ylab = "Density")
hist(sdo_samples_df$`mcmcRun[, "sdo"]`, main = "Posterior Distribution of sdo", xlab = "sdo", ylab = "Density")

