library(ggplot2)
library(data.table)
library(dplyr)
library(nimbleSMC)
library(coda)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
data_path <- "../data/combined_Data_imm.csv"
Combined_Data <- read.csv(data_path)
deviation <- 0.015 * Combined_Data$Population[length(Combined_Data$Population) - 1]

################################################################################################################

stateSpaceModelCode <- nimbleCode({
  # Priors
  phi ~ dunif(0, 0.1)   # phi (death rate)
  beta ~ dunif(0, 0.1)  # beta (birth rate)
  
  # Initial state
  S[1] ~ dpois(lambda = N_initial) # Initial surviving (total) population
  F[1] ~ dpois(lambda = F_initial) # Initial female population
  Y[1] ~ dlnorm(meanlog = log(S[1]), sdlog = log(sdo)) # Observation for initial state
  
  # Process model for subsequent years
  for (t in 2:T) {
    b[t] ~ dbin(size = F[t-1], prob = beta) # Binomial birth process based on female population
    d[t] ~ dbin(size = S[t-1], prob = phi)  # Binomial death process on total population
    S[t] ~ dpois(lambda = S[t-1] + b[t] - d[t] +I[t]) # Total surviving population
    F[t] ~ dpois(lambda = F[t-1] + b[t] * (F[t-1]/S[t-1]) - d[t] * (F[t-1]/S[t-1])) # Proportional female surviving population 
    
    # Observation process
    Y[t] ~ dlnorm(meanlog = log(S[t]), sdlog = log(sdo)) 
  }
})

data <- list(
  Y = Combined_Data$Population,  # Observed total population data
  I = Combined_Data$Immigration_Count
)

constants <- list(T = nrow(Combined_Data), N_initial = Combined_Data$Population[1], 
                  F_initial = Combined_Data$Female_pop[1], sdo = deviation)

# Initial values for the MCMC
phi_initial <- mean((Combined_Data$Death_Count / Combined_Data$Population))
beta_initial <- mean(Combined_Data$Birth_Count / Combined_Data$Female_pop) # Adjusted to use female population for birth rate
initial_values <- list(phi = phi_initial, beta = beta_initial)

inits <- list(
  phi = phi_initial, 
  beta = beta_initial,
  N_initial = Combined_Data$Population[1], # First year's total population
  F_initial = Combined_Data$Female_pop[1]  # First year's female population
)

# Build the model
stateSpaceModel <- nimbleModel(stateSpaceModelCode, data = data, constants = constants, inits = inits)

# Compile the model
compiledModel <- compileNimble(stateSpaceModel)

# build bootstrap filter for state estimation
bootstrapFilter <- buildBootstrapFilter(stateSpaceModel, nodes = "S",
                                        control = list(saveAll = TRUE, thresh = 0.9))
compiledFilter <- compileNimble(bootstrapFilter)

# Number of particles
parNum <- 200
# Run bootstrap filter to estimate state sequences
compiledFilter$run(parNum)

# Extract equally weighted posterior samples of state variables
posteriorSamples <- as.matrix(compiledFilter$mvEWSamples)

# Generate time series estimates of the population dynamic
timeSeriesEstimates <- apply(posteriorSamples, 2, mean)
# Add a Year column to timeSeriesEstimates that corresponds to the years in Population_Data
years = Combined_Data$Year  # assuming the years in Population_Data align with the indices of timeSeriesEstimates
timeSeriesEstimates_df = data.frame(Year = years, Estimated_Population = timeSeriesEstimates)

# Estimation as line and Obs as points
plot(timeSeriesEstimates_df$Year, timeSeriesEstimates_df$Estimated_Population, type = "l", col = "red",
     xlab = "Year", ylab = "Population", main = "SMC - Population Estimates vs Actual")
points(Combined_Data$Year, Combined_Data$Population, col = "blue", pch = 16)



