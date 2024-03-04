library(ggplot2)
library(data.table)
library(dplyr)
library(nimbleSMC)
library(coda)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data_path <- "../data/combined_Data_imm.csv"
Combined_Data <- read.csv(data_path)
deviation <- 0.015 * Combined_Data$Population[length(Combined_Data$Population) - 1]

################################################################################################################


stateSpaceModelCode <- nimbleCode({
  # Priors for time-varying rates
  for(t in 1:T) {
    phi[t] ~ dunif(0, 0.1)   # Time-varying phi (death rate) for each year
    beta[t] ~ dunif(0, 0.1)  # Time-varying beta (birth rate) for each year
  }
  
  # Initial state
  S[1] ~ dpois(lambda = N_initial) # Initial surviving population
  Y[1] ~ dlnorm(meanlog = log(S[1]), sdlog = log(sdo)) # Observation for initial state
  
  # Process model for subsequent years
  for (t in 2:T) {
    b[t] ~ dbin(size = S[t-1], prob = beta[t-1]) # Binomial birth process with time-varying rate
    d[t] ~ dbin(size = S[t-1], prob = phi[t-1])  # Binomial death process with time-varying rate
    S[t] ~ dpois(lambda = S[t-1] + b[t] - d[t] + I[t]) 
    
    # Observation process
    Y[t] ~ dlnorm(meanlog = log(S[t]), sdlog = log(sdo)) 
  }
})




data <- list(
  Y = Combined_Data$Population,  # Observed population data
  I = Combined_Data$Immigration_Count
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

# build bootstrap filter for state estimation
# build bootstrap filter for state estimation
# Specify all nodes you want to monitor, including 'S[t]' for all t, 'beta', and 'phi'
nodes_to_monitor <- c("beta", "phi")

bootstrapFilter <- buildBootstrapFilter(stateSpaceModel, nodes = "beta",
                                        control = list(saveAll = TRUE, thresh = 0.5))
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

length(unique(posteriorSamples[,183]))
hist(posteriorSamples[,183])


