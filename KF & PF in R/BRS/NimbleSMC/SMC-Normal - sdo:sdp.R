library(ggplot2)
library(data.table)
library(dplyr)
library(nimble)
library(nimbleSMC)
library(coda)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
data_path <- "../data/combined_Data.csv"
Combined_Data <- read.csv(data_path)
Combined_Data <- Combined_Data[,-1]

################################################################################################################
# Normal - Normal // Estimate sdp & sdo & known phi/beta // Bootstrap Filter //1839-2021 England & Wales

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

################################################################################################################

# Compile the model
compiledModel <- compileNimble(stateSpaceModel)

# build bootstrap filter for state estimation
bootstrapFilter <- buildBootstrapFilter(stateSpaceModel, nodes = "S",
                                        control = list(saveAll = TRUE, thresh = 0.9))
compiledFilter <- compileNimble(bootstrapFilter)

# Number of particles
parNum <- 5000
# Run bootstrap filter to estimate state sequences
compiledFilter$run(parNum)

# Extract equally weighted posterior samples of state variables
posteriorSamples <- as.matrix(compiledFilter$mvEWSamples)

# hist(posteriorSamples) # Final year population posterior 
# unique(posteriorSamples)

# Generate time series estimates of the population dynamic
timeSeriesEstimates <- apply(posteriorSamples, 2, mean)
# Add a Year column to timeSeriesEstimates that corresponds to the years in Population_Data
years = Combined_Data$Year  # assuming the years in Population_Data align with the indices of timeSeriesEstimates
timeSeriesEstimates_df = data.frame(Year = years[2:184], Estimated_Population = timeSeriesEstimates)

# Plot the estimated population data as a line
plot(timeSeriesEstimates_df$Year, timeSeriesEstimates_df$Estimated_Population, type = "l", col = "red", 
     xlab = "Year", ylab = "Population", main = "Exact - Population Estimates vs Actual")
# Add the actual population data as points
points(Combined_Data$Year, Combined_Data$Population, col = "blue", pch = 16)

