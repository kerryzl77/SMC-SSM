library(ggplot2)
library(data.table)
library(dplyr)
library(nimble)
library(nimbleSMC)
library(coda)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
data_path <- "../data/combined_Data_imm.csv"
Combined_Data <- read.csv(data_path)

deviation <- 0.015 * Combined_Data$Population[length(Combined_Data$Population) - 1]

################################################################################################################
# Normal - Normal // Estimate phi & beta & constant sdp/sdo // Bootstrap Filter //1839-2021 England & Wales

stateSpaceModelCode <- nimbleCode({
  # Priors 
  phi ~ dunif(0.9, 1)   # phi (survival rate)
  beta ~ dunif(0, 0.1)  # beta (birth rate)
  
  # Initial state
  S[1] ~ dnorm(N_initial, sd = sdp) # Initial surviving population
  Y[1] ~ dnorm((S[1]), sd = sdo) # Observation for initial state
  
  # Process model for subsequent years
  for (t in 2:T) {
    S[t] ~ dnorm(S[t-1] * (phi) * (1 + beta), sd = sdp) # Net surviving population
    
    # Log-normal observation process
    Y[t] ~ dnorm((S[t]), sd = sdo)
  }
})


data <- list(Y = Combined_Data$Population) # Observed population data

constants <- list(T = nrow(Combined_Data), N_initial = Combined_Data$Population[1], 
                  sdo = deviation, sdp = deviation)

# Initial values for the MCMC
phi_initial <- mean(1 - (Combined_Data$Death_Count / Combined_Data$Population))
beta_initial <- mean(Combined_Data$Birth_Count / Combined_Data$Population)

# Use var_S_empirical in your model
inits <- list(phi = phi_initial, beta = beta_initial)

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

# hist(posteriorSamples)
# unique(posteriorSamples)

# Generate time series estimates of the population dynamic
timeSeriesEstimates <- apply(posteriorSamples, 2, mean)
# Add a Year column to timeSeriesEstimates that corresponds to the years in Population_Data
years = Combined_Data$Year  # assuming the years in Population_Data align with the indices of timeSeriesEstimates
timeSeriesEstimates_df = data.frame(Year = years, Estimated_Population = timeSeriesEstimates)

# Estimation as line and Obs as points
plot(timeSeriesEstimates_df$Year, timeSeriesEstimates_df$Estimated_Population, type = "l", col = "red",
     xlab = "Year", ylab = "Population", main = "Var - Population Estimates vs Actual")
points(Combined_Data$Year, Combined_Data$Population, col = "blue", pch = 16)

############################## 10-Year Period Zoom In #################################################
plot_population_interval <- function(start_year, end_year) {
  # Filter the data for the specified interval
  interval_data <- timeSeriesEstimates_df[timeSeriesEstimates_df$Year >= start_year & timeSeriesEstimates_df$Year <= end_year, ]
  actual_data <- Combined_Data[Combined_Data$Year >= start_year & Combined_Data$Year <= end_year, ]
  
  # Plot the estimated population data for the interval
  plot(interval_data$Year, interval_data$Estimated_Population, type = "l", col = "red", 
       xlab = "Year", ylab = "Population", main = paste("Population Estimates vs Actual Data (", start_year, "-", end_year, ")", sep=""),
       xlim = c(start_year, end_year))
  
  # Add the actual population data as points for the same interval
  points(actual_data$Year, actual_data$Population, col = "blue", pch = 16)
}

# Determine the range of years in your dataset
min_year <- min(Combined_Data$Year)
max_year <- max(Combined_Data$Year)

# Loop through each 10-year interval and generate a plot
for (start_year in seq(min_year, max_year, by = 10)) {
  end_year <- start_year + 9
  if (end_year > max_year) { end_year <- max_year } # Adjust the end year if it exceeds the max year in your data
  plot_population_interval(start_year, end_year)
}

