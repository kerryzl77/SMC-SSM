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

# stateSpaceModelCode <- nimbleCode({
#   # Priors
#   phi ~ dunif(0, 0.1)   # phi (death rate)
#   beta ~ dunif(0, 0.1)  # beta (birth rate)
# 
#   # Initial state
#   S[1] ~ dpois(lambda = N_initial) # Initial surviving population
#   Y[1] ~ dlnorm(meanlog = log(S[1]), sdlog = log(sdo)) # Observation for initial state
# 
#   # Process model for subsequent years
#   for (t in 2:T) {
#     b[t] ~ dbin(size = S[t-1], prob = beta) # Binomial birth process
#     d[t] ~ dbin(size = S[t-1], prob = phi)  # Binomial death process
#     S[t] ~ dpois(lambda = S[t-1] + b[t] - d[t] + I[t])
# 
#     # Observation process
#     Y[t] ~ dlnorm(meanlog = log(S[t]), sdlog = log(sdo))
#   }
# 
# })

stateSpaceModelCode <- nimbleCode({
  # Priors
  phi ~ dunif(0, 0.1)   # phi (death rate)
  beta ~ dunif(0, 0.1)  # beta (birth rate)

  # Initial state
  S[1] ~ dpois(lambda = N_initial) # Initial population
  Y[1] ~ dlnorm(meanlog = log(S[1]), sdlog = log(sdo)) # Initial state observation

  # Process model
  for (t in 2:T) {
    # Binomial death process
    d[t] ~ dbin(size = S[t-1], prob = phi)

    # Surviving population after deaths but before births
    S_post_death[t] <- S[t-1] - d[t]

    # Poisson birth process based on the updated surviving population
    b[t] ~ dpois(lambda = beta * S_post_death[t])

    # Update total surviving population including new births
    S[t] ~ dpois(lambda = S_post_death[t] + b[t] + I[t])

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
bootstrapFilter <- buildBootstrapFilter(stateSpaceModel, nodes = "S",
                                        control = list(saveAll = TRUE, thresh = 0.7))
compiledFilter <- compileNimble(bootstrapFilter)

# Number of particles deafault 10000
parNum <- 100000
# Run bootstrap filter to estimate state sequences
compiledFilter$run(parNum)
ESS <- compiledFilter$returnESS()

# Extract equally weighted posterior samples of state variables
posteriorSamples <- as.matrix(compiledFilter$mvEWSamples)

# Generate time series estimates of the population dynamic
timeSeriesEstimates <- apply(posteriorSamples, 2, mean)
lowerCI <- apply(posteriorSamples, 2, quantile, probs = 0.025)
upperCI <- apply(posteriorSamples, 2, quantile, probs = 0.975)
# Add a Year column to timeSeriesEstimates that corresponds to the years in Population_Data
years = Combined_Data$Year  # assuming the years in Population_Data align with the indices of timeSeriesEstimates
timeSeriesEstimates_df <- data.frame(
    Year = years,
    Estimated_Population = timeSeriesEstimates,
    Lower_CI = lowerCI,
    Upper_CI = upperCI
)

ggplot(timeSeriesEstimates_df, aes(x = Year, y = Estimated_Population)) +
  geom_line(col = "red") +  
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.3, col = "grey") +  # Confidence interval band 
  geom_point(data = Combined_Data, aes(x = Year, y = Population), col = "blue", pch = 16) +
  xlab("Year") + ylab("Population") +
  ggtitle("SMC - Population Estimates vs Actual with 95% CI") +
  theme_minimal()

uniqueParticleCounts <- apply(posteriorSamples, 2, function(x) length(unique(x)))
plot(ESS)
plot(uniqueParticleCounts)
hist(posteriorSamples[,183])


