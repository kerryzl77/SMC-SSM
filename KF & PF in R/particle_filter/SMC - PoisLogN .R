library(ggplot2)
library(data.table)
library(dplyr)
library(nimbleSMC)
library(coda)
library(nimble)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data_path <- "../data/combined_Data_imm.csv"
Combined_Data <- read.csv(data_path)
deviation <- 0.015 * Combined_Data$Population[length(Combined_Data$Population) - 1]

################################################################################################################
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

# Initial values for the SMC
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
parNum <- 10000
# Run bootstrap filter to estimate state sequences
compiledFilter$run(parNum)

###################################
## create MCMC specification for the state space model
stateSpaceMCMCconf <- configureMCMC(stateSpaceModel, nodes = NULL)

## add a block pMCMC sampler for a, b, sigPN, and sigOE 
stateSpaceMCMCconf$addSampler(target = c('phi', 'beta'),
                              type = 'RW_PF_block', control = list(latents = 'S'))

## build and compile pMCMC sampler
stateSpaceMCMC <- buildMCMC(stateSpaceMCMCconf)
compiledList <- compileNimble(stateSpaceModel, stateSpaceMCMC, resetFunctions = TRUE)

compiledList$stateSpaceMCMC$run(5000)

par(mfrow = c(2,1))
posteriorSamps <- as.mcmc(as.matrix(compiledList$stateSpaceMCMC$mvSamples))
traceplot(posteriorSamps[,'phi'], ylab = 'phi')
traceplot(posteriorSamps[,'beta'], ylab = 'beta')
hist(posteriorSamps[,'phi'])
hist(posteriorSamps[,'beta'])
par(mfrow = c(1,1))

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
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.3, col = "grey", size = 2) +  # Confidence interval band
  geom_line(col = "red", size = 0.5) +
  geom_point(data = Combined_Data, aes(x = Year, y = Population), col = "blue", pch = 16) +
  xlab("Year") + ylab("Population") +
  ggtitle("SMC - Population Estimates vs Actual with 95% CI") +
  theme_minimal()



uniqueParticleCounts <- apply(posteriorSamples, 2, function(x) length(unique(x)))
plot(ESS)
plot(uniqueParticleCounts, main = 'Unique Particles over Iterations', ylab = 'Count', xlab = 'Index')
hist(posteriorSamples[,183], main = 'Posterior Distribution of 2021 Population Estimate', xlab = 'Count')
abline(v = tail(Combined_Data$Population, 1), col = "red", lwd = 2)
# 59660524 not in the region  

# Calculate residuals as the difference between observed and estimated populations
residuals <- Combined_Data$Population - timeSeriesEstimates_df$Estimated_Population
timeSeriesEstimates_df$Residuals <- residuals

plot(timeSeriesEstimates_df$Year, timeSeriesEstimates_df$Residuals,  ylab = "One-Step Ahead Residuals", xlab = "Year")
hist(residuals)


