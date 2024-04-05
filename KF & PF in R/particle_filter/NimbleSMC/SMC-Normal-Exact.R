library(ggplot2)
library(data.table)
library(readxl)
library(dplyr)
library(nimble)
library(nimbleSMC)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
data_path <- "../data/combined_Data_imm.csv"
Combined_Data <- read.csv(data_path)
sdo <- 0.015 * Combined_Data$Population[length(Combined_Data$Population) - 1]

#####################################
normalModelCode <- nimbleCode({
  # Initial state
  S[1] ~ dnorm(Population_initial, sd = 100000)
  # Process Equation using exact phi and beta
  for (t in 2:T) {
    S[t] ~ dnorm((S[t-1] - death[t-1] + birth[t-1]+ immigration[t-1]), sd = sqrt(S[t-1]))
  }
  # Measurement Equation 
  for (t in 1:T) {
    y[t] ~ dnorm(S[t], sd = sdo)
  }
})

# Data for the model
data <- list(
  y = Combined_Data$Population,
  death = Combined_Data$Death_Count,
  birth = Combined_Data$Birth_Count,
  immigration = Combined_Data$Immigration_Count
)
constants <- list(T = nrow(Combined_Data))

# Initial values 
initial_values <- list(
  Population_initial = Combined_Data$Population[1],
  sdo = sdo
)


# Build the model
normalModel <- nimbleModel(
  normalModelCode, 
  data = data, 
  constants = constants, 
  inits = initial_values
)
normalModel$initializeInfo()

# Compile the model
compiledModel <- compileNimble(normalModel)


# build bootstrap filter for state estimation
bootstrapFilter <- buildBootstrapFilter(normalModel, nodes = "S",
                                        control = list(saveAll = TRUE, thresh = 0.8))
compiledFilter <- compileNimble(bootstrapFilter)

# Number of particles
parNum <- 5000
# Run bootstrap filter to estimate state sequences
compiledFilter$run(parNum)

# Extract equally weighted posterior samples of state variables
posteriorSamples <- as.matrix(compiledFilter$mvEWSamples)

# Particle Degeneracy Check
hist(posteriorSamples[,105])
unique_values_lengths <- numeric(183)
for (i in 1:183) {
  unique_values_lengths[i] <- length(unique(posteriorSamples[,i]))
}
years = Combined_Data$Year
plot(x = years, y = unique_values_lengths)

# Generate time series estimates of the population dynamic
timeSeriesEstimates <- apply(posteriorSamples, 2, mean)
lowerCI <- apply(posteriorSamples, 2, quantile, probs = 0.025)
upperCI <- apply(posteriorSamples, 2, quantile, probs = 0.975)
# Add a Year column to timeSeriesEstimates that corresponds to the years in Population_Data
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
  ggtitle("Population Estimates w Exact Birth/Death") +
  theme_minimal()



