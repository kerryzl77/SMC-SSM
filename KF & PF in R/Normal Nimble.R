library(ggplot2)
library(data.table)
library(readxl)
library(dplyr)
library(nimble)
library(nimbleSMC)


# Tt <- matrix(phi + beta, nrow = 1, ncol = 1) # 1x1 Transition matrix
# start the particle filter work with this constant 
################################################################################################################
# Set working directory
setwd('/Users/liuzikai/Desktop/Population Model')

# Read birth data
Birth_Data <- read.csv('Births Count 1838-2022.csv', skip = 7) %>%
  select(Year = 1, Birth_Count = 2) %>%
  mutate(Type = "Birth")

# Read death data
Death_Data <- read_excel('Death 1838-2021.xlsx', range = "A5:B189", sheet = 1) %>%
  rename(Year = 1, Death_Count = 2) %>%
  mutate(Type = "Death")

Combined_Data <- inner_join(Birth_Data, Death_Data, by = "Year")
Combined_Data <- Combined_Data %>%
  select(Year, Birth_Count, Death_Count)
# head(Combined_Data)

# # Plotting Birth and Death counts on the same graph
# ggplot(data = Combined_Data, aes(x = Year)) +
#   geom_line(aes(y = Birth_Count, color = "Births")) +
#   geom_line(aes(y = Death_Count, color = "Deaths")) +
#   labs(x = "Year", y = "Count", color = "Event Type") +
#   theme_minimal()

Population_Data <- read_excel('Population-Statistic/Population 1871-2021.xlsx', 
                              sheet = 'Data', 
                              range = "B6:C156",
                              col_names = c('Year', 'Population')) %>%
  mutate(Year = as.numeric(Year)) 

Combined_Data <- inner_join(Combined_Data, Population_Data, by = "Year")

# Observation Error set to 1.5% of the Population Variance*
GGt_Default <- 0.015 * Combined_Data$Population[length(Combined_Data$Population) - 1]
################################################################################################################
#####################################
# S[t] represents the true (but unobserved) state of the population
normalModelCode <- nimbleCode({
  # Priors for the parameters
  phi ~ dunif(0, 1)   # survival rate
  beta ~ dnorm(0, sd = 1000)   # birth rate
  sigma_observation ~ dnorm(GGt_Default, 100)  # observation error - around the 1.5% estimate
  
  # Initial state
  S[1] ~ dnorm(Population_initial, sd = sqrt(Population_initial))
  
  # Process Equation
  for (t in 2:T) {
    # Process equation survival & birth process
    S[t] ~ dnorm(phi * S[t-1] + beta * S[t-1], sd = sqrt(S[t-1]))  # Variance proportional to S[t-1]
  }
  
  # Measurement Equation 
  for (t in 1:T) {
    y[t] ~ dnorm(S[t], sd = sigma_observation)  # Observation model
  }
})

# Data for the model
data <- list(
  y = Combined_Data$Population
)
constants <- list(T = nrow(Combined_Data))

# Initial values 
phi_initial <- mean(1 - (Combined_Data$Death_Count / Combined_Data$Population))
beta_initial <- mean(Combined_Data$Birth_Count / Combined_Data$Population)
initial_values <- list(
  phi = phi_initial, 
  beta = beta_initial,
  sigma_observation = 10000,  
  Population_initial = Combined_Data$Population[1] 
)


# Build the model
normalModel <- nimbleModel(
  normalModelCode, 
  data = data, 
  constants = constants, 
  inits = initial_values
)

# Compile the model
compiledModel <- compileNimble(normalModel)

# build bootstrap filter for state estimation
bootstrapFilter <- buildBootstrapFilter(normalModel, nodes = "S",
                   control = list(saveAll = TRUE, thresh = 0.9))
compiledFilter <- compileNimble(bootstrapFilter)

# Number of particles
parNum <- 5000
# Run bootstrap filter to estimate state sequences
compiledFilter$run(parNum)

# Extract equally weighted posterior samples of state variables
posteriorSamples <- as.matrix(compiledFilter$mvEWSamples)

hist(posteriorSamples)
unique(posteriorSamples)

# Generate time series estimates of the population dynamic
timeSeriesEstimates <- apply(posteriorSamples, 2, mean)
# Add a Year column to timeSeriesEstimates that corresponds to the years in Population_Data
years = Population_Data$Year  # assuming the years in Population_Data align with the indices of timeSeriesEstimates
timeSeriesEstimates_df = data.frame(Year = years, Estimated_Population = timeSeriesEstimates)

# Plot the estimated population data as a line
plot(timeSeriesEstimates_df$Year, timeSeriesEstimates_df$Estimated_Population, type = "l", col = "red", 
     xlab = "Year", ylab = "Population", main = "Population Estimates vs Actual Data")
# Add the actual population data as points
points(Population_Data$Year, Population_Data$Population, col = "blue", pch = 16)


################################################################################################################
# YOY actual population change
actual_growth <- c(NA, diff(Combined_Data$Population))
# Expected growth from birth and death data
expected_growth <- Combined_Data$Birth_Count - Combined_Data$Death_Count
# Immigrant gap
immigrant_gap <- actual_growth - expected_growth
Population_Growth_Data <- data.frame(
  Year = Combined_Data$Year,
  Expected_Growth = expected_growth,
  Actual_Growth = actual_growth,
  Immigrant_Gap = immigrant_gap
)
Population_Growth_Data <- Population_Growth_Data[-1, ]
# View the new data frame
head(Population_Growth_Data)
plot(Population_Growth_Data$Immigrant_Gap)

# 

