library(ggplot2)
library(data.table)
library(readxl)
library(dplyr)
library(nimble)
library(nimbleSMC)
library(coda)

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
head(Combined_Data)

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

################################################################################################################

stateSpaceModelCode <- nimbleCode({
  # Priors 
  phi ~ dunif(0, 1)   # phi (survival rate)
  beta ~ dnorm(0, sd = 1000)   # beta (birth rate)
  
  # Initial state
  S[1] ~ dnorm(phi * N_initial, sd = sqrt(phi * (1 - phi) * N_initial))
  b[1] ~ dnorm(beta * S[1], sd = sqrt(beta * S[1]))
  B[1] <- S[1] + b[1]
  
  # Process model for survival
  for (t in 2:(T)) {
    S[t] ~ dnorm(phi * N[t-1], sd = sqrt(phi * (1 - phi) * N[t-1]))
    
    # Process model for births
    b[t] ~ dnorm(beta * S[t], sd = sqrt(beta * S[t]))
    
    # Net population increase
    B[t] <- S[t] + b[t]
    
    # Likelihood for observed data 
    observedBirths[t] ~ dnorm(b[t], sd = sqrt(b[t]))
    observedDeaths[t] ~ dnorm(N[t-1] - S[t], sd = sqrt(N[t-1] - S[t]))
  }
  
})

data <- list(
  N = Combined_Data$Population, # Assuming this is the right data variable
  observedBirths = Combined_Data$Birth_Count,
  observedDeaths = Combined_Data$Death_Count
)

constants <- list(T = nrow(Combined_Data))

# Initial values for the MCMC
phi_initial <- mean(1 - (Combined_Data$Death_Count / Combined_Data$Population))
beta_initial <- mean(Combined_Data$Birth_Count / Combined_Data$Population)
initial_values <- list(phi = phi_initial, beta = beta_initial)
inits <- list(
  phi = phi_initial, 
  beta = beta_initial,
  N_initial = Combined_Data$Population[1] # First year's population
)

# Build the model
stateSpaceModel <- nimbleModel(stateSpaceModelCode, data = data, constants = constants, inits = inits)

# build bootstrap filter and compile model and filter
bootstrapFilter <- buildBootstrapFilter(stateSpaceModel, nodes = 'S')
compiledList <- compileNimble(stateSpaceModel, bootstrapFilter)

# Number of particles
parNum <- 5000
# Run bootstrap filter and returns estimate of model log-likelihood
compiledList$bootstrapFilter$run(parNum)

################################################################################################################

## extract equally weighted posterior samples of x[10]  and create a histogram
posteriorSamples <- as.matrix(compiledList$bootstrapFilter$mvEWSamples)
hist(posteriorSamples)

# Extract particle matrix for 'S'
boot_X <- as.matrix(compiledList$bootstrapFilter$mvEWSamples, 'S')
plot(boot_X)
