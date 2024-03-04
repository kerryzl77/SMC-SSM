library(tidyverse)
library(dplyr)
library(ggplot2)

# Read the data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data_path <- "../data/combined_Data_imm.csv"
data <- read.csv(data_path)

bayesian_particle_filter_adjusted <- function(data, num_particles = 10000, birth_rate_prior = c(0, 0.1), death_rate_prior = c(0, 0.1), thresh = 0.9) {
  # Initialize particles
  birth_rate_particles  <- rnorm(num_particles, birth_rate_prior[1], birth_rate_prior[2])  
  death_rate_particles <- rnorm(num_particles, death_rate_prior[1], death_rate_prior[2])
  
  # Ensure particles between (0,1)
  birth_rate_particles <- pmin(pmax(birth_rate_particles, 0), 1)
  death_rate_particles <- pmin(pmax(death_rate_particles, 0), 1)
  
  weights <- rep(1/num_particles, num_particles)
  
  estimated_parameters_and_population <- data.frame()
  
  for (year in 2:nrow(data)) {
    prev_population <- data[year - 1, 'Population']
    observed_population <- data[year, 'Population']

    # Simulate process model with changes:
    died <- rbinom(num_particles, prev_population, death_rate_particles)  # Binomial death
    survived <- prev_population - died                                   # Survivors
    births <- rpois(num_particles, lambda = survived * birth_rate_particles) # Poisson birth
    estimated_population <- survived + births + data[year, "Immigration_Count"] 
    
    # Add normal variation for observation error
    estimated_population <- estimated_population + rnorm(num_particles, 0, prev_population * 0.01)
    
    # Update weights
    weights <- weights * dpois(observed_population, estimated_population)
    weights <- weights + 1e-300
    weights <- weights / sum(weights)
    
    # Resample particles
    ess <- sum(weights^2) * num_particles # Calculate Effective Sample Size
    if (ess < num_particles * thresh) { # Resample if ESS below threshold
      indices <- sample(1:num_particles, size = num_particles, replace = TRUE, prob = weights)
      birth_rate_particles <- birth_rate_particles[indices]
      death_rate_particles <- death_rate_particles[indices]
      weights <- rep(1/num_particles, num_particles) # Reset weights after resampling
    }
    
    # Estimate parameters
    estimated_birth_rate <- mean(birth_rate_particles)
    estimated_death_rate <- mean(death_rate_particles)
    
    # Calculate Confidence Intervals
    lower_CI <- quantile(estimated_population, probs = 0.025) 
    upper_CI <- quantile(estimated_population, probs = 0.975) 
    
    estimated_parameters_and_population <- rbind(estimated_parameters_and_population, data.frame(
      Year = data[year, 'Year'],
      Estimated_birth_Rate = estimated_birth_rate,
      Estimated_Death_Rate = estimated_death_rate,
      Observed_Population = observed_population,
      Estimated_Population = mean(estimated_population),
      Lower_CI = lower_CI,
      Upper_CI = upper_CI,
      ESS = ess 
    ))
    
  }
  
  return(estimated_parameters_and_population)
}

# Applying the adjusted Bayesian particle filter
adjusted_bayesian_estimated_population_df <- bayesian_particle_filter_adjusted(data)

# Plotting the results
ggplot(adjusted_bayesian_estimated_population_df, aes(x = Year)) +
  geom_line(aes(y = Observed_Population, color = "Actual Population")) +
  geom_line(aes(y = Estimated_Population, color = "Estimated Population")) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.3, col = "grey") + 
  xlab('Year') +
  ylab('Population') +
  ggtitle('Particle Filtering Hard-Coded') +
  scale_color_manual(values = c('Actual Population' = 'blue', 'Estimated Population' = 'red')) +
  theme_minimal()

ggplot(adjusted_bayesian_estimated_population_df, aes(x = Year, y = ESS)) +
  geom_line() + 
  ggtitle("Evolution of Effective Sample Size (ESS)") +
  xlab("Year") +
  ylab("ESS") + 
  theme_minimal()


