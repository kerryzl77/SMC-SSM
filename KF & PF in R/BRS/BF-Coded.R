library(tidyverse)
library(dplyr)
library(ggplot2)

# Read the data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data_path <- "../data/combined_Data_imm.csv"
data <- read.csv(data_path)

sdp <- sqrt(1.6e11) # process error set to the optimization result for kalman filter

bayesian_particle_filter_adjusted <- function(data, num_particles = 10000, birth_rate_prior = c(0, 0.1), death_rate_prior = c(0, 0.1), thresh = 0.8) {
  
  # Initialize particles
  birth_rate_particles  <- runif(num_particles, birth_rate_prior[1], birth_rate_prior[2])  
  death_rate_particles <- runif(num_particles, death_rate_prior[1], death_rate_prior[2])
  
  weights <- rep(1/num_particles, num_particles)
  
  estimated_parameters_and_population <- data.frame()
  
  for (year in 2:nrow(data)) {
    prev_population <- data[year - 1, 'Population']
    observed_population <- data[year, 'Population']
    
    sdo <- sqrt(0.03 * prev_population) # Standard Deviation of Observation each Year

    # Simulate process model with changes:
    died <- rbinom(num_particles, prev_population, death_rate_particles)  # Binomial death
    survived <- prev_population - died                                   # Survivors
    births <- rpois(num_particles, lambda = survived * birth_rate_particles) # Poisson birth
    estimated_population <- survived + births + data[year, "Immigration_Count"] 
    
    # Add normal variation for process error
    estimated_population <- estimated_population + rnorm(num_particles, 0, sdp)
    
    # Update weights
    weights <- weights * dlnorm(observed_population, meanlog = log(estimated_population), sdlog = log((sdo)), log = FALSE)

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

ggplot(adjusted_bayesian_estimated_population_df, aes(x = Year)) +
  geom_line(aes(y = Estimated_Death_Rate, color = "Estimated Death Rate")) +
  geom_line(aes(y = Estimated_birth_Rate, color = "Estimated Birth Rate")) +
  xlab('Year') +
  ylab('Rate %') +
  ggtitle('Estimated Birth & Death Rate') +
  scale_color_manual(values = c('Estimated Death Rate' = 'blue', 'Estimated Birth Rate' = 'red')) +
  theme_minimal()
  
plot(adjusted_bayesian_estimated_population_df$Estimated_birth_Rate, main = 'Estimated Birth Rate', ylab = 'Birth Rate %', xlab = 'Year', type = "l")
plot(adjusted_bayesian_estimated_population_df$Estimated_Death_Rate, main = 'Estimated Death Rate', ylab = 'Death Rate %', xlab = 'Year', type = "l")

ggplot(adjusted_bayesian_estimated_population_df, aes(x = Year, y = ESS)) +
  geom_line() + 
  ggtitle("Evolution of Effective Sample Size (ESS)") +
  xlab("Year") +
  ylab("ESS") + 
  theme_minimal()


