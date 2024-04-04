library(tidyverse)
library(dplyr)
library(ggplot2)

# Read the data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data_path <- "../data/combined_Data_imm.csv"
data <- read.csv(data_path)

bayesian_particle_filter_adjusted <- function(data, num_particles = 10000, birth_rate_prior = c(0, 0.1), death_rate_prior = c(0, 0.1), thresh = 0.8) {
  
  set.seed(4599) 
  overall_log_likelihood <- 0 
  
  birth_rate_particles  <- runif(num_particles, birth_rate_prior[1], birth_rate_prior[2])  
  death_rate_particles <- runif(num_particles, death_rate_prior[1], death_rate_prior[2])
  
  weights <- rep(1/num_particles, num_particles)
  
  estimated_parameters_and_population <- data.frame()
  
  estimated_birth_rate <-  data$Birth_rate[1]
  estimated_death_rate <- data$Death_rate[1]
  
  final_birth_rate_particles <- NULL
  final_death_rate_particles <- NULL
  final_estimated_population <- NULL 
  
  for (year in 2:nrow(data)) {
    prev_population <- data[year - 1, 'Population']
    observed_population <- data[year, 'Population']
    
    sdo <- sqrt(0.03 * prev_population) # Standard Deviation of Observation each Year
    
    # Sample around the estimates from the previous year
    birth_rate_particles <- pmax(rnorm(num_particles, mean = estimated_birth_rate, sd = 0.01), 0)
    death_rate_particles <- pmax(rnorm(num_particles, mean = estimated_death_rate, sd = 0.01), 0)
    
    # Simulate process model with changes:x
    died <- rbinom(num_particles, prev_population, death_rate_particles)  # Binomial death
    survived <- prev_population - died                                   # Survivors
    births <- rpois(num_particles, lambda = survived * birth_rate_particles) # Poisson birth
    estimated_population <- survived + births + data[year, "Immigration_Count"] 
    
    # Add normal variation for process error
    estimated_population <- estimated_population
    
    # Update weights
    unnormalized_weights <- weights * dlnorm(observed_population, meanlog = log(estimated_population), sdlog = log((sdo)), log = FALSE)
    
    overall_log_likelihood <- overall_log_likelihood + log(sum((unnormalized_weights))/num_particles)
    
    weights <- unnormalized_weights / sum(unnormalized_weights)
    
    # Resample particles
    ess <- 1 / (sum(weights^2))
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
    
    # Print unique particle counts
    cat("Year:", data[year, 'Year'], "- Unique Particles:", length(unique(estimated_population)))
    
    estimated_parameters_and_population <- rbind(estimated_parameters_and_population, data.frame(
      Year = data[year, 'Year'],
      Estimated_birth_Rate = estimated_birth_rate,
      Estimated_Death_Rate = estimated_death_rate,
      Observed_Population = observed_population,
      Estimated_Population = mean(estimated_population),
      Lower_CI = lower_CI,
      Upper_CI = upper_CI,
      ESS = ess,
      Likelihood = overall_log_likelihood
    ))
    
  }
  # Update final particles at each iteration
  final_birth_rate_particles <- birth_rate_particles
  final_death_rate_particles <- death_rate_particles
  final_estimated_population <- estimated_population  
  
  return(list(
    Estimated_Parameters = estimated_parameters_and_population,
    Final_Birth_Rate_Particles = final_birth_rate_particles,
    Final_Death_Rate_Particles = final_death_rate_particles,
    Final_Estimated_Population = final_estimated_population
  ))}

# Run the Boostrap Model
results <- bayesian_particle_filter_adjusted(data)
final_birth_rates <- results$Final_Birth_Rate_Particles
final_death_rates <- results$Final_Death_Rate_Particles
final_estimated_populations <- results$Final_Estimated_Population

par(mfrow=c(1,3))
hist(final_estimated_populations, main = 'Estimated Populations')
abline(v = tail(data$Population, 1), col = "red", lwd = 2)

hist(final_birth_rates, main = 'Estimated Birth rate')
abline(v = tail(data$Birth_rate, 1), col = "red", lwd = 2)

hist(final_death_rates, main = 'Estimated Death rate')
abline(v = tail(data$Death_rate, 1), col = "red", lwd = 2)
par(mfrow=c(1,1))

adjusted_bayesian_estimated_population_df <- results$Estimated_Parameters
cat("Log Likelihood:", tail(adjusted_bayesian_estimated_population_df$Likelihood,1), "\n")


# Plot results
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

ggplot(adjusted_bayesian_estimated_population_df, aes(x = Year, y = ESS)) +
  geom_line() + 
  ggtitle("Evolution of Effective Sample Size (ESS)") +
  xlab("Year") +
  ylab("ESS") + 
  theme_minimal()

residuals <- adjusted_bayesian_estimated_population_df$Observed_Population - adjusted_bayesian_estimated_population_df$Estimated_Population
adjusted_bayesian_estimated_population_df$Residuals <- residuals

plot(adjusted_bayesian_estimated_population_df$Year, adjusted_bayesian_estimated_population_df$Residuals,  ylab = "One-Step Ahead Residuals", xlab = "Year")
hist(residuals)


