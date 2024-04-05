import pandas as pd
import numpy as np
from numpy.random import multinomial
from scipy.stats import poisson
import matplotlib.pyplot as plt

data = pd.read_csv('/Users/liuzikai/Desktop/Population Model/combined_Data.csv')
def bayesian_particle_filter_adjusted(data, num_particles=1000, survival_rate_prior=(0.97, 0.01), death_rate_prior=(0.015, 0.005)):
    """
    Adjusted Bayesian Particle Filtering for a state space model with binomial survival and death processes.

    Args:
    data (DataFrame): The dataset containing the population data.
    num_particles (int): The number of particles used in the filter.
    survival_rate_prior (tuple): The mean and std of the prior distribution for the survival rate.
    death_rate_prior (tuple): The mean and std of the prior distribution for the death rate.

    Returns:
    DataFrame: Estimated population and parameters for each year.
    """

    # Initialize particles: Sample from prior distributions for survival and death rates
    survival_rate_particles = np.random.normal(survival_rate_prior[0], survival_rate_prior[1], num_particles)
    death_rate_particles = np.random.normal(death_rate_prior[0], death_rate_prior[1], num_particles)
    
    # Ensure particles are within reasonable bounds
    survival_rate_particles = np.clip(survival_rate_particles, 0, 1)
    death_rate_particles = np.clip(death_rate_particles, 0, 1)

    weights = np.ones(num_particles) / num_particles

    estimated_parameters_and_population = []

    for year in range(1, len(data)):
        prev_population = data.iloc[year - 1]['Population']
        observed_population = data.iloc[year]['Population']

        # Simulate survival and death process
        survived = np.random.binomial(prev_population, survival_rate_particles)
        died = np.random.binomial(survived, death_rate_particles)
        estimated_population = survived - died

        # Add normal variation for observation error
        estimated_population = estimated_population + np.random.normal(0, prev_population * 0.01, num_particles)

        # Update weights based on the likelihood of the observed population
        weights *= poisson.pmf(observed_population, estimated_population)
        weights += 1e-300  # avoid division by zero
        weights /= sum(weights)  # normalize

        # Resample particles
        indices = np.random.choice(np.arange(num_particles), size=num_particles, p=weights)
        survival_rate_particles = survival_rate_particles[indices]
        death_rate_particles = death_rate_particles[indices]

        # Estimate parameters: mean of particles
        estimated_survival_rate = np.mean(survival_rate_particles)
        estimated_death_rate = np.mean(death_rate_particles)

        estimated_parameters_and_population.append([data.iloc[year]['Year'], estimated_survival_rate, estimated_death_rate, observed_population, estimated_population.mean()])

    estimated_parameters_and_population_df = pd.DataFrame(estimated_parameters_and_population, columns=['Year', 'Estimated_Survival_Rate', 'Estimated_Death_Rate', 'Observed_Population', 'Estimated_Population'])

    return estimated_parameters_and_population_df

# Applying the adjusted Bayesian particle filter to estimate the population and parameters
adjusted_bayesian_estimated_population_df = bayesian_particle_filter_adjusted(data)
plt.plot(adjusted_bayesian_estimated_population_df['Year'], adjusted_bayesian_estimated_population_df['Observed_Population'], label='Actual Population', color='blue')
plt.plot(adjusted_bayesian_estimated_population_df['Year'], adjusted_bayesian_estimated_population_df['Estimated_Population'], label='Estimated Population', color='red')
plt.xlabel('Year')
plt.ylabel('Population')
plt.title('Population Comparison')
plt.legend()
plt.show()
