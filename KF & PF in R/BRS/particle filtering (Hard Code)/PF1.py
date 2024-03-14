import pandas as pd
import numpy as np
from numpy.random import multinomial
from scipy.stats import poisson
import matplotlib.pyplot as plt

data = pd.read_csv('/Users/liuzikai/Desktop/Population Model/combined_Data.csv')

def particle_filter_population_estimation(data, num_particles=1000):
    """
    Particle Filtering (Sequential Monte Carlo) for estimating population each year, along with death and survival rates.

    Args:
    data (DataFrame): The dataset containing the population data.
    num_particles (int): The number of particles used in the filter.

    Returns:
    DataFrame: Estimated population and parameters for each year.
    """

    # Initialize particles
    # Assuming death rate and survival rate follow a uniform distribution initially
    particles = np.random.uniform(0, 1, (num_particles, 2))  # columns: [death_rate, survival_rate]

    # Weights of particles
    weights = np.ones(num_particles) / num_particles

    estimated_parameters_and_population = []

    for year in range(1, len(data)):
        prev_population = data.iloc[year - 1]['Population']
        births = data.iloc[year]['Birth_Count']
        deaths = data.iloc[year]['Death_Count']

        # Expected births and deaths for each particle
        expected_births = prev_population * particles[:, 1]
        expected_deaths = prev_population * particles[:, 0]

        # Update weights based on the likelihood of the actual births and deaths
        weights *= poisson.pmf(births, expected_births) * poisson.pmf(deaths, expected_deaths)
        weights += 1e-300  # avoid division by zero
        weights /= sum(weights)  # normalize

        # Resample particles
        indices = np.random.choice(np.arange(num_particles), size=num_particles, p=weights)
        particles = particles[indices, :]

        # Estimate parameters: mean of particles
        estimated_death_rate = np.mean(particles[:, 0])
        estimated_survival_rate = np.mean(particles[:, 1])

        # Estimate population for next year
        estimated_population_next_year = (prev_population + births) * estimated_survival_rate - deaths * estimated_death_rate
        estimated_parameters_and_population.append([data.iloc[year]['Year'], estimated_death_rate, estimated_survival_rate, estimated_population_next_year])

    estimated_parameters_and_population_df = pd.DataFrame(estimated_parameters_and_population, columns=['Year', 'Estimated_Death_Rate', 'Estimated_Survival_Rate', 'Estimated_Population'])

    return estimated_parameters_and_population_df

# Applying the modified particle filter to estimate the population
estimated_population_df = particle_filter_population_estimation(data)
# estimated_population_df.to_csv('/Users/liuzikai/Desktop/Population Model/estimated_population_df.csv')

def particle_filter_corrected(data, num_particles=1000):
    """
    Particle Filtering (Sequential Monte Carlo) for estimating death and survival rates, with corrected resampling.

    Args:
    data (DataFrame): The dataset containing the population data.
    num_particles (int): The number of particles used in the filter.

    Returns:
    DataFrame: Estimated parameters for each year.
    """

    # Initialize particles
    # Assuming death rate and survival rate follow a uniform distribution initially
    particles = np.random.uniform(0, 1, (num_particles, 2))  # columns: [death_rate, survival_rate]

    # Weights of particles
    weights = np.ones(num_particles) / num_particles

    estimated_parameters = []

    for year in range(1, len(data)):
        # Predict next state for each particle
        prev_population = data.iloc[year - 1]['Population']
        births = data.iloc[year]['Birth_Count']
        deaths = data.iloc[year]['Death_Count']

        # Expected births and deaths for each particle
        expected_births = prev_population * particles[:, 1]
        expected_deaths = prev_population * particles[:, 0]

        # Update weights based on the likelihood of the actual births and deaths
        weights *= poisson.pmf(births, expected_births) * poisson.pmf(deaths, expected_deaths)
        weights += 1e-300  # avoid division by zero
        weights /= sum(weights)  # normalize

        # Resample particles
        indices = np.random.choice(np.arange(num_particles), size=num_particles, p=weights)
        particles = particles[indices, :]

        # Estimate parameters: mean of particles
        estimated_death_rate = np.mean(particles[:, 0])
        estimated_survival_rate = np.mean(particles[:, 1])
        estimated_parameters.append([data.iloc[year]['Year'], estimated_death_rate, estimated_survival_rate])

    estimated_parameters_df = pd.DataFrame(estimated_parameters, columns=['Year', 'Estimated_Death_Rate', 'Estimated_Survival_Rate'])

    return estimated_parameters_df

corrected_estimated_params_df = particle_filter_corrected(data)
# corrected_estimated_params_df.to_csv('/Users/liuzikai/Desktop/Population Model/corrected_estimated_params_df.csv')

def bayesian_particle_filter(data, num_particles=1000, birth_rate_prior=(0.01, 0.02), death_rate_prior=(0.01, 0.02)):
    """
    Bayesian Particle Filtering for a state space model with binomial birth and death processes.

    Args:
    data (DataFrame): The dataset containing the population data.
    num_particles (int): The number of particles used in the filter.
    birth_rate_prior (tuple): The mean and std of the prior distribution for the birth rate.
    death_rate_prior (tuple): The mean and std of the prior distribution for the death rate.

    Returns:
    DataFrame: Estimated population and parameters for each year.
    """

    # Initialize particles: Sample from prior distributions for birth and death rates
    birth_rate_particles = np.random.normal(birth_rate_prior[0], birth_rate_prior[1], num_particles)
    death_rate_particles = np.random.normal(death_rate_prior[0], death_rate_prior[1], num_particles)
    
    # Ensure particles are within reasonable bounds
    birth_rate_particles = np.clip(birth_rate_particles, 0, 1)
    death_rate_particles = np.clip(death_rate_particles, 0, 1)

    weights = np.ones(num_particles) / num_particles

    estimated_parameters_and_population = []

    for year in range(1, len(data)):
        prev_population = data.iloc[year - 1]['Population']
        observed_population = data.iloc[year]['Population']

        # Simulate next year's population for each particle
        births = np.random.binomial(prev_population, birth_rate_particles)
        deaths = np.random.binomial(prev_population, death_rate_particles)
        estimated_population = prev_population + births - deaths

        # Update weights based on the likelihood of the observed population
        weights *= poisson.pmf(observed_population, estimated_population)
        weights += 1e-300  # avoid division by zero
        weights /= sum(weights)  # normalize

        # Resample particles
        indices = np.random.choice(np.arange(num_particles), size=num_particles, p=weights)
        birth_rate_particles = birth_rate_particles[indices]
        death_rate_particles = death_rate_particles[indices]

        # Estimate parameters: mean of particles
        estimated_birth_rate = np.mean(birth_rate_particles)
        estimated_death_rate = np.mean(death_rate_particles)

        estimated_parameters_and_population.append([data.iloc[year]['Year'], estimated_birth_rate, estimated_death_rate, observed_population, estimated_population.mean()])

    estimated_parameters_and_population_df = pd.DataFrame(estimated_parameters_and_population, columns=['Year', 'Estimated_Birth_Rate', 'Estimated_Death_Rate', 'Observed_Population', 'Estimated_Population'])

    return estimated_parameters_and_population_df

# Applying the Bayesian particle filter to estimate the population and parameters
bayesian_estimated_population_df = bayesian_particle_filter(data)
# Merging the estimated data with the actual data for comparison
comparison_data = pd.merge(data[['Year', 'Population']], 
                           estimated_population_df, 
                           on='Year', 
                           how='inner')
plt.plot(comparison_data['Year'], comparison_data['Population'], label='Actual Population', color='blue')
plt.plot(comparison_data['Year'], comparison_data['Estimated_Population'], label='Estimated Population', color='red')
plt.xlabel('Year')
plt.ylabel('Population')
plt.title('Population Comparison')
plt.legend()
plt.show()