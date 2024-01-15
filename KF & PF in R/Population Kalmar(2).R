library(ggplot2)
library(data.table)
library(FKF)
library(readxl)
library(dplyr)

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

# Plotting Birth and Death counts on the same graph
ggplot(data = Combined_Data, aes(x = Year)) +
  geom_line(aes(y = Birth_Count, color = "Births")) +
  geom_line(aes(y = Death_Count, color = "Deaths")) +
  labs(x = "Year", y = "Count", color = "Event Type") +
  theme_minimal()

Population_Data <- read_excel('Population-Statistic/Population 1871-2021.xlsx', 
                              sheet = 'Data', 
                              range = "B6:C156",
                              col_names = c('Year', 'Population')) %>%
  mutate(Year = as.numeric(Year)) 

Combined_Data <- inner_join(Combined_Data, Population_Data, by = "Year")


################################################################################################################
# Univariate process
# Define the parameters (phi, beta)
phi <- mean(1 - (Combined_Data$Death_Count / Combined_Data$Population))
beta <- mean(Combined_Data$Birth_Count / Combined_Data$Population)
plot(Combined_Data$Year,Combined_Data$Population)
# Initial state
N_initial <- Combined_Data$Population[!is.na(Combined_Data$Population)][1]

# Initialize vectors
years <- Combined_Data$Year

# Constant Factor 1x1 matrix in measurement equation
Zt <- matrix(1, nrow = 1, ncol = 1)

# Initial state and covariance matrix
a0 <- N_initial # Initial state
P0 <- matrix(1) # Initial variance

# Define the state transition matrix (Tt)
Tt <- matrix(phi + beta, nrow = 1, ncol = 1) # 1x1 Transition matrix

# Observation noise (Set as 1.5% relative standard error of the previous year's population)
GGt_Default <- matrix(0.03 * Combined_Data$Population[length(Combined_Data$Population) - 1], nrow = 1, ncol = 1)
# Calculate GGt as 1.5% relative standard error of the previous year's population
# GGt_Default <- sapply(Combined_Data$Population[-length(Combined_Data$Population)], 
#                       function(pop) matrix(0.015 * pop, nrow = 1, ncol = 1))

########################
# Observed population
y <- Combined_Data$Population
y <- rbind(y)
y[is.na(y)] <- 0 # Replace NA with 0 for Kalman filtering

dt <- ct <- matrix(0) # intercepts of the transition and measurement equations

# Process noise
# Approximate with the variance of the yearly differences in population
yearly_differences <- diff(na.omit(Combined_Data$Population))
process_noise_estimate <- sqrt(var(yearly_differences))
# GGt_Default <- process_noise_estimate
HHt_Default <- matrix(process_noise_estimate)
measurement_noise_estimate <- 0
# # MLE Step: Estimate Measurement Noise Variance using Maximum Likelihood Estimation
# optim_results <- optim(c(measurement_noise = HHt_Default), 
#                        fn = function(params) {
#                          HHt <- matrix(params[1])
#                          GGt <- GGt_Default # Use the predefined GGt
#                          -fkf(HHt = HHt, GGt = GGt, Tt = Tt, Zt = Zt, a0 = a0, P0 = P0, yt = y, dt = dt, ct = ct)$logLik
#                        })
# # Extract the estimated measurement noise variance
# HHt_est <- matrix(optim_results$par[1])

# Running the Kalman filter
kf <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, Zt = Zt, HHt = HHt_Default, GGt = GGt_Default, yt = y)

# Estimated states (filtered values)
# kf$att # these are the filtered estimates of the state (population)

# Extract the filtered state estimates
filtered_estimates <- as.numeric(kf$att)

# Convert the Kalman Filter estimates to a data frame for plotting
kf_df <- data.frame(Year = years, Estimated_Population = filtered_estimates)
# Combine with the true population data
Combined_Plot_Data <- left_join(Combined_Data, kf_df, by = "Year")

################################################################################################################
# Running the Kalman filter
kf <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, Zt = Zt, HHt = HHt_Default, GGt = GGt_Default, yt = y)

# Extract the filtered state estimates
filtered_estimates <- as.numeric(kf$att)

# Calculate standard errors (square root of diagonal elements of Ptt)
standard_errors <- list()

# Loop through each time point in the array
for (i in 1:dim(kf$Ptt)[3]) {
  # Extract the standard error (square root of the variance)
  se <- (kf$Ptt[1, 1, i])
  # Append to the list
  standard_errors[i] <- as.numeric(se)
}
standard_errors_vector <- unlist(standard_errors)

# Compute 95% Confidence Interval
ci_upper <- filtered_estimates + (1.96 * standard_errors_vector)
ci_lower <- filtered_estimates - (1.96 * standard_errors_vector)

# Convert the Kalman Filter estimates and Confidence Interval to a data frame for plotting
kf_df <- data.frame(Year = years, Estimated_Population = filtered_estimates, 
                    CI_Upper = ci_upper, CI_Lower = ci_lower)

# Combine with the true population data
Combined_Plot_Data <- left_join(Combined_Data, kf_df, by = "Year")

# Plotting with Confidence Interval
ggplot() +
  geom_point(data = Combined_Plot_Data, aes(x = Year, y = Population), color = "blue", size = 2) +
  geom_line(data = Combined_Plot_Data, aes(x = Year, y = Estimated_Population), color = "red") +
  geom_ribbon(data = Combined_Plot_Data, aes(x = Year, ymin = ci_lower, ymax = ci_upper), alpha = 0.5, fill = "grey") +
  theme_minimal()

# 
# 
# # State Space Model Hard Coded 
# N <- S <- b <- B <- rep(NA, length(years))
# N[1] <- N_initial
# # Process model loop
# for(t in 1:(length(years) - 1)) {
#   # Survival Dynamics
#   S[t + 1] <- rnorm(1, mean = N[t] * phi, sd = sqrt(N[t] * phi * (1 - phi)))
#   # Birth Dynamics
#   b[t + 1] <- rnorm(1, mean = beta * S[t + 1], sd = sqrt(beta * S[t + 1]))
#   B[t + 1] <- S[t + 1] + b[t + 1]
#   # Population for next year
#   N[t + 1] <- B[t + 1] 
# }
# SSM_df <- data.frame(Year = years, SSM_Population = N)
# Combined_Plot_Data <- left_join(Combined_Plot_Data, SSM_df, by = "Year")
# ggplot() +
#   geom_point(data = Combined_Plot_Data, aes(x = Year, y = Population), color = "blue", size = 2) +
#   geom_line(data = Combined_Plot_Data, aes(x = Year, y = Estimated_Population), color = "red") +
#   geom_line(data = Combined_Plot_Data, aes(x = Year, y = SSM_Population), color = "green") +
#   theme_minimal()
