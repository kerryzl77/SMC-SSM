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

Population_Data <- read_excel('Population 1801 - 2021.xlsx', range = "A9:C32")
Population_Data <- Population_Data %>%
  select(Year = 1, Population = 3)

Combined_Data <- left_join(Combined_Data, Population_Data, by = "Year")
Combined_Data$Population <- as.numeric(as.character(Combined_Data$Population))


################################################################################################################
# Define the parameters (phi, beta)
valid_data <- Combined_Data[!is.na(Combined_Data$Population), ]

valid_survival_data <- valid_data[!is.na(valid_data$Death_Count), ]
phi <- mean(1 - (valid_survival_data$Death_Count / valid_survival_data$Population))

valid_birth_data <- valid_data[!is.na(valid_data$Birth_Count), ]
beta <- mean(valid_birth_data$Birth_Count / valid_birth_data$Population)

# Initial state
N_initial <- Combined_Data$Population[!is.na(Combined_Data$Population)][1]

# Initialize vectors
years <- Combined_Data$Year

# Constant Factor 1x1 matrix in measurement equation
Zt <- matrix(1, nrow = 1, ncol = 1)

# Initial state and covariance matrix
a0 <- N_initial # Initial state
P0 <- matrix(1) # Initial variance 

################################################################################################################

# Define the state transition matrix (Tt)
Tt <- matrix(phi + beta, nrow = 1, ncol = 1) 

# Measurement noise (Assume zero - Comparison)
HHt_Default <- (0) 
# Process noise (Assume zero - Comparison)
GGt_Default <- phi * (1 - phi)

########################
y <- Combined_Data$Population
y <- rbind(y)
y[is.na(y)] <- 0 # Replace NA with 0 for Kalman filtering purposes

dt <- ct <- matrix(0) # intercepts of the transition and measurement equations

# MLE Step: Estimate Variances using Maximum Likelihood Estimation
optim_results <- optim(c(measurement_noise = HHt_Default, process_noise = GGt_Default), 
                       fn = function(params) {
                           HHt <- matrix(params[1])
                           GGt <- matrix(params[2])
                           -fkf(HHt = HHt, GGt = GGt, Tt = Tt, Zt = Zt, a0 = a0, P0 = P0, yt = y, dt = dt, ct = ct)$logLik
                       })

# Extract the estimated variances
HHt_est <- matrix(optim_results$par[1])
GGt_est <- matrix(optim_results$par[2])

# Running the Kalman filter with estimated noise variances
kf <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, Zt = Zt, HHt = HHt_est, GGt = GGt_est, yt = y)


# Estimated states (filtered values)
# kf$att # these are the filtered estimates of the state (population)

# Extract the filtered state estimates
filtered_estimates <- as.numeric(kf$att)

# Convert the Kalman Filter estimates to a data frame for plotting
kf_df <- data.frame(Year = years, Estimated_Population = filtered_estimates)
# Combine with the true population data
Combined_Plot_Data <- left_join(Combined_Data, kf_df, by = "Year")

################################################################################################################
# State Space Model Hard Coded 
N <- S <- b <- B <- rep(NA, length(years))
N[1] <- N_initial
# Process model loop
for(t in 1:(length(years) - 1)) {
  # Survival Dynamics
  S[t + 1] <- rnorm(1, mean = N[t] * phi, sd = sqrt(N[t] * phi * (1 - phi)))
  # Birth Dynamics
  b[t + 1] <- rnorm(1, mean = beta * S[t + 1], sd = sqrt(beta * S[t + 1]))
  B[t + 1] <- S[t + 1] + b[t + 1]
  # Population for next year
  N[t + 1] <- B[t + 1] 
}
SSM_df <- data.frame(Year = years, SSM_Population = N)
Combined_Plot_Data <- left_join(Combined_Plot_Data, SSM_df, by = "Year")
ggplot() +
  geom_point(data = Combined_Plot_Data, aes(x = Year, y = Population), color = "blue", size = 2) +
  geom_line(data = Combined_Plot_Data, aes(x = Year, y = Estimated_Population), color = "red") +
  geom_line(data = Combined_Plot_Data, aes(x = Year, y = SSM_Population), color = "green") +
  theme_minimal()
