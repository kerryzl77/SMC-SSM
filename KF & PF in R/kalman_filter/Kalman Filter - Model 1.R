library(ggplot2)
library(data.table)
library(FKF)
library(readxl)
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data_path <- "../data/combined_Data_imm.csv"
Combined_Data <- read.csv(data_path)

## <--------------------------------------------------------------------------->
## Kalman Filter: Univariate SSMs for Population Growth
## <--------------------------------------------------------------------------->
## Transition equation:
## alpha[t+1] = Growth rate * alpha[t] + immigration[t], eta[t] ~ N(0, HHt)          
## Measurement equation:
## y[t] = alpha[t] + eps[t], eps[t] ~  N(0, GGt)

# Define the parameters (phi, beta)
phi <- mean(1 - (Combined_Data$Death_Count / Combined_Data$Population))
beta <- mean(Combined_Data$Birth_Count / Combined_Data$Population)

years <- Combined_Data$Year

# Constant factor of 1 in measurement equation
Zt <- matrix(1, nrow = 1, ncol = 1)

# Initial state and variance
a0 <- as.double(Combined_Data$Population[!is.na(Combined_Data$Population)][1]) # Initial state
P0 <- matrix(1000) # Initial variance

# Constant factor of state transition matrix
Tt <- matrix(phi + beta, nrow = 1, ncol = 1) 

# Observation noise (Set as 3% of the mean population count)
GGt_Default <- matrix(mean(0.03 * Combined_Data$Population), nrow = 1, ncol = 1)

########################
# Observed population
y <- Combined_Data$Population
yt <- rbind(y)
yt[is.na(yt)] <- 0 # Replace NA with 0 for Kalman filtering

dt <- matrix(as.double(Combined_Data$Immigration_Count),nrow = 1)
ct <- matrix(0) # intercepts of the transition and measurement equations

# Process noise
# Approximate with the variance of the yearly differences in population
yearly_differences <- diff(na.omit(Combined_Data$Population))
process_noise_estimate <- (var(yearly_differences))
process_noise_squared <- (sqrt(var(yearly_differences)))
HHt_Default <- matrix(process_noise_estimate)
HHt_squared <- matrix(process_noise_squared)

################################################################################################################
# Fit Kalman Filter and Plot
run_kalman_and_plot <- function(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt, years, data, plot_title) {
  kf <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt, yt = yt)
  
  log_likelihood <- kf$logLik
  print(paste("Log-likelihood:", log_likelihood)) 

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
  Combined_Plot_Data <- left_join(data, kf_df, by = "Year")
  
  # Plotting with Confidence Interval
  ggplot() +
    geom_point(data = Combined_Plot_Data, aes(x = Year, y = Population), color = "blue", size = 2) +
    geom_line(data = Combined_Plot_Data, aes(x = Year, y = Estimated_Population), color = "red") +
    geom_ribbon(data = Combined_Plot_Data, aes(x = Year, ymin = ci_lower, ymax = ci_upper), alpha = 0.5, fill = "grey") +
    theme_minimal()+
    labs(title = plot_title) # Adding title to the plot
  
}

source("optimizeKalmanParams.R")
optimized_HHt <- optimize_kalman_params(yt, Tt, Zt, a0, P0, dt, ct, HHt_Default, GGt_Default)

run_kalman_and_plot(a0, P0, dt, ct, Tt, Zt, HHt_Default, GGt_Default, yt, years, Combined_Data, 'Default process noise')
run_kalman_and_plot(a0, P0, dt, ct, Tt, Zt, optimized_HHt, GGt_Default, yt, years, Combined_Data, 'Optimized process noise')
run_kalman_and_plot(a0, P0, dt, ct, Tt, Zt, HHt_squared, GGt_Default, yt, years, Combined_Data, 'Squared process noise')







