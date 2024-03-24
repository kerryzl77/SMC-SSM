



Combined_Data <- read.csv('/Users/liuzikai/Desktop/Population Model/combined_Data.csv')
Combined_Data <- Combined_Data[,-1]
################################################################################################################

N_initial <- Combined_Data$Population[!is.na(Combined_Data$Population)][1]

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
Combined_Plot_Data <- left_join(Combined_Data, SSM_df, by = "Year")
ggplot() +
  geom_point(data = Combined_Plot_Data, aes(x = Year, y = Population), color = "blue", size = 2) +
  geom_line(data = Combined_Plot_Data, aes(x = Year, y = SSM_Population), color = "red") +
  theme_minimal()