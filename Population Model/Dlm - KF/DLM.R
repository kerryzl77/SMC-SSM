library(dlm)

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
# Step 1: Define the DLM model function
m2pDlm <- function(theta) {
  dlm(m0 = Population_Data$Population[1], C0 = 100, GG = 1, W = exp(theta[1]), FF = 1, V = exp(theta[2]))
}

# Step 2: Fit the model by maximizing the likelihood
# Assuming sdp and sdo are defined. If not, you need to initialize them
sdp <- 1  # Initial guess or calculated value for process noise
sdo <- 1  # Initial guess or calculated value for observation noise
f2pDlm <- dlmMLE(y = Population_Data$Population, parm = c(log(sdp^2), log(sdo^2)), build = m2pDlm)
sqrt(exp(f2pDlm$par))

# Create the DLM model object with the estimated parameters
mMle2pDlm <- m2pDlm(f2pDlm$par)

# Step 3: Filtering to get state estimates
filter2pDlm <- dlmFilter(y = Population_Data$Population, mod = mMle2pDlm)

# Extract filtered estimates
filtered_estimates_dlm <- filter2pDlm$m
filtered_estimates_dlm <- filtered_estimates_dlm[2:152]

# Compute standard errors and Confidence Interval
fvar2pDlm <- unlist(dlmSvd2var(filter2pDlm$U.C,
                               filter2pDlm$D.C))
fvar2pDlm <- fvar2pDlm[2:152]
sd <- sqrt(fvar2pDlm)
ci_upper_dlm <- filtered_estimates_dlm + (1.96 * sd)
ci_lower_dlm <- filtered_estimates_dlm - (1.96 * sd)
# Convert to data frame for plotting
kf_df_dlm <- data.frame(Year = Population_Data$Year, Estimated_Population = filtered_estimates_dlm, 
                        CI_Upper = ci_upper_dlm, CI_Lower = ci_lower_dlm)

# Combine with the true population data
Combined_Plot_Data_DLM <- left_join(Combined_Data, kf_df_dlm, by = "Year")

# Plotting with Confidence Interval using DLM estimates
ggplot(Combined_Plot_Data_DLM) +
  geom_point(aes(x = Year, y = Population), color = "blue", size = 2) +
  geom_line(aes(x = Year, y = Estimated_Population), color = "red") +
  geom_ribbon(aes(x = Year, ymin = CI_Lower, ymax = CI_Upper), alpha = 0.5, fill = "grey") +
  theme_minimal()
