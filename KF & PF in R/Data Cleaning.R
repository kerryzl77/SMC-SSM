library(dplyr)
library(readxl)
################################################################################################################
# Set working directory
folder_path <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "/data/Data - England & Wales")
setwd(folder_path)

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

Population_Data <- read.csv('England & Wales by Sex 1838-2022.csv', skip = 1) %>%
  select(Year = 1, Population = 2, Male_Pop = 3, Female_pop = 4) %>%
  mutate(
    Year = as.numeric(Year),
    Population = as.numeric(gsub(",", "", Population)),  
    Male_Pop = as.numeric(gsub(",", "", Male_Pop)),      
    Female_pop = as.numeric(gsub(",", "", Female_pop))  
  )

Migration_Data <- read_excel('Migration 1964-2023.xlsx', range = "A1:B60", sheet = 3, col_names = FALSE) %>%
  setNames(c("Year", "Immigration_Count")) %>%
  mutate(
    Year = as.numeric(Year),
    Immigration_Count = as.numeric(Immigration_Count)
  )

Combined_Data <- inner_join(Combined_Data, Population_Data, by = "Year")
Combined_Data <- left_join(Combined_Data, Migration_Data, by = "Year") %>%
  mutate(Immigration_Count = ifelse(is.na(Immigration_Count), 0, Immigration_Count))
Combined_Data$Birth_rate <- round(Combined_Data$Birth_Count / Combined_Data$Population,3)
Combined_Data$Death_rate <- round(Combined_Data$Death_Count / Combined_Data$Population,3)
head(Combined_Data)

################################################################################################################
# Reproduce Figure 1 
ggplot(data = Combined_Data, aes(x = Year)) +
  geom_line(aes(y = Birth_Count, color = "Births")) +  
  geom_line(aes(y = Death_Count, color = "Deaths")) +  
  labs(x = "Year", y = "Count", color = "Event Type") +
  scale_color_manual(values = c(Births = 'red', Deaths = 'blue')) +  
  ggtitle('Birth and Death Count 1838-2021') +
  theme_minimal()

# Reproduce Figure 2
ggplot(data = Combined_Data, aes(x = Year)) +
  geom_line(aes(y = Birth_rate, color = "Births")) +  
  geom_line(aes(y = Death_rate, color = "Deaths")) +  
  labs(x = "Year", y = "Rate %", color = "Event Type") +
  scale_color_manual(values = c(Births = 'red', Deaths = 'blue')) +  
  ggtitle('Birth and Death Rate 1838-2021') +
  theme_minimal()

# Reproduce Figure 3
actual_growth <- c(NA, diff(Combined_Data$Population)) # YOY actual population change
expected_growth <- Combined_Data$Birth_Count - Combined_Data$Death_Count
immigrant_gap <- actual_growth - expected_growth # Gap (actual - expected)
Population_Growth_Data <- data.frame(
  Year = Combined_Data$Year,
  Expected_Growth = expected_growth,
  Actual_Growth = actual_growth,
  Immigrant_Gap = immigrant_gap
)
Population_Growth_Data <- Population_Growth_Data[-1, ]
plot(Population_Growth_Data$Year, Population_Growth_Data$Immigrant_Gap, type = "l", main = "Yearly Immigration Gap", xlab = "Year", ylab = "Immigration Gap")

