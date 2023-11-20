# Parameters
N_initial <- 1000
phi <- 0.98
beta <- 0.03
mu <- 50
sigma <- 5
epsilon <- 0.01
years <- 1871:1875

# Initialize vectors
N <- numeric(length(years))
S <- numeric(length(years))
b <- numeric(length(years))
B <- numeric(length(years))
i <- numeric(length(years))
I <- numeric(length(years))
E <- numeric(length(years))

N[1] <- N_initial

for(t in 1:(length(years)-1)){
  
  # Survival Dynamics
  S[t+1] <- rnorm(1, N[t]*phi, sqrt(N[t]*phi*(1-phi)))
  
  # Birth Dynamics
  b[t+1] <- rnorm(1, beta*S[t+1], sqrt(beta*S[t+1]))
  B[t+1] <- S[t+1] + b[t+1]
  
  # Immigration Dynamics
  i[t+1] <- rnorm(1, mu, sigma)
  I[t+1] <- B[t+1] + i[t+1]
  
  # Emigration Dynamics
  E[t+1] <- rnorm(1, I[t+1]*epsilon, sqrt(I[t+1]*epsilon*(1-epsilon)))
  
  # T+1 Population
  N[t+1] <- I[t+1] - E[t+1]
}

data.frame(Year = years, Population = N)
