library(nimble)

# Define the model
nimbleCode <- nimbleCode({
  for(t in 1:T) {
    # Process equation
    z[t] ~ dnorm(beta * z[t-1], sd = sigma_p) # State transition
    
    # Observation equation
    y[t] ~ dnorm(alpha * z[t], sd = sigma_o)   # Observation model
    
    # Prior for the initial state
    z[1] ~ dnorm(z_init, sd = sigma_p)
  }
})

# Initial values for the state
z_init <- 0 # This would be your initial guess or prior for the first state

# Constants and data (replace these with actual values or data)
constants <- list(T = 100, alpha = 1, beta = 0.9, sigma_p = 0.1, sigma_o = 0.1)
data <- list(y = rnorm(constants$T, 0, 1)) # replace with actual data

# Model configuration
rModel <- nimbleModel(code = nimbleCode, data = data, constants = constants, inits = list(z = rep(z_init, constants$T)))

# Compile the model
cModel <- compileNimble(rModel)

# Configure the SMC algorithm
smcConf <- configureMCMC(rModel, monitors = 'z')

# Use a particle filter (sequential Monte Carlo)
smcAlgorithm <- buildBootstrapFilter(rModel, smcConf)

# Compile the SMC algorithm
cSMC <- compileNimble(smcAlgorithm, project = rModel)

# Run the SMC
cSMC$run(niter = 1000)

# Get the estimates
estimates <- cSMC$calculate(c('z'))
