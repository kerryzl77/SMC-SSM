# Examples
## For illustration only.
exampleCode <- nimbleCode({
  x0 ~ dnorm(0, var = 1)
  x[1] ~ dnorm(.8 * x0, var = 1)
  y[1] ~ dnorm(x[1], var = .5)
  for(t in 2:10){
    x[t] ~ dnorm(.8 * x[t-1], var = 1)
    y[t] ~ dnorm(x[t], var = .5)
  }
})
model <- nimbleModel(code = exampleCode, data = list(y = rnorm(10)),
                     inits = list(x0 = 0, x = rnorm(10)))
my_BootF <- buildBootstrapFilter(model, 'x',
                                 control = list(saveAll = TRUE, thresh = 1))
## Now compile and run, e.g.,
Cmodel <- compileNimble(model)
Cmy_BootF <- compileNimble(my_BootF, project = model)
logLik <- Cmy_BootF$run(m = 1000)
ESS <- Cmy_BootF$returnESS()
boot_X <- as.matrix(Cmy_BootF$mvEWSamples, 'x')