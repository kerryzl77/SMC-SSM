#####################################################################################################
# Optimizing Kalman Filter Parameters process error (HHt) 
# by minimizing the negative log-likelihood, keeping observation noise (GGt) default using L-BFGS-B 
#####################################################################################################

optimize_kalman_params <- function(yt, Tt, Zt, a0, P0, dt, ct, HHt_Default, GGt_Default) {
  # Define the negative log-likelihood function
  kalman_negllik <- function(params, yt, Tt, Zt, a0, P0, dt, ct, GGt_Default) {
    HHt <- matrix(exp(params))  # Only optimizing for HHt
    GGt <- GGt_Default          # Using the default GGt value
    
    logLik <- -fkf(HHt = HHt, GGt = GGt, Tt = Tt, Zt = Zt, a0 = a0, P0 = P0, yt = yt, dt = dt, ct = ct)$logLik
    return(logLik)
  }
  
  # Perform the optimization
  optim_results <- optim(c(HHt = log(HHt_Default)),  # Initial param only for HHt
                         fn = kalman_negllik,
                         yt = yt,
                         Tt = Tt,
                         Zt = Zt,
                         a0 = a0,
                         P0 = P0,
                         dt = dt,
                         ct = ct,
                         GGt_Default = GGt_Default,  # Pass the default GGt
                         control = list(trace = 1),
                         method = 'BFGS',
                         hessian = FALSE)
  
  # Extract and return optimized HHt
  optimized_HHt <- matrix(exp(optim_results$par))
  return(optimized_HHt)
}

