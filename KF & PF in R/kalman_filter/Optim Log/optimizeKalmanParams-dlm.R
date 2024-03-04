
library(ggplot2)
library(data.table)
library(dplyr)
library(dlm)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data_path <- "../data/combined_Data.csv"
Combined_Data <- read.csv(data_path)
Combined_Data <- Combined_Data[,-1] # Remove Index

###################################################################################################
# Optimizing Kalman Filter Parameters  process error (HHt) and observation noise (GGt) using dlmMLE
###################################################################################################

# Initial state
N_initial <- Combined_Data$Population[!is.na(Combined_Data$Population)][1]

# Initialize vectors
years <- Combined_Data$Year
TT <- length(years)
# Constant Factor 1x1 matrix in measurement equation
Zt <- matrix(1, nrow = 1, ncol = 1)

# Initial state and covariance matrix
a0 <- as.double(N_initial) # Ensure initial state is of type 'double'
P0 <- matrix(1000) # Initial variance

# Define the state transition matrix (Tt)
phi <- mean(1 - (Combined_Data$Death_Count / Combined_Data$Population))
beta <- mean(Combined_Data$Birth_Count / Combined_Data$Population)
Tt <- matrix(phi + beta, nrow = 1, ncol = 1) # 1x1 Transition matrix

# Process error
GGt_Default<- matrix(0.03 * Combined_Data$Population[length(Combined_Data$Population) - 1], nrow = 1, ncol = 1)

# Observation noise
# Approximate with the variance of the yearly differences in population
yearly_differences <- diff(na.omit(Combined_Data$Population))
process_noise_estimate <- sqrt(var(yearly_differences))
HHt_Default <- matrix(process_noise_estimate)

########################
# Observed population
y <- Combined_Data$Population
yt <- rbind(y)
yt[is.na(yt)] <- 0 # Replace NA with 0 for Kalman filtering

dt <- ct <- matrix(0) # intercepts of the transition and measurement equations

#################################################################################################################
m2pDlm <- function(theta) {
  #   m0 = initial value of the state, z[0]
  #   C0 = variance of the initial value
  #   GG = matrix multiplier in the State Eq: z[t] = GG*z[t-1] + epsilon[t]
  #   W  = state process variance
  #   FF = matrix multiplier in the Obs. Eq: y[t] = FF*z[t] + eta[t] 
  #   V = observation error variance
  
  print(paste0("Iteration values, theta: ", theta))
  
  dlm(m0=a0, C0=1000, GG=Tt, W=exp(theta[1]), FF=1, V=exp(theta[2]))
}
f2pDlm <- dlmMLE(y=y, parm=c(log(GGt_Default^2), log(HHt_Default^2)),
                 build=m2pDlm)
f2pDlm$value
f2pDlm$convergence
mMle2pDlm <- m2pDlm(f2pDlm$par)

optimized_GGt <- matrix((exp(f2pDlm$par))[1])
optimized_HHt <- matrix((exp(f2pDlm$par))[2])
optimized_HHt; optimized_GGt

################################################################################################################
# Use filtering to get state estimates
filter2pDlm <- dlmFilter(y=y, mod=mMle2pDlm)
# Use smoothing to get state estimates
smooth2pDlm <- dlmSmooth(y=y, mod=mMle2pDlm)

# Variance of filtered values
fvar2pDlm <- unlist(dlmSvd2var(filter2pDlm$U.C,
                               filter2pDlm$D.C))
# CI of filtered estimates
zfCIl2pDlm <- filter2pDlm$m + qnorm(0.025, sd = sqrt(fvar2pDlm))
zfCIu2pDlm <- filter2pDlm$m + qnorm(0.975, sd = sqrt(fvar2pDlm))
# Variance of smoothed values
svar2pDlm <- unlist(dlmSvd2var(smooth2pDlm$U.S,
                               smooth2pDlm$D.S))
# CI of smoothed estimates
zsCIl2pDlm <- smooth2pDlm$s + qnorm(0.025, sd = sqrt(svar2pDlm))
zsCIu2pDlm <- smooth2pDlm$s + qnorm(0.975, sd = sqrt(svar2pDlm))

plot(1:TT, y,
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(z[t]),
     xlim = c(0,TT), ylim = c(min(y), max(y)+max(y)/5),
     las = 1)
polygon(c(0:TT, TT:0),
        c(zfCIl2pDlm,rev(zfCIu2pDlm)),
        col=rgb(0,0,0,0.2), border=FALSE)
lines(0:TT, filter2pDlm$m,
      lwd=2)
polygon(c(0:TT, TT:0),
        c(zsCIl2pDlm,rev(zsCIu2pDlm)),
        col=rgb(1,0.7,0.4,0.3), border=FALSE)
lines(0:TT, smooth2pDlm$s,
      col="darkgoldenrod1")
legend("top",
       legend = c("Obs.", "Filter. states", "Smooth. states"),
       pch = c(3, 19, NA, NA),
       col = c("blue", "black", "darkgoldenrod1"),
       lwd = c(1, 1, 2, 1.5), lty = c(3, 1, 1, 1),
       horiz=TRUE, bty="n", cex=0.9)