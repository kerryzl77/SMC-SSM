library(ggplot2)
library(data.table)
library(FKF)

# visualize
data <- data.table(x = seq(from = 1912, to = 1971, by = 1),
                   y = as.numeric(nhtemp))
ggplot(data, aes(x = x, y = y)) +
  geom_line() +
  xlab("year") + xlab("temperature") +
  ggtitle("Measured Yearly Average Temperature in New Haven") +
  theme_bw()

# let a0 be our initial guess for  θ - temperature in the first period y[1]. 
# We let P0 be our best guess for the variance of  θ
y <- nhtemp
a0 <- y[1]
P0 <- matrix(1)
dt <- ct <- matrix(0)
Zt <- Tt <- matrix(1)
# Estimate parameters:
fit.fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
                   GGt = var(y, na.rm = TRUE) * .5),
                 fn = function(par, ...)
                   -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                 yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
                 Zt = Zt, Tt = Tt)


# recover values
HHt <- as.numeric(fit.fkf$par[1])
GGt <- as.numeric(fit.fkf$par[2])
HHt; GGt

y_fkf <- fkf(a0, P0, dt, ct, Tt, Zt,
             HHt = matrix(HHt), GGt = matrix(GGt),
             yt = rbind(y))

data <- data.table(x = seq(from = 1912, to = 1971, by = 1),
                   y = as.numeric(nhtemp),
                   y_kalman = as.numeric(y_fkf$att))


ggplot(data, aes(x = x, y = y)) +
  geom_line() +
  geom_line(data = data, aes(x = x, y = y_kalman), col = "blue") +
  xlab("year") + xlab("temperature") +
  ggtitle("Measured Yearly Average Temperature in New Haven (with Kalman Filter)") +
  theme_bw()
