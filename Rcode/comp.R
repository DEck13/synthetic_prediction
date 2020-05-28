# require
require('tikzDevice')
# set seed
set.seed(2020)
# parameter setup
n = 40; T = 250; mu.alpha = 5; sigma.alpha = 1; sigma.X = 1; sigma = 1

# T* is the shock time point
Tstar <- floor(T / 2)
# random effects
eta <- rnorm(n + 1)
# with normalization
phi <- runif(n + 1, 0, .1)
theta <- rnorm(n + 1)
theta <- theta / as.vector(sqrt(crossprod(theta)))
beta <- rnorm(n + 1)
beta <- beta / as.vector(sqrt(crossprod(beta)))
beta.X <- rnorm(n + 1)
beta.X <- beta.X / as.vector(sqrt(crossprod(beta.X)))
alpha <- rnorm(n + 1, mean = mu.alpha, sd = sigma.alpha)
# T: time points; n: different time series
X <- matrix(0, nrow = T + 1, ncol = n + 1)
X[1, ] <- rnorm(n + 1)
Y <- matrix(0, nrow = T + 1, ncol = n + 1)
Y[1,] <- rnorm(n + 1)

# i represents for time -- T
for(i in 2:(T + 1)){
  X[i,] <- beta.X * X[i - 1,] + rnorm(n + 1, sd = sigma.X)
  if (i - 2 == Tstar) {
    Y[i, ] <- eta + phi * Y[i - 1,] + beta * X[i,]  + rnorm(n + 1, sd = sigma) + alpha
  } else {
    Y[i, ] <- eta + phi * Y[i - 1,] + beta * X[i,]  + rnorm(n + 1, sd = sigma)
    } 
}

# compute mu alpha hat
mu.alpha.hat <- rep(0, n)
se <- rep(0, n)
for (j in 2:(n + 1)){
  # OLS Regression
  m1 <- lm(Y[2:(T + 1), j] ~ 1 + Y[1:T, j] + X[2:(T + 1), j] + ifelse((2:(T + 1)) == Tstar + 2, 1, 0))
  # find alpha hat
  mu.alpha.hat[j - 1] <- coefficients(m1)[4]
  se[j - 1] <- summary(m1)$coef[4, 2]
}

alpha.hat <- mean(mu.alpha.hat)

# sum(mu.alpha.hat / se ^ 2) /  (sum(1 / se ^ 2))
# SCM method
scm <- function(X, Tstar) {
  # X is a list of covariates for disparate time series
  # X[[1]] should be the covariate of the time series to predict
  # X[[p]] for p = 2,...,n+1 are covariates for time series pool
  
  # T^* is a vector of shock-effects time points
  # shock effect point must be > 2
  
  # package for constrained optiminization
  require('Rsolnp')
  
  # number of time series for pool
  n <- length(X) - 1
  
  # covariate for time series for prediction
  X1 <- X[[1]][c(Tstar[1] - 1, Tstar[1]), , drop = FALSE]
  
  # covariates for time series pool
  X0 <- c()
  for (i in 1:n) {
    X0[[i]] <- X[[i + 1]][c(Tstar[i + 1] - 1, Tstar[i + 1]), ,drop = FALSE]
  }
  
  # objective function
  weightedX0 <- function(W) {
    # W is a vector of weight of the same length of X0
    n <- length(W)
    p <- ncol(X1)
    XW <- matrix(0, nrow = 2, ncol = p)
    for (i in 1:n) {
      XW <- XW + W[i] * X0[[i]]
    }
    norm <- as.numeric(crossprod(matrix(X1 - XW)))
    return(norm)
  }
  
  # constraint for W
  Wcons <- function(W) sum(W) - 1
  
  
  # optimization
  solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, eqB = 0, LB = rep(0, n), UB = rep(1, n))
}
# Xs
# Xs <- c()
# for (i in 1:(n + 1)) {
#  Xs[[i]] <- X[, i, drop = FALSE]
# }
# Wstar <- scm(X = Xs, Tstar = rep(127, n + 1))$par

# SCM Estimator
# alpha.hat <- sum(Wstar * mu.alpha.hat)
# alpha.hat
# sd(mu.alpha.hat) / sqrt(length(mu.alpha.hat))

# compute estimate for Tstar without post shock
m1 <- lm(Y[2:(Tstar + 1), 1] ~ 1 + Y[1:Tstar, 1] + X[2:(Tstar + 1), 1])
beta.hat <- coefficients(m1)
Ypred.original <- as.vector(cbind(1, Y[1:(Tstar + 1), 1], X[2:(Tstar + 2), 1]) %*% beta.hat)
Ypred.adjusted <- Ypred.original + alpha.hat * as.matrix(ifelse(2:(Tstar + 2) == Tstar + 2, 1, 0))

# set working directory
setwd('/Users/mac/Desktop/Research/Post-Shock Prediction/')
# figure setting
tikz('comp.tex', standAlone = TRUE, width = 8, height = 4)
# plot
par(mar = c(4, 4, 1, 2))
# original
matplot(2:(Tstar + 2), cbind(Y[2:(Tstar + 2), 1], Ypred.original, Ypred.adjusted), 
        xlab = 'Time $T$', lwd = 1.5,
        ylab = 'Response', lty = c(1, 2, 3),
        col = c('black', 'indianred1', 'magenta'),
        xlim = c(0, 165),
        type = 'l', ylim = c(min(mu.alpha.hat, Y), max(c(Ypred.adjusted, mu.alpha.hat))))
# add shock
segments(x0 = Tstar + 1, y0 = min(mu.alpha.hat, Y), y1 = max(c(Ypred.adjusted, mu.alpha.hat)), 
         col = 'deepskyblue', lty = 2, lwd = 2)
points(x = rep(Tstar, n), y = mu.alpha.hat, col = 'magenta')
# add arrows
arrows(x0 = Tstar + 7, y0 = Ypred.original[Tstar + 1], col = 'magenta',
       y1 = Ypred.adjusted[Tstar + 1],  code = 3, length = 0.07)
arrows(x0 = Tstar + 4, y0 = Y[Tstar + 2] - alpha[1],
       y1 = Y[Tstar + 2],  code = 3, length = 0.07)
# \hat{\alpha}
text(x = 150, y = mean(c(Ypred.original[Tstar + 1], max(Ypred.adjusted[Tstar + 1]))),
     col = 'magenta',
     labels = paste0('$\\displaystyle \\hat{\\alpha} = \\frac{1}{n}\\sum_{i=1}^{n} \\hat{\\alpha}_i= \\;$',
                     round(alpha.hat, digits = 2)))
text(x = 142, y = mean(Y[Tstar + 2] - alpha[1], Y[Tstar + 2]) + 1, col = 'black',
     labels = paste0('$\\alpha = ', round(alpha[1], digits = 2), '$'))
# shock time points
text(x = Tstar + 25, y = -4, 'shock time point $T^*_1 + 1$')
# forecasts
text(x = Tstar - 20, y = Y[Tstar + 2], '$y_{T_1^*+1}$')
text(x = Tstar - 30.5, y = Ypred.adjusted[Tstar + 1], '$\\hat{y}^{2}_{T_1^*+1}=\\hat{y}^{1}_{T_1^*+1} + \\hat{\\alpha}$',  col = 'magenta')
arrows(x0 = Tstar - 15, x1 = Tstar + 2, col = 'magenta',
       y0 = Ypred.adjusted[Tstar + 1], code = 2, length = 0.07)
arrows(x0 = Tstar - 15, x1 = Tstar + 2, 
       y0 = Y[Tstar + 2], code = 2, length = 0.07)
# legends
legend(x = 0, y = max(Ypred.adjusted), 
       col = c('indianred1', 'magenta'),
       lwd = 2, lty = c(2, 3),
       c('$\\hat{y}_{t}^{1}$ without shock effects', '$\\hat{y}_{t}^{2}$ with shock effects'))
# output
dev.off()



