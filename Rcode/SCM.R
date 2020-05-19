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

# this function returns the W^* estimated by synthetic control method (SCM)
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
  
  # optimization
  solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, LB = rep(0, n), UB = rep(1, n))
}



########### Simulation Example

# generate a list of X; n = 40
# different Tis
X <- c()
Tis <- c()
for (i in 1:41) {
 Ti <- sample(300:500, size = 1)
 Tis <- c(Tis, Ti)
 X[[i]] <- cbind(rnorm(Ti), rnorm(Ti))
}
# generate a vector of T^*
Tstar <- c()
for (i in 1:41) {
  Tstar <- c(Tstar, sample(3:Tis[i], size = 1))
}

# covariate for time series for prediction
X1 <- X[[1]][c(Tstar[1] - 1, Tstar[1]),]

# number of time series for pool
n <- length(X) - 1
# covariates for time series pool
X0 <- c()
for (i in 1:n) {
  X0[[i]] <- X[[i + 1]][c(Tstar[i + 1] - 1, Tstar[i + 1]),]
}

# output
scm(X = X, Tstar = Tstar)$par

