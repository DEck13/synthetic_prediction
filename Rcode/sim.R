# Bootstrap Simulation
set.seed(2014)
# parameter setup
n <- 10 # pool size
T <- sample(50:100, size = n + 1) # Time Length
Tstar <- c() # Shock Time Points
for (t in T) {
 Tstar <- c(Tstar, sample(3:t, size = 1))
}
phi <- round(runif(n + 1, 0, 1), 3) # autoregressive parameters

# parameter 
mu.alpha <- 5; sigma.alpha <- 0.1; sigma <- 1

# construction of design matrix and shock effects
X <- c()
alpha <- c()
delta <- c()
gamma <- c()
for (i in 1:(n + 1)) {
  Ti <- T[i]
  Tstari <- Tstar[i]
  X[[i]] <- #rnorm(Ti + 1, sd = 10) 
    rgamma(Ti + 1, shape = 1, scale = 10)
  # parameter setup
  delta[i] <- rnorm(1, mean = 5, sd = 0.1)
  gamma[i] <- rnorm(1, mean = 5, sd = 0.1)
  epsilontildei <- rnorm(n = 1, sd = sigma.alpha)
  # alpha
  alpha <- c(alpha, mu.alpha + delta[i] * X[[i]][Tstari + 1] + 
               gamma[i] * X[[i]][Tstari] + epsilontildei)
}

mu.alpha + delta[1] * X[[i]][Tstar[1] + 1] + gamma[1] * X[[i]][Tstar[1]]
mu.alpha + delta[1] * X[[1]][Tstar[1] + 1] + gamma[1] * X[[1]][Tstar[1]]
mu.alpha + delta[5] * X[[5]][Tstar[5] + 1] + gamma[5] * X[[5]][Tstar[5]]


# generation of yit
Y <- c()
for (i in 1:(n + 1)) {
  
  # initial value
  yi0 <- rnorm(1)
  
  # setup
  Tstari <- Tstar[i]
  alphai <- alpha[i]
  phii <- phi[i]
  xi <- X[[i]]
  
  # parameter setup
  thetai <- rnorm(1)
  betai <- rnorm(1)
  etai <- rnorm(1)
  
  yi <- yi0
  for (t in 2:(T[i] + 1)) {
    epsilonit <- rnorm(n = 1, sd = sigma)
    yi <- c(yi, etai + alphai * ifelse(t == Tstari + 2, yes = 1, no = 0) +
              phii * yi[t - 1] + thetai * xi[t] + betai * xi[t - 1] + epsilonit)
  }
  
  Y[[i]] <- yi
}

# OLS
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
  
  # covariate for time series for prediction
  X1 <- X[[1]][c(Tstar[1], Tstar[1] + 1), , drop = FALSE]
  
  # covariates for time series pool
  X0 <- c()
  for (i in 1:n) {
    X0[[i]] <- X[[i + 1]][c(Tstar[i + 1], Tstar[i + 1] + 1), ,drop = FALSE]
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
  solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, 
        eqB = 0, LB = rep(0, n), UB = rep(1, n), control = list(trace = 0))
}
ols.est.alphahat <- function(Tstar, Y, X) {
  
  n <- length(Y) - 1
  T <- c()
  for (i in 1:(n + 1)) {
    T <- c(T, length(Y[[i]]) - 1)
  }
  
  # empty
  alphahat <- c()
  se <- c()
  res <- c()
  lmod <- c()
  
  
  for (i in 2:(n + 1)) {
    
    # set up
    Ti <- T[i]
    Tstari <- Tstar[i]
    yi <- Y[[i]][-1]
    xi <- X[[i]][-1]
    
    # lag
    yilag <- Y[[i]][-(Ti + 1)]
    xilag <- X[[i]][-(Ti + 1)]
    
    # OLS
    lmodi <- lm(yi ~ 1 + ifelse(1:Ti == Tstari + 1, yes = 1, no = 0) + 
                  yilag + xi + xilag)
    
    # find estimates
    alphahat <- c(alphahat, coef(lmodi)[2])
    se <- c(se, summary(lmodi)$coef[2, 2])
    res[[i - 1]] <- residuals(lmodi)
    lmod[[i - 1]] <- lmodi
  }
  
  # uname
  names(alphahat) <- names(se) <- NULL
  
  # adjustment estimator
  alphahat.adj <- mean(alphahat)
  
  # weighted adjustment estimator
  if (is.matrix(X[[1]]) == FALSE) {
    for (i in 1:(n + 1)) {
      X[[i]] <- as.matrix(X[[i]])
    }
  }
  # Weights
  W <- round(scm(X = X, Tstar = Tstar)$par, digits = 3)
  # Computation
  alphahat.wadj <- as.numeric(crossprod(W, alphahat))
  
  # Inverse-Variance Estimator
  alphahat.IVW <- sum(alphahat / se ^ 2) /  (sum(1 / se ^ 2))
  
  est <- c(alphahat.adj, alphahat.wadj, alphahat.IVW)
  names(est) <- c('adj', 'wadj', 'IVW')
  
  # output
  return(list(alphahat = alphahat, est = est, Tstar = Tstar, X = X, Y = Y, lmod = lmod, res = res))
  
}

# estimates
alphahat <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)$alphahat
est <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)$est

# ols object
ols <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)
ols$est

# Bootstrap
B <- 10000

# function based on ols object
boots <- function(B, ols) {
  # extract
  est <- ols$est
  res <- ols$res
  Tstar <- ols$Tstar
  Y <- ols$Y
  X <- ols$X
  lmod <- ols$lmod
  
  alphahat.adj.B <- c()
  alphahat.wadj.B <- c()
  alphahat.IVW.B <- c()
  
  for (b in 1:B) {
    Yb <- c(); Yb[[1]] <- Y[[1]]
    for (i in 2:(n + 1)) {
      
      # setup
      Ti <- length(Y[[i]]) - 1
      coefitb <- coef(lmod[[i - 1]])
      resit <- base::sample(res[[i - 1]], replace = TRUE)
      yitb <- Y[[i]][1]
      Tstari <- Tstar[i]
      lmodi <- lmod[[i - 1]]
      xi <- X[[i]][-1]
      xilag <- X[[i]][-(Ti + 1)]
      
      # parameter
      etahati <- coef(lmodi)[1]
      alphahati <- coef(lmodi)[2]
      phihati <- coef(lmodi)[3]
      thetahati <- coef(lmodi)[4]
      betahati <- coef(lmodi)[5]
      
      for (t in 1: Ti) {
        yitb <- c(yitb, etahati + alphahati * ifelse(t == Tstari + 1, yes = 1, no = 0) +
          phihati * yitb[t] + thetahati * xi[t] + betahati * xilag[t] + resit[t])
      }
      Yb[[i]] <- yitb
    }
    olsb <- ols.est.alphahat(Tstar = Tstar, Y = Yb, X = X)
    alphahat.adj.B <- c(alphahat.adj.B, olsb$est['adj'])
    alphahat.wadj.B <- c(alphahat.wadj.B, olsb$est['wadj'])
    alphahat.IVW.B <- c(alphahat.IVW.B, olsb$est['IVW'])
  }
  return(list(adj = alphahat.adj.B, wadj = alphahat.wadj.B, IVW = alphahat.IVW.B))
}


# bootstrap samples
bootsamp <- boots(B = 500, ols = ols)

par(mfrow = c(2,2))
hist(bootsamp[[1]])
summary(bootsamp[[1]])

hist(bootsamp[[2]])
summary(bootsamp[[2]])

hist(bootsamp[[3]])
summary(bootsamp[[3]])


# experiment 
Bs <- c(100, 1000, 10000)
means <- matrix(0, nrow = 3, ncol = 3)
vars <- matrix(0, nrow = 3, ncol = 3)
for (i in 1:3) {
  bootsamp <- boots(B = Bs[i], ols = ols)
  for (j in 1:3) {
    means[i, j] <- mean(bootsamp[[j]])
    vars[i, j] <- var(bootsamp[[j]])
  }
}
colnames(means) <- colnames(vars) <- c('adj', 'wadj', 'IVW')
rownames(means) <- rownames(vars) <- c('B = 100', 'B = 1000', 'B = 10000')


# risk-reduction conditions evaluation
risk.reduction <- function(means, vars) {
  rr.adj <- (means['wadj']) ^ 2 - vars['adj'] - (means['adj'] - means['wadj']) ^ 2
  rr.wadj <- (means['wadj']) ^ 2 - vars['wadj']
  rr.IVW <- (means['wadj']) ^ 2 - vars['IVW'] - (means['IVW'] - means['wadj']) ^ 2
  est <- c(rr.adj, rr.wadj, rr.IVW)
  names(est) <- c('adj', 'wadj', 'IVW')
  return(list(usable = which(est > 0), best = which.max(est)))
}
risk.reduction(means = means[3,], vars = vars[3,])
# in this case adj is the best
