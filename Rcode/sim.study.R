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
  solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, eqB = 0, 
        LB = rep(0, n), UB = rep(1, n), control = list(trace = 0))
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
    xi <- X[[i]][-1, ]
    
    # lag
    yilag <- Y[[i]][-(Ti + 1)]
    xilag <- X[[i]][-(Ti + 1), ]
    
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
  return(list(alphahat = alphahat, est = est, Tstar = Tstar, X = X, Y = Y, 
              lmod = lmod, res = res, Wstar = W, se = se))
  
}
# function based on ols object
boots <- function(B, ols) {
  # extract
  est <- ols$est
  res <- ols$res
  Tstar <- ols$Tstar
  Y <- ols$Y
  X <- ols$X
  lmod <- ols$lmod
  p <- ncol(X[[1]])
  n <- length(ols$Y) - 1
  
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
      xi <- X[[i]][-1,]
      xilag <- X[[i]][-(Ti + 1), ]
      
      # parameter
      etahati <- coef(lmodi)[1]
      alphahati <- coef(lmodi)[2]
      phihati <- coef(lmodi)[3]
      thetahati <- coef(lmodi)[4:(4 + p - 1)]
      betahati <- coef(lmodi)[(4 + p):(3 + 2 * p)]
      
      for (t in 1: Ti) {
        yitb <- c(yitb, etahati + alphahati * ifelse(t == Tstari + 1, yes = 1, no = 0) +
                    phihati * yitb[t] + thetahati %*% xi[t, ] + betahati %*% xilag[t, ] + resit[t])
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
# function based on ols object
boots.mix <- function(B, ols) {
  # extract
  est <- ols$est
  res <- ols$res
  Tstar <- ols$Tstar
  Y <- ols$Y
  X <- ols$X
  lmod <- ols$lmod
  p <- ncol(X[[1]])
  n <- length(ols$Y) - 1
  
  alphahat.adj.B <- c()
  alphahat.wadj.B <- c()
  alphahat.IVW.B <- c()
  
  for (b in 1:B) {
    Yb <- Xb <- Tstarb <- c();
    Yb[[1]] <- Y[[1]]; Xb[[1]] <- X[[1]]; Tstarb <- Tstar[1]
    
    # main difference from boots function
    pool <- sample(2:(n + 1), size = n, replace = TRUE)
    for (i in 1:n) {
      
      # setup
      index <- pool[i]
      Ti <- length(Y[[index]]) - 1
      coefitb <- coef(lmod[[index - 1]])
      resit <- base::sample(res[[index - 1]], replace = TRUE)
      yitb <- Y[[index]][1]
      Tstari <- Tstar[index]
      lmodi <- lmod[[index - 1]]
      xi <- X[[index]][-1,]
      xilag <- X[[index]][-(Ti + 1), ]
      
      # parameter
      etahati <- coef(lmodi)[1]
      alphahati <- coef(lmodi)[2]
      phihati <- coef(lmodi)[3]
      thetahati <- coef(lmodi)[4:(4 + p - 1)]
      betahati <- coef(lmodi)[(4 + p):(3 + 2 * p)]
      
      for (t in 1:Ti) {
        yitb <- c(yitb, etahati + alphahati * ifelse(t == Tstari + 1, yes = 1, no = 0) +
                    phihati * yitb[t] + thetahati %*% xi[t, ] + betahati %*% xilag[t, ] + resit[t])
      }
      Yb[[i + 1]] <- yitb
      Xb[[i + 1]] <- X[[index]]
      Tstarb[[i + 1]] <- Tstari
    }
    olsb <- ols.est.alphahat(Tstar = Tstarb, Y = Yb, X = Xb)
    alphahat.adj.B <- c(alphahat.adj.B, olsb$est['adj'])
    alphahat.wadj.B <- c(alphahat.wadj.B, olsb$est['wadj'])
    alphahat.IVW.B <- c(alphahat.IVW.B, olsb$est['IVW'])
  }
  return(list(adj = alphahat.adj.B, wadj = alphahat.wadj.B, IVW = alphahat.IVW.B))
}
# risk condition
risk.reduction <- function(means, vars) {
  rr.adj <- (means['wadj']) ^ 2 - vars['adj'] - (means['adj'] - means['wadj']) ^ 2
  rr.wadj <- (means['wadj']) ^ 2 - vars['wadj']
  rr.IVW <- (means['wadj']) ^ 2 - vars['IVW'] - (means['IVW'] - means['wadj']) ^ 2
  est <- c(rr.adj, rr.wadj, rr.IVW)
  names(est) <- c('adj', 'wadj', 'IVW')
  return(list(usable = ifelse(est > 0, yes = 1, no = 0), best = which.max(est)))
}
risk.reduction2 <- function(est, vars) {
  rr.adj <- (est['wadj']) ^ 2 - vars['adj'] - (est['adj'] - est['wadj']) ^ 2
  rr.wadj <- (est['wadj']) ^ 2 - vars['wadj']
  rr.IVW <- (est['wadj']) ^ 2 - vars['IVW'] - (est['IVW'] - est['wadj']) ^ 2
  rest <- c(rr.adj, rr.wadj, rr.IVW)
  names(rest) <- c('adj', 'wadj', 'IVW')
  return(list(usable = ifelse(rest > 0, yes = 1, no = 0), best = which.max(rest), rr = rest))
}
# nowcasting 
nowcast.alpha <- function(X, Y, Tstar, best = c('adj', 'wadj', 'IVW')) {
  # set up
  T1 <- length(Y[[1]]) - 1
  Tstar1 <- Tstar[1]
  y1 <- Y[[1]][-1]
  x1 <- X[[1]][-1, ]
  
  # lag
  y1lag <- Y[[1]][-(T1 + 1)]
  x1lag <- X[[1]][-(T1 + 1), ]
  
  # OLS
  lmod1 <- lm(y1 ~ 1 + y1lag + x1 + x1lag, subset = 1:Tstar1)
  
  # ols object
  ols <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)
  
  # design matrix
  design <- rbind(c(1, y1lag[Tstar1 + 1], x1[Tstar1 + 1, ], x1lag[Tstar1 + 1, ]))
  
  # forecast 1
  yhat1 <- design %*% coef(lmod1)
  
  # forecast 2
  yhat2 <- yhat1 + ols$est[best]
  
  # output
  return(list(yhat2 = yhat2, yhat1 = yhat1, alpha1est = ols$est[best]))
}

# double bootstrap
## main function for double bootstrap
ddboots <- function(B, K, ols, np) {
  
  # extract
  est <- ols$est
  res <- ols$res
  Tstar <- ols$Tstar
  Y <- ols$Y
  X <- ols$X
  lmod <- ols$lmod
  p <- ncol(X[[1]])
  n <- length(ols$Y) - 1
  
  # Distribution of Risk-Reduction
  RR.adj <- c()
  RR.wadj <- c()
  RR.IVW <- c()
  
  # FIRST BOOTSTRAP
  
  for (b in 1:B) {
    
    # setup
    Yb <- Xb <- Tstarb <- c();
    Yb[[1]] <- Y[[1]]; Xb[[1]] <- X[[1]]; Tstarb <- Tstar[1]
    
    # non-parametric bootstrap or not of the donor pool
    if (np == TRUE) {
      pool <- sample(2:(n + 1), size = n, replace = TRUE)
    } else {
      pool <- 2:(n + 1)
    }
    
    for (i in 1:n) {
      
      # setup
      index <- pool[i]
      Ti <- length(Y[[index]]) - 1
      coefitb <- coef(lmod[[index - 1]])
      resit <- base::sample(res[[index - 1]], replace = TRUE)
      yitb <- Y[[index]][1]
      Tstari <- Tstar[index]
      lmodi <- lmod[[index - 1]]
      xi <- X[[index]][-1,]
      xilag <- X[[index]][-(Ti + 1), ]
      
      # parameter
      etahati <- coef(lmodi)[1]
      alphahati <- coef(lmodi)[2]
      phihati <- coef(lmodi)[3]
      thetahati <- coef(lmodi)[4:(4 + p - 1)]
      betahati <- coef(lmodi)[(4 + p):(3 + 2 * p)]
      
      for (t in 1:Ti) {
        yitb <- c(yitb, etahati + alphahati * ifelse(t == Tstari + 1, yes = 1, no = 0) +
                    phihati * yitb[t] + thetahati %*% xi[t, ] + betahati %*% xilag[t, ] + resit[t])
      }
      Yb[[i + 1]] <- yitb
      Xb[[i + 1]] <- X[[index]]
      Tstarb[[i + 1]] <- Tstari
    }
    
    # OLS objects for b 
    # if np = FALSE, Xb = X
    olsb <- ols.est.alphahat(Tstar = Tstarb, Y = Yb, X = Xb)
    
    # SECOND BOOTSTRAP
    samp.b <- boots(B = K, ols = olsb)
    
    # estimate variance
    vars <- c()
    for (j in 1:3) {
      vars <- c(vars, var(samp.b[[j]], na.rm = TRUE))
    }
    names(vars) <- c('adj', 'wadj', 'IVW')
    
    # risk-reduction quantity
    RR.adj <- c(RR.adj, risk.reduction2(est = olsb$est, vars = vars)$rr['adj'])
    RR.wadj <- c(RR.wadj, risk.reduction2(est = olsb$est, vars = vars)$rr['wadj'])
    RR.IVW <- c(RR.IVW, risk.reduction2(est = olsb$est, vars = vars)$rr['IVW'])
  }
  
  return(list(adj = RR.adj, wadj = RR.wadj, IVW = RR.IVW))
}

# simulation study normal
sim.study.normal.gammaX <- function(mu.gamma.delta, mu.alpha, sigma, 
                                    sigma.alpha, sigma.delta.gamma, 
                                    p, B, n, scale, np = c(TRUE, FALSE)) {
  T <- round(rgamma(n = n + 1, shape = 15, scale = 10)) # Time Length
  T[which(T < 90)] <- 90
  Tstar <- c() # Shock Time Points
  for (t in T) {
    Tstar <- c(Tstar, sample((2 * p + 3 + 1):(t - 1), size = 1))
  }
  phi <- round(runif(n + 1, 0, 1), 3) # autoregressive parameters

  
  # construction of design matrix and shock effects
  X <- c()
  alpha <- c()
  delta <- c()
  gamma <- c()
  for (i in 1:(n + 1)) {
    Ti <- T[i]
    Tstari <- Tstar[i]
    X[[i]] <- matrix(rgamma(n = p * (Ti + 1), shape = 1, scale = scale), ncol = p, byrow = T) 
    # matrix(rnorm(n = (Ti + 1) * p), ncol = p, byrow = T)
    # parameter setup
    delta[[i]] <- matrix(rnorm(p, mean = mu.gamma.delta, sd = sigma.delta.gamma), nrow = 1)
    gamma[[i]] <- matrix(rnorm(p, mean = mu.gamma.delta, sd = sigma.delta.gamma), nrow = 1)
    epsilontildei <- rnorm(n = 1, sd = sigma.alpha)
    # alpha
    alpha <- c(alpha, mu.alpha + delta[[i]] %*% X[[i]][Tstari + 1, ] + 
                 gamma[[i]] %*% X[[i]][Tstari, ] + epsilontildei)
  }
  
  # E(alpha1)
  Ealpha1 <- mu.alpha + matrix(mu.gamma.delta, nrow = 1, ncol = p) %*% X[[1]][Tstar[1] + 1, ] + 
    matrix(mu.gamma.delta, nrow = 1, ncol = p) %*% X[[1]][Tstar[1], ]
  
  # E(alpha_adj)
  alphas.2.np1 <- c()
  for (i in 2:(n + 1)) {
    alphas.2.np1 <- c(alphas.2.np1, mu.alpha + matrix(mu.gamma.delta, nrow = 1, ncol = p) %*% X[[i]][Tstar[i] + 1, ] + 
                        matrix(mu.gamma.delta, nrow = 1, ncol = p) %*% X[[i]][Tstar[i], ])
  }
  Ealpha.adj <- mean(alphas.2.np1)
  
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
    thetai <- matrix(rnorm(p), nrow = 1)
    betai <- matrix(rnorm(p), nrow = 1)
    etai <- rnorm(1)
    
    yi <- yi0
    for (t in 2:(T[i] + 1)) {
      epsilonit <- rnorm(n = 1, sd = sigma)
      yi <- c(yi, etai + alphai * ifelse(t == Tstari + 2, yes = 1, no = 0) +
                phii * yi[t - 1] + thetai %*% xi[t, ] + betai %*% xi[t - 1, ] + epsilonit)
    }
    
    Y[[i]] <- yi
  }
  
  # estimates
  alphahat <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)$alphahat
  est <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)$est
  
  # ols object
  ols <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)
  
  # bootstrap samples
  if (np == TRUE) {
    bootsamp <- boots.mix(B = B, ols = ols)
  } else {
    bootsamp <- boots(B = B, ols = ols)
  }
  
  # Distance
  dist <- abs(est - alpha[1])
  
  # risk conditions
  means <- c()
  vars <- c()
  for (j in 1:3) {
    means <- c(means, mean(bootsamp[[j]], na.rm = TRUE))
    vars <- c(vars, var(bootsamp[[j]], na.rm = TRUE))
  }
  names(means) <- names(vars) <- c('adj', 'wadj', 'IVW')
  
  # yhats
  yhat2s <- sapply(c('adj', 'wadj', 'IVW'), function(d) nowcast.alpha(X = X, Y = Y, Tstar = Tstar, best = d)[[1]])
  yhat1 <- nowcast.alpha(X = X, Y = Y, Tstar = Tstar, best = 'wadj')[[2]]
  names(yhat1) <- 'no'
  names(yhat2s) <- c('adj', 'wadj', 'IVW')
  
  # truth
  truth <- ifelse(as.numeric(abs(Y[[1]][Tstar[1] + 2] - yhat1)) - as.numeric(abs(Y[[1]][Tstar[1] + 2] - yhat2s)) > 0,
                  yes = 1, no = 0) 
  names(truth) <- c('adj', 'wadj', 'IVW')
  
  # risk
  rmse <- sqrt(abs(Y[[1]][Tstar[1] + 2] - c(yhat1, yhat2s)) ^ 2)
  
  # empty result
  result <- c()
  
  for (z in 1:2) {
    if (z == 1) {
      # Bias
      bias <- c(abs(mean(bootsamp[[1]], na.rm = TRUE) - Ealpha.adj), 
                abs(mean(bootsamp[[2]], na.rm = TRUE) - Ealpha1))
      names(bias) <- c('E.adj', 'E.wadj')
      
      # comparison
      guess <- risk.reduction(means = means, vars = vars)$usable
      consistency <- ifelse(truth == guess, yes = 1, no = 0)
      
      # best consistency
      best.consistency <- ifelse(which.min(abs(Y[[1]][Tstar[1] + 2] - yhat2s)) == 
                                   risk.reduction(means = means, vars = vars)$best, yes = 1, no = 0)
      
      result[[1]] <- rbind(c(bias, dist, guess, consistency, best.consistency, rmse))
    } else {
      # Bias
      bias <- c(abs(est['adj'] - Ealpha.adj), 
                abs(est['wadj'] - Ealpha1))
      names(bias) <- c('E.adj', 'E.wadj')
      
      # comparison
      guess <- risk.reduction2(est = est, vars = vars)$usable
      consistency <- ifelse(truth == guess, yes = 1, no = 0)
      
      # best consistency
      best.consistency <- ifelse(which.min(abs(Y[[1]][Tstar[1] + 2] - yhat2s)) == 
                                   risk.reduction2(est = est, vars = vars)$best, yes = 1, no = 0)
      
      result[[2]] <- rbind(c(bias, dist, guess, consistency, best.consistency, rmse))
    }
  }
  names(result) <- c('boot', 'samp')
  return(result)
}

require('DescTools')
# simulation study normal
sim.study.normal.gammaX.ddboots <- function(mu.gamma.delta, mu.alpha, sigma, 
                                    sigma.alpha, sigma.delta.gamma, 
                                    p, B, K, n, scale, np = c(TRUE, FALSE)) {
  T <- round(rgamma(n = n + 1, shape = 15, scale = 10)) # Time Length
  T[which(T < 90)] <- 90
  Tstar <- c() # Shock Time Points
  for (t in T) {
    Tstar <- c(Tstar, sample((2 * p + 3 + 1):(t - 1), size = 1))
  }
  phi <- round(runif(n + 1, 0, 1), 3) # autoregressive parameters
  
  
  # construction of design matrix and shock effects
  X <- c()
  alpha <- c()
  delta <- c()
  gamma <- c()
  for (i in 1:(n + 1)) {
    Ti <- T[i]
    Tstari <- Tstar[i]
    X[[i]] <- matrix(rgamma(n = p * (Ti + 1), shape = 1, scale = scale), ncol = p, byrow = T) 
    # matrix(rnorm(n = (Ti + 1) * p), ncol = p, byrow = T)
    # parameter setup
    delta[[i]] <- matrix(rnorm(p, mean = mu.gamma.delta, sd = sigma.delta.gamma), nrow = 1)
    gamma[[i]] <- matrix(rnorm(p, mean = mu.gamma.delta, sd = sigma.delta.gamma), nrow = 1)
    epsilontildei <- rnorm(n = 1, sd = sigma.alpha)
    # alpha
    alpha <- c(alpha, mu.alpha + delta[[i]] %*% X[[i]][Tstari + 1, ] + 
                 gamma[[i]] %*% X[[i]][Tstari, ] + epsilontildei)
  }
  
  # E(alpha1)
  Ealpha1 <- mu.alpha + matrix(mu.gamma.delta, nrow = 1, ncol = p) %*% X[[1]][Tstar[1] + 1, ] + 
    matrix(mu.gamma.delta, nrow = 1, ncol = p) %*% X[[1]][Tstar[1], ]
  
  # E(alpha_adj)
  alphas.2.np1 <- c()
  for (i in 2:(n + 1)) {
    alphas.2.np1 <- c(alphas.2.np1, mu.alpha + matrix(mu.gamma.delta, nrow = 1, ncol = p) %*% X[[i]][Tstar[i] + 1, ] + 
                        matrix(mu.gamma.delta, nrow = 1, ncol = p) %*% X[[i]][Tstar[i], ])
  }
  Ealpha.adj <- mean(alphas.2.np1)
  
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
    thetai <- matrix(rnorm(p), nrow = 1)
    betai <- matrix(rnorm(p), nrow = 1)
    etai <- rnorm(1)
    
    yi <- yi0
    for (t in 2:(T[i] + 1)) {
      epsilonit <- rnorm(n = 1, sd = sigma)
      yi <- c(yi, etai + alphai * ifelse(t == Tstari + 2, yes = 1, no = 0) +
                phii * yi[t - 1] + thetai %*% xi[t, ] + betai %*% xi[t - 1, ] + epsilonit)
    }
    
    Y[[i]] <- yi
  }
  
  # estimates
  alphahat <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)$alphahat
  est <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)$est
  
  # ols object
  ols <- ols.est.alphahat(Tstar = Tstar, Y = Y, X = X)
  
  # double bootstrap
  bootsamp <- ddboots(B = B, K = K, ols = ols, np = np)

  # Distance
  dist <- abs(est - alpha[1])
  
  # yhats
  yhat2s <- sapply(c('adj', 'wadj', 'IVW'), function(d) nowcast.alpha(X = X, Y = Y, Tstar = Tstar, best = d)[[1]])
  yhat1 <- nowcast.alpha(X = X, Y = Y, Tstar = Tstar, best = 'wadj')[[2]]
  names(yhat1) <- 'no'
  names(yhat2s) <- c('adj', 'wadj', 'IVW')
  
  # risk
  rmse <- sqrt(abs(Y[[1]][Tstar[1] + 2] - c(yhat1, yhat2s)) ^ 2)
  
  # empty result
  result <- c()
  
  # Bias
  bias <- c(abs(est['adj'] - Ealpha.adj), 
            abs(est['wadj'] - Ealpha1))
  names(bias) <- c('E.adj', 'E.wadj')
  
  # truth
  truth <- ifelse(as.numeric(abs(Y[[1]][Tstar[1] + 2] - yhat1)) - as.numeric(abs(Y[[1]][Tstar[1] + 2] - yhat2s)) > 0,
                  yes = 1, no = 0) 
  names(truth) <- c('adj', 'wadj', 'IVW')
  # guess
  guess <- c()
  LBs <- c()
  UBs <- c()
  for (k in 1:3) {
    LBs[k] <- quantile(bootsamp[[k]], probs = 0.025)
    UBs[k] <- quantile(bootsamp[[k]], probs = 0.975)
    if (LBs[k] >= 0) {
      guess <- c(guess, 1)
    } else {
      guess <- c(guess, 0)
    }
  }
  names(LBs) <- names(LBs) <- names(guess) <- c('adj', 'wadj', 'IVW')
  consistency <- ifelse(truth == guess, yes = 1, no = 0)
  
  # best consistency
  overlaps <- c()
  for (l in 1:2) {
    for (f in (l + 1):3) {
      overlaps <- c(overlaps, c(LBs[l], UBs[l]) %overlaps% c(LBs[f], UBs[f]))
    }
  }
  # if there is any overlaps, the ranking is not effective
  best <- c()
  if (TRUE %in% overlaps) {
    best <- NA
  } else {
    best <- which.max(LBs)
  }
  # best consistency
  if (is.na(best) == TRUE) {
    best.consistency <- NA
  } else {
    best.consistency <- ifelse(which.min(abs(Y[[1]][Tstar[1] + 2] - yhat2s)) == 
                                 best, yes = 1, no = 0)
  }
  
  # output result
  result <- rbind(c(bias, dist, guess, consistency, best.consistency, rmse))
  return(result)
}



# collect results
frame <- c()
for (z in 1:2) {
  result <- c()
  for (i in 1:nrow(sim_params)) {
  table <- output[[i]]
  # means and sds
  means <- apply(table[0:199 * 2 + z, ], 2, function(x) mean(x))
  sds <- apply(table[0:199 * 2 + z, ], 2, function(x) sd(x))
  result.i <- c()
  for (j in 1:16) {
    result.i <- cbind(result.i, paste0(round(means[j], digits = 3), 
                                         ' (', round(sds[j] / sqrt(50), 
                                                     digits = 3), ')'))
    }
  result <- rbind(result, result.i)
  }
  result <- cbind(sim_params[, c(2,1)], result)
  rownames(result) <- 1:20
  frame[[z]] <- result
}
frame <- do.call('rbind', frame)
# results
riskprop <- frame[, c(1:4, 8:14)]
pred <- frame[, -c(3:4, 8:14)]
# xtable
require('xtable')
xtable(riskprop)
xtable(pred[1:20,])


# MC
library("parallel")
library("doParallel")
library("foreach")
# 8 cores -- use 7
ncores <- detectCores() - 1
registerDoParallel(cores = ncores)
set.seed(2020)
RNGkind("L'Ecuyer-CMRG")
nsim <- 200
# parameter setup
sigma <- c(1, 5, 10, 25, 100)
sigma.alphas <- c(1, 5, 10, 25, 100)
sim_params <- expand.grid(list(sigma.alphas = sigma.alphas, sigma = sigma))

# simulation time
system.time(
  output <- lapply(1:nrow(sim_params), FUN = function(j) {
    # parameters
    sigma.alpha <- sim_params[j, 1]
    sigma <- sim_params[j, 2]
    # %do% evaluates sequentially
    # %dopar% evaluates in parallel
    # .combine results
    out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
      # result
      study <- sim.study.normal.gammaX(mu.gamma.delta = 2, 
                                       mu.alpha = 10, sigma = sigma, 
                                       sigma.alpha = sigma.alpha,
                                       sigma.delta.gamma = 1, 
                                       p = 13, B = 500, scale = 10, 
                                       n = 10, np = FALSE)
      result <- rbind(study$boot, study$samp)
      return(result)
    }
    # return results
    out
  })
)

# store results
# load packages
require('readxl')
require('writexl')
setwd('/Users/mac/Desktop/Research/Post-Shock Prediction/')
write_xlsx(lapply(output, as.data.frame), 'parametricssigma.xlsx')



# collect results
frame <- c()
for (z in 1:2) {
  result <- c()
  for (i in 1:nrow(sim_params)) {
    table <- output[[i]]
    # means and sds
    means <- apply(table[0:199 * 2 + z, ], 2, function(x) mean(x))
    sds <- apply(table[0:199 * 2 + z, ], 2, function(x) sd(x))
    result.i <- c()
    for (j in 1:16) {
      result.i <- cbind(result.i, paste0(round(means[j], digits = 3), 
                                         ' (', round(sds[j] / sqrt(50), 
                                                     digits = 3), ')'))
    }
    result <- rbind(result, result.i)
  }
  result <- cbind(sim_params[, c(2,1)], result)
  rownames(result) <- 1:25
  frame[[z]] <- result
}
frame <- do.call('rbind', frame)
# results
riskprop <- frame[, c(1:4, 8:14)]
pred <- frame[, -c(3:4, 8:14)]
# xtable
require('xtable')
xtable(riskprop[-(1:25),])
xtable(pred[-(1:25),])





# MC
library("parallel")
library("doParallel")
library("foreach")
# 8 cores -- use 7
ncores <- detectCores() - 1
registerDoParallel(cores = ncores)
set.seed(2020)
RNGkind("L'Ecuyer-CMRG")
nsim <- 50


# parameter setup
ns <- c(5, 10, 15, 25)
sigma.alphas <- c(1, 5, 10, 25, 100)
sim_params <- expand.grid(list(sigma.alphas = sigma.alphas, ns = ns))

# simulation time
system.time(
  output <- lapply(1:nrow(sim_params), FUN = function(j) {
    # parameters
    sigma.alpha <- sim_params[j, 1]
    n <- sim_params[j, 2]
    # %do% evaluates sequentially
    # %dopar% evaluates in parallel
    # .combine results
    out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
      # result
      study <- sim.study.normal.gammaX.ddboots(mu.gamma.delta = 2, 
                                               mu.alpha = 10, sigma = 1, 
                                               sigma.alpha = sigma.alpha, 
                                               sigma.delta.gamma = 1, 
                                               p = 13, B = 50, scale = 10, 
                                               n = n, K = 50, np = FALSE)
      return(study)
    }
    # return results
    out
  })
)
# store results
# load packages
require('readxl')
require('writexl')
setwd('/Users/mac/Desktop/Research/Post-Shock Prediction/')
write_xlsx(lapply(output, as.data.frame), 'ddparametricssigma.xlsx')

result <- c()
for (i in 1:nrow(sim_params)) {
  table <- output[[i]]
  # means and sds
  means <- apply(table, 2, function(x) mean(x))
  sds <- apply(table, 2, function(x) sd(x))
  result.i <- c()
  for (j in 1:16) {
    result.i <- cbind(result.i, paste0(round(means[j], digits = 3), 
                                       ' (', round(sds[j] / sqrt(50), 
                                                   digits = 3), ')'))
  }
  result <- rbind(result, result.i)
}
result <- cbind(sim_params[, c(2,1)], result)
# results
riskprop <- result[, c(1:4, 8:13)]
pred <- result[, -c(3:4, 8:14)]
require('xtable')
xtable(riskprop)


# MC
library("parallel")
library("doParallel")
library("foreach")
# 8 cores -- use 7
ncores <- detectCores() - 1
registerDoParallel(cores = ncores)
set.seed(2020)
RNGkind("L'Ecuyer-CMRG")
nsim <- 50

# parameter setup
sigma <- c(1, 5, 10, 25, 100)
sigma.alphas <- c(1, 5, 10, 25, 100)
sim_params <- expand.grid(list(sigma.alphas = sigma.alphas, sigma = sigma))

# simulation time
system.time(
  output <- lapply(1:nrow(sim_params), FUN = function(j) {
    # parameters
    sigma.alpha <- sim_params[j, 1]
    sigma <- sim_params[j, 2]
    # %do% evaluates sequentially
    # %dopar% evaluates in parallel
    # .combine results
    out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
      # result
      study <- sim.study.normal.gammaX.ddboots(mu.gamma.delta = 2, 
                                               mu.alpha = 10, sigma = sigma, 
                                               sigma.alpha = sigma.alpha, 
                                               sigma.delta.gamma = 1, 
                                               p = 13, B = 50, scale = 10, 
                                               n = 10, K = 50, np = FALSE)
      return(study)
    }
    # return results
    out
  })
)
# store results
# load packages
require('readxl')
require('writexl')
setwd('/Users/mac/Desktop/Research/Post-Shock Prediction/')
write_xlsx(lapply(output, as.data.frame), 'ddparametricssigma.xlsx')

result <- c()
for (i in 1:nrow(sim_params)) {
  table <- output[[i]]
  # means and sds
  means <- apply(table, 2, function(x) mean(x))
  sds <- apply(table, 2, function(x) sd(x))
  result.i <- c()
  for (j in 1:16) {
    result.i <- cbind(result.i, paste0(round(means[j], digits = 3), 
                                       ' (', round(sds[j] / sqrt(50), 
                                                   digits = 3), ')'))
  }
  result <- rbind(result, result.i)
}
result <- cbind(sim_params[, c(2,1)], result)
# results
riskprop <- result[, c(1:4, 8:13)]
pred <- result[, -c(3:4, 8:14)]
require('xtable')
xtable(riskprop)















# MC
library("parallel")
library("doParallel")
library("foreach")
# 8 cores -- use 7
ncores <- detectCores() - 1
registerDoParallel(cores = ncores)
set.seed(2020)
RNGkind("L'Ecuyer-CMRG")
nsim <- 100


# parameter setup
ns <- c(5, 10, 15, 25)
sigma.alphas <- c(5, 10, 25, 50, 100)
sim_params <- expand.grid(list(sigma.alphas = sigma.alphas, ns = ns))

# simulation time
system.time(
  output <- lapply(1:nrow(sim_params), FUN = function(j) {
    # parameters
    sigma.alpha <- sim_params[j, 1]
    n <- sim_params[j, 2]
    # %do% evaluates sequentially
    # %dopar% evaluates in parallel
    # .combine results
    out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
      # result
      study <- sim.study.normal.gammaX(mu.gamma.delta = 1, 
                                       mu.alpha = 2, sigma = 10, 
                                       sigma.alpha = sigma.alpha,
                                       sigma.delta.gamma = 0.5, 
                                       p = 13, B = 200, scale = 2, 
                                       n = n, np = FALSE)
      result <- study$samp
      return(result)
    }
    # return results
    out
  })
)

# store results
# load packages
require('readxl')
require('writexl')
setwd('/Users/mac/Desktop/Research/Post-Shock Prediction/')
write_xlsx(lapply(output, as.data.frame), 'parametricnsigma.xlsx')

result <- c()
for (i in 1:nrow(sim_params)) {
  table <- output[[i]]
  # means and sds
  means <- apply(table, 2, function(x) mean(x))
  sds <- apply(table, 2, function(x) sd(x))
  result.i <- c()
  for (j in 1:16) {
    result.i <- cbind(result.i, paste0(round(means[j], digits = 3), 
                                       ' (', round(sds[j] / sqrt(100), 
                                                   digits = 3), ')'))
  }
  result <- rbind(result, result.i)
}
result <- cbind(sim_params[, c(2,1)], result)
# results
riskprop <- result[, c(1:4, 8:14)]
pred <- result[, -c(3:4, 8:14)]
require('xtable')
xtable(riskprop)
xtable(pred)
