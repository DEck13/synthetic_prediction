################## methodology section
# load packages
library("data.table")
library("dplyr")
library("tseries")
library("quantmod")
library('Rsolnp')
library('msos')
library('tikzDevice')
library('xtable')

# load packages
require('forecast')

# this function returns the W^* estimated by synthetic control method (SCM)
scm <- function(X, Tstar, scale = FALSE) {
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
  X1 <- X[[1]][Tstar[1] + 1, , drop = FALSE]
  
  # covariates for time series pool
  X0 <- c()
  for (i in 1:n) {
    X0[[i]] <- X[[i + 1]][Tstar[i + 1] + 1, , drop = FALSE]
  }
  
  if (scale == TRUE) {
    dat <- rbind(X1, do.call('rbind', X0))
    dat <- apply(dat, 2, function(x) scale(x, center = TRUE, scale = TRUE))
    X1 <- dat[1, , drop = FALSE]
    X0 <- c()
    for (i in 1:n) {
      X0[[i]] <- dat[i + 1, , drop = FALSE]
    }
  }
  
  
  # objective function
  weightedX0 <- function(W) {
    # W is a vector of weight of the same length of X0
    n <- length(W)
    p <- ncol(X1)
    XW <- matrix(0, nrow = 1, ncol = p)
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

# squared error loss
sel <- function(yhat, y) {
  (y - yhat) ^ 2
}

# taSPA test for comparison of two models
taSPA.mhfc <- function(ell, d, B, bw = 4) {
  # ell is the block length
  # d is the matrix of loss differential
  # B is the bootstrap size
  # bw is the bandwidth parameter for computing HAC estimator
  
  # compute dbar
  d.bar <- mean(matrix(apply(d, 1, mean)))
  
  # moving-block bootstrap
  Tt <- nrow(d); H <- ncol(d)
  # T = ell * K
  K <- Tt / ell
  # Bootstrap replication
  t.aSPA.B <- c()
  for (b in 1:B) {
    # uniform draw from 1,...,T - ell +1
    Iks <- sample(1:(Tt - ell + 1), K, replace = TRUE)
    tau <- matrix(sapply(Iks, function(x) x:(x + ell - 1)))
    tau <- as.numeric(tau)
    # moving-block bootstrap
    d.b <- d[tau, ]
    d.b.bar <- mean(matrix(apply(d.b, 1, mean)))
    # compute moving block bootstrap standard error
    d.b.w <- apply(d.b, 1, mean)
    xi.b <- 0
    for (k in 1:K) {
      xi.b <- xi.b + (sum(d.b.w[(k - 1) * ell + 1:ell]) - mean(d.b.w)) ^ 2 / ell
    }
    xi.b <- xi.b / K
    # compute t.aSPA statistic
    t.aSPA.B[b] <- sqrt(Tt) * (d.b.bar - d.bar) / xi.b
  }
  # compute t.aSPA statistic
  ## compute Heteroskedasticity and autocorrelation Consistent (HAC) estimator
  require('tsapp')
  Omega <- HAC(d, method = 'Quadratic Spectral', bw = bw)
  w <- matrix(1 / H, nrow = H)
  xi.hat <- sqrt(t(w) %*% Omega %*% w)
  t.aSPA <- as.numeric(sqrt(Tt) * d.bar / xi.hat)
  p <- mean(t.aSPA < t.aSPA.B)
  # return output
  return(p)
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
    yilag <- Y[[i]][-Ti]
    
    # OLS
    lmodi <- lm(yi ~ 1 + ifelse(1:Ti == Tstari + 1, yes = 1, no = 0) + 
                  yilag + xi)
    
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

# functions that return alpha.hat and synthetic weights
ps.indic.W.permanent <- function(Tstar, Y, X, K, H, Ts, ell, B, bw, sig.levl = .05, 
                                 sc.x = NA,
                                 retro = TRUE, selfW = NA, scale = FALSE) {
  
  n <- length(Y) - 1
  
  # empty
  ps <- c()
  
  if (retro == TRUE) k <- 1 else 2
  
  for (i in k:(n + 1)) {
    
    # set up
    Ti <- Ts[i]
    Ki <- K[i]
    Tstari <- Tstar[i]
    
    if (i > 2) {
      sym <- ols.est.alphahat(Tstar = Tstar[c(i, 2:(i - 1))], X = X[c(i, 2:(i - 1))], Y = Y[c(i, 2:(i - 1))])
      alpha.wadj <- sym$est[2]
    }
    
    m1.L.i <- matrix(NA, nrow = Ti - Ki - H - 1, ncol = H)
    m2.L.i <- matrix(NA, nrow = Ti - Ki - H - 1, ncol = H)
    
    # compute losses
    for (h in 1:H) {
      for (t in 1:(Ti - Ki - H - 1)) {
        
        yi <- Y[[i]][-1][(t + H - h + 1):(t + Ki + H - h)]
        xi <- X[[i]][-1, ][(t + H - h + 1):(t + Ki + H - h),]
        
        # lag
        yilag <- Y[[i]][-(Ti + 1)][(t + H - h + 1):(t + Ki + H - h)]
        
        # special case that needs covariates without lags
        sc <- c()
        
        if (length(sc.x) == 1) {
          sc <- 0
        } else {
          if (i %in% sc.x[[1]]) {
            sc <- 1
          } else {
            sc <- 0
          }
        }
        
        if (sc == 0) {
          # OLS
          lmodi.adj <- lm(yi ~ 1 + ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0) + 
                            yilag + xi)
          lmodi.unadj <- lm(yi ~ 1 + yilag + xi)
          
          # beta.hat
          beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
          beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
          beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
          
          yhat.adj <- tail(yi, 1)
          yhat.unadj <- tail(yi, 1)
          
          for (j in 1:h) {
            yhat.adj.h <- matrix(c(1, ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0),
                                   yhat.adj[j], X[[i]][-1, ][t + Ki + H - h + j, ]))
            yhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, yhat.unadj[j], X[[i]][-1, ][t + Ki + H - h + j, ]))
            
            if (t + Ki + H - h == Tstari & i >= 3 & h == 1) {
              yhat.adj.h <- yhat.adj.h + alpha.wadj
            }
            yhat.adj <- c(yhat.adj, yhat.adj.h)
            yhat.unadj <- c(yhat.unadj, yhat.unadj.h)
          }
          
        } else {
          sc.xx <- sc.x[[2]][[which(sc.x[[1]] == i)]][(t + H - h + 1):(t + Ki + H - h)]
          
          # OLS
          lmodi.adj <- lm(yi ~ 1 + ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0) + 
                            yilag + xi + sc.xx)
          lmodi.unadj <- lm(yi ~ 1 + yilag + xi + sc.xx)
          
          # beta.hat
          beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
          beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
          beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
          
          yhat.adj <- tail(yi, 1)
          yhat.unadj <- tail(yi, 1)
          
          for (j in 1:h) {
            yhat.adj.h <- matrix(c(1, ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0),
                                   yhat.adj[j], X[[i]][-1, ][t + Ki + H - h + j, ],
                                   sc.xx[[2]][[which(sc.xx[[1]] == i)]][t + Ki + H - h + j - k]))
            yhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, yhat.unadj[j], X[[i]][-1, ][t + Ki + H - h + j, ],
                                                        sc.xx[[2]][[which(sc.xx[[1]] == i)]][t + Ki + H - h + j - k]))
            
            if (t + Ki + H - h == Tstari & i >= 3 & h == 1) {
              yhat.adj.h <- yhat.adj.h + alpha.wadj
            }
            
            yhat.adj <- c(yhat.adj, yhat.adj.h)
            yhat.unadj <- c(yhat.unadj, yhat.unadj.h)
          }
        }
        
        
        # losses
        m1.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.adj[h + 1])
        m2.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.unadj[h + 1])
      }
    }
    # loss differential
    d <- m2.L.i - m1.L.i
    ps[i] <- taSPA.mhfc(ell = ell, d = d, B = B, bw = bw)
  }
  
  # weights
  if (is.matrix(X[[1]]) == FALSE) {
    for (i in 1:(n + 1)) {
      X[[i]] <- as.matrix(X[[i]])
    }
  }
  # Weights
  W <- round(scm(X = X, Tstar = Tstar, scale = scale)$par, digits = 3)
  
  if (is.na(selfW)[1] == FALSE) {
    W <- selfW
  }
  
  # test
  Is <- ifelse(ps <= .05, yes = 1, no = 0)
  
  # output
  return(list(ps = ps, W = W, Is = Is))
}

# functions that return alpha.hat and synthetic weights for decaying shock effects
ps.indic.W.dynamic <- function(Tstar, Y, X, K, H, Ts,
                               q1, q2, 
                               ell, B, bw, sig.levl = .05, retro = TRUE,
                               nolag.i.x = NA, selfW = NA, scale = FALSE) {
  
  n <- length(Y) - 1
  
  # empty
  ps <- c()
  
  if (retro == TRUE) k <- 1 else 2
  
  for (i in k:(n + 1)) {
    
    # set up
    Ti <- Ts[i]
    Ki <- K[i]
    Tstari <- Tstar[i]
    
    if (i > 2) {
      sym <- ols.est.alphahat(Tstar = Tstar[c(i, 2:(i - 1))], X = X[c(i, 2:(i - 1))], Y = Y[c(i, 2:(i - 1))])
      alpha.wadj <- sym$est[2]
    }
    
    m1.L.i <- matrix(NA, nrow = Ti - Ki - H - max(q1, q2 - 1), ncol = H)
    m2.L.i <- matrix(NA, nrow = Ti - Ki - H - max(q1, q2 - 1), ncol = H)
    
    # compute losses
    for (h in 1:H) {
      for (t in 1:(Ti - Ki - H - max(q1, q2 - 1))) {
        TL <- length(Y[[i]])
        yi <- Y[[i]][-(1:max(q1, q2 - 1))][(t + H - h + 1):(t + Ki + H - h)]
        xi <- X[[i]][-(1:max(q1, q2 - 1)), ][(t + H - h + 1):(t + Ki + H - h),]
        
        # yi lags
        yilags <- c()
        for (j in 1:q1) {
          yijlag <- Y[[i]][-c(1:(q1 - j), (TL - (q1 - 1)):(TL - (q1 - j)))][(t + H - h + 1):(t + Ki + H - h)]
          yilags <- cbind(yilags, yijlag)
        }
        
        # x and x lags
        x.xlags <- xi
        if (q2 > 1) {
          
          for (j in 1:(q2 - 1)) {
            xj.xlags <- X[[i]][-c(1:(q2 - 1 - j), (TL - (q2 - 1 - 1)):(TL - (q2 - 1 - j))), ][(t + H - h + 1):(t + Ki + H - h), ]
            x.xlags <- cbind(x.xlags, xj.xlags)
          }
        }
        
        # functional
        D.i.t <- ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0)
        yilags.D.i.t <- yilags
        x.xlags.D.i.t <- x.xlags
        for (d in 1:nrow(yilags)) {
          yilags.D.i.t[d, ] <- yilags[d, ] * D.i.t[d]
          x.xlags.D.i.t[d, ] <- x.xlags[d, ] * D.i.t[d]
        }
        
        # OLS
        
        # special case that needs covariates without lags
        sc <- c()
        
        if (length(nolag.i.x) == 1) {
          sc <- 0
        } else {
          if (i %in% nolag.i.x[[1]]) {
            sc <- 1
          } else {
            sc <- 0
          }
        }
        
        if (sc == 0) {
          lmodi.adj <- lm(yi ~ 1 + yilags + x.xlags + yilags.D.i.t + x.xlags.D.i.t + D.i.t)
          lmodi.unadj <- lm(yi ~ 1 + yilags + x.xlags)
          
          # beta.hat
          beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
          beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
          beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
          beta.hat.unadj[which(is.na(beta.hat.unadj) == TRUE)] <- 0
          
          
          yhat.adj <- Y[[i]][(t + Ki + H - h + 1 - q1):(t + Ki + H - h)]
          yhat.unadj <- Y[[i]][(t + Ki + H - h + 1 - q1):(t + Ki + H - h)]
          
          for (j in 1:h) {
            
            x.xlags.for.pred <- X[[i]][-(1:q1), ][t + Ki + H - h + j, ]
            if (q2 > 1) {
              
              for (k in 1:(q2 - 1)) {
                x.xlags.for.pred <- cbind(x.xlags.for.pred, 
                                          X[[i]][-1, ][t + Ki + H - h + j - k, ])
              }
            }
            D.i.t <- ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0) 
            x.xlags.for.pred.D.i.t <- x.xlags.for.pred * D.i.t
            yilags.for.pred.D.i.t <- yhat.adj[j:(j + q1 - 1)] * D.i.t
            
            yhat.adj.h <- beta.hat.adj %*% matrix(c(1, yhat.adj[j:(j + q1 - 1)], 
                                                    x.xlags.for.pred, 
                                                    yilags.for.pred.D.i.t,
                                                    x.xlags.for.pred.D.i.t,
                                                    D.i.t))
            yhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, yhat.adj[j:(j + q1 - 1)], 
                                                        x.xlags.for.pred))
            
            if (t + Ki + H - h == Tstari & i >= 3 & h == 1) {
              yhat.adj.h <- yhat.adj.h + alpha.wadj
            }
            
            yhat.adj <- c(yhat.adj, yhat.adj.h)
            yhat.unadj <- c(yhat.unadj, yhat.unadj.h)
          }
        } else {
          
          no.lag.x <- nolag.i.x[[2]][[which(nolag.i.x[[1]] == i)]][-(1:max(q1, q2 - 1))][(t + H - h + 1):(t + Ki + H - h)]
          
          lmodi.adj <- lm(yi ~ 1 + yilags + x.xlags + yilags.D.i.t + x.xlags.D.i.t + D.i.t + no.lag.x)
          lmodi.unadj <- lm(yi ~ 1 + yilags + x.xlags + no.lag.x)
          
          # beta.hat
          beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
          beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
          beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
          beta.hat.unadj[which(is.na(beta.hat.unadj) == TRUE)] <- 0
          
          
          yhat.adj <- Y[[i]][(t + Ki + H - h + 1 - q1):(t + Ki + H - h)]
          yhat.unadj <- Y[[i]][(t + Ki + H - h + 1 - q1):(t + Ki + H - h)]
          
          for (j in 1:h) {
            
            x.xlags.for.pred <- X[[i]][-(1:q1), ][t + Ki + H - h + j, ]
            if (q2 > 1) {
              
              for (k in 1:(q2 - 1)) {
                x.xlags.for.pred <- cbind(x.xlags.for.pred, 
                                          X[[i]][-1, ][t + Ki + H - h + j - k, ])
              }
            }
            D.i.t <- ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0) 
            x.xlags.for.pred.D.i.t <- x.xlags.for.pred * D.i.t
            yilags.for.pred.D.i.t <- yhat.adj[j:(j + q1 - 1)] * D.i.t
            
            yhat.adj.h <- beta.hat.adj %*% matrix(c(1, yhat.adj[j:(j + q1 - 1)], 
                                                    x.xlags.for.pred, 
                                                    yilags.for.pred.D.i.t,
                                                    x.xlags.for.pred.D.i.t,
                                                    D.i.t, nolag.i.x[[2]][[which(nolag.i.x[[1]] == i)]][-(1:max(q1, q2 - 1))][t + Ki + H - h + j - k]))
            yhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, yhat.adj[j:(j + q1 - 1)], 
                                                        x.xlags.for.pred, nolag.i.x[[2]][[which(nolag.i.x[[1]] == i)]][-(1:max(q1, q2 - 1))][t + Ki + H - h + j - k]))
            
            if (t + Ki + H - h == Tstari & i >= 3 & h == 1) {
              yhat.adj.h <- yhat.adj.h + alpha.wadj
            }
            
            yhat.adj <- c(yhat.adj, yhat.adj.h)
            yhat.unadj <- c(yhat.unadj, yhat.unadj.h)
          }
        }
        
        
        
        # losses
        m1.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.adj[h + q1])
        m2.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.unadj[h + q1])
      }
    }
    # loss differential
    d <- m2.L.i - m1.L.i
    ps[i] <- taSPA.mhfc(ell = ell, d = d, B = B, bw = bw)
  }
  
  
  # weights
  if (is.matrix(X[[1]]) == FALSE) {
    for (i in 1:(n + 1)) {
      X[[i]] <- as.matrix(X[[i]])
    }
  }
  # Weights
  W <- round(scm(X = X, Tstar = Tstar, scale = scale)$par, digits = 3)
  
  if (is.na(selfW)[1] == FALSE) {
    W <- selfW
  }
  
  # test
  Is <- ifelse(ps <= .05, yes = 1, no = 0)
  
  # output
  return(list(ps = ps, W = W, Is = Is))
}


## 'PAYEMS' is the Federal Reserve of St. Louis' name for the monthly employment figure,
## (i.e. the number of persons employed in the US, in thousands)
## a time series that goes back to the 1950s.

### In this file, we produce two simple kinds output:

# (1) a set of shock times
# (2) a large matrix of data indexed by month, where the index
# includes the shock times as a proper subset

#significant digits
options(scipen = 7)

# load packages (not all of these may be necessary)
library("data.table")
library("dplyr")
library("tseries")
library("quantmod")
library('Rsolnp')
library('msos')
library('tikzDevice')
library('xtable')
# load packages
require('forecast')

## Time series under study: # persons on nonfarm payrolls
getSymbols("PAYEMS", src = 'FRED')

###          COVARIATES            ###
###           START                ###

## https://fred.stlouisfed.org/series/W825RC1: transfer receipts
getSymbols("W825RC1", src = 'FRED')

## https://fred.stlouisfed.org/series/RPI: real personal income
getSymbols("RPI", src = 'FRED')

## https://fred.stlouisfed.org/series/PCE: Personal Consumption Expenditures
getSymbols("PCE", src = 'FRED')

## https://fred.stlouisfed.org/series/INDPRO: Industrial Production
getSymbols("INDPRO", src = 'FRED') 

## https://fred.stlouisfed.org/series/CPIAUCSL: Consumer Price Index
getSymbols("CPIAUCSL", src = 'FRED')

## https://fred.stlouisfed.org/series/FEDFUNDS: Federal Funds Index
getSymbols("FEDFUNDS", src = 'FRED')

## https://fred.stlouisfed.org/series/LNS12000031: Black Employment Count
getSymbols("LNS12000031", src = 'FRED') ######## NOTE: IF THIS is not super important as variable, let's drop it

###          COVARIATES          ###
###             END              ###



### CONSTRUCTION OF THE DONOR POOL ###
###           START                ###

#https://www.federalreservehistory.org/essays/oil-shock-of-1973-74
#https://www.federalreservehistory.org/essays/recession-of-1981-82
#credit-control program initiated in March 1980 by the Carter administration
#https://www.history.com/news/us-economic-recessions-timeline
#https://fred.stlouisfed.org/series/FEDFUNDS#0

# As a rule, when the shock occurs mid-month, we take that to be the shock-time, even though
# the shock effect is distributed across that month as well as the following month(s)
shock_time_vec <- c('1957-04-01', ## Flu hits US
                    '1958-08-01', ## Fed Funds rate increases from .68 to 1.53 in one month
                    '1973-10-01', ## OAPEC oil embargo begins
                    '1980-03-01', ## program announced on March 14th by Carter
                    '2001-09-01', ## 9/11 attacks occur 1/3 of way into month
                    '2020-03-01') ## COVID shutdown - included in the time series of interest

                    # These are all T* points not T*+1.

### CONSTRUCTION OF THE DONOR POOL ###
###           END                  ###


## Now, we take the outcome variable and covariates and smash them into a 
## large matrix, including lags of the covariates.
datasets <- list(W825RC1, RPI, PCE, INDPRO, CPIAUCSL, FEDFUNDS, LNS12000031)
df <- PAYEMS
for (i in 1:length(datasets)) {df <- merge(df, datasets[[i]])}

# Hit every column with the differenced-log transformation
difflog_df <- data.frame(diff(as.matrix(log(df))))

# We create lags of the covariates
# https://stackoverflow.com/questions/38119225/debugging-function-to-create-multiple-lags-for-multiple-columns-dplyr

difflog_df.lag <- shift(difflog_df, n=1:2, give.names = T)  ##column indexes of columns to be lagged as "[,startcol:endcol]", "n=1:3" specifies the number of lags (lag1, lag2 and lag3 in this case)

# We we combine and original series and the lags
merged <- bind_cols(difflog_df, difflog_df.lag)

#Now add the row names
row.names(merged) <- row.names(difflog_df)

# We do not need any rows prior to 1954-07-01
merged <- merged[row.names(merged) >= '1954-07-01', ]

#Now, due to the release date for each of these monthly series,
#we unfortunately cannot use the same-month data points for some of 
#these columns.  We now drop them...

merged <- subset(merged, select = -c(W825RC1, 
                                    RPI, 
                                    PCE, 
                                    INDPRO, 
                                    CPIAUCSL,
                                    LNS12000031))


#Finally, we have missing data in the late 1950s, so we are faced with a choice:

#(1) We can drop all rows with NA entries, which will remove the 1950s donors.
complete_cases_merged <- merged[complete.cases(merged),]
paste('This is a dataset of dimension', dim(complete_cases_merged)[1], 'by', dim(complete_cases_merged)[2])

#(2) We can drop all columns with NA values, which will drop some covariates, but keep all 5 donors.
no_NA_cols_merged <- merged[ , colSums(is.na(merged)) == 0]
paste('This is a dataset of dimension', dim(no_NA_cols_merged)[1], 'by', dim(no_NA_cols_merged)[2])


# And we're done!  Either of the two datasets above will work.  It's a question of whether we 
# want 5 donors with fewer covariates, or fewer donors with more covariates.


dat <- no_NA_cols_merged
date <- rownames(dat)
# setup
H <- 12
K <- 24
L <- 12

# Time Series 2
shock.t2 <- which(date == '1958-08-01') + 1
shock.adj <- which(date == '1957-04-01') + 1
TS2 <- dat[(shock.t2 - (1 + 12 + K + H)):(shock.t2 + L),]
TS2$shock.t2 <- as.numeric((shock.t2 - (1 + 12 + K + H)):(shock.t2 + L) %in% shock.t2:(shock.t2 + L))
TS2$shock.adj <- as.numeric((shock.t2 - (1 + 12 + K + H)):(shock.t2 + L) == shock.adj)
TS2$shock.adj.2 <- as.numeric((shock.t2 - (1 + 12 + K + H)):(shock.t2 + L) %in% shock.adj:(shock.t2 + L))
## time series 2 modeling

# use AIC to choose between models
mod.t2.null <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1, data = TS2)
mod.t2.permanent.m1 <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1 + 
                                                  shock.adj + shock.t2, data = TS2)
mod.t2.permanent.m2 <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1 + 
                            shock.adj.2 + shock.t2, data = TS2)

AIC(mod.t2.null);AIC(mod.t2.permanent.m1);AIC(mod.t2.permanent.m2)
# choose the second model
mod.t2.permanent <- mod.t2.permanent.m2

## dynamic
co.D <- c()
for (d in 1:nrow(TS2)) {
  co.D <- rbind(co.D, TS2[d, c(2, 4, 3, 5, 6, 7)] * TS2$shock.t2[d])
}
co.D <- as.matrix(co.D)
## modeling
# choose between models
mod.t2.dynamic.null <- lm(TS2$PAYEMS ~ as.matrix(TS2[, c(2, 4, 3, 5, 6, 7)]) + co.D + TS2$shock.t2)
mod.t2.dynamic.m1 <- lm(TS2$PAYEMS ~ as.matrix(TS2[, c(2, 4, 3, 5, 6, 7)]) + co.D + TS2$shock.t2 + TS2$shock.adj)
mod.t2.dynamic.m2 <- lm(TS2$PAYEMS ~ as.matrix(TS2[, c(2, 4, 3, 5, 6, 7)]) + co.D + TS2$shock.t2 + TS2$shock.adj.2)
AIC(mod.t2.dynamic.null); AIC(mod.t2.dynamic.m1); AIC(mod.t2.dynamic.m2)
# choose m2
mod.t2.dynamic <- lm(TS2$PAYEMS ~ as.matrix(TS2[, c(2, 4, 3, 5, 6, 7)]) + co.D + TS2$shock.t2 +  TS2$shock.adj.2)
# AICs
AICs.19580801 <- c(AIC(mod.t2.null), AIC(mod.t2.permanent), AIC(mod.t2.dynamic))

# Time Series 3
shock.t3 <- which(date == '1973-10-01') + 1
TS3 <- dat[(shock.t3 - (1 + 12 + K + H)):(shock.t3 + L),]
TS3$shock.t3 <- as.numeric((shock.t3 - (1 + 12 + K + H)):(shock.t3 + L) %in% shock.t3:(shock.t3 + L))
## time series 2 modeling
mod.t3.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1 + shock.t3, data = TS3)
mod.t3.null <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1, data = TS3)
## dynamic
co.D <- c()
for (d in 1:nrow(TS3)) {
  co.D <- rbind(co.D, TS3[d, c(2, 4, 3, 5, 6, 7)] * TS3$shock.t3[d])
}
co.D <- as.matrix(co.D)
## modeling
mod.t3.dynamic <- lm(TS3$PAYEMS ~ as.matrix(TS3[, c(2, 4, 3, 5, 6, 7)]) + co.D + TS3$shock.t3)
# AICs
AICs.19731001 <- c(AIC(mod.t3.null), AIC(mod.t3.permanent), AIC(mod.t3.dynamic))

# Time Series 4
shock.t4 <- which(date == '1980-03-01') + 1
TS4 <- dat[(shock.t4 - (1 + 12 + K + H)):(shock.t4 + L),]
TS4$shock.t4 <- as.numeric((shock.t4 - (1 + 12 + K + H)):(shock.t4 + L) %in% shock.t4:(shock.t4 + L))
## time series 2 modeling
mod.t4.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1 + shock.t4, data = TS4)
mod.t4.null <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1, data = TS4)
## dynamic
co.D <- c()
for (d in 1:nrow(TS4)) {
  co.D <- rbind(co.D, TS4[d, c(2, 4, 3, 5, 6, 7)] * TS4$shock.t4[d])
}
co.D <- as.matrix(co.D)
## modeling
mod.t4.dynamic <- lm(TS4$PAYEMS ~ as.matrix(TS4[, c(2, 4, 3, 5, 6, 7)]) + co.D + TS4$shock.t4)
# AICs
AICs.19800301 <- c(AIC(mod.t4.null), AIC(mod.t4.permanent), AIC(mod.t4.dynamic))


# Time Series 5
shock.t5 <- which(date == '2001-09-01') + 1
TS5 <- dat[(shock.t5 - (1 + 12 + K + H)):(shock.t5 + L),]
TS5$shock.t5 <- as.numeric((shock.t5 - (1 + 12 + K + H)):(shock.t5 + L) %in% shock.t5:(shock.t5 + L))
## time series 2 modeling
mod.t5.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1 + shock.t5, data = TS5)
mod.t5.null <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1, data = TS5)
## dynamic
co.D <- c()
for (d in 1:nrow(TS5)) {
  co.D <- rbind(co.D, TS5[d, c(2, 4, 3, 5, 6, 7)] * TS5$shock.t5[d])
}
co.D <- as.matrix(co.D)
## modeling
mod.t5.dynamic <- lm(TS5$PAYEMS ~ as.matrix(TS5[, c(2, 4, 3, 5, 6, 7)]) + co.D + TS5$shock.t5)
# AICs
AICs.20010901 <- c(AIC(mod.t5.null), AIC(mod.t5.permanent), AIC(mod.t5.dynamic))

# Time Series 1
shock.t1 <- which(date == '2020-03-01') + 1
TS1 <- dat[(shock.t1 - (1 + 12 + K + H)):(shock.t1 + L),]
TS1$shock.t1 <- as.numeric((shock.t1 - (1 + 12 + K + H)):(shock.t1 + L) %in% shock.t1:(shock.t1 + L))
## time series 2 modeling
mod.t1.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1 + shock.t1, data = TS1)
mod.t1.null <- lm(PAYEMS ~ PAYEMS_lag_1 + INDPRO_lag_1 + CPIAUCSL_lag_1, data = TS1)
## dynamic
co.D <- c()
for (d in 1:nrow(TS1)) {
  co.D <- rbind(co.D, TS1[d, c(2, 4, 3, 5, 6, 7)] * TS1$shock.t1[d])
}
co.D <- as.matrix(co.D)
## modeling
mod.t1.dynamic <- lm(TS1$PAYEMS ~ as.matrix(TS1[, c(2, 4, 3, 5, 6, 7)]) + co.D + TS1$shock.t1)
# AICs
AICs.20200301 <- c(AIC(mod.t1.null), AIC(mod.t1.permanent), AIC(mod.t1.dynamic))

# Ti
Ts <- rep(62, 5)
# Tstar
Tstar <- rep(49, 5)
# Y
Y <- c()
X <- c()
TS <- list(TS1, TS2, TS3, TS4, TS5)
for (i in 1:5) {
  Y[[i]] <- TS[[i]]$PAYEMS
  X[[i]] <- as.matrix(TS[[i]][, c(4, 6)])
}

# additional covariates that are needed to be adjusted in time series 2
sc.x <- nolag.i.x <- list(2, list(TS2$shock.adj.2))

# testing 
res1 <- ps.indic.W.permanent(Tstar = Tstar, Y = Y, X = X, K = rep(K, 5), H = H,
                             Ts = Ts, ell = 4, B = 200, bw = 4, 
                             sc.x = sc.x,
                             scale = TRUE)

res2 <- ps.indic.W.dynamic(Tstar = Tstar
                           , Y = Y, X = X, K = rep(K, 5), 
                           q1 = 2, q2 = 2, nolag.i.x = nolag.i.x,
                           H = H, Ts = Ts, ell = 4, B = 200, bw = 4,
                           scale = TRUE)

Wstar <- res1$W


AICs <- rbind(AICs.19580801, AICs.19731001, AICs.19800301, AICs.20010901)


# plot shock transience
setwd('~/Desktop/Research/synthetic prediction/')

# data
ts <- dat[(shock.t1 - (1 + 15 + K + H)):nrow(dat),]

# load package
require('tikzDevice')
# tex package specification
options(tikzLatexPackages 
        = c(getOption( "tikzLatexPackages" ),
            "\\usepackage{amsmath,amsfonts,amsthm, palatino, mathpazo}"))
# file setting
tikz(file = 'UEMtransience.tex', width = 5.5, height = 4, standAlone = TRUE)
# plot
plot(x = as.Date(rownames(ts)), y = ts$PAYEMS, type = 'l',
     xlab = 'Time $t$', ylab = 'Persons on nonfarm payrolls $y_{1,t}$',
     main = 'Shock transience of persons on nonfarm payrolls $y_{1,t}$',
     col = 'deepskyblue')
# segments
segments(x0 = as.Date("2020-04-01"), y0 = -.12, y1 = .02, lty = 3, col = 'magenta')
segments(x0 = as.Date("2020-03-01"), y0 = -.12, y1 = .02, lty = 3, col = 'magenta')
segments(x0 = as.Date("2020-03-01"), x1 = as.Date("2021-12-01"),
         y0 = ts$PAYEMS[which(rownames(ts) == as.Date("2020-03-01"))],
         col = 'indianred1', lty = 4)
# output
dev.off()

