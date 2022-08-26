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

# tuSPA test for comparison of two models
tuSPA.mhfc <- function(ell, d, B, bw = 4) {
  # ell is the block length
  # d is the matrix of loss differential
  # B is the bootstrap size
  # bw is the bandwidth parameter for computing HAC estimator
  
  # compute dbar
  d.bar <- matrix(apply(d, 2, mean))
  
  # moving-block bootstrap
  Tt <- nrow(d); H <- ncol(d)
  # T = ell * K
  K <- Tt / ell
  # Bootstrap replication
  t.uSPA.B <- c()
  for (b in 1:B) { 
    # uniform draw from 1,...,T - ell +1
    Iks <- sample(1:(Tt - ell + 1), K, replace = TRUE)
    tau <- matrix(sapply(Iks, function(x) x:(x + ell - 1)))
    tau <- as.numeric(tau)
    t.uSPA.b.H <- c()
    for (h in 1:H) { 
      # moving-block bootstrap
      d.b.h <- d[tau, h]
      d.b.bar.h <- mean(d.b.h)
      # compute moving block bootstrap standard error
      xi.b.h <- 0
      for (k in 1:K) {
        xi.b.h <- xi.b.h + (sum(d.b.h[(k - 1) * ell + 1:ell] - d.b.bar.h)) ^ 2 / ell
      }
      xi.b.h <- xi.b.h / K
      # compute t.uSPA statistic
      t.uSPA.b.H[h] <- sqrt(Tt) * (d.b.bar.h - d.bar[h]) / sqrt(xi.b.h)
    }
    t.uSPA.B[b] <- min(t.uSPA.b.H)
  }
  # compute t.aSPA statistic
  ## compute Heteroskedasticity and autocorrelation Consistent (HAC) estimator
  require('tsapp')
  Omega <- HAC(d, method = 'Quadratic Spectral', bw = bw)
  t.uSPA <- min(as.numeric(sqrt(Tt) * d.bar / sqrt(diag(Omega))))
  p <- mean(t.uSPA < t.uSPA.B)
  # return output
  return(p)
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
      xi.b <- xi.b + (sum(d.b.w[(k - 1) * ell + 1:ell] - mean(d.b.w))) ^ 2 / ell
    }
    xi.b <- xi.b / K
    # compute t.aSPA statistic
    t.aSPA.B[b] <- sqrt(Tt) * (d.b.bar - d.bar) / sqrt(xi.b)
  }
  # compute t.aSPA statistic
  ## compute Heteroskedasticity and autocorrelation Consistent (HAC) estimator
  require('tsapp')
  Omega <- HAC(d, method = 'Quadratic Spectral', bw = bw)
  w <- matrix(1 / H, nrow = H)
  xi.hat <- sqrt(t(w) %*% Omega %*% w)
  t.aSPA <- as.numeric(sqrt(Tt) * d.bar / sqrt(xi.hat))
  p <- mean(t.aSPA < t.aSPA.B)
  # return output
  return(p)
}

# voting
vote <- function(ps, sig = 0.05, weight = FALSE, W) {
  indic.p <- ifelse(ps < 0.05, yes = 1, no = 0)
  if (weight == FALSE) {
    mean <- mean(indic.p)
  } else {
    mean <- crossprod(indic.p, W)
  }
  if (mean > 0.5) output <- 1
  else if (mean == 0.5) output <- sample(c(1, 0), size = 1)
  else output <- 0
  return(output)
}

# functions that return alpha.hat and synthetic weights
# functions that return alpha.hat and synthetic weights
ps.indic.W.permanent <- function(Tstar, Y, X, K, H, Ts, ell, B, bw, 
                                 sig.levl = .05, q1, q2, subset = FALSE, 
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
    
    m1.L.i <- matrix(NA, nrow = Ti - Ki - H - q1, ncol = H)
    m2.L.i <- matrix(NA, nrow = Ti - Ki - H - q1, ncol = H)
    
    # compute losses
    for (h in 1:H) {
      for (t in 1:(Ti - Ki - H - q1)) {
        TL <- length(Y[[i]])
        yi <- Y[[i]][-(1:q1)][(t + H - h + 1):(t + Ki + H - h)]
        xi <- X[[i]][-(1:q1), ][(t + H - h + 1):(t + Ki + H - h),]
        
        # yi lags
        yilags <- c()
        for (j in 1:q1) {
          yijlag <- dplyr::lag(Y[[i]][-(1:q1)], n = j)[(t + H - h + 1):(t + Ki + H - h)]
          yilags <- cbind(yilags, yijlag)
        }
        
        # x and x lags
        x.xlags <- xi
        if (q2 > 1) {
          
          for (j in 1:(q2 - 1)) {
            xj.xlags <- dplyr::lag(X[[i]][-(1:max(q1, q2 - 1)), ], n = j)[(t + H - h + 1):(t + Ki + H - h), ]
            x.xlags <- cbind(x.xlags, xj.xlags)
          }
        }
        
        # OLS
        lmodi.adj <- lm(yi ~ 1 + ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0) + 
                          yilags + x.xlags)
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
                                        X[[i]][-(1:q1), ][t + Ki + H - h + j - k, ])
            }
          }
          yhat.adj.h <- beta.hat.adj %*% matrix(c(1, ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0),
                                                  yhat.adj[j:(j + q1 - 1)], x.xlags.for.pred))
          yhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, yhat.adj[j:(j + q1 - 1)], x.xlags.for.pred))
          
          yhat.adj <- c(yhat.adj, yhat.adj.h)
          yhat.unadj <- c(yhat.unadj, yhat.unadj.h)
        }
        
        # losses
        m1.L.i[t, h] <- sel(y = Y[[i]][-(1:q1)][t + Ki + H], yhat = yhat.adj[h + q1])
        m2.L.i[t, h] <- sel(y = Y[[i]][-(1:q1)][t + Ki + H], yhat = yhat.unadj[h + q1])
      }
    }
    # loss differential
    d <- m2.L.i - m1.L.i
    if (subset == FALSE) {
      ps[i] <- taSPA.mhfc(ell = ell, d = d, B = B, bw = bw) 
    } else {
      sub.index <- (Ti - Ki - max(q1, q2 - 1) - (Ti - Tstari) + ceiling(1.5 * sqrt(Ti - Tstari))):(Ti - Ki - H - max(q1, q2 - 1))
      d.subset <- d[sub.index,  ]
      ps[i] <- taSPA.mhfc(ell = ell, d = d.subset, B = B, bw = bw)
    }
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
  return(list(ps = ps, W = W, Is = Is, dim = dim(d.subset)))
}

# functions that return alpha.hat and synthetic weights for decaying shock effects
ps.indic.W.dynamic <- function(Tstar, Y, X, K, H, Ts,
                               q1, q2, subset = TRUE, 
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
          yijlag <- dplyr::lag(Y[[i]][-(1:max(q1, q2 - 1))], n = j)[(t + H - h + 1):(t + Ki + H - h)]
          yilags <- cbind(yilags, yijlag)
        }
        
        # x and x lags
        x.xlags <- xi
        if (q2 > 1) {
          
          for (j in 1:(q2 - 1)) {
            xj.xlags <- dplyr::lag(X[[i]][-(1:max(q1, q2 - 1)), ], n = j)[(t + H - h + 1):(t + Ki + H - h), ]
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
                                          X[[i]][-(1:q1), ][t + Ki + H - h + j - k, ])
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
            yhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, yhat.unadj[j:(j + q1 - 1)], 
                                                        x.xlags.for.pred))
            
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
                                          X[[i]][-(1:q1), ][t + Ki + H - h + j - k, ])
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
            yhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, yhat.unadj[j:(j + q1 - 1)], 
                                                        x.xlags.for.pred, nolag.i.x[[2]][[which(nolag.i.x[[1]] == i)]][-(1:max(q1, q2 - 1))][t + Ki + H - h + j - k]))
            
            yhat.adj <- c(yhat.adj, yhat.adj.h)
            yhat.unadj <- c(yhat.unadj, yhat.unadj.h)
          }
        }
        
        
        
        # losses
        m1.L.i[t, h] <- sel(y = Y[[i]][-(1:max(q1, q2 - 1))][t + Ki + H], yhat = yhat.adj[h + q1])
        m2.L.i[t, h] <- sel(y = Y[[i]][-(1:max(q1, q2 - 1))][t + Ki + H], yhat = yhat.unadj[h + q1])
      }
    }
    # loss differential
    d <- m2.L.i - m1.L.i
    if (subset == FALSE) {
      ps[i] <- taSPA.mhfc(ell = ell, d = d, B = B, bw = bw) 
    } else {
      sub.index <- (Ti - Ki - max(q1, q2 - 1) - (Ti - Tstari) + ceiling(1.5 * sqrt(Ti - Tstari))):(Ti - Ki - H - max(q1, q2 - 1))
      d.subset <- d[sub.index,  ]
      ps[i] <- taSPA.mhfc(ell = ell, d = d.subset, B = B, bw = bw)
    }
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
  return(list(ps = ps, W = W, Is = Is, dim = dim(d.subset)))
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

difflog_df.lag <- shift(difflog_df, n=1:1, give.names = T)  ##column indexes of columns to be lagged as "[,startcol:endcol]", "n=1:3" specifies the number of lags (lag1, lag2 and lag3 in this case)

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


dat <- complete_cases_merged
dat <- dat[, c(1, 3, 2, 4, 6, 8, 10)]
date <- rownames(dat)
# setup
K <- 24
# number of horizons
Hs <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
# tab
tab <- data.frame()
L <- 15

# result
for (H in Hs) {

# Time Series 2
shock.t2 <- which(date == '1980-03-01') + 1
TS2 <- dat[(shock.t2 - (1 + 12 + K + H)):(shock.t2 + L),]
TS2$shock.t2 <- as.numeric((shock.t2 - (1 + 12 + K + H)):(shock.t2 + L) %in% shock.t2:(shock.t2 + L))
## time series 2 modeling
mod.t2.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                         PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1 + shock.t2, data = TS2)
mod.t2.null <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                    PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1, data = TS2)
## dynamic
co.D <- c()
for (d in 1:nrow(TS2)) {
  co.D <- rbind(co.D, TS2[d, 2:7] * TS2$shock.t2[d])
}
co.D <- as.matrix(co.D)
## modeling
mod.t2.dynamic <- lm(TS2$PAYEMS ~ as.matrix(TS2[, 2:7]) + co.D + TS2$shock.t2)
# AICs
AICs.19800301 <- c(AIC(mod.t2.null), AIC(mod.t2.permanent), AIC(mod.t2.dynamic))
I.AIC.2 <- ifelse(AIC(mod.t2.permanent) > AIC(mod.t2.dynamic), yes = 1, no = 0)

# Time Series 3
shock.t3 <- which(date == '2001-09-01') + 1
TS3 <- dat[(shock.t3 - (1 + 12 + K + H)):(shock.t3 + L),]
TS3$shock.t3 <- as.numeric((shock.t3 - (1 + 12 + K + H)):(shock.t3 + L) %in% shock.t3:(shock.t3 + L))
## time series 2 modeling
mod.t3.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                         PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1 + shock.t3, data = TS3)
mod.t3.null <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                    PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1, data = TS3)
## dynamic
co.D <- c()
for (d in 1:nrow(TS3)) {
  co.D <- rbind(co.D, TS3[d, 2:7] * TS3$shock.t3[d])
}
co.D <- as.matrix(co.D)
## modeling
mod.t3.dynamic <- lm(TS3$PAYEMS ~ as.matrix(TS3[, 2:7]) + co.D + TS3$shock.t3)
# AICs
AICs.20010901 <- c(AIC(mod.t3.null), AIC(mod.t3.permanent), AIC(mod.t3.dynamic))
I.AIC.3 <- ifelse(AIC(mod.t3.permanent) > AIC(mod.t3.dynamic), yes = 1, no = 0)

# Time Series 4
shock.t4 <- which(date == '2008-03-01') + 1
TS4 <- dat[(shock.t4 - (1 + 12 + K + H)):(shock.t4 + L),]
TS4$shock.t4 <- as.numeric((shock.t4 - (1 + 12 + K + H)):(shock.t4 + L) %in% shock.t4:(shock.t4 + L))
## time series 2 modeling
mod.t4.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                         PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1 + shock.t4, data = TS4)
mod.t4.null <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                    PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1, data = TS4)
## dynamic
co.D <- c()
for (d in 1:nrow(TS4)) {
  co.D <- rbind(co.D, TS4[d, 2:7] * TS4$shock.t4[d])
}
co.D <- as.matrix(co.D)
## modeling
mod.t4.dynamic <- lm(TS4$PAYEMS ~ as.matrix(TS4[, 2:7]) + co.D + TS4$shock.t4)
# AICs
AICs.20080301 <- c(AIC(mod.t4.null), AIC(mod.t4.permanent), AIC(mod.t4.dynamic))
I.AIC.4 <- ifelse(AIC(mod.t4.permanent) > AIC(mod.t4.dynamic), yes = 1, no = 0)


# Time Series 1
shock.t1 <- which(date == '2020-03-01') + 1
TS1 <- dat[(shock.t1 - (1 + 12 + K + H)):(shock.t1 + L),]
TS1$shock.t1 <- as.numeric((shock.t1 - (1 + 12 + K + H)):(shock.t1 + L) %in% shock.t1:(shock.t1 + L))
## time series 2 modeling
mod.t1.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                         PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1 + shock.t1, data = TS1)
mod.t1.null <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                    PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1, data = TS1)
## dynamic
co.D <- c()
for (d in 1:nrow(TS1)) {
  co.D <- rbind(co.D, TS1[d, 2:7] * TS1$shock.t1[d])
}
co.D <- as.matrix(co.D)
## modeling
mod.t1.dynamic <- lm(TS1$PAYEMS ~ as.matrix(TS1[, 2:7]) + co.D + TS1$shock.t1)
# AICs
AICs.20200301 <- c(AIC(mod.t1.null), AIC(mod.t1.permanent), AIC(mod.t1.dynamic))

# Ti
Ts <- rep(nrow(TS1), 5)
# Tstar
Tstar <- rep(39, 5)
# Y
Y <- c()
X <- c()
TS <- list(TS1, TS2, TS3, TS4)
for (i in 1:4) {
  Y[[i]] <- TS[[i]]$PAYEMS
  X[[i]] <- as.matrix(TS[[i]][, c(3, 4, 5, 6, 7)])
}


# testing 
res1 <- ps.indic.W.permanent(Tstar = Tstar, q1 = 1, q2 = 1, subset = TRUE, 
                             Y = Y, X = X, K = rep(K, 4), H = H,
                             Ts = Ts, ell = 4, B = 500, bw = 4, 
                             scale = TRUE)

res2 <- ps.indic.W.dynamic(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), 
                           q1 = 1, q2 = 1, subset = TRUE,  
                           H = H, Ts = Ts, ell = 4, B = 500, bw = 4,
                           scale = TRUE)

# p-values
res1.ps <- res1$ps
res2.ps <- res2$ps
# phats
res1.phat <- sum(res1$ps[-1] * res1$W)
res2.phat <- sum(res2$ps[-1] * res2$W)
# p-value distance
res1.pd <- abs(res1$ps[1] - res1.phat)
res2.pd <- abs(res2$ps[1] - res2.phat)
# weighted indicators
res1.wi <- sum(res1$Is[-1] * res1$W)
res2.wi <- sum(res2$Is[-1] * res2$W)
# voting
res1.vote <- ifelse(mean(res1$Is[-1]) >= 0.5, yes = 1, no = 0)
res2.vote <- ifelse(mean(res2$Is[-1]) >= 0.5, yes = 1, no = 0)
# weighted voting
res1.wvote <- ifelse(res1.wi >= 0.5, yes = 1, no = 0)
res2.wvote <- ifelse(res2.wi >= 0.5, yes = 1, no = 0)
# correctness
res1.vote.c <-
  ifelse(res1.vote == res1$Is[1], yes = 'Yes', no = 'No')
res2.vote.c <-
  ifelse(res2.vote == res2$Is[1], yes = 'Yes', no = 'No')
res1.wvote.c <-
  ifelse(res1.wvote == res1$Is[1], yes = 'Yes', no = 'No')
res2.wvote.c <-
  ifelse(res2.wvote == res2$Is[1], yes = 'Yes', no = 'No')
# correctness output
res1.vote.output <- paste0(res1.vote, ' (', res1.vote.c, ')')
res2.vote.output <- paste0(res2.vote, ' (', res2.vote.c, ')')
res1.wvote.output <- paste0(res1.wvote, ' (', res1.wvote.c, ')')
res2.wvote.output <- paste0(res2.wvote, ' (', res2.wvote.c, ')')
# table
tab.i1 <-
  data.frame(
    H = H,
    model = 'M0',
    p1 = res1$ps[1],
    phat = res1.phat,
    pd = res1.pd,
    vote.output = res1.vote.output,
    wvote.output = res1.wvote.output
  )
tab.i2 <-
  data.frame(
    H = H,
    model = 'M1',
    p1 = res2$ps[1],
    phat = res2.phat,
    pd = res2.pd,
    vote.output = res2.vote.output,
    wvote.output = res2.wvote.output
  )
tab.i <- rbind(tab.i1, tab.i2)
tab <- rbind(tab, tab.i)
}

res <- c()
for (t in H:(nrow(dat) - shock.t1)) {
   TS1 <- dat[(shock.t1 - (1 + 12 + K + t - H)):(shock.t1 + t - H),]
   TS1$shock.t1 <- as.numeric((shock.t1 - (1 + 12 + K + t - H)):(shock.t1 + t - H) %in% shock.t1:(shock.t1 + t - H))
   mod.t1.permanent <- lm(PAYEMS ~ PAYEMS_lag_1 + FEDFUNDS + W825RC1_lag_1 + 
                               +                            PCE_lag_1 + CPIAUCSL_lag_1 + LNS12000031_lag_1 + shock.t1, data = TS1)
   res <- rbind(res, coef(summary(mod.t1.permanent))[8, c(1, 2, 4)])
}


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
     xlab = 'Time $t$', ylab = 'Change in Log Nonfarm Payrolls $y_{1,t}$',
     main = 'Shock transience of $y_{1,t}$',
     col = 'deepskyblue')
# segments
segments(x0 = as.Date("2020-04-01"), y0 = -.12, y1 = .02, lty = 3, col = 'magenta')
segments(x0 = as.Date("2020-03-01"), y0 = -.12, y1 = .02, lty = 3, col = 'magenta')
segments(x0 = as.Date("2020-03-01"), x1 = as.Date("2021-12-01"),
         y0 = ts$PAYEMS[which(rownames(ts) == as.Date("2020-03-01"))],
         col = 'indianred1', lty = 4)
# output
dev.off()

