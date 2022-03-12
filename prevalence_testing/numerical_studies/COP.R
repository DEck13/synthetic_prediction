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
ps.indic.W.permanent <- function(Tstar, Y, X, K, H, Ts, ell, B, bw, 
                                 sig.levl = .05, q1, 
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
          yijlag <- lag(Y[[i]][-(1:q1)], n = j)[(t + H - h + 1):(t + Ki + H - h)]
          yilags <- cbind(yilags, yijlag)
        }
        
        # OLS
        lmodi.adj <- lm(yi ~ 1 + ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0) + 
                          yilags + xi)
        lmodi.unadj <- lm(yi ~ 1 + yilags + xi)
        
        # beta.hat
        beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
        beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
        beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
        
        yhat.adj <- Y[[i]][(t + Ki + H - h + 1 - q1):(t + Ki + H - h)]
        yhat.unadj <- Y[[i]][(t + Ki + H - h + 1 - q1):(t + Ki + H - h)]
        
        for (j in 1:h) {
          yhat.adj.h <- matrix(c(1, ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0),
                                 yhat.adj[j:(j + q1 - 1)], X[[i]][-(1:q1), ][t + Ki + H - h + j, ]))
          yhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, yhat.adj[j:(j + q1 - 1)], X[[i]][-(1:q1), ][t + Ki + H - h + j, ]))
          
          if (t + Ki + H - h == Tstari & i >= 3 & h == 1) {
            yhat.adj.h <- yhat.adj.h + alpha.wadj
          }
          
          yhat.adj <- c(yhat.adj, yhat.adj.h)
          yhat.unadj <- c(yhat.unadj, yhat.unadj.h)
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
          yijlag <- lag(Y[[i]][-(1:max(q1, q2 - 1))], n = j)[(t + H - h + 1):(t + Ki + H - h)]
          yilags <- cbind(yilags, yijlag)
        }
        
        # x and x lags
        x.xlags <- xi
        if (q2 > 1) {
          
          for (j in 1:(q2 - 1)) {
            xj.xlags <- lag(X[[i]][-(1:max(q1, q2 - 1)), ], n = j)[(t + H - h + 1):(t + Ki + H - h), ]
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

# set working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Research/Post-Shock Prediction/")


## Conoco Phillips
getSymbols('COP', from = "2000-01-01")
COP <- as.data.frame(COP)
COP <- COP %>% mutate(Date = rownames(COP))

## S&P 500
getSymbols('^GSPC', from = "1970-01-01")
GSPC <- as.data.frame(GSPC)
GSPC <- GSPC %>% mutate(Date = rownames(GSPC))

## Brent Crude prices
Brent_Crude <- read.csv("https://pkgstore.datahub.io/core/oil-prices/brent-daily_csv/data/d93216330ab2c27aae3d177b2f0f0921/brent-daily_csv.csv") %>%
  rename(Oil_Close = Price)

## WTI Crude prices
WTI_Crude <- read.csv("https://pkgstore.datahub.io/core/oil-prices/wti-daily_csv/data/c414c9d375ec3c8f9d7c276d866fb2a4/wti-daily_csv.csv") %>%
  rename(WTI_Close = Price)

## Gold Price
getSymbols('GC=F', from = "2000-01-01")
gold <- as.data.frame(`GC=F`)
gold <- na.omit(gold)
gold <- gold %>% mutate(Date = rownames(gold))

## Dollar Index
getSymbols('DX-Y.NYB', from = "2000-01-01")
USD <- as.data.frame(`DX-Y.NYB`)
USD <- na.omit(USD)
USD <- USD %>% mutate(Date = rownames(USD))

## 13-Week T-Bill
# ~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Research/Post-Shock Prediction/
TB <- read.csv('^IRX.csv', na.strings = 'null')
TB <- na.omit(TB)

## Volatility Index
getSymbols('^VIX', from = "2000-01-01")
VIX <- as.data.frame(VIX)
VIX <- na.omit(VIX)
VIX <- VIX %>% mutate(Date = rownames(VIX))


## inflation adjustment
getSymbols("CPIAUCSL", src = 'FRED')
avg.cpi <- apply.yearly(CPIAUCSL, mean)
inflation_adj <- as.numeric(avg.cpi['2020'])/avg.cpi
inflation_adj <- as.data.frame(inflation_adj)
colnames(inflation_adj) <- c("dollars_2020")
inflation_adj <- inflation_adj %>% mutate(year = 1947:2022)


# Data Preparation
COP_close <- COP %>% dplyr::select(COP.Close, Date) %>% rename(COP_Close = COP.Close)
GSPC_close <- GSPC %>% dplyr::select(GSPC.Close, Date) %>% rename(GSPC_Close = GSPC.Close)
USD_close <- USD %>% dplyr::select(`DX-Y.NYB.Close`, Date) %>% rename(USD_Close = `DX-Y.NYB.Close`)
TB_close <- TB %>% dplyr::select(Close, Date) %>% rename(TB_Close = Close)
VIX_close <- VIX %>% dplyr::select(VIX.Close, Date) %>% rename(VIX_Close = VIX.Close)

tom <- list(GSPC_close, WTI_Crude, USD_close, TB_close, VIX_close)
for (i in 1:length(tom)) {
  COP_close <- merge(COP_close, tom[[i]])
}


# response
Y <- COP_close$COP_Close[-1]

# data frame
COP_close <- data.frame(COP_close[-nrow(COP_close), ], Y)

COP_close <- data.frame(COP_close[-1, ], COP_close[-nrow(COP_close), -c(1,8)])

colnames(COP_close)[9:14] <- c('ylag2', paste0(colnames(COP_close)[3:7], '.lag'))


# number of horizons
H <- 10
# training sample size
K <- 30
# number of days after shock date
L <- 30

#### Monday, March 17th, 2008

## March 17th; 1 day nowcast

# shock effect date
start <- which(COP_close$Date == "2008-03-14")
start_day_20080317 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
COP_close <- COP_close %>% mutate(start_day_20080317 = start_day_20080317)
TS2 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
# inflation adjustment
TS2[, 2:8] <- TS2[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008] 


m_COP_3_17.null <- lm(Y ~ COP_Close + ylag2 + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close + 
                        GSPC_Close.lag + WTI_Close.lag + USD_Close.lag + TB_Close.lag + VIX_Close.lag, 
                      data = TS2)
m_COP_3_17 <- lm(Y ~ COP_Close + ylag2 + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close + 
                   GSPC_Close.lag + WTI_Close.lag + USD_Close.lag + TB_Close.lag + VIX_Close.lag + 
                   start_day_20080317, 
                 data = TS2)

# dynamic
co.D <- c()
for (d in 1:nrow(TS2)) {
  co.D <- rbind(co.D, TS2[d, c(2:7, 9:14)] * TS2$start_day_20080317[d])
}
co.D <- as.matrix(co.D)

m_COP.dynamic <- lm(TS2$Y ~ as.matrix(TS2[, c(2:7, 9:14)]) + co.D + TS2$start_day_20080317)

# AICs
AIC.20080314.null <- AIC(m_COP_3_17.null)
AIC.20080314.permanent <- AIC(m_COP_3_17)
AIC.20080314.dynamic <- AIC(m_COP.dynamic)

AIC.20080314 <- c(AIC.20080314.null, AIC.20080314.permanent, 
                  AIC.20080314.dynamic)
I.AIC.2 <- ifelse(AIC(m_COP_3_17) > AIC(m_COP.dynamic), yes = 1, no = 0)

#### 2008 shock effects

# shock effect date
start <- which(COP_close$Date == "2008-09-08")
start_day_20080908 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
COP_close <- COP_close %>% mutate(start_day_20080908 = start_day_20080908)
TS3 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
# adjust for inflation
TS3[, 2:8] <- TS3[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008] 

start2 <- which(TS3$Date == "2008-09-08")

# models
m_COP_null <- lm(Y ~ COP_Close + ylag2 + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close + 
                   GSPC_Close.lag + WTI_Close.lag + USD_Close.lag + TB_Close.lag + VIX_Close.lag, 
                 data = TS3)
m_COP_Sept_08 <- lm(Y ~ COP_Close + ylag2 + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close + 
                      GSPC_Close.lag + WTI_Close.lag + USD_Close.lag + TB_Close.lag + VIX_Close.lag +
                      start_day_20080908, data = TS3)
# dynamic
co.D <- c()
for (d in 1:nrow(TS3)) {
  co.D <- rbind(co.D, TS3[d, c(2:7, 9:14)] * TS3$start_day_20080908[d])
}
co.D <- as.matrix(co.D)

m_COP.dynamic <- lm(TS3$Y ~ as.matrix(TS3[, c(2:7, 9:14)]) + co.D + TS3$start_day_20080908)

# AICs
AIC.20080908.null <- AIC(m_COP_null) 
AIC.20080908.permanent <- AIC(m_COP_Sept_08) 
AIC.20080908.dynamic <- AIC(m_COP.dynamic) 

AIC.20080908 <- c(AIC.20080908.null, AIC.20080908.permanent,
                  AIC.20080908.dynamic)
I.AIC.3 <- ifelse(AIC(m_COP_Sept_08) > AIC(m_COP.dynamic), yes = 1, no = 0)

#### Thursday, November 27, 2014

# shock effect date
start <- which(COP_close$Date == "2014-11-26")
start_day_20141127 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
COP_close <- COP_close %>% mutate(start_day_20141127 = start_day_20141127)
# time window
TS4 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
# adjust for inflation
TS4[, 2:8] <- TS4[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2014] 

# models
m_COP_11_27_14.null <- lm(Y ~ COP_Close + ylag2 + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close + 
                            GSPC_Close.lag + WTI_Close.lag + USD_Close.lag + TB_Close.lag + VIX_Close.lag, data = TS4)
m_COP_11_27_14 <- lm(Y ~ COP_Close + ylag2 + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close + 
                       GSPC_Close.lag + WTI_Close.lag + USD_Close.lag + TB_Close.lag + VIX_Close.lag + start_day_20141127, 
                     data = TS4)
# dynamic
co.D <- c()
for (d in 1:nrow(TS4)) {
  co.D <- rbind(co.D, TS4[d, c(2:7, 9:14)] * TS4$start_day_20141127[d])
}
co.D <- as.matrix(co.D)

m_COP.dynamic <- lm(TS4$Y ~ as.matrix(TS4[, c(2:7, 9:14)]) + co.D + TS4$start_day_20141127)

# AICs of different models
AIC.20141127.null <- AIC(m_COP_11_27_14.null) 
AIC.20141127.permanent <- AIC(m_COP_11_27_14) 
AIC.20141127.dynamic <- AIC(m_COP.dynamic) 

AIC.20141127 <- c(AIC.20141127.null, AIC.20141127.permanent, 
                  AIC.20141127.dynamic)

I.AIC.4 <- ifelse(AIC(m_COP_11_27_14) > AIC(m_COP.dynamic), yes = 1, no = 0)
#### The March 9th, 2020 shock effect:

# shock effect date
start <- which(COP_close$Date == "2020-03-06")
start_day_20200309 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
COP_close <- COP_close %>% mutate(start_day_20200309 = start_day_20200309)
# time window
TS1 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]

m_COP_03_09_20.null <- lm(Y ~ COP_Close + ylag2 + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close + 
                            GSPC_Close.lag + WTI_Close.lag + USD_Close.lag + TB_Close.lag + VIX_Close.lag,  
                     data = TS1)
m_COP_03_09_20 <- lm(Y ~ COP_Close + ylag2 + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close + 
                       GSPC_Close.lag + WTI_Close.lag + USD_Close.lag + TB_Close.lag + VIX_Close.lag + start_day_20200309,  
                     data = TS1)
# dynamic
co.D <- c()
for (d in 1:nrow(TS4)) {
  co.D <- rbind(co.D, TS1[d, c(2:7, 9:14)] * TS1$start_day_20200309[d])
}
co.D <- as.matrix(co.D)

m_COP.dynamic <- lm(TS1$Y ~ as.matrix(TS1[, c(2:7, 9:14)]) + co.D + TS1$start_day_20200309)

AIC.20200309.null <- AIC(m_COP_03_09_20.null)
AIC.20200309.permanent <- AIC(m_COP_03_09_20)
AIC.20200309.dynamic <- AIC(m_COP.dynamic)

AIC.20200309 <- c(AIC.20200309.null, AIC.20200309.permanent, 
                  AIC.20200309.dynamic)


# Ti
Ts <- c(nrow(TS1), nrow(TS2), nrow(TS3), nrow(TS4))
# Tstar
Tstar <- c(which(TS1$Date == "2020-03-05"),
           which(TS2$Date == "2008-03-13"),
           which(TS3$Date == "2008-09-05"), 
           which(TS4$Date == "2014-11-25"))
# Y
Y <- c()
X <- c()
TS <- list(TS1, TS2, TS3, TS4)
for (i in 1:4) {
  Y[[i]] <- TS[[i]]$Y
  X[[i]] <- as.matrix(TS[[i]][, c(3:7, 10:14)])
}

res1 <- ps.indic.W.permanent(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), H = H,
                             q1 = 2, 
                             Ts = Ts, ell = 4, B = 200, bw = 4, scale = TRUE)
res2 <- ps.indic.W.dynamic(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), 
                              q1 = 2, q2 = 2, 
                              H = H, Ts = Ts, ell = 4, B = 200, bw = 4,
                           scale = TRUE)

res1$ps
sum(res1$ps[-1] * res1$W)
sum(res1$Is[-1] * res1$W)

res2$ps
sum(res2$ps[-1] * res2$W)
sum(res2$Is[-1] * res2$W)

AICs <- rbind(AIC.20080314, AIC.20080908,  AIC.20141127)

rownames(AICs) <- c('TS2', 'TS3', 'TS4')
colnames(AICs) <- c('null', 'Permanent',  'dynamic')

apply(AICs, 1, which.min)
AICs <- rbind(AICs, apply(AICs, 2, mean))
rownames(AICs)[4] <- 'Mean'


wAICs <- c()
for (j in 1:3) {
  wAICs <- c(wAICs, sum(AICs[1:3, j] *  res1$W))
}

AICs <- rbind(AICs, wAICs, AIC.20200309)
rownames(AICs)[c(5, 6)] <- c('wAIC', 'TS1' )
xtable(AICs)

# weighted comparison
I.AICs <- c(I.AIC.2, I.AIC.3, I.AIC.4)
matrix(res1$W, nrow = 1) %*% I.AICs

# plot shock transience
setwd('~/Desktop/Research/synthetic prediction/')

# data
ts <- COP_close[(start - (1 + 15 + K + H)):nrow(COP_close),]

# load package
require('tikzDevice')
# tex package specification
options(tikzLatexPackages 
        = c(getOption( "tikzLatexPackages" ),
            "\\usepackage{amsmath,amsfonts,amsthm, palatino, mathpazo}"))
# file setting
tikz(file = 'COPtransience.tex', width = 5.5, height = 4, standAlone = TRUE)
# plot
plot(x = as.Date(ts$Date), y = ts$COP_Close, type = 'l',
     xlab = 'Date $t$', ylab = 'Price $y_{1,t}$',
     main = 'Shock transience of Conoco Phillips Stock Price $y_{1,t}$',
     col = 'deepskyblue')
# segments
segments(x0 = as.Date("2020-03-09"), y0 = 30, y1 = 60, lty = 3, col = 'magenta')
segments(x0 = as.Date("2020-03-06"), y0 = 30, y1 = 60, lty = 3, col = 'magenta')
segments(x0 = as.Date("2020-03-06"), x1 = as.Date("2020-05-15"),
         y0 = ts$COP_Close[which(ts$Date == as.Date("2020-03-06"))],
         col = 'indianred1', lty = 4)
# output
dev.off()
