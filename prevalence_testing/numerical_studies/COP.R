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
  X1 <- X[[1]][Tstar[1] + 1, , drop = FALSE]
  
  # covariates for time series pool
  X0 <- c()
  for (i in 1:n) {
    X0[[i]] <- X[[i + 1]][Tstar[i + 1] + 1, , drop = FALSE]
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

# functions that return alpha.hat and synthetic weights
ps.indic.W.permanent <- function(Tstar, Y, X, K, H, Ts, ell, B, bw, sig.levl = .05) {
  
  n <- length(Y) - 1
  
  # empty
  ps <- c()
  
  AICs <- c()
  BICs <- c()
  for (i in 1:(n + 1)) {
    
    # set up
    Ti <- Ts[i]
    Ki <- K[i]
    Tstari <- Tstar[i]
    
    m1.L.i <- matrix(NA, nrow = Ti, ncol = H)
    m2.L.i <- matrix(NA, nrow = Ti, ncol = H)
    
    AICs.i <- c()
    BICs.i <- c()
    # compute losses
    for (h in 1:H) {
      AICs.h.i <- c()
      BICs.h.i <- c()
      for (t in 1:Ti) {
        
        yi <- Y[[i]][-1][(t + H - h + 1):(t + Ki + H - h)]
        xi <- X[[i]][-1, ][(t + H - h + 1):(t + Ki + H - h),]
        
        # lag
        yilag <- Y[[i]][-(Ti + 1)][(t + H - h + 1):(t + Ki + H - h)]
        
        # OLS
        lmodi.adj <- lm(yi ~ 1 + ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0) + 
                          yilag + xi)
        lmodi.unadj <- lm(yi ~ 1 + yilag + xi)
        AICs.h.i <- c(AICs.h.i, AIC(lmodi.adj))
        BICs.h.i <- c(BICs.h.i, BIC(lmodi.adj))
        
        # beta.hat
        beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
        beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
        beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
        
        yhat.adj <- tail(yi, 1)
        yhat.unadj <- tail(yi, 1)
        for (j in 1:h) {
          yhat.adj <- c(yhat.adj, beta.hat.adj %*% matrix(c(1, ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0),
                                                            yhat.adj[j], X[[i]][-1, ][t + Ki + H - h + j, ])))
          yhat.unadj <- c(yhat.unadj, beta.hat.unadj %*% matrix(c(1, yhat.unadj[j], X[[i]][-1, ][t + Ki + H - h + j, ])))
        }
        
        # losses
        m1.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.adj[h + 1])
        m2.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.unadj[h + 1])
      }
    }
    AICs.i <- c(AICs.i, mean(AICs.h.i))
    BICs.h.i <- c(BICs.h.i, BIC(lmodi.adj))
    # loss differential
    d <- m2.L.i - m1.L.i
    ps[i] <- taSPA.mhfc(ell = ell, d = d, B = B, bw = bw)
  }
  
  AIC <- mean(AICs.i)
  BIC <- mean(BICs.i)
  # weights
  if (is.matrix(X[[1]]) == FALSE) {
    for (i in 1:(n + 1)) {
      X[[i]] <- as.matrix(X[[i]])
    }
  }
  # Weights
  W <- round(scm(X = X, Tstar = Tstar)$par, digits = 3)
  
  # test
  Is <- ifelse(ps <= .05, yes = 1, no = 0)
  
  # output
  return(list(ps = ps, W = W, Is = Is, AIC = AIC, BIC = BIC))
}

ps.indic.W.decay <- function(Tstar, Y, X, K, H, Ts, ell, B, bw, sig.levl = .05) {
  
  n <- length(Y) - 1
  
  # empty
  ps <- c()
  
  AICs <- c()
  BICs <- c()
  for (i in 1:(n + 1)) {
    
    # set up
    Ti <- Ts[i]
    Ki <- K[i]
    Tstari <- Tstar[i]
    
    m1.L.i <- matrix(NA, nrow = Ti, ncol = H)
    m2.L.i <- matrix(NA, nrow = Ti, ncol = H)
    
    AICs.i <- c()
    BICs.i <- c()
    # compute losses
    for (h in 1:H) {
      AICs.h.i <- c()
      BICs.h.i <- c()
      for (t in 1:Ti) {
        
        yi <- Y[[i]][-1][(t + H - h + 1):(t + Ki + H - h)]
        xi <- X[[i]][-1, ][(t + H - h + 1):(t + Ki + H - h),]
        
        # lag
        yilag <- Y[[i]][-(Ti + 1)][(t + H - h + 1):(t + Ki + H - h)]
        
        # decay
        decay <- exp(- ((t + H - h + 1):(t + Ki + H - h) - Tstari - 1)) * ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0)
        
        # OLS
        lmodi.adj <- lm(yi ~ 1 + decay + yilag + xi)
        lmodi.unadj <- lm(yi ~ 1 + yilag + xi)
        AICs.h.i <- c(AICs.h.i, AIC(lmodi.adj))
        BICs.h.i <- c(BICs.h.i, BIC(lmodi.adj))
        
        # beta.hat
        beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
        beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
        beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
        
        yhat.adj <- tail(yi, 1)
        yhat.unadj <- tail(yi, 1)
        for (j in 1:h) {
          yhat.adj <- c(yhat.adj, beta.hat.adj %*% matrix(c(1, exp(-(t + Ki + H - h + j - Tstari - 1)) * 
                                                              ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0),
                                                            yhat.adj[j], X[[i]][-1, ][t + Ki + H - h + j, ])))
          yhat.unadj <- c(yhat.unadj, beta.hat.unadj %*% matrix(c(1, yhat.unadj[j], X[[i]][-1, ][t + Ki + H - h + j, ])))
        }
        
        # losses
        m1.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.adj[h + 1])
        m2.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.unadj[h + 1])
      }
      AICs.i <- c(AICs.i, mean(AICs.h.i))
      BICs.i <- c(BICs.i, mean(BICs.h.i))
    }
    # loss differential
    d <- m2.L.i - m1.L.i
    ps[i] <- taSPA.mhfc(ell = ell, d = d, B = B, bw = bw)
  }
  
  AIC <- mean(AICs.i)
  BIC <- mean(BICs.i)
  # weights
  if (is.matrix(X[[1]]) == FALSE) {
    for (i in 1:(n + 1)) {
      X[[i]] <- as.matrix(X[[i]])
    }
  }
  # Weights
  W <- round(scm(X = X, Tstar = Tstar)$par, digits = 3)
  
  # test
  Is <- ifelse(ps <= .05, yes = 1, no = 0)
  
  # output
  return(list(ps = ps, W = W, Is = Is, AIC = AIC, BIC = BIC))
}

# functions that return alpha.hat and synthetic weights for decaying shock effects
ps.indic.W.complicate <- function(Tstar, Y, X, K, H, Ts,
                                  q1, q2, 
                                  ell, B, bw, sig.levl = .05) {
  
  n <- length(Y) - 1
  
  # empty
  ps <- c()
  
  AICs <- c()
  BICs <- c()
  for (i in 1:(n + 1)) {
    
    # set up
    Ti <- Ts[i]
    Ki <- K[i]
    Tstari <- Tstar[i]
    
    m1.L.i <- matrix(NA, nrow = Ti, ncol = H)
    m2.L.i <- matrix(NA, nrow = Ti, ncol = H)
    
    AICs.i <- c()
    BICs.i <- c()
    # compute losses
    for (h in 1:H) {
      AICs.h.i <- c()
      BICs.h.i <- c()
      for (t in 1:Ti) {
        TL <- length(Y[[i]])
        yi <- Y[[i]][-(1:q1)][(t + H - h + 1):(t + Ki + H - h)]
        xi <- X[[i]][-(1:q1), ][(t + H - h + 1):(t + Ki + H - h),]
        
        # yi lags
        yilags <- c()
        for (j in 1:q1) {
          yilags <- cbind(yilags, Y[[i]][(t + H - h + 2 - j):(t + Ki + H - h - j + 1)])
        }
        
        # x and x lags
        x.xlags <- xi
        if (q2 > 1) {
          
          for (j in 1:(q2 - 1)) {
            x.xlags <- cbind(x.xlags, X[[i]][(t + H - h + 2 - j):(t + Ki + H - h - j + 1),])
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
        lmodi.adj <- lm(yi ~ 1 + yilags + x.xlags + yilags.D.i.t + x.xlags.D.i.t + D.i.t)
        lmodi.unadj <- lm(yi ~ 1 + yilags + x.xlags)
        AICs.h.i <- c(AICs.h.i, AIC(lmodi.adj))
        BICs.h.i <- c(BICs.h.i, BIC(lmodi.adj))
        
        # beta.hat
        beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
        beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
        beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
        
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
          
          
          yhat.adj <- c(yhat.adj, beta.hat.adj %*% matrix(c(1,   yhat.adj[j:(j + q1 - 1)], 
                                                            x.xlags.for.pred, 
                                                            yilags.for.pred.D.i.t,
                                                            x.xlags.for.pred.D.i.t,
                                                            D.i.t)))
          yhat.unadj <- c(yhat.unadj, beta.hat.unadj %*% matrix(c(1, yhat.adj[j:(j + q1 - 1)], 
                                                                  x.xlags.for.pred)))
        }
        
        # losses
        m1.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.adj[h + q1])
        m2.L.i[t, h] <- sel(y = Y[[i]][-1][t + Ki + H], yhat = yhat.unadj[h + q1])
      }
      AICs.i <- c(AICs.i, mean(AICs.h.i))
      BICs.i <- c(BICs.i, mean(BICs.h.i))
    }
    # loss differential
    d <- m2.L.i - m1.L.i
    ps[i] <- taSPA.mhfc(ell = ell, d = d, B = B, bw = bw)
  }
  
  AIC <- mean(AICs.i)
  BIC <- mean(BICs.i)
  
  # weights
  if (is.matrix(X[[1]]) == FALSE) {
    for (i in 1:(n + 1)) {
      X[[i]] <- as.matrix(X[[i]])
    }
  }
  # Weights
  W <- round(scm(X = X, Tstar = Tstar)$par, digits = 3)
  
  # test
  Is <- ifelse(ps <= .05, yes = 1, no = 0)
  
  # output
  return(list(ps = ps, W = W, Is = Is, AIC = AIC, BIC = BIC))
}

ps.indic.W.decay.nls <- function(Tstar, Y, X, K, H, Ts, ell, B, bw, sig.levl = .05) {
  
  n <- length(Y) - 1
  
  # empty
  ps <- c()
  
  for (i in 1:(n + 1)) {
    
    # set up
    Ti <- Ts[i]
    Ki <- K[i]
    Tstari <- Tstar[i]
    
    m1.L.i <- matrix(NA, nrow = Ti, ncol = H)
    m2.L.i <- matrix(NA, nrow = Ti, ncol = H)
    
    # compute losses
    for (h in 1:H) {
      for (t in 1:Ti) {
        
        yi <- Y[[i]][-1][(t + H - h + 1):(t + Ki + H - h)]
        xi <- X[[i]][-1, ][(t + H - h + 1):(t + Ki + H - h),]
        p <- ncol(xi)
        
        # lag
        yilag <- Y[[i]][-(Ti + 1)][(t + H - h + 1):(t + Ki + H - h)]
        
        # decay
        decay.indic <-  ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0)
        # OLS
        lmodi.adj <- nls(yi ~ eta + exp(- ((t + H - h + 1):(t + Ki + H - h) - Tstari - 1) * c) * alpha + yilag * phi + xi %*% beta,
                         start = list(eta = rnorm(1), c = rnorm(1), alpha = rnorm(1), phi = runif(1), beta = rnorm(p)))
        lmodi.unadj <- lm(yi ~ 1 + yilag + xi)
        
        # beta.hat
        beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
        beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
        beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
        
        yhat.adj <- tail(yi, 1)
        yhat.unadj <- tail(yi, 1)
        for (j in 1:h) {
          yhat.adj <- c(yhat.adj, 
                        beta.hat.adj[1] + exp(-(t + Ki + H - h + j - Tstari - 1) * beta.hat.adj[2]) * beta.hat.adj[3] + 
                          yhat.adj[j] * beta.hat.adj[4] + 
                          X[[i]][-1, ][t + Ki + H - h + j, ] %*% matrix(beta.hat.adj[-(1:4)]))
          yhat.unadj <- c(yhat.unadj, beta.hat.unadj %*% matrix(c(1, yhat.unadj[j], X[[i]][-1, ][t + Ki + H - h + j, ])))
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
  W <- round(scm(X = X, Tstar = Tstar)$par, digits = 3)
  
  # test
  Is <- ifelse(ps <= .05, yes = 1, no = 0)
  
  # output
  return(list(ps = ps, W = W, Is = Is))
}

# set working directory
setwd("/Users/mac/Desktop/Research/Post-Shock Prediction/")

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
inflation_adj <- inflation_adj %>% mutate(year = 1947:2021)


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
L <- 7

#### Monday, March 17th, 2008

## March 17th; 1 day nowcast

# shock effect date
start <- which(COP_close$Date == "2008-03-14")
start_day_20080317 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
COP_close <- COP_close %>% mutate(start_day_20080317 = start_day_20080317)
TS2 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
# inflation adjustment
TS2[, 2:8] <- TS2[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008] 

start2 <- which(TS2$Date == "2008-03-14")
expterm <- exp(-(1:nrow(TS2) - start2)) * ifelse(1:nrow(TS2) >= start2, yes = 1, no = 0)

m_COP_3_17.null <- lm(Y ~ COP_Close + GSPC_Close + WTI_Close + 
                   USD_Close + TB_Close, 
                 data = TS2)
m_COP_3_17 <- lm(Y ~ COP_Close + start_day_20080317 + GSPC_Close + WTI_Close + 
                   USD_Close + TB_Close, 
                 data = TS2)
m_COP.exp <- lm(Y ~ COP_Close + expterm + GSPC_Close + WTI_Close + 
                       USD_Close + TB_Close,  
                     data = TS2)



# decay

# NLS
m_COP.decay <- nls(Y ~ eta + exp(- (1:nrow(TS2) - (start2 + 1)) * c) * alpha + COP_Close * phi + 
                     cbind(GSPC_Close, WTI_Close, USD_Close, TB_Close) %*% beta, data = TS2,
                 start = list(eta = rnorm(1), c = rnorm(1), alpha = rnorm(1), phi = runif(1), beta = rnorm(4)))
# complicate
co.D <- c()
for (d in 1:nrow(TS2)) {
  co.D <- rbind(co.D, TS2[d, c(2:7, 9:14)] * TS2$start_day_20080317[d])
}
co.D <- as.matrix(co.D)

m_COP.complicate <- lm(TS2$Y ~ as.matrix(TS2[, c(2:7, 9:14)]) + co.D + TS2$start_day_20080317)

# AICs
AIC.20080314.null <- AIC(m_COP_3_17.null)
AIC.20080314.permanent <- AIC(m_COP_3_17)
AIC.20080314.decay1 <- AIC(m_COP.exp)
AIC.20080314.decay2 <- AIC(m_COP.decay)
AIC.20080314.complicate <- AIC(m_COP.complicate)

AIC.20080314 <- c(AIC.20080314.null, AIC.20080314.permanent, 
                  AIC.20080314.decay1, AIC.20080314.decay2,
                  AIC.20080314.complicate)

#### 2008 shock effects

# shock effect date
start <- which(COP_close$Date == "2008-09-08")
start_day_20080908 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
COP_close <- COP_close %>% mutate(start_day_20080908 = start_day_20080908)
TS3 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
# adjust for inflation
TS3[, 2:8] <- TS3[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008] 

start2 <- which(TS3$Date == "2008-09-08")
expterm <- exp(-(1:nrow(TS3) - start2)) * ifelse(1:nrow(TS3) >= start2, yes = 1, no = 0)

# models
m_COP_null <- lm(Y ~ COP_Close + GSPC_Close + WTI_Close + USD_Close + TB_Close, data = TS3)
m_COP_Sept_08 <- lm(Y ~ COP_Close + start_day_20080908 + GSPC_Close + WTI_Close + USD_Close + TB_Close, data = TS3)
m_COP_exp <- lm(Y ~ COP_Close + expterm + GSPC_Close + WTI_Close + 
                  USD_Close + TB_Close,  data = TS3)
m_COP_decay <- nls(Y ~ eta + exp(- (1:nrow(TS3) - (start2 + 1)) * c) * alpha + COP_Close * phi + 
                     cbind(GSPC_Close, WTI_Close, USD_Close, TB_Close) %*% beta, data = TS3,
                   start = list(eta = rnorm(1), c = rnorm(1), alpha = rnorm(1), phi = runif(1, min = -1), beta = rnorm(4)))
# complicate
co.D <- c()
for (d in 1:nrow(TS3)) {
  co.D <- rbind(co.D, TS3[d, c(2:7, 9:14)] * TS3$start_day_20080908[d])
}
co.D <- as.matrix(co.D)

m_COP.complicate <- lm(TS3$Y ~ as.matrix(TS3[, c(2:7, 9:14)]) + co.D + TS3$start_day_20080908)


# AICs
AIC.20080908.null <- AIC(m_COP_null) 
AIC.20080908.permanent <- AIC(m_COP_Sept_08) 
AIC.20080908.decay1 <- AIC(m_COP_exp) 
AIC.20080908.decay2 <- AIC(m_COP_decay) 
AIC.20080908.complicate <- AIC(m_COP.complicate) 

AIC.20080908 <- c(AIC.20080908.null, AIC.20080908.permanent, 
                  AIC.20080908.decay1, AIC.20080908.decay2,
                  AIC.20080908.complicate)


#### Thursday, November 27, 2014

# shock effect date
start <- which(COP_close$Date == "2014-11-26")
start_day_20141127 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
COP_close <- COP_close %>% mutate(start_day_20141127 = start_day_20141127)
# time window
TS4 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
# adjust for inflation
TS4[, 2:8] <- TS4[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2014] 

start2 <- which(TS4$Date == "2014-11-26")
expterm <- exp(-(1:nrow(TS4) - start2)) * ifelse(1:nrow(TS4) >= start2, yes = 1, no = 0)

# models
m_COP_11_27_14.null <- lm(Y ~ COP_Close + GSPC_Close + WTI_Close + USD_Close + TB_Close, data = TS4)
m_COP_11_27_14 <- lm(Y ~ COP_Close + start_day_20141127 + GSPC_Close + WTI_Close + USD_Close + TB_Close, 
                     data = TS4)
m_COP_exp <- lm(Y ~ COP_Close + expterm + GSPC_Close + WTI_Close + USD_Close + TB_Close,  data = TS4)
m_COP_decay <- nls(Y ~ eta + exp(- (1:nrow(TS4) - (start2 + 1)) * c) * alpha + COP_Close * phi + 
                     cbind(GSPC_Close, WTI_Close, USD_Close, TB_Close) %*% beta, data = TS4,
                   start = list(eta = rnorm(1), c = rnorm(1), alpha = rnorm(1), phi = runif(1, min = -1), beta = rnorm(4)))
# complicate
co.D <- c()
for (d in 1:nrow(TS4)) {
  co.D <- rbind(co.D, TS4[d, c(2:7, 9:14)] * TS4$start_day_20141127[d])
}
co.D <- as.matrix(co.D)

m_COP.complicate <- lm(TS4$Y ~ as.matrix(TS4[, c(2:7, 9:14)]) + co.D + TS4$start_day_20141127)


# AICs of different models
AIC.20141127.null <- AIC(m_COP_11_27_14.null) 
AIC.20141127.permanent <- AIC(m_COP_11_27_14) 
AIC.20141127.decay1 <- AIC(m_COP_exp) 
AIC.20141127.decay2 <- AIC(m_COP_decay) 
AIC.20141127.complicate <- AIC(m_COP.complicate) 

AIC.20141127 <- c(AIC.20141127.null, AIC.20141127.permanent, 
                  AIC.20141127.decay1, AIC.20141127.decay2,
                  AIC.20141127.complicate)


#### The March 9th, 2020 shock effect:

# shock effect date
start <- which(COP_close$Date == "2020-03-06")
start_day_20200309 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
COP_close <- COP_close %>% mutate(start_day_20200309 = start_day_20200309)
# time window
TS1 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]

start2 <- which(TS1$Date == "2020-03-06")
expterm <- exp(-(1:nrow(TS1) - start2)) * ifelse(1:nrow(TS1) >= start2, yes = 1, no = 0)

m_COP_03_09_20.null <- lm(Y ~ COP_Close + GSPC_Close + WTI_Close + 
                       USD_Close + TB_Close,  
                     data = TS1)
m_COP_03_09_20 <- lm(Y ~ COP_Close + start_day_20200309 + GSPC_Close + WTI_Close + 
                       USD_Close + TB_Close,  
                     data = TS1)
m_COP_exp <- lm(Y ~ COP_Close + expterm + GSPC_Close + WTI_Close + USD_Close + TB_Close,  data = TS1)
m_COP_decay <- nls(Y ~ eta + exp(- (1:nrow(TS1) - (start2 + 1)) * c) * alpha + COP_Close * phi + 
                     cbind(GSPC_Close, WTI_Close, USD_Close, TB_Close) %*% beta, data = TS1,
                   start = list(eta = rnorm(1), c = rnorm(1), alpha = rnorm(1), phi = runif(1, min = -1), beta = rnorm(4)))

# complicate
co.D <- c()
for (d in 1:nrow(TS4)) {
  co.D <- rbind(co.D, TS1[d, c(2:7, 9:14)] * TS1$start_day_20200309[d])
}
co.D <- as.matrix(co.D)

m_COP.complicate <- lm(TS1$Y ~ as.matrix(TS1[, c(2:7, 9:14)]) + co.D + TS1$start_day_20200309)


AIC.20200309.null <- AIC(m_COP_03_09_20.null)
AIC.20200309.permanent <- AIC(m_COP_03_09_20)
AIC.20200309.decay1 <- AIC(m_COP_exp)
AIC.20200309.decay2 <- AIC(m_COP_decay)
AIC.20200309.complicate <- AIC(m_COP.complicate)

AIC.20200309 <- c(AIC.20200309.null, AIC.20200309.permanent, 
                  AIC.20200309.decay1, AIC.20200309.decay2,
                  AIC.20200309.complicate)

AICs <- rbind(AIC.20080314, AIC.20080908, AIC.20141127)
rownames(AICs) <- c('TS2', 'TS3', 'TS4')
colnames(AICs) <- c('null', 'Permanent', 'Decay without scale', 'Decay with scale', 'Complicate')

apply(AICs, 1, which.min)
AICs <- rbind(AICs, apply(AICs, 2, mean))
rownames(AICs)[4] <- 'Mean'
AICs <- rbind(AICs, AIC.20200309)
rownames(AICs)[5] <- 'TS1' 
xtable(AICs)



# Ti
Ts <- c(nrow(TS1), nrow(TS2), nrow(TS3), nrow(TS4))
# Tstar
Tstar <- c(which(TS1$start_day_20200309 == 1),
           which(TS2$start_day_20080317 == 1),
           which(TS3$start_day_20080908 == 1),
           which(TS4$start_day_20141127 == 1))
# Y
Y <- c()
X <- c()
TS <- list(TS1, TS2, TS3, TS4)
for (i in 1:4) {
  Y[[i]] <- TS[[i]]$Y
  X[[i]] <- as.matrix(TS[[i]][, 3:7])
}

res1 <- ps.indic.W.permanent(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), H = H, Ts = Ts, ell = 4, B = 200, bw = 4)
res2 <- ps.indic.W.decay(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), H = H, Ts = Ts, ell = 4, B = 200, bw = 4)
res3 <- ps.indic.W.complicate(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), 
                              q1 = 2, q2 = 2, 
                              H = H, Ts = Ts, ell = 4, B = 200, bw = 4)
res2 <- ps.indic.W.decay.nls(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), H = H, Ts = Ts, ell = 4, B = 200, bw = 4)

ps1 <- res1$ps
ps2 <- res2$ps
ps3 <- res3$ps

fisher <- function(ps) {
  abs(1 - pchisq( -2 * sum(log(ps[-1] + 1e-10)) , df = 2 * length(ps[-1])) - ps[1])
}
Pearson <- function(ps) {
  abs(1 - pchisq( -2 * sum(log(1 - ps[-1] + 1e-10)) , df = 2 * length(ps[-1])) - ps[1])
}

fisher(ps1)
Pearson(ps1)
