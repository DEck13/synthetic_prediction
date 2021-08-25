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
ps.indic.W <- function(Tstar, Y, X, K, H, Ts, ell, B, bw, sig.levl = .05) {
  
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
        
        # lag
        yilag <- Y[[i]][-(Ti + 1)][(t + H - h + 1):(t + Ki + H - h)]
        
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
          yhat.adj <- c(yhat.adj, beta.hat.adj %*% matrix(c(1, ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0),
                                                            yhat.adj[j], X[[i]][-1, ][t + Ki + H - h + j, ])))
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

ps.indic.W.decay <- function(Tstar, Y, X, K, H, Ts, ell, B, bw, sig.levl = .05) {
  
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
        
        # lag
        yilag <- Y[[i]][-(Ti + 1)][(t + H - h + 1):(t + Ki + H - h)]
        
        # decay
        decay <- exp(- ((t + H - h + 1):(t + Ki + H - h) - Tstari - 1)) * ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0)
        
        # OLS
        lmodi.adj <- lm(yi ~ 1 + decay + yilag + xi)
        lmodi.unadj <- lm(yi ~ 1 + yilag + xi)
        
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
start_day_20080317 <- as.numeric(1:nrow(COP_close) == start)
COP_close <- COP_close %>% mutate(start_day_20080317 = start_day_20080317)
TS2 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]


# inflation adjustment
TS2[, 2:8] <- TS2[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008] 

#### 2008 shock effects

# shock effect date
start <- which(COP_close$Date == "2008-09-08")
start_day_20080908 <- as.numeric(1:nrow(COP_close) %in% start)
COP_close <- COP_close %>% mutate(start_day_20080908 = start_day_20080908)
TS3 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
# adjust for inflation
TS3[, 2:8] <- TS3[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008] 

#### Thursday, November 27, 2014

# shock effect date
start <- which(COP_close$Date == "2014-11-26")
start_day_20141127 <- as.numeric(1:nrow(COP_close) == start)
COP_close <- COP_close %>% mutate(start_day_20141127 = start_day_20141127)
# time window
TS4 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
# adjust for inflation
TS4[, 2:8] <- TS4[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2014] 



#### The March 9th, 2020 shock effect:

# shock effect date
start <- which(COP_close$Date == "2020-03-06")
start_day_20200309 <- as.numeric(1:nrow(COP_close) == start)
COP_close <- COP_close %>% mutate(start_day_20200309 = start_day_20200309)
# time window
TS1 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]

# Ti
Ts <- rep(15 + L, 4)
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

res1 <- ps.indic.W(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), H = H, Ts = Ts, ell = 4, B = 200, bw = 4)
res2 <- ps.indic.W.decay(Tstar = Tstar, Y = Y, X = X, K = rep(K, 4), H = H, Ts = Ts, ell = 4, B = 200, bw = 4)

ps1 <- res1$ps
ps2 <- res2$ps
abs(1 - pchisq( -2 * sum(log(ps1[-1])) , df = 2 * length(ps1[-1])) - ps1[1])

abs(1 - pchisq( -2 * sum(log(ps2[-1] + 1e-10)) , df = 2 * length(ps2[-1])) - ps2[1])

abs(1 - pchisq( -2 * sum(log(1 - ps1[-1] + 1e-10)) , df = 2 * length(ps1[-1])) - ps2[1])

abs(1 - pchisq( -2 * sum(log(1 - ps2[-1])) , df = 2 * length(ps2[-1])) - ps2[1])
