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

# function for dual shock
# this assumes two shocks occur at the same time point in the aggregated series
ps.indic.W.dynamic.additive <- function(Tstar, 
                                        Y, X, 
                                        Y0, X0,
                                        K, H, Ts,
                                        q1, q2, X1, 
                                        ell, B, bw, sig.levl = .05,  
                                        selfW = NA, scale = FALSE) {
  
  n <- length(Y) 
  
  # empty
  ps <- c()
  
  for (i in 1:n) {
    
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
        zi <- Y[[i]][-(1:max(q1, q2 - 1))][(t + H - h + 1):(t + Ki + H - h)] + 
          Y0[-(1:max(q1, q2 - 1))][(t + H - h + 1):(t + Ki + H - h)]
        xi <- X[[i]][-(1:max(q1, q2 - 1)), ][(t + H - h + 1):(t + Ki + H - h),]
        x0 <- X0[-(1:max(q1, q2 - 1)), ][(t + H - h + 1):(t + Ki + H - h),]
        
        # yi lags
        zilags <- c()
        for (j in 1:q1) {
          zijlag <- lag(Y[[i]][-(1:max(q1, q2 - 1))] + Y0[-(1:max(q1, q2 - 1))],
                        n = j)[(t + H - h + 1):(t + Ki + H - h)]
          zilags <- cbind(zilags, zijlag)
        }
        
        # x and x lags
        x.xlags <- xi
        if (q2 > 1) {
          for (j in 1:(q2 - 1)) {
            xj.xlags <- lag(X[[i]][-(1:max(q1, q2 - 1)), ], n = j)[(t + H - h + 1):(t + Ki + H - h), ]
            x.xlags <- cbind(x.xlags, xj.xlags)
          }
        }
        x.xlags.0 <- x0
        if (q2 > 1) {
          for (j in 1:(q2 - 1)) {
            xj.xlags.0 <- lag(X0[-(1:max(q1, q2 - 1)), ], n = j)[(t + H - h + 1):(t + Ki + H - h), ]
            x.xlags.0 <- cbind(x.xlags.0, xj.xlags.0)
          }
        }
        
        # functional
        D.i.t <- ifelse((t + H - h + 1):(t + Ki + H - h) >= Tstari + 1, yes = 1, no = 0)
        zilags.D.i.t <- zilags
        x.xlags.D.i.t <- x.xlags
        x.xlags0.D.i.t <- x.xlags.0
        
        for (d in 1:nrow(zilags)) {
          zilags.D.i.t[d, ] <- zilags[d, ] * D.i.t[d]
          x.xlags.D.i.t[d, ] <- x.xlags[d, ] * D.i.t[d]
          x.xlags0.D.i.t[d, ] <- x.xlags.0[d, ] * D.i.t[d]
        }
        
        # OLS
        lmodi.adj <- lm(zi ~ 1 + zilags + x.xlags + x.xlags.0 + 
                          zilags.D.i.t + x.xlags.D.i.t + 
                          x.xlags0.D.i.t + D.i.t)
        lmodi.unadj <- lm(zi ~ 1 + zilags + x.xlags + x.xlags.0)
        
        # beta.hat
        beta.hat.adj <- matrix(coef(lmodi.adj), nrow = 1)
        beta.hat.adj[which(is.na(beta.hat.adj) == TRUE)] <- 0
        beta.hat.unadj <- matrix(coef(lmodi.unadj), nrow = 1)
        beta.hat.unadj[which(is.na(beta.hat.unadj) == TRUE)] <- 0
        
        zhat.adj <- Y[[i]][(t + Ki + H - h + 1 - q1):(t + Ki + H - h)] + 
          Y0[(t + Ki + H - h + 1 - q1):(t + Ki + H - h)]
        zhat.unadj <- Y[[i]][(t + Ki + H - h + 1 - q1):(t + Ki + H - h)] + 
          Y0[(t + Ki + H - h + 1 - q1):(t + Ki + H - h)]
        
        for (j in 1:h) {
          x.xlags.for.pred <- X[[i]][-(1:q1),][t + Ki + H - h + j,]
          x.xlags0.for.pred <- X0[-(1:q1),][t + Ki + H - h + j,]
          if (q2 > 1) {
            for (k in 1:(q2 - 1)) {
              x.xlags.for.pred <- cbind(x.xlags.for.pred,
                                        X[[i]][-1,][t + Ki + H - h + j - k,])
              x.xlags0.for.pred <- cbind(x.xlags0.for.pred,
                                        X0[-1,][t + Ki + H - h + j - k,])
            }
          }
        
          D.i.t <- ifelse(t + Ki + H - h + j >= Tstari + 1, yes = 1, no = 0)
          x.xlags.for.pred.D.i.t <- x.xlags.for.pred * D.i.t
          zilags.for.pred.D.i.t <- zhat.adj[j:(j + q1 - 1)] * D.i.t
          
          zhat.adj.h <- beta.hat.adj %*% matrix(c(1, zhat.adj[j:(j + q1 - 1)], 
                                                  x.xlags.for.pred, x.xlags0.for.pred,
                                                  zilags.for.pred.D.i.t, x.xlags.for.pred.D.i.t, 
                                                  x.xlags0.for.pred, D.i.t))
          zhat.unadj.h <- beta.hat.unadj %*% matrix(c(1, zhat.adj[j:(j + q1 - 1)], 
                                                      x.xlags.for.pred, x.xlags0.for.pred))
          
          zhat.adj <- c(zhat.adj, zhat.adj.h)
          zhat.unadj <- c(zhat.unadj, zhat.unadj.h)
        }
          
          # losses
          m1.L.i[t, h] <- sel(y = Y[[i]][-(1:max(q1, q2 - 1))][t + Ki + H] + 
                                Y0[-(1:max(q1, q2 - 1))][t + Ki + H],
                              yhat = zhat.adj[h + q1])
          m2.L.i[t, h] <- sel(y = Y[[i]][-(1:max(q1, q2 - 1))][t + Ki + H] + 
                                Y0[-(1:max(q1, q2 - 1))][t + Ki + H], 
                              yhat = zhat.unadj[h + q1])
        }
    }
    # loss differential
    d <- m2.L.i - m1.L.i
    
    sub.index <- (Ti - Ki - max(q1, q2 - 1) - (Ti - Tstari) + 1):(Ti - Ki - H - max(q1, q2 - 1))
    d.subset <- d[sub.index,  ]
    ps[i] <- taSPA.mhfc(ell = ell, d = d.subset, B = B, bw = bw)
  }
  
  X1 <- as.matrix(X1)
  X0 <- as.matrix(X0)
  X.scm <- c()
  X.scm[[1]] <- as.matrix(cbind(X1, X1))
  for (i in 2:(n + 1)) {
    X.scm[[i]] <- cbind(as.matrix(X[[i - 1]]), X0)
  }
  
  # Weights
  W <- round(scm(X = X.scm, Tstar = Tstar, scale = scale)$par, digits = 3)
  
  if (is.na(selfW)[1] == FALSE) {
    W <- selfW
  }
  
  # test
  Is <- ifelse(ps <= .05, yes = 1, no = 0)
  
  # output
  return(list(ps = ps, W = W, Is = Is, dim = dim(d.subset)))
}

# set working directory
setwd("~/Desktop/Research/Post-Shock Prediction/")

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

Hs <- c(2, 4, 6, 8)
tab <- c()
for (H in Hs) {
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
  
  #### 2008 shock effects
  
  # shock effect date
  start <- which(COP_close$Date == "2008-09-08")
  start_day_20080908 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
  COP_close <- COP_close %>% mutate(start_day_20080908 = start_day_20080908)
  TS3 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
  # adjust for inflation
  TS3[, 2:8] <- TS3[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008] 
  
  
  #### Thursday, November 27, 2014
  
  # shock effect date
  start <- which(COP_close$Date == "2014-11-26")
  start_day_20141127 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
  COP_close <- COP_close %>% mutate(start_day_20141127 = start_day_20141127)
  # time window
  TS4 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
  # adjust for inflation
  TS4[, 2:8] <- TS4[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2014] 
  
  
  #### The March 9th, 2020 shock effect:
  
  # shock effect date
  start <- which(COP_close$Date == "2020-03-06")
  start_day_20200309 <- as.numeric(1:nrow(COP_close) %in% start:(start + L))
  COP_close <- COP_close %>% mutate(start_day_20200309 = start_day_20200309)
  # time window
  TS1 <- COP_close[(start - (1 + 15 + K + H)):(start + L),]
  
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
    X[[i]] <- as.matrix(TS[[i]][, c(3:7)])
  }
  
  res <- ps.indic.W.dynamic.additive(Tstar = Tstar, Y = Y[c(2, 3)], X = X[c(2, 3)], 
                                     Y0 = TS4$Y, X0 = as.matrix(TS4[, c(3:7)]),
                                     X1 = as.matrix(TS1[, c(3:7)]), q1 = 1, q2 = 1,
                                     Ts = Ts[c(2, 3)], K = rep(K, 2), H = H,
                                     ell = 4, B = 200, bw = 4, scale = TRUE)
  # p-values
  res.ps <- res$ps
  # phats
  res.phat <- sum(res$ps[-1] * res$W)
  # weighted indicators
  res.wi <- sum(res$Is[-1] * res$W)
  # voting
  res.vote <- ifelse(mean(res$Is[-1]) >= 0.5, yes = 1, no = 0)
  # weighted voting
  res.wvote <- ifelse(res.wi >= 0.5, yes = 1, no = 0)
  tab.i <- data.frame(H = H, phat = res.phat, vote.output = res.vote, 
                      wvote.output = res.wvote)
  tab <- rbind(tab, tab.i)
}

