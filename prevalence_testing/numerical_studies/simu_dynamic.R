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
      sub.index <- (Ti - Ki - max(q1, q2 - 1) - (Ti - Tstari) + 1):(Ti - Ki - H - max(q1, q2 - 1))
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

# simulation study normal for decaying shock effects
sim <- function(mu.delta = 1, mu.alpha, sigma, 
                sigma.alpha, sigma.delta = 0.1, 
                p, B, n, H, c = 2, ell = ell, bw = 4,
                q1, q2, lambda1, lambda2,
                L, tshape, tscale) {
  # training sample size
  Ks <- ceiling(rgamma(n = n + 1, shape = tshape, scale = tscale))
  # shock time point
  Tstar <- c()
  for (i in 1:(n + 1)) {
    Tstar <- c(Tstar,  Ks[i] + max(Ks[i] + 1, ceiling(rgamma(n = 1, shape = tshape, scale = tscale))))
  }
  # Time length
  Ts <- Tstar + L + H
  # autoregressive parameters
  phis <- matrix(runif(q1 * (n + 1), 0, 1), ncol = q1)
  thetas <- matrix(runif(q2 * (n + 1), 0, 1), ncol = q2)
  
  # construction of design matrix and shock effects
  X <- c()
  alpha <- c()
  delta <- c()
  phis.tilde <- c()
  thetas.tilde <- c()
  thetas <- c()
  a <- rexp(n = 1, rate = lambda1)
  b <- rexp(n = 1, rate = lambda2)
  for (i in 1:(n + 1)) {
    Ti <- Ts[i]
    Tstari <- Tstar[i]
    Ki <- Ks[i]
    X[[i]] <- matrix(rgamma(n = p * (Ti + q2 - 1), shape = 1, scale = c), 
                     ncol = p, byrow = TRUE) 
    # matrix(rnorm(n = (Ti + 1) * p), ncol = p, byrow = T)
    # parameter setup
    delta[[i]] <- matrix(rnorm(p, mean = mu.delta, sd = sigma.delta), nrow = 1)
    epsilontildei <- rnorm(n = 1, sd = sigma.alpha)
    # alpha
    alpha <- c(alpha, mu.alpha + delta[[i]] %*% X[[i]][Tstari + 1 + q2 - 1, ] + epsilontildei)
    norm.i <- sqrt(crossprod(X[[i]][Tstari + 1 + q2 - 1, ]))
    phis.tilde <- rbind(phis.tilde, 
                        rbeta(n = q1, shape1 = a * norm.i + b, shape2 = 1))
    thetas[[i]] <- matrix(runif(q2 * p, 0, 1), ncol = p)
    thetas.tilde[[i]] <- matrix(rbeta(n = q2 * p, 
                                      shape1 = a * norm.i + b, 
                                      shape2 = a * norm.i + b),
                                ncol = p)
  }
  
  # autoregressive parameters
  
  
  # generation of yit
  Y <- c()
  for (i in 1:(n + 1)) {
    
    # initial value
    yi0 <- rnorm(q1)
    
    # setup
    Tstari <- Tstar[i]
    alphai <- alpha[i]
    xi <- X[[i]]
    yi <- yi0
    
    # parameter setup
    etai <- rnorm(1)
    
    for (t in (max(q1, q2 - 1) + 1):(Ts[i] + max(q1, q2 - 1))) {
      epsilonit <- rnorm(n = 1, sd = sigma)
      D.it <- ifelse(t >= Tstari + 1 + q1, yes = 1, no = 0)
      
      theta.t.sum <- 0
      theta.t.sum.tilde <- 0
      for (k in 1:q2) {
        theta.t.sum <- theta.t.sum + 
          xi[t - (k - 1) - abs(q1 - q2 + 1), ] %*% matrix(thetas[[i]][k, ])
        theta.t.sum.tilde <- theta.t.sum.tilde +
          xi[t - (k - 1) - abs(q1 - q2 + 1), ] %*% matrix(thetas.tilde[[i]][k, ])
      }
      
      
      yi <- c(yi, etai + (1 - D.it) * 
                (phis[i, , drop = FALSE] %*% matrix(yi[(t - q1):(t - 1)]) + 
                   theta.t.sum) + 
                D.it * (phis.tilde[i, , drop = FALSE] %*% matrix(yi[(t - q1):(t - 1)]) + 
                          theta.t.sum.tilde) + epsilonit)
    }
    
    Y[[i]] <- yi[-(1:q1)]
    if (q2 - 1 > 0) {
      X[[i]] <- xi[-(1:(q2 - 1)),]
    } else {
      X[[i]] <- xi
    }
  }
  
  # estimates
  est <- ps.indic.W.dynamic(Tstar = Tstar,
                            Y = Y, X = X, 
                            K = Ks, H = H, Ts = Ts, 
                            q1 = q1, q2 = q2, scale = TRUE,
                            ell = ell, B = B, bw = bw)
  ps <- est$ps
  W <- est$W
  Is <- est$Is
  # compute forecasts
  votes <- vote(ps = ps[-1], sig = 0.05)
  wvotes <- vote(ps = ps[-1], sig = 0.05, weight = TRUE, W = W)
  phat <- sum(W * ps[-1])
  # difference
  p.diff <- abs(phat - ps[1])
  indic.diff <- abs(votes - Is[1])
  indic.diff.W <- abs(wvotes - Is[1])
  
  return(c(p.diff = p.diff, indic.diff = indic.diff, indic.diff.W))
}


# MC
library("parallel")
library("doParallel")
library("foreach")
# 8 cores -- use 7
ncores <- detectCores() - 2
registerDoParallel(cores = ncores)
set.seed(2022)
RNGkind("L'Ecuyer-CMRG")
nsim <- 200


# parameter setup
paras <- c(1 / 5 ^ 2, 1 / 5, 1, 5, 5 ^ 2)

# simulation time
system.time(
  output <- lapply(1:length(paras), FUN = function(j) {
    # %do% evaluates sequentially
    # %dopar% evaluates in parallel
    # .combine results
    lambda2.j <- paras[j]
    out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
      # result
      study <- sim(mu.delta = 1, mu.alpha = 10, sigma = 1, 
                   sigma.alpha = 2, sigma.delta = 0.05, 
                   p = 4, B = 200, H = 10, ell = 4, bw = 4,
                   q1 = 2, q2 = 2, c = 2, lambda1 = 50, lambda2 = lambda2.j,
                   L = 30, n = 10, tshape = 15, tscale = 10)
      return(study)
    }
    # return results
    out
  })
)
