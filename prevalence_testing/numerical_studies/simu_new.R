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

vote <- function(ps, sig = 0.05) {
  indic.p <- ifelse(ps < 0.05, yes = 1, no = 0)
  mean <- mean(indic.p)
  if (mean > 0.5) output <- 1
  else if (mean == 0.5) output <- sample(c(1, 0), size = 1)
  else output <- 0
  return(output)
}

Wp <- function(ps, W, type = c('Fisher', 'Pearson')) {
  if (type == 'Fisher') {
    Mis <- -log(ps + 1e-10)
  } else if (type == 'Pearson') {
    Mis <- -log(1 - ps + 1e-10)
  }
  M <- Mis * W
  p <- 1 - pgamma(2 * M, shape = 1 / (sum(W ^ 2)), scale = 2)
  return(p)
}


# functions that return alpha.hat and synthetic weights for decaying shock effects
ps.indic.W.decay <- function(Tstar, Y, X, K, H, Ts,
                             q1, q2, 
                             ell, B, bw, sig.levl = .05) {
  
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


# simulation study normal for decaying shock effects
sim.normal.gammaX <- function(mu.delta = 1, mu.alpha, sigma, 
                              sigma.alpha, sigma.delta = 0.1, 
                              B, n, H, scale, ell = ell, bw = 4,
                              Kscale = 1 / 2, Tscale = 1 / 2, 
                              Kshape, Tshape, p, q1, q2) {
  K <- ceiling(rgamma(n + 1, scale = Kscale, shape = Kshape)) # training sample size
  Ts <- ceiling(rgamma(n + 1, scale = Tscale, shape = Tshape)) # Time Length
  Tstar <- c()
  for (i in 1:(n + 1)) {
    Tstar <- c(Tstar,  max(Ts[i] + 1, ceiling(0.5 * (Ts[i] + K[i] + H))))
  }
  phi <- round(matrix(runif((n + 1) * q1, 0, 0.01), nrow = n + 1), 3) # autoregressive parameters
  
  X <- c()
  alpha <- c()
  # construction of design matrix and shock effects
  delta <- c()
  theta <- c()
  theta.tilde <- c()
  phi.tilde <- c()
  for (i in 1:(n + 1)) {
    Ti <- Ts[i]
    Tstari <- Tstar[i]
    Ki <- K[i]
    # generation of covariates
    X[[i]] <- matrix(rgamma(n = p * (Ti + Ki + H + q1 + q2), shape = 1, scale = scale), ncol = p, byrow = TRUE) 
    
    # generation of phi
    x.star <- X[[i]][Tstari + 1, ]
    phi.tilde[[i]] <- rbeta(q1, shape1 = 1, shape2 = sqrt(crossprod(x.star)))
    # matrix(rnorm(n = (Ti + 1) * p), ncol = p, byrow = T)
    
    # generation of theta 
    theta.i <- c()
    theta.i.tilde <- c()
    for (j in 1:q2) {
      theta.i.j <- matrix(rnorm(p), nrow = 1)
      theta.i[[j]] <- theta.i.j
      
      # generation of gamma
      beta.0 <- matrix(rnorm(p + 1, mean = 0, sd = 0.1))
      gamma.i.j <- matrix(rnorm(p, mean = matrix(c(1, x.star), nrow = 1) %*% matrix(beta.0), sd = 1), nrow = 1)
      
      # generation of theta.tilde
      theta.i.tilde.j <- theta.i.j + gamma.i.j
      theta.i.tilde[[j]] <- theta.i.tilde.j
    }
    theta[[i]] <- theta.i
    theta.tilde[[i]] <- theta.i.tilde
    
    # generation of delta
    delta[[i]] <- matrix(rnorm(p, mean = mu.delta, sd = sigma.delta), nrow = 1)
    
    # generation of noise 
    epsilontildei <- rnorm(n = 1, sd = sigma.alpha)
    
    # generation of alpha
    alpha.i <- mu.alpha + delta[[i]] %*% X[[i]][Tstari + 1, ] + epsilontildei
    alpha <- c(alpha, alpha.i)
  }
  
  
  # generation of yit
  Y <- c()
  for (i in 1:(n + 1)) {
    
    # initial value
    yi0 <- rnorm(q1)
    
    # setup
    Tstari <- Tstar[i]
    alphai <- alpha[i]
    phii <- phi[i, , drop = FALSE]
    xi <- X[[i]]
    yi <- yi0
    
    # parameter setup
    etai <- rnorm(1)
    
    for (t in 1:(K[i] + Ts[i] + H + q1)) {
      epsilonit <- rnorm(n = 1, sd = sigma)
      
      product.x <- 0
      # product of x and theta
      for (j in 0:(q2 - 1)) {
        product.x <- 0 + xi[t - j + q2, , drop = FALSE] %*% matrix(theta[[i]][[j + 1]])
      }
      
      product.x.tilde <- 0
      # product of x and theta
      for (j in 0:(q2 - 1)) {
        product.x.tilde <- product.x.tilde + xi[t - j + q2, , drop = FALSE] %*% matrix(theta.tilde[[i]][[j + 1]])
      }
      
      # functional of alpha
      f.alpha.i <- (alpha[i] + matrix(phi.tilde[[i]], nrow = 1) %*% matrix(yi[t: (t + q1 - 1)]) +
                      product.x.tilde) * ifelse(t >= Tstari + 1, yes = 1, no = 0)
      
      yi <- c(yi, etai + phii %*% matrix(yi[t: (t + q1 - 1)]) + product.x + f.alpha.i + epsilonit)
    }
    
    Y[[i]] <- yi
  }
  
  # estimates
  est <- ps.indic.W.decay(Tstar = Tstar, Y = Y, X = X, K = K, H = H, 
                          q1 = q1, q2 = q2,
                          Ts = Ts, ell = ell, B = B, bw = bw)
  ps <- est$ps
  W <- est$W
  Is <- est$Is
  # compute forecasts
  votes <- vote(ps = ps[-1], sig = 0.05)
  phat <- sum(W * ps[-1])
  p.mean <- mean(ps[-1])
  # difference
  p.diff <- abs(phat - ps[1])
  p.mean.diff <- abs(p.mean - ps[1])
  indic.diff <- abs(votes - Is[1])
  
  # Fisher
  p.fisher <- 1 - pchisq( -2 * sum(log(ps[-1] + 1e-10)) , df = 2 * length(ps[-1]))
  p.fisher.diff <- abs(p.fisher - ps[1])
  
  # Pearson
  p.pearson <- 1 - pchisq( -2 * sum(log(1 - ps[-1]  + 1e-10)) , df = 2 * length(ps[-1]))
  p.pearson.diff <- abs(p.pearson - ps[1])
  
  return(c(p.diff = p.diff, p.mean.diff = p.mean.diff,
           p.fisher.diff = p.fisher.diff, p.pearson.diff = p.pearson.diff,
           indic.diff = indic.diff))
}

system.time(result <- sim.normal.gammaX(mu.delta = 0.2, 
                                        mu.alpha = -0.2, sigma = 0.1, 
                                        sigma.alpha = 0.05, 
                                        sigma.delta = 0.1, 
                                        p = 2, B = 200, scale = 2, 
                                        n = 20, H = 8, ell = 4,
                                        Kshape = 100, Tshape = 100,
                                        q1 = 2, q2 = 3))


# MC
library("parallel")
library("doParallel")
library("foreach")
# 8 cores -- use 7
ncores <- detectCores() - 2
registerDoParallel(cores = ncores)
set.seed(2020)
RNGkind("L'Ecuyer-CMRG")
nsim <- 50


# parameter setup
shape.K.Ts <- c(200, 400, 600, 800)
ns <- c(5, 10, 20, 30)
sim_params <- expand.grid(list(shape.K.Ts = shape.K.Ts, ns = ns))

# simulation time
system.time(
  output <- lapply(1:nrow(sim_params), FUN = function(j) {
    # parameters
    shape.K.T <- sim_params[j, 1]
    n <- sim_params[j, 2]
    # %do% evaluates sequentially
    # %dopar% evaluates in parallel
    # .combine results
    out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
      # result
      study <- sim.normal.gammaX.decay(mu.gamma.delta = 2, 
                                       mu.alpha = 10, sigma = 0.1, 
                                       sigma.alpha = 0.05, 
                                       sigma.delta.gamma = 0.1, 
                                       p = 13, B = 200, scale = 2, 
                                       n = n, H = 8, ell = 4,
                                       Kshape = shape.K.T, Tshape = shape.K.T)
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
write_xlsx(lapply(output, as.data.frame), 'ntrial2.xlsx')

result <- c()
for (i in 1:nrow(sim_params)) {
  table <- output[[i]]
  # means and sds
  means <- apply(table, 2, function(x) mean(x))
  sds <- apply(table, 2, function(x) sd(x))
  result.i <- c()
  for (j in 1:5) {
    result.i <- cbind(result.i, paste0(round(means[j], digits = 3), 
                                       ' (', round(sds[j] / sqrt(100), 
                                                   digits = 3), ')'))
  }
  result <- rbind(result, result.i)
}
result <- cbind(sim_params[, c(1,2)], result)
require('xtable')
xtable(result)

