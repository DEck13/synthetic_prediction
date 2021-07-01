# load packages
require('forecast')
ncores <- 20


# this function returns the W^* estimated by synthetic control method (SCM)
scm <- function(X, Tstar) {
	# X is a list of covariates for disparate time series
	# X[[1]] should be the covariate of the time series to predict
	# X[[p]] for p = 2,...,n+1 are covariates for time series pool
	
	# T^* is a vector of shock-effects time points
	# shock effect point must be > 2
	
	# package for constrained optimization
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


# simulation study normal
sim.normal.gammaX <- function(mu.gamma.delta = 1, mu.alpha, sigma, 
															sigma.alpha, sigma.delta.gamma = 0.5, 
															p, B, n, H, scale, ell = ell, bw = 4) {
	K <- round(rgamma(n = n + 1, shape = 15, scale = 10)) # training sample size
	Ts <- round(rgamma(n = n + 1, shape = 15, scale = 10)) # Time Length
	Ts[which(Ts < 90)] <- 90
	K[which(K < 90)] <- 90
	Tstar <- c() # Shock Time Points
	for (i in 1:(n + 1)) {
		Tstar <- c(Tstar, sample((ceiling(1 / 4 * Ts[i]) + 1):(ceiling(3 / 4 * Ts[i]) + K[i] + H), size = 1))
	}
	phi <- round(runif(n + 1, 0, 1), 3) # autoregressive parameters
	
	X <- c()
	alpha <- c()
	# construction of design matrix and shock effects
	gamma <- c()
	for (i in 1:(n + 1)) {
		Ti <- Ts[i]
		Tstari <- Tstar[i]
		Ki <- K[i]
		X[[i]] <- matrix(rgamma(n = p * (Ti + Ki + H + 1), shape = 1, scale = scale), ncol = p, byrow = TRUE) 
		# matrix(rnorm(n = (Ti + 1) * p), ncol = p, byrow = T)
		# parameter setup
		gamma[[i]] <- matrix(rnorm(p, mean = mu.gamma.delta, sd = sigma.delta.gamma), nrow = 1)
		epsilontildei <- rnorm(n = 1, sd = sigma.alpha)
		# alpha
		alpha <- c(alpha, mu.alpha + gamma[[i]] %*% X[[i]][Tstari + 1, ] + epsilontildei)
	}
	
	
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
		yi <- yi0
		
		# parameter setup
		thetai <- matrix(rnorm(p), nrow = 1)
		etai <- rnorm(1)
		
		for (t in 2:(K[i] + Ts[i] + H + 1)) {
			epsilonit <- rnorm(n = 1, sd = sigma)
			yi <- c(yi, etai + alphai * ifelse(t >= Tstari + 2, yes = 1, no = 0) +
								phii * yi[t - 1] + thetai %*% xi[t, ] + epsilonit)
		}
		
		Y[[i]] <- yi
	}
	
	# estimates
	est <- ps.indic.W(Tstar = Tstar, Y = Y, X = X, K = K, H = H, Ts = Ts, ell = ell, B = B, bw = bw)
	ps <- est$ps
	W <- est$W
	Is <- est$Is
	# compute forecasts
	phat <- sum(W * ps[-1])
	Ihat <- sum(W * Is[-1])
	p.diff <- abs(phat - ps[1])
	I.diff <- abs(Ihat - Is[1])
	return(c(p.diff = p.diff, I.diff = I.diff))
}
