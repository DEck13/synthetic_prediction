


library(data.table) 


# squared error loss
sel <- function(yhat, y) {
	(y - yhat) ^ 2
}


data_setup <- function(T, n, H, lag, shock, ar){

	y <- arima.sim(n = n + H + T, 
								 model = list(order = c(lag, 0, 0), 
								 						 ar = rep(ar, lag)))
	
	shift <- rep(0, length(y))
	shift[1:(T+n+H) > (T+n)/2] <- shock
	y <- y + shift
	
	df <- data.frame(t=y)
	setDT(df)[, paste("t", 1:lag, sep = "") := shift(t, 1:lag)]
	df$shock <- as.factor(as.numeric(1:nrow(df) > (T+n)/2))
  
	output = list(df = df, shift = shift)
	output

}



diff <- function(dat, T, H, lag){
	df <- dat$df
	shift <- dat$shift
	
	m <- lm(t ~ ., data = df)
	shock_hat <- coef(m)[lag + 2]
	shift_hat <- rep(0, nrow(df))
	shift_hat[1:(T+n+H) > (T+n)/2] <- shock_hat
	
	m1.L <- matrix(NA, nrow = T, ncol = H)
	m2.L <- matrix(NA, nrow = T, ncol = H)
	
	# compute losses
	y <- df$t
	for (h in 1:H) {
		for (t in 1:T) {
			# AR(lag) with shift removed
			m1.t.h <- arima(y[(t + H - h + 1):(t + n + H - h)] - 
												shift[(t + H - h + 1):(t + n + H - h)], 
											order = c(lag, 0, 0))
			# AR(lag)
			m2.t.h <- arima(y[(t + H - h + 1):(t + n + H - h)], 
											order = c(lag, 0, 0))
			m1.L[t, h] <- sel(y = y[t + n + H], yhat = predict(m1.t.h, h)$pred[h])
			m2.L[t, h] <- sel(y = y[t + n + H], yhat = predict(m2.t.h, h)$pred[h])
		}
	}
	d <- m1.L - m2.L
	d
	
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
	T <- nrow(d); H <- ncol(d)
	# T = ell * K
	K <- T / ell
	
	# Bootstrap replication
	t.aSPA.B <- c()
	for (b in 1:B) {
		# uniform draw from 1,...,T - ell + 1
		Iks <- sample(1:(T - ell + 1), K, replace = TRUE)
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
		t.aSPA.B[b] <- sqrt(T) * (d.b.bar - d.bar) / xi.b
	}
	# compute t.aSPA statistic
	## compute Heteroskedasticity and autocorrelation Consistent (HAC) estimator
	require('tsapp')
	Omega <- HAC(d, method = 'Quadratic Spectral', bw = bw)
	w <- matrix(1 / H, nrow = H)
	xi.hat <- sqrt(t(w) %*% Omega %*% w)
	t.aSPA <- as.numeric(sqrt(T) * d.bar / xi.hat)
	p <- mean(t.aSPA < t.aSPA.B)
	# return output
	return(p)
}



# decision
# h = 1 to 24
# T = 120 length of series to evaluate
# y[145:246] is of interest
# n: training sample
set.seed(13)
T <- 150
n <- 50
H <- 20
lag <- 5
shock <- 1
ar <- 0.10

dat <- data_setup(T = T, n = n, H = H, lag = lag, shock = shock, ar = ar)
d <- diff(dat = dat, T = T, H = H, lag = lag)
taSPA.mhfc(ell = ceiling(T^(1/2-0.01)), d = d, B = 3000)

