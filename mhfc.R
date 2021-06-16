# squared error loss
sel <- function(yhat, y) {
  (y - yhat) ^ 2
}

# h = 1 to 24
# T = 120 length of series to evaluate
# y[145:246] is of interest
# n: training sample
T <- 240
n <- 120
H <- 12
m1.L <- matrix(NA, nrow = T, ncol = H)
m2.L <- matrix(NA, nrow = T, ncol = H)

# simulate a time series
y <- arima.sim(n = n + H + T, model = list(order = c(8, 0, 0), ar = rep(0.1, 8)))

# compute losses
for (h in 1:H) {
  for (t in 1:T) {
    # AR(3)
    m1.t.h <- arima(y[(t + H - h + 1):(t + n + H - h)], order = c(1, 0, 0))
    # ARMA(3,3)
    m2.t.h <- arima(y[(t + H - h + 1):(t + n + H - h)], order = c(8, 0, 0))
    m1.L[t, h] <- sel(y = y[t + n + H], yhat = predict(m1.t.h, h)$pred[h])
    m2.L[t, h] <- sel(y = y[t + n + H], yhat = predict(m2.t.h, h)$pred[h])
  }
}
d <- m1.L - m2.L

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
  B <- 500
  t.aSPA.B <- c()
  for (b in 1:B) {
    # uniform draw from 1,...,T - ell +1
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
# fail to reject
taSPA.mhfc(ell = 3, d = d, B = 1000)
