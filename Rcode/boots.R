# require
require('tikzDevice')
# set seed
set.seed(2020)
# parameter setup
n = 100; T = 250; mu.alpha = 5; sigma.alpha = 1; sigma.X = 1; sigma = 1

# T* is the shock time point
Tstar <- floor(T / 2)
# random effects
eta <- rnorm(n + 1)
# with normalization
phi <- runif(n + 1, 0, 1)
theta <- rnorm(n + 1)
theta <- theta / as.vector(sqrt(crossprod(theta)))
beta <- rnorm(n + 1)
beta <- beta / as.vector(sqrt(crossprod(beta)))
beta.X <- rnorm(n + 1)
beta.X <- beta.X / as.vector(sqrt(crossprod(beta.X)))
alpha <- rnorm(n + 1, mean = mu.alpha, sd = sigma.alpha)
# T: time points; n: different time series
X <- matrix(0, nrow = T + 1, ncol = n + 1)
X[1, ] <- rnorm(n + 1)
Y <- matrix(0, nrow = T + 1, ncol = n + 1)
Y[1,] <- rnorm(n + 1)

# i represents for time -- T
for(i in 2:(T + 1)){
  X[i,] <- beta.X * X[i - 1,] + rnorm(n + 1, sd = sigma.X)
  if (i - 2 == Tstar) {
    Y[i, ] <- eta + phi * Y[i - 1,] + beta * X[i,]  + rnorm(n + 1, sd = sigma) + alpha
  } else {
    Y[i, ] <- eta + phi * Y[i - 1,] + beta * X[i,]  + rnorm(n + 1, sd = sigma)
  } 
}

# compute mu alpha hat
mu.alpha.hat <- rep(0, n)
se <- rep(0, n)
for (j in 2:(n + 1)){
  # OLS Regression
  m1 <- lm(Y[2:(T + 1), j] ~ 1 + Y[1:T, j] + X[2:(T + 1), j] + ifelse((2:(T + 1)) == Tstar + 2, 1, 0))
  # find alpha hat
  mu.alpha.hat[j - 1] <- coefficients(m1)[4]
  se[j - 1] <- summary(m1)$coef[4, 2]
}

# Non-parametric Bootstrap of Shock-Effects
B <- 50000
IVWB <- c()
for (b in 1:B) {
  index <- base::sample(1:100, size = 100, replace = T)
  alphahatb <- mu.alpha.hat[index]
  seb <- se[index]
  IVWB <- c(IVWB, sum(alphahatb / seb ^ 2) /  (sum(1 / seb ^ 2)))
}
# histogram
hist(IVWB)



