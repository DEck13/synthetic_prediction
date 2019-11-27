


simstudy <- function(J = 20, T = 20, mu.alpha = 4, 
  sigma.alpha = 1, sigma.X = 1, sigma = 1){

  Tstar <- floor(T/2)
  eta <- rnorm(J+1)
  theta <- rnorm(J+1)
  theta <- theta / as.vector(sqrt(crossprod(theta)))
  beta <- rnorm(J+1)
  beta <- beta / as.vector(sqrt(crossprod(beta)))
  beta.X <- rnorm(J+1)
  beta.X <- beta.X / as.vector(sqrt(crossprod(beta.X)))
  alpha <- rnorm(J+1, mean = mu.alpha, sd = sigma.alpha)

  X <- matrix(0, nrow = T+1, ncol = J+1)
  X[1, ] <- rnorm(J+1)
  Y <- matrix(0, nrow = T+1, ncol = J+1)
  Y[1, ] <- rnorm(J+1)

  for(i in 2:(T+1)){
    X[i, ] <- beta.X * X[i-1, ] + rnorm(J+1, sd = sigma.X)
    Y[i, ] <- eta + theta * Y[i-1, ] + beta * X[i, ] + alpha * ifelse(i-1 > Tstar, 1, 0) + rnorm(J+1, sd = sigma) 
  }

  mu.alpha.hat <- rep(0, J)
  for(j in 2:(J+1)){
    m1 <- lm(Y[2:(T+1), j] ~ 1 + Y[1:T, j] + X[2:(T+1), j] + ifelse(1:T > Tstar,1,0))
    mu.alpha.hat[j-1] <- coefficients(m1)[4]
  }

  alpha.hat <- mean(mu.alpha.hat)
  #alpha.hat
  #sd(mu.alpha.hat) / sqrt(length(mu.alpha.hat))

  m1 <- lm(Y[2:(Tstar+1), 1] ~ 1 + Y[1:Tstar, 1] + X[2:(Tstar+1), 1])
  beta.hat <- coefficients(m1)
  Ypred.original <- as.vector(beta.hat %*% c(1, Y[Tstar+1,1], X[Tstar+2,1]))
  Ypred.adjusted <- Ypred.original + alpha.hat

  R.original <- (Ypred.original - Y[Tstar+2,1])^2
  R.adjusted <- (Ypred.adjusted - Y[Tstar+2,1])^2

  c(R.original, R.adjusted)
}


library(matrixStats)
set.seed(13)

n.sim <- 100
test <- sapply(1:n.sim, FUN = function(j) 
  simstudy(J = 100, T = 10, sigma.alpha = 1, sigma = 1))
rowMeans(test) # average forecast error for each forecast
rowSds(test) / sqrt(n.sim) # Monte Carlo standard errors corresponding to the above averages

mean(apply(test, 2, which.min) - 1) # proportion that adjusting wins
diffs <- apply(test, 2, diff)

# average reduction achieved from adjustment
# theory says this value should be close to -16 under the current function defualts
mean(diffs)

# Monte Carlo standard error corresponding to the above average
sd(diffs) / sqrt(n.sim) 



test <- sapply(1:n.sim, FUN = function(j) 
  simstudy(J = 100, T = 100, sigma.alpha = 1, sigma = 1))
rowMeans(test) # average forecast error for each forecast
rowSds(test) / sqrt(n.sim) # Monte Carlo standard errors corresponding to the above averages

mean(apply(test, 2, which.min) - 1) # proportion that adjusting wins
diffs <- apply(test, 2, diff)

# average reduction achieved from adjustment
# theory says this value should be close to -16 under the current function defualts
mean(diffs)

# Monte Carlo standard error corresponding to the above average
sd(diffs) / sqrt(n.sim) 
