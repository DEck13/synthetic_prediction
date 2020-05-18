# require
require('tikzDevice')
# set seed
set.seed(2017)
# parameter setup
n = 40; T = 250; mu.alpha = 200; sigma.alpha = 100; sigma.X = 100; sigma = 100

# T* is the shock time point
Tstar <- floor(T / 2)
# random effects
eta <- rnorm(n + 1)
# with normalization
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
  Y[i, ] <- eta + theta * Y[i - 1,] + beta * X[i,] + 
    alpha * ifelse(i - 1 > Tstar, 1, 0) + rnorm(n + 1, sd = sigma) 
}

# compute mu alpha hat
mu.alpha.hat <- rep(0, n)
for (j in 2:(n + 1)){
  # OLS Regression
  m1 <- lm(Y[2:(T + 1), j] ~ 1 + Y[1:T, j] + X[2:(T + 1), j] + ifelse(1:T > Tstar, 1, 0))
  # find alpha hat
  mu.alpha.hat[j - 1] <- coefficients(m1)[4]
}

alpha.hat <- mean(mu.alpha.hat)
# alpha.hat
# sd(mu.alpha.hat) / sqrt(length(mu.alpha.hat))

# compute estimate for Tstar without post shock
m1 <- lm(Y[2:(Tstar + 1), 1] ~ 1 + Y[1:Tstar, 1] + X[2:(Tstar + 1), 1])
beta.hat <- coefficients(m1)
Ypred.original <- as.vector(cbind(1, Y[1:T, 1], X[2:(T+1), 1]) %*% beta.hat)
Ypred.adjusted <- Ypred.original + alpha.hat * as.matrix(ifelse(2:(T + 1) > Tstar, 1, 0))

# set working directory
setwd('/Users/mac/Desktop/Research/Post-Shock Prediction/')
# figure setting
tikz('comp.tex', standAlone = TRUE, width = 5.5, height = 3)
# plot
par(mar = c(4, 4, 1, 2))
# original
matplot(2:(T + 1), cbind(Y[2:(T + 1), 1], Ypred.original, Ypred.adjusted), 
        xlab = 'Time $T$', lwd = 1.5,
        ylab = 'Response', lty = c(1, 2, 3),
        col = c('black', 'indianred1', 'magenta'),
        type = 'l')
# add shock
segments(x0 = Tstar, y0 = -400, y1 = 600, col = 'deepskyblue', lty = 2, lwd = 2)
text(x = Tstar + 40, y = -350, 'shock time point $T^*_1$')
# legends
legend(x = 0, y = 600, 
       col = c('indianred1', 'magenta'),
       lwd = 2, lty = c(2, 3),
       c('$\\hat{Y}$ without shock effects', '$\\hat{Y}$ with shock effects'))
# output
dev.off()