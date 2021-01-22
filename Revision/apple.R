# load package
library("quantmod")
library('dplyr')
library('Rsolnp')
library('ggplot2')
library('tikzDevice')
# set working directory
setwd("/Users/mac/Desktop/Research/Post-Shock Prediction/")

# Apple 
getSymbols('AAPL', from = "2000-01-01")
AAPL <- as.data.frame(AAPL)
AAPL <- AAPL %>% mutate(Date = rownames(AAPL))


## S&P 500
getSymbols('^GSPC', from = "1970-01-01")
GSPC <- as.data.frame(GSPC)
GSPC <- GSPC %>% mutate(Date = rownames(GSPC))

## 13-Week T-Bill
TB <- read.csv('^IRX.csv', na.strings = 'null')
TB <- na.omit(TB)


## inflation adjustment
getSymbols("CPIAUCSL", src = 'FRED')
avg.cpi <- apply.yearly(CPIAUCSL, mean)
inflation_adj <- as.numeric(avg.cpi['2020'])/avg.cpi
inflation_adj <- as.data.frame(inflation_adj)
colnames(inflation_adj) <- c("dollars_2020")
inflation_adj <- inflation_adj %>% mutate(year = 1947:2020)


# Data Preparation
AAPL_close <- AAPL %>% dplyr::select(AAPL.Close, Date) %>% rename(AAPL_Close = AAPL.Close)
GSPC_close <- GSPC %>% dplyr::select(GSPC.Close, Date) %>% rename(GSPC_Close = GSPC.Close)
TB_close <- TB %>% dplyr::select(Close, Date) %>% rename(TB_Close = Close)

tom <- list(AAPL_close, GSPC_close, TB_close)
for (i in 1:length(tom)) {
  AAPL_close <- merge(AAPL_close, tom[[i]])
}

# response
Y <- AAPL_close$AAPL_Close[-1]

# data frame
AAPL_close <- data.frame(AAPL_close[-nrow(AAPL_close), ], Y)


# shock effect date
start <- which(AAPL_close$Date %in% c("2017-06-06", "2017-06-07"))
start_day_201706 <- as.numeric(1:nrow(AAPL_close) %in% start)
AAPL_close <- AAPL_close %>% mutate(start_day_201706 = start_day_201706)
TS2 <- AAPL_close[(start[1] - 5):(start[1] + 10), ]
# inflation adjustment
TS2[, 2:5] <- TS2[, 2:5] * inflation_adj$dollars_2020[inflation_adj$year == 2017] 
m_AAPL_201706 <- lm(Y ~ AAPL_Close + start_day_201706 + GSPC_Close + TB_Close, 
                 data = TS2)
alpha_201706 <- summary(m_AAPL_201706)$coef[3,1:2] 
# shock-effects
alpha_201706


# shock effect date
start <- which(AAPL_close$Date %in% c("2018-10-29", "2018-10-30"))
start_day_201810 <- as.numeric(1:nrow(AAPL_close) %in% start)
AAPL_close <- AAPL_close %>% mutate(start_day_201810 = start_day_201810)
TS3 <- AAPL_close[(start[1] - 5):(start[1] + 10), ]
# inflation adjustment
TS3[, 2:5] <- TS3[, 2:5] * inflation_adj$dollars_2020[inflation_adj$year == 2018] 
m_AAPL_201810 <- lm(Y ~ AAPL_Close + start_day_201810 + GSPC_Close + TB_Close, 
                    data = TS3)
alpha_201810 <- summary(m_AAPL_201810)$coef[3,1:2] 
# shock-effects
alpha_201810

# shock effect date
start <- which(AAPL_close$Date %in% c("2019-06-03", "2019-06-04"))
start_day_201906 <- as.numeric(1:nrow(AAPL_close) %in% start)
AAPL_close <- AAPL_close %>% mutate(start_day_201906 = start_day_201906)
TS4 <- AAPL_close[(start[1] - 5):(start[1] + 10), ]
# inflation adjustment
TS4[, 2:5] <- TS4[, 2:5] * inflation_adj$dollars_2020[inflation_adj$year == 2019] 
m_AAPL_201906 <- lm(Y ~ AAPL_Close + start_day_201906 + GSPC_Close + TB_Close, 
                    data = TS4)
alpha_201906 <- summary(m_AAPL_201906)$coef[3,1:2] 
# shock-effects
alpha_201906

# IVW Estimator
estimates <- rbind(alpha_201706, alpha_201810, alpha_201906)
estimates[, 2] <- estimates[, 2] ^ 2
colnames(estimates) <- c("alpha_hat", "var")

# adjustment estimator
alpha_adj <- mean(estimates[, 1])

# IVW estimator
weights <- (1 / estimates[,2]) / sum(1 / estimates[, 2])
alpha_IVW <- sum(weights * estimates[, 1])


# TS1 
start <- which(AAPL_close$Date == "2020-11-06")
TS1 <- AAPL_close[(start - 15):start, ]

# ar(TS1$GSPC_Close) chooses 3
sp500.ar <- arima(x = TS1$GSPC_Close, order = c(3, 0, 0))
# ar(TS1$TB_Close) chooses 1
TB.ar <- arima(TS1$TB_Close,  order = c(1, 0, 0))


# predicted covariates
X1 <- cbind(predict(sp500.ar, n.ahead = 1)$pred,
            predict(TB.ar, n.ahead = 1)$pred)
X1 <- as.matrix(X1)
X1 <- rbind(as.matrix(TS1[nrow(TS1), c(3, 4)]), X1)
X1p <- X1

# weighted adjustment estimator
Tstar.Date <- c("2017-06-05", "2018-10-26", "2019-05-31")
Tstar <- sapply(Tstar.Date, function(x) which(AAPL_close$Date == x))
# X1 <- as.matrix(COP_close[c(Tstar[1], Tstar[1] + 1), c(3, 4)])
# X0
X0 <- c()
for (i in 1:3) {
  X0[[i]] <- as.matrix(AAPL_close[c(Tstar[i] + 1, Tstar[i] + 2), 3:4])
}

# SCM
dat <- scale(rbind(X1, do.call('rbind', X0)), center = T, scale = T)
X1 <- dat[1:2, ]
X0 <- c()
for (i in 1:3) {
  X0[[i]] <- dat[c(2 * i + 1, 2 * i + 2), , drop = FALSE]
}


scmm <- function(X1, X0) {
  weightedX0 <- function(W) {
    # W is a vector of weight of the same length of X0
    n <- length(W)
    p <- ncol(X1)
    XW <- matrix(0, nrow = 2, ncol = p)
    for (i in 1:n) {
      XW <- XW + W[i] * X0[[i]]
    }
    norm <- as.numeric(crossprod(matrix(X1 - XW)))
    return(norm)
  }
  # constraint for W
  Wcons <- function(W) sum(W) - 1
  n <- length(X0)
  # optimization
  outs <- solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, eqB = 0, LB = rep(0, n), UB = rep(1, n))
  
  # output weights
  Wstar <- outs$pars
  
  return(Wstar)
}

# objective function is not 0; the fit may not be good
Wstar <- scmm(X1 = X1, X0 = X0)

weightedX0 <- function(W) {
  # W is a vector of weight of the same length of X0
  n <- length(W)
  p <- ncol(X1)
  XW <- matrix(0, nrow = 2, ncol = p)
  for (i in 1:n) {
    XW <- XW + W[i] * X0[[i]]
  }
  norm <- as.numeric(crossprod(matrix(X1 - XW)))
  return(norm)
}
# constraint for W
Wcons <- function(W) sum(W) - 1
n <- length(X0)
# optimization
outs <- solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, eqB = 0, LB = rep(0, n), UB = rep(1, n))

# output weights
Wstar <- outs$pars
weightedX0(round(outs$pars, digits = 3))

# weighted adjustment estimator
alpha_wadj <- sum(Wstar * estimates[,1])



# Parametric Bootstrap Estimation

# Set a seed
set.seed(2020)
# Bootstrap replications
B <- 1000
# List of linear models
lmod <- list(m_AAPL_201706, m_AAPL_201810, m_AAPL_201906)
# List of Data
TS <- list(TS2, TS3, TS4)
# List of T*
Tstar.Date <- c('2020-11-06', "2017-06-05", "2018-10-26", "2019-05-31")
# Empty List for storation
alphas <- vector(mode = 'list', length = 3)
# Loop begins
for (b in 1:B) {
  
  # Vector for storing alpha hats
  alphahatsb <- c()
  # Weights for IVW Estimator
  weights <- c()
  
  for (i in 1:3) {
    # preparation
    res <- residuals(lmod[[i]])
    dat <- TS[[i]]
    Ti <- nrow(dat)
    coef <- matrix(coef(lmod[[i]]), nrow = 1)
    
    # BOOTSTRAP
    resb <- sample(res, size = Ti, replace = TRUE)
    
    # New response
    yi0 <- dat$AAPL_Close[1]
    yib <- yi0
    Tstari <- which(dat$Date == Tstar.Date[i + 1])
      
      for (t in 1:Ti) {
        datt <- matrix(c(1, yib[t], ifelse(t %in% c(Tstari + 1, Tstari + 2), 
                                           yes = 1, no = 0), 
                         as.numeric(dat[t, c('GSPC_Close', 'TB_Close')])))
        yib <- c(yib, resb[t] + coef %*% datt)
      }
      
      # Prepare for new data
      yb <- yib[-1]; yblag <- yib[-(Ti + 1)]
      datbi <- data.frame(yblag, ifelse(1:Ti %in% c(Tstari + 1, Tstari + 2), yes = 1, no = 0),
                          dat[, c('GSPC_Close', 'TB_Close')])
      
      # New colnames
      colnames(datbi) <- c('yblag', 'shock', c('GSPC_Close', 'TB_Close'))
      
      # New Linear Model
      lmodbi <- lm(yb ~ 1 + yblag + shock + GSPC_Close + TB_Close, dat = datbi)
      
      # Shock Effects
      alphahatsb <- c(alphahatsb, coef(lmodbi)[3])
      
      # Weights
      weights <- c(weights, 1 / summary(lmodbi)$coef[3, 2] ^ 2)
    }
  
  
  # Store Computed Shock-Effects Estimators
  alphas[[1]] <- c(alphas[[1]], mean(alphahatsb))
  alphas[[2]] <- c(alphas[[2]], sum(alphahatsb * weights / sum(weights)))
  alphas[[3]] <- c(alphas[[3]], sum(Wstar * alphahatsb))
}

# Parameters
means <- c()
vars <- c()
for (j in 1:3) {
  means <- c(means, mean(alphas[[j]]))
  vars <- c(vars, var(alphas[[j]]))
}
names(means) <- names(vars) <- c('adj', 'wadj', 'IVW')

# sample version
risk.reduction2 <- function(est, vars) {
  rr.adj <- (est['wadj']) ^ 2 - vars['adj'] - (est['adj'] - est['wadj']) ^ 2
  rr.wadj <- (est['wadj']) ^ 2 - vars['wadj']
  rr.IVW <- (est['wadj']) ^ 2 - vars['IVW'] - (est['IVW'] - est['wadj']) ^ 2
  rest <- c(rr.adj, rr.wadj, rr.IVW)
  names(rest) <- c('adj', 'wadj', 'IVW')
  return(list(usable = ifelse(rest > 0, yes = 1, no = 0), best = which.max(rest), rr = rest))
}
est <- c(alpha_adj, alpha_wadj, alpha_IVW)
names(est) <- c('adj', 'wadj', 'IVW')
risk.reduction2(est = est, vars = vars) 


## Post-shock forecasts
m_AAPL_2020_11 <- lm(Y ~ AAPL_Close + GSPC_Close + TB_Close, data = TS1[-nrow(TS1)])
# Yhat 1
y1 <- coef(m_AAPL_2020_11) %*% t(as.matrix(cbind(1, AAPL_Close = TS1$AAPL_Close[nrow(TS1)], 
                                                    GSPC_Close = X1p[1, 1],
                                                    TB_Close = X1p[1, 2])))
Yhat_nothing <- coef(m_AAPL_2020_11) %*% t(as.matrix(cbind(1, AAPL_Close = y1, 
                                                 GSPC_Close = X1p[2, 1],
                                                 TB_Close = X1p[2, 2])))

Yhat_wadj <- Yhat_nothing + alpha_wadj

y <- AAPL_close[which(AAPL_close$Date == '2020-11-11'), 2]

abs(Yhat_nothing - y)

Yhat_adj <- c(adj = alpha_adj + Yhat_nothing, IVW = alpha_IVW + Yhat_nothing, wadj = alpha_wadj + Yhat_nothing)

abs(Yhat_adj - y)



# plot data
start <- which(AAPL_close$Date == "2020-11-09")
TS1 <- AAPL_close[(start - 15):(start + 2), ]
TS1$id <- 1:nrow(TS1)
mat <- cbind(TS1$id[nrow(TS1)], c(Yhat_adj))
colnames(mat) <- c("id", "Yhat_adj")
dat <- as.data.frame(mat)
colnames(Yhat_nothing) <- "y"
# set working directory
setwd('/Users/mac/Desktop/Research/Post-Shock Prediction/')
# plot setting
tikz('applepsp.tex', standAlone = TRUE, width = 7, height = 5)
# plot
ggplot(TS1, mapping = aes(x = id, y = AAPL_Close)) + 
  labs(title = "Apple Stock Forecasting (2020 October 20th to 2020 November 11th)", 
       x = "Day", y = "Closing Stock price (in USD)") +
  geom_point() + 
  geom_point(data = dat, aes(x = id, y = Yhat_adj), 
             col = c("magenta", "deepskyblue", "indianred2"), 
             pch = 2:4, cex = 1.5) + 
  geom_point(data = data.frame(x = unique(dat$id), y = Yhat_nothing), 
             aes(x = x, y = y), col = "violet", cex = 1.5) + 
  geom_line(aes(x = id, y = c(m_AAPL_2020_11$fitted.values, y1, Yhat_nothing)), 
            col = "violet") + 
  annotate("text", x = 1, y = seq(from = 110, to = 104, length.out = 4), 
           label = c("$\\hat{y}_{1, T_1^* + 2}^{1}$",
                     "$\\hat{y}_{1, T_1^* + 2}^{1} + \\hat{\\alpha}_{\\rm adj}$",
                     "$\\hat{y}_{1, T_1^* + 2}^{1} + \\hat{\\alpha}_{\\rm IVW}$",
                     "$\\hat{y}_{1, T_1^* + 2}^{1} + \\hat{\\alpha}_{\\rm wadj}$"), 
           hjust = 0, size = 5) + 
  annotate("point", x = 0, y = seq(from = 110, to = 104, length.out = 4), 
           pch = c(16, 2:4), 
           color = c("violet", "magenta", "deepskyblue", "indianred2"),
           size = 2) +
  # add margin
  theme(plot.margin = unit(c(.5, .3, .3, .5), "cm")) + 
  # no grid
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  # minimal
  theme_minimal() +
  # center title
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
# output
dev.off()

