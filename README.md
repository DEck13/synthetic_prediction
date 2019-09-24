# synthetic_prediction

## Ideas for synthetic prediction

We seek to develop a forecasting methodology for time series data that is 
thought to have undergone a shock which has origins that have not been 
previously observed.  We still can provide credible forecasts for a time 
series in the presence of such systematic shocks by drawing from disparate 
time series that have undergone similar shocks for which post-shock 
outcome data is recorded.  These disparate time series are assumed to have 
mechanistic similarities to the time series under study but are otherwise 
independent (Granger noncausal).  The inferential goal of our forecasting 
methodology is to supplement observed time series data with post-shock 
data from the disparate time series in order to minimize average forecast 
risk. 

## Avenues of research (so far)

1) Dynamic panel models with shock effects that are generated from some 
common process.

2) Synthetic control methods with shock effects that are generated from 
some common process and are correlated with observed covariates.

3) Bayesian risk minimization under general random effects model

4) Minimax prediction from the frequentist perspective when shock effects are unknown scalars
