

# MC
library("parallel")
library("doParallel")
library("foreach")
source("simu_source.R")

#ncores <- detectCores() - 1
registerDoParallel(cores = ncores)
set.seed(15)
RNGkind("L'Ecuyer-CMRG")
nsim <- 100

# parameter setup
ns <- 15
sigma.alphas <- c(1, 5, 10, 25, 100)
sim_params <- expand.grid(list(sigma.alphas = sigma.alphas, ns = ns))



# simulation time
system.time(
	output_ns5 <- lapply(1:nrow(sim_params), FUN = function(j) {
		# parameters
		sigma.alpha <- sim_params[j, 1]
		n <- sim_params[j, 2]
		# %do% evaluates sequentially
		# %dopar% evaluates in parallel
		# .combine results
		out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
			# result
			study <- sim.normal.gammaX(mu.gamma.delta = 2, 
																 mu.alpha = 10, sigma = 1, 
																 sigma.alpha = sigma.alpha, 
																 sigma.delta.gamma = 1, 
																 p = 13, B = 200, scale = 10, 
																 n = n, H = 12, ell = 3)
			return(study)
		}
		# return results
		out
	})
)

save(output_ns15, file = "output_ns15.RData")


