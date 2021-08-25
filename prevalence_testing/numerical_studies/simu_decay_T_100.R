

# MC
library("parallel")
library("doParallel")
library("foreach")
source("simu_source.R")

registerDoParallel(cores = ncores)
set.seed(13)
RNGkind("L'Ecuyer-CMRG")
nsim <- 200


# parameter setup
#shape.K.Ts <- c(200, 400, 600, 800)
shape.K.Ts <- 100
ns <- c(50, 100, 200, 300)
sim_params <- expand.grid(list(shape.K.Ts = shape.K.Ts, ns = ns))

# simulation time
system.time(
	output_T_100 <- lapply(1:nrow(sim_params), FUN = function(j) {
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
																			 p = 13, B = 1000, scale = 2, 
																			 n = n, H = 8, ell = 4,
																			 Kshape = shape.K.T, Tshape = shape.K.T)
			return(study)
		}
		# return results
		out
	})
)


save(output_T_100, file = "output_decay_T_100.RData")
