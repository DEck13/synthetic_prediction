

# MC
library("parallel")
library("doParallel")
library("foreach")
source("simu_source.R")

registerDoParallel(cores = ncores)
set.seed(14)
RNGkind("L'Ecuyer-CMRG")
nsim <- 150


# parameter setup
#Hs <- c(2, 4, 8, 16)
Hs <- 4
mu.alphas <- c(-1, 0, 1, 2)

sim_params <- expand.grid(list(mu.alphas = mu.alphas, Hs = Hs))

# simulation time
system.time(suppressMessages(
	output_mu_H4 <- lapply(1:nrow(sim_params), FUN = function(j) {
		# parameters
		mu.alpha <- sim_params[j, 1]
		H <- sim_params[j, 2]
		# %do% evaluates sequentially
		# %dopar% evaluates in parallel
		# .combine results
		out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
			# result
			study <- sim.normal.gammaX.decay(mu.gamma.delta = 0, 
																			 mu.alpha = mu.alpha, sigma = 1, 
																			 sigma.alpha = 1, 
																			 sigma.delta.gamma = 1, 
																			 p = 13, B = 2000, scale = 4, 
																			 n = 10, H = H, ell = 3)
			return(study)
		}
		# return results
		out
	})
))


# parameter setup
ells <- c(2, 4, 8, 16)
sim_params <- expand.grid(list(ells = ells, Hs = Hs))

# simulation time
system.time(suppressMessages(
	output_ell_H4 <- lapply(1:nrow(sim_params), FUN = function(j) {
		# parameters
		ell <- sim_params[j, 1]
		H <- sim_params[j, 2]
		# %do% evaluates sequentially
		# %dopar% evaluates in parallel
		# .combine results
		out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
			# result
			study <- sim.normal.gammaX.decay(mu.gamma.delta = 0, 
																			 mu.alpha = 2, sigma = 1, 
																			 sigma.alpha = 1, 
																			 sigma.delta.gamma = 1, 
																			 p = 13, B = 2000, scale = 4, 
																			 n = 10, H = H, ell = ell)
			return(study)
		}
		# return results
		out
	})
))

save(output_mu_H4, output_ell_H4, file = "output_decay_H4.RData")

