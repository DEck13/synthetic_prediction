

rm(list = ls())

## load in output
load("output_ns5.RData")
load("output_ns10.RData")
load("output_ns15.RData")
load("output_ns25.RData")

output <- c(output_ns5, output_ns10, output_ns15, 
						output_ns25)

## create simulation parameters
ns <- c(5, 10, 15, 25)
sigma.alphas <- c(1, 5, 10, 25, 100)
sim_params <- expand.grid(list(sigma.alphas = sigma.alphas, ns = ns))
nsim <- 100

# store results
# load packages
require('readxl')
require('writexl')
write_xlsx(lapply(output, as.data.frame), 'trial.xlsx')

result <- c()
for (i in 1:nrow(sim_params)) {
	table <- output[[i]]
	# means and sds
	means <- apply(table, 2, function(x) mean(x))
	sds <- apply(table, 2, function(x) sd(x))
	result.i <- c()
	for (j in 1:2) {
		result.i <- cbind(result.i, paste0(round(means[j], digits = 3), 
																			 ' (', round(sds[j] / sqrt(nsim), 
																			 						digits = 3), ')'))
	}
	result <- rbind(result, result.i)
}
result <- cbind(sim_params[, c(2,1)], result)
require('xtable')
xtable(result)
