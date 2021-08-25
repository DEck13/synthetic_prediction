
rm(list = ls())

## load in output
load("output_decay_T_100.RData")
load("output_decay_T_200.RData")
load("output_decay_T_300.RData")
load("output_decay_T_400.RData")

output <- c(output_T_100, output_T_200, output_T_300, 
						output_T_400)

## create simulation parameters
shape.K.Ts <- c(100, 200, 300, 400)
ns <- c(50, 100, 200, 300)

sim_params <- expand.grid(list(shape.K.Ts = shape.K.Ts, ns = ns))
sim_params <- expand.grid(list(ns = ns, shape.K.Ts = shape.K.Ts))
nsim <- 200

# store results
# load packages
require('readxl')
require('writexl')
write_xlsx(lapply(output, as.data.frame), 'ntrial.xlsx')

result <- c()
for (i in 1:nrow(sim_params)) {
	table <- output[[i]]
	# means and sds
	means <- apply(table, 2, function(x) mean(x))
	sds <- apply(table, 2, function(x) sd(x))
	result.i <- c()
	for (j in 1:7) {
		result.i <- cbind(result.i, paste0(round(means[j], digits = 3), 
																			 ' (', round(sds[j] / sqrt(nsim), 
																			 						digits = 3), ')'))
	}
	result <- rbind(result, result.i)
}
result <- cbind(sim_params[, c(1,2)], result)
require('xtable')
xtable(result, digits = c(0,0,0,3,3,3,3,3,3,3))

