#!/bin/bash
#
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name=simu_T300
#SBATCH --partition=stat
#SBATCH --output=myjob.o%j
##SBATCH --error=myjob.e%j
##SBATCH --mail-user=dje13@illinois.edu
##SBATCH --mail-type=BEGIN,END
#
# End of embedded SBATCH options
#	
module load R/4.0.3_sandybridge

export R_LIBS=~/Rlibs

R CMD BATCH --vanilla simu_decay_T_300.R
