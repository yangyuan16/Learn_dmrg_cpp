#!/bin/bash
#SBATCH -p G1Part_vip12
#SBATHC -N 1
#SBATCH -n 1
#SBATCH -c 12
source /es01/paratera/parasoft/module.sh
module load mpi/intel/20.0.4
module load gsl/2.4
srun runDMRG
