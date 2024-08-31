#!/bin/bash
#SBATCH -J DMRG
#SBATCH -p intel
#SBATHC -N 1
#SBATCH -n 12
#SBATCH --ntasks-per-node=12
#SBATCH -o out.%j
#SBATCH -e err.%j

export PATH=/public/software/mpi/openmpi/3.0.0/intel/bin:$PATH
export LD_LIBRARY_PATH=/public/software/mpi/openmpi/3.0.0/intel/lib:$LD_LIBRARY_PATH
export INCLUDE=/public/software/mpi/openmpi/3.0.0/intel/include:$INCLUDE

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/gsl/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/mkl/mkl/2024.2/lib

module purge

module load compiler/intel/intel-8458p

srun runDMRG

rm -r -f mid

rm -r -f truncated_density_eigenvector

rm -r -f truncated_wave_function
