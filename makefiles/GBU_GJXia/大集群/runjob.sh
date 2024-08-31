#!/bin/bash
#提交单个作业
##SBATCH -J cputest     
##SBATCH --partition=amd
#SBATCH --partition=intel
#SBATCH --nodes=1         
#SBATCH --cpus-per-task=32     
#SBATCH --mem=200GB
##SBATCH -n 1      
##SBATCH --ntasks-per-node=16  
##SBATCH --nodelist=XHPC-6150-20
##SBATCH --exclude=XHPC-9654-[01-03] 
##SBATCH --time=dd-hh:mm:ss    
##SBATCH --output=file_name   
##SBATCH --error=file_name    




#For module 
module purge
module load compiler/intel/intel-8458p apps/gsl/gsl-2.8 
#echo $PATH > path.log

#source /public/software/mkl/setvars.sh intel64
make clean
make slurm_cpus_per_node=$SLURM_JOB_CPUS_PER_NODE > Make_out.log
#make > Make_out.log

./runDMRG > DMRG.log

rm -r -f mid

rm -r -f truncated_density_eigenvector

rm -r -f truncated_wave_function
