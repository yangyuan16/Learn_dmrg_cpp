#!/bin/bash
#提交单个作业
#SBATCH -J 2orbital     
#SBATCH --partition=CPU-Gong
##SBATCH --partition=CPU-AMD
#SBATCH --nodes=1         
#SBATCH --cpus-per-task=12     
##Suggestion: CPU-INT use 36 cores; CPU-AMD use 32 64 or 96 cores;
#SBATCH --mem=100GB
##SBATCH -n 1      
##SBATCH --ntasks-per-node=16  
##SBATCH --nodelist=XHPC-6150-20
#SBATCH --exclude=XHPC-9654-[01-03] 
##IF partition=CPU-AMD and use 96 cores use --exclude=XHPC-7B13-[01-05]
##SBATCH --time=dd-hh:mm:ss    
##SBATCH --output=file_name   
##SBATCH --error=file_name    


#export PATH=/public/software/mpi/openmpi/3.0.0/intel/bin:$PATH
#export LD_LIBRARY_PATH=/public/software/mpi/openmpi/3.0.0/intel/lib:$LD_LIBRARY_PATH
#export INCLUDE=/public/software/mpi/openmpi/3.0.0/intel/include:$INCLUDE

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/gsl/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/mkl/mkl/2024.2/lib


#For module 
export MODULEPATH=/app/modulefiles
module purge
module load intel gsl-2.8
#echo $PATH > path.log

#source /public/software/mkl/setvars.sh intel64
make clean
make slurm_cpus_per_node=$SLURM_JOB_CPUS_PER_NODE > Make_out.log
#make > Make_out.log

./runDMRG 

#rm -r -f mid

#rm -r -f truncated_density_eigenvector

#rm -r -f truncated_wave_function
