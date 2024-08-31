#!/bin/sh

module purge

module load compiler/intel/intel-8458p

source /public/software/mkl/setvars.sh intel64

make
