#!/bin/bash
#SBATCH -p G1Part_vip12
#SBATHC -N 1
#SBATCH -n 1
#SBATCH -c 12

mv measurement_density_correlation.dat measurement_density_correlation_ref44.dat 
mv measurement_spin_correlation.dat measurement_spin_correlation_ref44.dat
mv measurement_green_function.dat measurement_green_function_ref44.dat 
mv measurement_pairing.dat measurement_pairing_zz_ref44.dat
mv measurement_pairing_yy.dat measurement_pairing_yy_ref44.dat
