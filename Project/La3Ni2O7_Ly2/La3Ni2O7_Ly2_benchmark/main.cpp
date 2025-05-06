//=====================================================================================
// 	   		       DMRG program for t-J model
//			The algorithm includes the SU(2) symmetry.
//				     SHOU-SHU GONG
//=====================================================================================
#include"omp.h"
#include"mkl.h"
#include<iostream>
#include<time.h>
#include<math.h>
#include<stdio.h>
#include<sys/stat.h>
using namespace std;

#include"parameter.h"
#include"sweep.h"

//------Main Program---------
int main() {

	//------Choose a Job
     	char Sign = 'N'; 

	//------OMP settings
        int thread = 12;
        omp_set_num_threads(thread);

	//------DMRG parameters
	int Spin = 1;			//Spin=1 for spin-1/2.
	int Total_h = 96;		//Quantum number of hole.
	int Total_J = 0;		//Quantum number of spin, Total_J = 2 * J. 对角化的时候
	int Number_hole = Total_h;		//Number of holes in the system.
	int StateNoKept = 1000;		//DMRG optimal state number for infinite step.

	//------Lattice geometry and size
	int N_u = 2;
	int N_x = 48;
	int N_y = 2;
	int Total_N = N_u * N_x * N_y;

	//------Model parameters
	double T_n = 3.0;  // in plane tj
        double J_n = 1.0; // in plane tj	

	double T_nn = 0.0; // between two layers
        double J_nn = 4.0; // between two layers

	double N_n = -0.25 * J_n;
	double N_nn = 0.0 * J_nn;
        
	int Desiredzeros = 0;
	int Seedvalue = 42;

	//------Open files
        mkdir( "./e", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./entanglement", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./entanglement/spectrum", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./entanglement/entropy", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./spin_correlation", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./wavefunction", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./mid", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./mid/space", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./mid/operator", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./mid/operator/H", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		mkdir( "./mid/operator/S_Dia", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		mkdir( "./mid/operator/S_M_Dia", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		mkdir( "./mid/operator/T", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		mkdir( "./mid/operator/NN", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./6j_factor", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./6j_factor/S_Dia_old", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./6j_factor/S_Dia_n", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./6j_factor/S_M_Dia_old", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./6j_factor/S_M_Dia_n", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		mkdir( "./6j_factor/T_old", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		mkdir( "./6j_factor/T_n", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./6j_factor/H", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./new_block", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./new_block_T", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./truncated_wave_function", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        mkdir( "./truncated_density_eigenvector", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

	//------Sweep
        Parameter para( Sign, Spin, Total_h, Total_J, Number_hole, StateNoKept, N_u, N_x, N_y, Total_N, T_n, T_nn, J_n, J_nn, N_n, N_nn, Desiredzeros, Seedvalue );

	Sweep sweep( para );

	cout << endl;
        cout << "finished!!!" << endl;
        
        return 0;

}
