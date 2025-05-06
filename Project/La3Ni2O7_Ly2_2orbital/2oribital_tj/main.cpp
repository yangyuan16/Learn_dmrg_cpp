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
        int thread = slurm_cpus_per_node;
        omp_set_num_threads(thread);

	//------DMRG parameters
	int Spin = 1;			//Spin=1 for spin-1/2.
	int Total_h = 48;		//Quantum number of hole.
	int Total_J = 0;		//Quantum number of spin, Total_J = 2 * J. 对角化的时候
	int Number_hole = Total_h;		//Number of holes in the system.
	int StateNoKept = 1000;		//DMRG optimal state number for infinite step.

	//------Lattice geometry and size
	int N_u = 2;
	int N_x = 48;
	int N_y = 2;
	int Total_N = N_u * N_x * N_y;

	//------Model parameters
        // layer a, xx orbital,  
	double T_a_xx = 1.0;  // in plane tj
        double J_a_xx = 0.5; // in plane tj	
        double N_a_xx = -0.25 * J_a_xx;
	// layer a, zz orbital,
	double T_a_zz = 0.27;
	double J_a_zz = 0.0;
	double N_a_zz = -0.25 * J_a_zz;
	// layer a, xz hopping,
	double T_a_xz = -0.491;
	double J_a_xz = 0.0;
	double N_a_xz = -0.25 * J_a_xz;
	//---------------------------
	// layer b, xx orbital,
	double T_b_xx = 1.0;
	double J_b_xx = 0.5;
	double N_b_xx = -0.25 * J_b_xx;
	// layer b, zz orbital,
	double T_b_zz = 0.27;
	double J_b_zz = 0.0;
	double N_b_zz = -0.25 * J_b_zz;
	// layer b, xz hopping,
	double T_b_xz = -0.491;
	double J_b_xz = 0.0;
	double N_b_xz = -0.25 * J_b_xz;
	//---------------------------
	// between layer a and b, xx orbital,
	double T_ab_xx = 0.0;
	double J_ab_xx = 0.0;
	double N_ab_xx = -0.25 * J_ab_xx;
	// between layer a and b, zz orbital,
	double T_ab_zz = 0.942;
	double J_ab_zz = 0.444;
	double N_ab_zz = -0.25 * J_ab_zz;
        //
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
        Parameter para( Sign, Spin, Total_h, Total_J, Number_hole, StateNoKept, N_u, N_x, N_y, Total_N, T_a_xx, J_a_xx, N_a_xx,T_a_zz,J_a_zz,N_a_zz,
		       	T_a_xz,J_a_xz,N_a_xz, T_b_xx, J_b_xx, N_b_xx, T_b_zz, J_b_zz, N_b_zz, T_b_xz, J_b_xz, N_b_xz,
		        T_ab_xx,J_ab_xx,N_ab_xx, T_ab_zz, J_ab_zz, N_ab_zz,Desiredzeros, Seedvalue );

	Sweep sweep( para );

	cout << endl;
        cout << "finished!!!" << endl;
        
        return 0;

}
