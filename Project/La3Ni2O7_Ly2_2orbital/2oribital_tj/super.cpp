#include<iostream>
using namespace std;
#include<math.h>
#include<time.h>
#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<mkl.h>
#include<gsl/gsl_sf_coupling.h>

#include"sub.h"
#include"super.h"
#include"common.h"

#define constant sqrt(6.0)*0.5
#define hoping -sqrt(2.0)

//================================BLAS ROUTINES============================
extern "C" {

void daxpy_( const int *n, const double *alpha, double *x, const int *incx, double *y, const int *incy );

void dsymm_( char *side, char *uplo, const int *m, const int *n, const double *alpha, double *a, const int *lda, double *b, const int *ldb, const double *beta, double *c, const int *ldc );

void dgemm_( char *transa, char *transb, const int *m, const int *n, const int *k, const double *alpha, double *a, const int *lda, double *b, const int *ldb, const double *beta, double *c, const int *ldc );

}

//====Constructe the subspace with given quantum number for Infinite sweep===
Super::Super( Parameter &para, Sub *sys_space1, Sub *env_space1, Sub *sysnew_space1, Sub *envnew_space1 ) {

	//----
        destruct = 's';   //"s" stands for creating space for block

        sys_space = sys_space1;
        env_space = env_space1;
        sysnew_space = sysnew_space1;
        envnew_space = envnew_space1;

	//----
	int total_length = sys_space -> TotSiteNo + env_space -> TotSiteNo + 2 - 1;

	QuantumNumber_hole = para.QuantumNumber_hole[ total_length ];

	QuantumNumber_J = para.QuantumNumber_J[ total_length ];

	cout<<"\n hole number="<<QuantumNumber_hole<<"\t angular momentum number="<<QuantumNumber_J<<endl;

        AllocateBlockNumber(para);

}

//==============Allocate Block Number for TargetBlock=======================
//             Find the Hibert space for the TargetBlock
//  Configuration 1: (sys+ns & env+ne); Configuration 2: (sys+ne & env+ns)
//==========================================================================
inline void Super::AllocateBlockNumber( Parameter &para ) {

	//----Find the total block number of the wavefunction matrix for a given targetspin
	int J_min, J_max, J_num, oldh_sys, oldh_env, oldJ_sys, oldJ_env;

	Dim = 0;  // total dimension of the diagonalized space

	BlockNumber_for_TargetBlock = 0;  // number of the subblocks with different quantum numbers

	//----
	for( int n_1 = 0; n_1 < sysnew_space -> Block_Num_hole; n_1++ )
	for( int n_2 = 0; n_2 < envnew_space -> Block_Num_hole; n_2++ )
	if( sysnew_space -> Num_hole_block[ n_1 ] + envnew_space -> Num_hole_block[ n_2 ] == QuantumNumber_hole ) {

	        for( int j_sysnew = 0; j_sysnew < sysnew_space -> Num_J_hole_block[ n_1 ]; j_sysnew++ )
        	for( int j_envnew = 0; j_envnew < envnew_space -> Num_J_hole_block[ n_2 ]; j_envnew++ ) {

			//----
                	J_min = abs( sysnew_space -> Value_J_block[ n_1 ][ j_sysnew ] - envnew_space -> Value_J_block[ n_2 ][ j_envnew ] );

	                J_max = ( sysnew_space -> Value_J_block[ n_1 ][ j_sysnew ] + envnew_space -> Value_J_block[ n_2 ][ j_envnew ] );

        	        J_num = ( J_max - J_min ) / 2 + 1;

			//----
                	for( int n = 0; n < J_num; n++ )
	                if( QuantumNumber_J == J_min + 2 * n ) {

                                for( int j_n = 0; j_n < 3; j_n++ )
                                for( int j_e = 0; j_e < 3; j_e++ )
                                if( sysnew_space -> Hole_blockOld[ n_1 ][ j_n ][ j_sysnew ] != -1  &&  envnew_space -> Hole_blockOld[ n_2 ][ j_e ][ j_envnew ] != -1  &&  sysnew_space -> J_blockOld[ n_1 ][ j_n ][ j_sysnew ] != -1  &&  envnew_space -> J_blockOld[ n_2 ][ j_e ][ j_envnew ] != -1 ) {

                                        BlockNumber_for_TargetBlock++;

				}

         	       }

        	}

	}

//	cout<<"\n BlockNumber_for_TargetBlock="<<BlockNumber_for_TargetBlock<<endl;

  //------Create Space for sys, env, sysnew, envnew, Dim_block-------------------------------------------
        H_sys = new int [ BlockNumber_for_TargetBlock ];
        H_env = new int [ BlockNumber_for_TargetBlock ];
        H_sysnew = new int [ BlockNumber_for_TargetBlock ];
        H_envnew = new int [ BlockNumber_for_TargetBlock ];

        J_sys = new int [ BlockNumber_for_TargetBlock ];
        J_env = new int [ BlockNumber_for_TargetBlock ];
        J_sysnew = new int [ BlockNumber_for_TargetBlock ];
        J_envnew = new int [ BlockNumber_for_TargetBlock ];

        Dim_block = new int [ BlockNumber_for_TargetBlock ];

        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) {

		H_sys[ i ] = 0;     H_env[ i ] = 0;     
		H_sysnew[ i ] = 0;  H_envnew[ i ] = 0;
                J_sys[ i ] = 0;     J_env[ i ] = 0;     
		J_sysnew[ i ] = 0;  J_envnew[ i ] = 0;  
		Dim_block[ i ] = 0;

        }

  //----Initialize the values of the above arraies
	index = 0;

        for( int n_1 = 0; n_1 < sysnew_space -> Block_Num_hole; n_1++ )
        for( int n_2 = 0; n_2 < envnew_space -> Block_Num_hole; n_2++ )
        if( sysnew_space -> Num_hole_block[ n_1 ] + envnew_space -> Num_hole_block[ n_2 ] == QuantumNumber_hole ) {

		//----
                for( int j_sysnew = 0; j_sysnew < sysnew_space -> Num_J_hole_block[ n_1 ]; j_sysnew++ )
                for( int j_envnew = 0; j_envnew < envnew_space -> Num_J_hole_block[ n_2 ]; j_envnew++ ) {

			//----
                        J_min = abs( sysnew_space -> Value_J_block[ n_1 ][ j_sysnew ] - envnew_space -> Value_J_block[ n_2 ][ j_envnew ] );

                        J_max = ( sysnew_space -> Value_J_block[ n_1 ][ j_sysnew ] + envnew_space -> Value_J_block[ n_2 ][ j_envnew ] );

                        J_num = ( J_max - J_min ) / 2 + 1;

			//----
                        for( int n = 0; n < J_num; n++ )
                        if( QuantumNumber_J == J_min + 2 * n ) {

                                for( int j_n = 0; j_n < 3; j_n++ )
                                for( int j_e = 0; j_e < 3; j_e++ )
                                if( ( oldh_sys = sysnew_space -> Hole_blockOld[ n_1 ][ j_n ][ j_sysnew ] ) != -1  &&  ( oldh_env = envnew_space -> Hole_blockOld[ n_2 ][ j_e ][ j_envnew ] ) != -1  &&  ( oldJ_sys = sysnew_space -> J_blockOld[ n_1 ][ j_n ][ j_sysnew ] ) != -1  &&  ( oldJ_env = envnew_space -> J_blockOld[ n_2 ][ j_e ][ j_envnew ] ) != -1 ) {

					//----
					H_sys[ index ] = oldh_sys;  H_env[ index ] = oldh_env;

					H_sysnew[ index ] = n_1;  H_envnew[ index ] = n_2;

					J_sys[ index ] = oldJ_sys;  J_env[ index ] = oldJ_env;

					J_sysnew[ index ] = j_sysnew;  J_envnew[ index ] = j_envnew;

					Dim_block[ index ] = sys_space -> Dim_J_block[ oldh_sys ][ oldJ_sys ] * env_space -> Dim_J_block[ oldh_env ][ oldJ_env ];

					index++;

                                }

				break;

                        }

                }

        }

        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ )	//find total dimension for wave function

                Dim += Dim_block[ i ];

        cout<<"\n Dimension of the diagonalizing subspace="<<Dim;

  //----Initialize the wavefunction

	//----
        WaveFunction = new double [ Dim ];

        for( int i = 0; i < Dim; i++ )

                WaveFunction[ i ] = 0.0; 

	//----
        WaveFunction_block = new double * [ BlockNumber_for_TargetBlock ];

        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) {

                WaveFunction_block[ i ] = new double [ Dim_block[ i ] ];

                for( int j = 0; j < Dim_block[ i ]; j++ )

                        WaveFunction_block[ i ][ j ] = 0.0;

        }

}

//=================Set the variables for diagonalization process=====================
Super::Super( char &sign, Parameter &para, Sub *sys1, Sub *env1, Sub *sysnew1, Sub *envnew1, Super &space ) {

	//----
        destruct = 'd';   //'d' stands for creating space for diagonalizaiton

        sys = sys1;  
	env = env1;

	sysnew = sysnew1;  
	envnew = envnew1;

  	//Initialize the parameters for BLAS subroutines
	inc = 1;
        N = para.Total_N; 
        side_L = 'L';     side_R = 'R';     uplo = 'U';
        trans_N = 'N';    trans_T = 'T';    
        beta = 1.0;       alpha_p = 1.0;    beta_p = 0.0;

	int N_square = N * N;

	QuantumNumber_hole = space.QuantumNumber_hole;
	QuantumNumber_J = space.QuantumNumber_J;

        //geometry (path) dependent
        if( sign == 'F' ) {  // "StartSite" is the first site of the env block, 'F' is for finite-DMRG

		StartSite = N - 1;

	}

        else if( sign == 'I' ) {  // 'I' is for infinite-DMRG

//		StartSite = 2 * sys -> TotSiteNo + 1;   //For 1D chain

                int totalsite = sys->TotSiteNo + env->TotSiteNo + 2;

                for( int i = 1; i < para.N_x; i++ ) {

                        if( ( i + 1 ) * para.N_u * para.N_y >= totalsite ) {

                                StartSite = ( i + 1 ) * para.N_u * para.N_y - 1;

                                break;

                        }

                }

        }

  	//Allocate Interaction Table
        Table_T = new int [ N_square ];  
	Table_J = new int [ N_square ];  
	Table_N = new int [ N_square ];

        Interaction_T = new double [ N_square ];  
	Interaction_J = new double [ N_square ];  
	Interaction_N = new double [ N_square ];

        for( int i = 0; i < N_square; i++ ) {

                Table_T[i] = para.Table_T[i];  
		Table_J[i] = para.Table_J[i];  
		Table_N[i] = para.Table_N[i];

                Interaction_T[i] = para.Interaction_T[i];
		Interaction_J[i] = para.Interaction_J[i];
		Interaction_N[i] = para.Interaction_N[i];

        }

	//Find the sites with interactions outside the block in system and environment
        operator_number_T_sys = para.Table_T_sys[sys->TotSiteNo-1];
        operator_number_T_env = para.Table_T_env[env->TotSiteNo-1];

        operator_number_J_sys = para.Table_J_sys[sys->TotSiteNo-1];
        operator_number_J_env = para.Table_J_env[env->TotSiteNo-1];

        operator_number_N_sys = para.Table_N_sys[sys->TotSiteNo-1];
        operator_number_N_env = para.Table_N_env[env->TotSiteNo-1];

        Table_T_sys = new int [operator_number_T_sys];
        for( int i=0; i<operator_number_T_sys; i++ )
                Table_T_sys[i]=para.Table_T_sys_site[sys->TotSiteNo-1][i];

	Table_T_env = new int [operator_number_T_env];
        for( int i=0; i<operator_number_T_env; i++ )
                Table_T_env[i]=para.Table_T_env_site[env->TotSiteNo-1][i];

        Table_J_sys = new int [operator_number_J_sys];
        for( int i=0; i<operator_number_J_sys; i++ )
                Table_J_sys[i]=para.Table_J_sys_site[sys->TotSiteNo-1][i];

	Table_J_env = new int [operator_number_J_env];
        for( int i=0; i<operator_number_J_env; i++ )
                Table_J_env[i]=para.Table_J_env_site[env->TotSiteNo-1][i];

        Table_N_sys = new int [operator_number_N_sys];
        for( int i=0; i<operator_number_N_sys; i++ )
                Table_N_sys[i]=para.Table_N_sys_site[sys->TotSiteNo-1][i];

	Table_N_env = new int [operator_number_N_env];
        for( int i=0; i<operator_number_N_env; i++ )
                Table_N_env[i]=para.Table_N_env_site[env->TotSiteNo-1][i];

  //------Read space variables from subroutine space
        BlockNumber_for_TargetBlock = space.BlockNumber_for_TargetBlock;
        Dim = space.Dim;

        H_sys = new int [BlockNumber_for_TargetBlock];
        H_env = new int [BlockNumber_for_TargetBlock];
        H_sysnew = new int [BlockNumber_for_TargetBlock];
        H_envnew = new int [BlockNumber_for_TargetBlock];

        J_sys = new int [BlockNumber_for_TargetBlock];
        J_env = new int [BlockNumber_for_TargetBlock];
        J_sysnew = new int [BlockNumber_for_TargetBlock];
        J_envnew = new int [BlockNumber_for_TargetBlock];
        Dim_block = new int [BlockNumber_for_TargetBlock];

        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) {

		H_sys[i] = space.H_sys[i];        
		H_env[i] = space.H_env[i];     

		H_sysnew[i] = space.H_sysnew[i];
		H_envnew[i] = space.H_envnew[i];  

		J_sys[i] = space.J_sys[i];     
		J_env[i] = space.J_env[i];     
		
		J_sysnew[i] = space.J_sysnew[i];  
		J_envnew[i] = space.J_envnew[i];  

		Dim_block[i] = space.Dim_block[i];

        }

	//Allocate Tables denoting the index between *f and **f for diagonalization
    
	//----create Space
        Table_1to2_Num = new int [Dim];
        Table_1to2_Site = new int [Dim];

        for( int i=0; i<Dim; i++ ) {

                Table_1to2_Num[i] = 0;    

		Table_1to2_Site[i] = 0;

        }

	//----
        Table_2to1 = new int * [ BlockNumber_for_TargetBlock ];

        for( int i=0; i<BlockNumber_for_TargetBlock; i++ ) {

                Table_2to1[i] = new int [ Dim_block[i] ];

                for( int j=0; j<Dim_block[i]; j++ )

                        Table_2to1[i][j] = 0;

        }

	//Initialize the values for the above two Tables      
      
	//----For diagonalization with column-storage vectors
        index = 0;        //"index" denotes the "site" here!!!, "site_s" denotes "site_block"!!!
	site_s = 0;

        for(int i=0; i<BlockNumber_for_TargetBlock; i++)
        for(int j=0; j<env->Dim_J_block[H_env[i]][J_env[i]]; j++)  //These two lines indicate the
        for(int k=0; k<sys->Dim_J_block[H_sys[i]][J_sys[i]]; k++) {//column-main storage!!!

                site_s = j * sys->Dim_J_block[H_sys[i]][J_sys[i]] + k;

                Table_1to2_Num[index] = i;                //index==site

                Table_1to2_Site[index] = site_s;          //index==site

                Table_2to1[i][site_s] = index;            //index==site

                index++;                                //index==site

        }

	//Allocate f1, f2, g1 and g2 for Matrix-Vector multiplication
       
	//----
	f1 = new double * [BlockNumber_for_TargetBlock];
        g1 = new double * [BlockNumber_for_TargetBlock];
        f2 = new double * [BlockNumber_for_TargetBlock];
        g2 = new double * [BlockNumber_for_TargetBlock];

        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) {

                f1[i] = new double [ Dim_block[i] ];        
		g1[i] = new double [ Dim_block[i] ];

                f2[i] = new double [ Dim_block[i] ];        
		g2[i] = new double [ Dim_block[i] ];

                for( int j = 0; j < Dim_block[i]; j++ ) {

                        f1[i][j] = 0.0;           
			g1[i][j] = 0.0;
                        
			f2[i][j] = 0.0;           
			g2[i][j] = 0.0;

                }

        }
 
	//create space for operators

	//----
	T_sys = new double ** [ sys -> Block_Num_hole - 1 ];

        for( int i = 0; i < sys -> Block_Num_hole - 1; i++ ) {

		T_sys[i] = new double * [ sys -> Num_block_T[i] ];

	        for( int j = 0; j < sys -> Num_block_T[i]; j++ ) {

                	dimension = sys -> Dim_J_block[i][ sys -> J_block_T_bra[i][j] ] * sys -> Dim_J_block[ i + 1 ][ sys -> J_block_T_ket[i][j] ];

			T_sys[i][j] = new double [ dimension ];

                        for( int k = 0; k < dimension; k++ )

                        	T_sys[i][j][k] = 0.0;

                }

	}

	//----
	T_env = new double ** [ env -> Block_Num_hole - 1 ];

        for( int i = 0; i < env -> Block_Num_hole - 1; i++ ) {

		T_env[i] = new double * [ env -> Num_block_T[i] ];

	        for( int j = 0; j < env -> Num_block_T[i]; j++ ) {

                	dimension = env -> Dim_J_block[i][ env -> J_block_T_bra[i][j] ] * env -> Dim_J_block[ i + 1 ][ env -> J_block_T_ket[i][j] ];

			T_env[i][j] = new double [ dimension ];

                        for( int k = 0; k < dimension; k++ )

                        	T_env[i][j][k] = 0.0;

                }

	}

	//----
        S_Dia_sys = new double ** [ sys -> Block_Num_hole ];

        for( int i = 0; i < sys -> Block_Num_hole; i++ ) {

        	S_Dia_sys[i] = new double * [ sys -> Num_J_hole_block[i] ];

                for( int j = 0; j < sys -> Num_J_hole_block[i]; j++ ) {

                	dimension = sys -> Dim_J_block[i][j] * sys -> Dim_J_block[i][j];

                        S_Dia_sys[i][j] = new double [ dimension ];

                        for( int k =0; k < dimension; k++ )

                        	S_Dia_sys[i][j][k] = (double) 0;

                }

	}

	//----
        S_M_Dia_sys = new double ** [ sys -> Block_Num_hole ];

        for( int i = 0; i < sys -> Block_Num_hole; i++ ) {

        	S_M_Dia_sys[i] = new double * [ sys -> Num_J_hole_block[i] - 1 ];

                for( int j = 0; j < sys -> Num_J_hole_block[i] - 1; j++ ) {

                	dimension = sys -> Dim_J_block[i][j] * sys -> Dim_J_block[i][j+1];

                        S_M_Dia_sys[i][j] = new double [ dimension ];

                        for( int k =0; k < dimension; k++ )

                        	S_M_Dia_sys[i][j][k] = (double) 0;

                }

	}

	//----
        S_Dia_env = new double ** [ env -> Block_Num_hole ];

        for( int i = 0; i < env -> Block_Num_hole; i++ ) {

        	S_Dia_env[i] = new double * [ env -> Num_J_hole_block[i] ];

                for( int j = 0; j < env -> Num_J_hole_block[i]; j++ ) {

                	dimension = env -> Dim_J_block[i][j] * env -> Dim_J_block[i][j];

                        S_Dia_env[i][j] = new double [ dimension ];

                        for( int k =0; k < dimension; k++ )

	                        S_Dia_env[i][j][k] = (double) 0;

                }

	}

	//----
        S_M_Dia_env = new double ** [ env -> Block_Num_hole ];

        for( int i = 0; i < env -> Block_Num_hole; i++ ) {

        	S_M_Dia_env[i] = new double * [ env -> Num_J_hole_block[i] - 1 ];

                for( int j = 0; j < env -> Num_J_hole_block[i] - 1; j++ ) {

                	dimension = env -> Dim_J_block[i][j] * env -> Dim_J_block[i][j+1];

                        S_M_Dia_env[i][j] = new double [ dimension ];

                        for( int k =0; k < dimension; k++ )

                        	S_M_Dia_env[i][j][k] = (double) 0;

                }

	}

	//----
        NN_sys = new double ** [ sys -> Block_Num_hole ];

        for( int i = 0; i < sys -> Block_Num_hole; i++ ) {

        	NN_sys[i] = new double * [ sys -> Num_J_hole_block[i] ];

                for( int j = 0; j < sys -> Num_J_hole_block[i]; j++ ) {

                	dimension = sys -> Dim_J_block[i][j] * sys -> Dim_J_block[i][j];

                        NN_sys[i][j] = new double [ dimension ];

                        for( int k =0; k < dimension; k++ )

                        	NN_sys[i][j][k] = (double) 0;

                }

	}

	//----
        NN_env = new double ** [ env -> Block_Num_hole ];

        for( int i = 0; i < env -> Block_Num_hole; i++ ) {

        	NN_env[i] = new double * [ env -> Num_J_hole_block[i] ];

                for( int j = 0; j < env -> Num_J_hole_block[i]; j++ ) {

                	dimension = env -> Dim_J_block[i][j] * env -> Dim_J_block[i][j];

                        NN_env[i][j] = new double [ dimension ];

                        for( int k =0; k < dimension; k++ )

                        	NN_env[i][j][k] = (double) 0;

                }

	}

	//6j and 9j coefficients for configurations 1, and 2
	Allocate_6j_config_1_2( space );

	Allocate_9j_config_2( para, space );

	//Allocate the block information for the wave function of configure 3
	AllocateBlockNumber_config_3( para );

        Allocate_6j_config_3();

        Allocate_9j_config_3( para, space );

}

//=====================================================================================
//                    Super for the calculations of measurements
//=====================================================================================
Super::Super( Sub *sys_space1, Sub *env_space1, Sub *sysnew_space1, Sub *envnew_space1, Parameter &para ) {

        destruct = 'm';   

        sys_space = sys_space1;
        env_space = env_space1;
        sysnew_space = sysnew_space1;
        envnew_space = envnew_space1;

	QuantumNumber_hole = para.Total_h;
	QuantumNumber_J = para.Total_J; 

  //----create space for quantum numbers
	int J_min, J_max, J_num, oldh_sys, oldh_env, oldJ_sys, oldJ_env;

        BlockNumber_for_TargetBlock = 0;   Dim = 0;

	for( int n_1=0; n_1<sysnew_space->Block_Num_hole; n_1++ )
	for( int n_2=0; n_2<envnew_space->Block_Num_hole; n_2++ )
	if( sysnew_space->Num_hole_block[n_1]+envnew_space->Num_hole_block[n_2]== QuantumNumber_hole ) {

	        for( int j_sysnew=0; j_sysnew<sysnew_space->Num_J_hole_block[n_1]; j_sysnew++ )
        	for( int j_envnew=0; j_envnew<envnew_space->Num_J_hole_block[n_2]; j_envnew++ ) {

                	J_min = abs(sysnew_space->Value_J_block[n_1][j_sysnew]-envnew_space->Value_J_block[n_2][j_envnew]);
	                J_max = (sysnew_space->Value_J_block[n_1][j_sysnew]+envnew_space->Value_J_block[n_2][j_envnew]);
        	        J_num = (J_max-J_min)/2+1;

                	for( int n=0; n<J_num; n++ )
	                if( QuantumNumber_J ==(J_min+2*n) ) {

                                for( int j_n=0; j_n<3; j_n++ )
                                for( int j_e=0; j_e<3; j_e++ )
                                if( sysnew_space->Hole_blockOld[n_1][j_n][j_sysnew]!=-1 && envnew_space->Hole_blockOld[n_2][j_e][j_envnew]!=-1 && sysnew_space->J_blockOld[n_1][j_n][j_sysnew]!=-1 && envnew_space->J_blockOld[n_2][j_e][j_envnew]!=-1 ) {
                                        BlockNumber_for_TargetBlock++;
                                }

         	       }

        	}

	}

  //----Create Space for sys, env, sysnew, envnew, Dim_block-------------------------------------------
        H_sys = new int [BlockNumber_for_TargetBlock];
        H_env = new int [BlockNumber_for_TargetBlock];
        H_sysnew = new int [BlockNumber_for_TargetBlock];
        H_envnew = new int [BlockNumber_for_TargetBlock];

        J_sys = new int [BlockNumber_for_TargetBlock];
        J_env = new int [BlockNumber_for_TargetBlock];
        J_sysnew = new int [BlockNumber_for_TargetBlock];
        J_envnew = new int [BlockNumber_for_TargetBlock];

        Dim_block = new int [BlockNumber_for_TargetBlock];

        for( int i=0; i<BlockNumber_for_TargetBlock; i++ ) {

		H_sys[i]=0;     H_env[i]=0;     H_sysnew[i]=0;  H_envnew[i]=0;
                J_sys[i]=0;     J_env[i]=0;     J_sysnew[i]=0;  J_envnew[i]=0;  Dim_block[i]=0;

        }

  //----Initialize the values of the above arraies-------------------------------------------------------------
	index = 0;

        for( int n_1=0; n_1<sysnew_space->Block_Num_hole; n_1++ )
        for( int n_2=0; n_2<envnew_space->Block_Num_hole; n_2++ )
        if( sysnew_space->Num_hole_block[n_1]+envnew_space->Num_hole_block[n_2]== QuantumNumber_hole ) {

                for( int j_sysnew=0; j_sysnew<sysnew_space->Num_J_hole_block[n_1]; j_sysnew++ )
                for( int j_envnew=0; j_envnew<envnew_space->Num_J_hole_block[n_2]; j_envnew++ ) {

                        J_min = abs(sysnew_space->Value_J_block[n_1][j_sysnew]-envnew_space->Value_J_block[n_2][j_envnew]);
                        J_max = (sysnew_space->Value_J_block[n_1][j_sysnew]+envnew_space->Value_J_block[n_2][j_envnew]);
                        J_num = (J_max-J_min)/2+1;

                        for( int n=0; n<J_num; n++ )
                        if( QuantumNumber_J ==(J_min+2*n) ) {

                                for( int j_n=0; j_n<3; j_n++ )
                                for( int j_e=0; j_e<3; j_e++ )
                                if( (oldh_sys=sysnew_space->Hole_blockOld[n_1][j_n][j_sysnew])!=-1 && (oldh_env=envnew_space->Hole_blockOld[n_2][j_e][j_envnew])!=-1 && (oldJ_sys=sysnew_space->J_blockOld[n_1][j_n][j_sysnew])!=-1 && (oldJ_env=envnew_space->J_blockOld[n_2][j_e][j_envnew])!=-1 ) {

					H_sys[index] = oldh_sys;  H_env[index] = oldh_env;
					H_sysnew[index] = n_1;  H_envnew[index] = n_2;
					J_sys[index] = oldJ_sys;  J_env[index] = oldJ_env;
					J_sysnew[index] = j_sysnew;  J_envnew[index] = j_envnew;
					Dim_block[index++] = sys_space->Dim_J_block[oldh_sys][oldJ_sys] * env_space->Dim_J_block[oldh_env][oldJ_env];

                                }

				break;

                        }

                }

        }

	//---- for the configuration 3
        BlockNumber_for_TargetBlock_config_3 = 0;

	for(int n_1=0; n_1<sys_space->Block_Num_hole; n_1++)
	for(int n_2=0; n_2<env_space->Block_Num_hole; n_2++)
	for(int n_n=0; n_n<2; n_n++)
	for(int n_e=0; n_e<2; n_e++)
	if( sys_space->Num_hole_block[n_1] + env_space->Num_hole_block[n_2] + n_n + n_e == QuantumNumber_hole ) {

	        for(int jsys=0; jsys<sys_space->Num_J_hole_block[n_1]; jsys++)
        	for(int jenv=0; jenv<env_space->Num_J_hole_block[n_2]; jenv++) {

                	J_min=abs(sys_space->Value_J_block[n_1][jsys] - env_space->Value_J_block[n_2][jenv]);
                	J_max=(sys_space->Value_J_block[n_1][jsys] + env_space->Value_J_block[n_2][jenv]);
                	J_num=(J_max-J_min)/2+1;

                	for(int n=0; n<J_num; n++) {

                        	int J_mid=(J_min+2*n);

				int J_ne_min=abs(n_n-n_e);
				int J_ne_max=(2-n_n-n_e);
				int J_ne_num=(J_ne_max-J_ne_min)/2+1;
		
				for(int m=0; m<J_ne_num; m++) {

					int J_ne=(J_ne_min+2*m);

					int J_total_min=abs(J_mid-J_ne);
					int J_total_max=(J_mid+J_ne);
					int J_total_num=(J_total_max-J_total_min)/2+1;
				
					for(int l=0; l<J_total_num; l++) {
						int J_total=(J_total_min+2*l);
						if(J_total== QuantumNumber_J ) 
							BlockNumber_for_TargetBlock_config_3++;
					}

                        	}

                	}

        	}

	}

  //----Create Space for sys_3, env_3, sysnew_3, envnew_3, Dim_block_config_3
        H_sys_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        H_env_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        H_ns_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        H_ne_config_3=new int [BlockNumber_for_TargetBlock_config_3];

        J_sys_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        J_env_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        J_sysnew_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        J_envnew_config_3=new int [BlockNumber_for_TargetBlock_config_3];

        Dim_block_config_3=new int [BlockNumber_for_TargetBlock_config_3];

        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) {

		H_sys_config_3[i]=0;     H_env_config_3[i]=0;     
		H_ns_config_3[i]=0;	H_ne_config_3[i]=0;
                J_sys_config_3[i]=0;     J_env_config_3[i]=0;     
		J_sysnew_config_3[i]=0;  J_envnew_config_3[i]=0;  
		Dim_block_config_3[i]=0;

        }

  //----Initialize the quantum number arraies
	index=0;

	for(int n_1=0; n_1<sys_space->Block_Num_hole; n_1++)
	for(int n_2=0; n_2<env_space->Block_Num_hole; n_2++)
	for(int n_n=0; n_n<2; n_n++)
	for(int n_e=0; n_e<2; n_e++)
	if(sys_space->Num_hole_block[n_1] + env_space->Num_hole_block[n_2]+n_n+n_e== QuantumNumber_hole ) {

	        for(int jsys=0; jsys<sys_space->Num_J_hole_block[n_1]; jsys++)
        	for(int jenv=0; jenv<env_space->Num_J_hole_block[n_2]; jenv++) {

                	J_min=abs(sys_space->Value_J_block[n_1][jsys] - env_space->Value_J_block[n_2][jenv]);
                	J_max=(sys_space->Value_J_block[n_1][jsys] + env_space->Value_J_block[n_2][jenv]);
                	J_num=(J_max-J_min)/2+1;

                	for(int n=0; n<J_num; n++) {

                        	int J_mid=(J_min+2*n);

				int J_ne_min=abs(n_n-n_e);
				int J_ne_max=(2-n_n-n_e);
				int J_ne_num=(J_ne_max-J_ne_min)/2+1;
		
				for(int m=0; m<J_ne_num; m++) {

					int J_ne=(J_ne_min+2*m);

					int J_total_min=abs(J_mid-J_ne);
					int J_total_max=(J_mid+J_ne);
					int J_total_num=(J_total_max-J_total_min)/2+1;
				
					for(int l=0; l<J_total_num; l++) {
						int J_total=(J_total_min+2*l);
						if(J_total== QuantumNumber_J ) {
							H_sys_config_3[index]=n_1;	H_env_config_3[index]=n_2;
							H_ns_config_3[index]=n_n;	H_ne_config_3[index]=n_e;
							J_sys_config_3[index]=jsys;	J_env_config_3[index]=jenv;
							J_sysnew_config_3[index]=J_mid;
							J_envnew_config_3[index]=J_ne;
							Dim_block_config_3[index++]=sys_space->Dim_J_block[n_1][jsys]*env_space->Dim_J_block[n_2][jenv];
						}
					}

                        	}

                	}

        	}

	}

	//---- 9-j coefficient
        nine_j_config_3=new double * [BlockNumber_for_TargetBlock_config_3];

        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) {

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock; j++)
		if(H_sys[j]==H_sys_config_3[i] && H_env[j]==H_env_config_3[i] && J_sys[j]==J_sys_config_3[i] && J_env[j]==J_env_config_3[i] && H_ns_config_3[i]==(sysnew_space->Num_hole_block[H_sysnew[j]] - sys_space->Num_hole_block[H_sys[j]]) && H_ne_config_3[i]==(envnew_space->Num_hole_block[H_envnew[j]]-env_space->Num_hole_block[H_env[j]]))
			index++;

                nine_j_config_3[i]=new double [index];

                for(int n=0; n<index; n++)
                        nine_j_config_3[i][n]=(double) 0;

        }

        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) {

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock; j++)
		if(H_sys[j]==H_sys_config_3[i] && H_env[j]==H_env_config_3[i] && J_sys[j]==J_sys_config_3[i] && J_env[j]==J_env_config_3[i] && H_ns_config_3[i]==(sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]]) && H_ne_config_3[i]==(envnew_space->Num_hole_block[H_envnew[j]]-env_space->Num_hole_block[H_env[j]])) {

			alpha=1.0;

			if(H_ns_config_3[i]==0 && (env_space->TotSiteNo - env_space->Num_hole_block[H_env[j]])%2==1)	
				alpha=-1.0;
		
			nine_j_config_3[i][index++] = alpha * gsl_sf_coupling_9j( sys_space->Value_J_block[H_sys[j]][J_sys[j]], env_space->Value_J_block[H_env[j]][J_env[j]], J_sysnew_config_3[i], 1-H_ns_config_3[i], 1-H_ne_config_3[i], J_envnew_config_3[i], sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]], envnew_space->Value_J_block[H_envnew[j]][J_envnew[j]], QuantumNumber_J) * sqrt((J_sysnew_config_3[i]+1.0) * (J_envnew_config_3[i]+1.0) * (sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]]+1.0) * (envnew_space->Value_J_block[H_envnew[j]][J_envnew[j]]+1.0));
                        
		}

        }

	//9j for config-2
        nine_j_config_2=new double * [BlockNumber_for_TargetBlock];

        for(int i=0; i<BlockNumber_for_TargetBlock; i++) {      //"i" denote the basis of configuration II

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock; j++)
                if(H_sys[i]==H_sys[j] && H_env[i]==H_env[j] && J_sys[i]==J_sys[j] && J_env[i]==J_env[j] && (sysnew_space->Num_hole_block[H_sysnew[i]]-sys_space->Num_hole_block[H_sys[i]])==(envnew_space->Num_hole_block[H_envnew[j]]-env_space->Num_hole_block[H_env[j]]) && (sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]])==(envnew_space->Num_hole_block[H_envnew[i]]-env_space->Num_hole_block[H_env[i]]))
                        index++;

                nine_j_config_2[i]=new double [index];

                for(int n=0; n<index; n++) {
                        nine_j_config_2[i][n]=0.0;

                }

        }

        for(int i=0; i<BlockNumber_for_TargetBlock; i++) {

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock; j++)
                if(H_sys[i]==H_sys[j] && H_env[i]==H_env[j] && J_sys[i]==J_sys[j] && J_env[i]==J_env[j] && (sysnew_space->Num_hole_block[H_sysnew[i]]-sys_space->Num_hole_block[H_sys[i]])==(envnew_space->Num_hole_block[H_envnew[j]]-env_space->Num_hole_block[H_env[j]]) && (sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]])==(envnew_space->Num_hole_block[H_envnew[i]]-env_space->Num_hole_block[H_env[i]])) {

                        alpha=sqrt((sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]]+1.0)*(envnew_space->Value_J_block[H_envnew[i]][J_envnew[i]]+1.0)*(sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]]+1.0)*(envnew_space->Value_J_block[H_envnew[j]][J_envnew[j]]+1.0))*(QuantumNumber_J+1.0)*pow(-1.0,(-sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]]+envnew_space->Value_J_block[H_envnew[i]][J_envnew[i]]-sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]]+envnew_space->Value_J_block[H_envnew[j]][J_envnew[j]]-2*(sys_space->Value_J_block[H_sys[i]][J_sys[i]]+env_space->Value_J_block[H_env[i]][J_env[i]])+2*(2-sysnew_space->Num_hole_block[H_sysnew[i]]+sys_space->Num_hole_block[H_sys[i]]-sysnew_space->Num_hole_block[H_sysnew[j]]+sys_space->Num_hole_block[H_sys[j]])-4*QuantumNumber_J)/2);

                        if(((sysnew_space->Num_hole_block[H_sysnew[i]]-sys_space->Num_hole_block[H_sys[i]])==0 && (sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]])==0) || ((sysnew_space->Num_hole_block[H_sysnew[i]]-sys_space->Num_hole_block[H_sys[i]])==0 && (sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]])==1 && (env_space->TotSiteNo-env_space->Num_hole_block[H_env[i]])%2==1) || ((sysnew_space->Num_hole_block[H_sysnew[i]]-sys_space->Num_hole_block[H_sys[i]])==1 && (sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]])==0 && (env_space->TotSiteNo-env_space->Num_hole_block[H_env[i]])%2==1))
                               alpha=-alpha;

                        for(int m_n=0; m_n<(2-sysnew_space->Num_hole_block[H_sysnew[j]]+sys_space->Num_hole_block[H_sys[j]]); m_n++)
                        for(int m_e=0; m_e<(2-sysnew_space->Num_hole_block[H_sysnew[i]]+sys_space->Num_hole_block[H_sys[i]]); m_e++)
                        for(int m_1=0; m_1<sys_space->Value_J_block[H_sys[i]][J_sys[i]]+1; m_1++)
                        for(int m_2=0; m_2<env_space->Value_J_block[H_env[i]][J_env[i]]+1; m_2++)
                        if((sysnew_space->Num_hole_block[H_sysnew[i]]-sys_space->Num_hole_block[H_sys[i]]+sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]]+2*(m_n+m_e+m_1+m_2-1)-sys_space->Value_J_block[H_sys[i]][J_sys[i]]-env_space->Value_J_block[H_env[i]][J_env[i]])==QuantumNumber_J)

		nine_j_config_2[i][index] += alpha * gsl_sf_coupling_3j(sys_space->Value_J_block[H_sys[i]][J_sys[i]], 1-sysnew_space->Num_hole_block[H_sysnew[i]]+sys_space->Num_hole_block[H_sys[i]], sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]], -sys_space->Value_J_block[H_sys[i]][J_sys[i]]+2*m_1, sysnew_space->Num_hole_block[H_sysnew[i]]-sys_space->Num_hole_block[H_sys[i]]-1+2*m_e, -sysnew_space->Num_hole_block[H_sysnew[i]]+sys_space->Num_hole_block[H_sys[i]]+1-2*(m_e+m_1)+sys_space->Value_J_block[H_sys[j]][J_sys[j]])*gsl_sf_coupling_3j(env_space->Value_J_block[H_env[j]][J_env[j]], 1-envnew_space->Num_hole_block[H_envnew[i]]+env_space->Num_hole_block[H_env[i]], envnew_space->Value_J_block[H_envnew[i]][J_envnew[i]], -env_space->Value_J_block[H_env[j]][J_env[j]]+2*m_2, envnew_space->Num_hole_block[H_envnew[i]]-env_space->Num_hole_block[H_env[i]]-1+2*m_n, -envnew_space->Num_hole_block[H_envnew[i]]+env_space->Num_hole_block[H_env[i]]+1-2*(m_n+m_2)+env_space->Value_J_block[H_env[j]][J_env[j]])*gsl_sf_coupling_3j(sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]], envnew_space->Value_J_block[H_envnew[j]][J_envnew[j]], QuantumNumber_J, sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]]-1+2*(m_n+m_1)-sys_space->Value_J_block[H_sys[j]][J_sys[j]], envnew_space->Num_hole_block[H_envnew[j]]-env_space->Num_hole_block[H_env[j]]-1+2*(m_e+m_2)-env_space->Value_J_block[H_env[j]][J_env[j]], -QuantumNumber_J)*gsl_sf_coupling_3j(sys_space->Value_J_block[H_sys[j]][J_sys[j]], 1-sysnew_space->Num_hole_block[H_sysnew[j]]+sys_space->Num_hole_block[H_sys[j]], sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]], -sys_space->Value_J_block[H_sys[j]][J_sys[j]]+2*m_1, sysnew_space->Num_hole_block[H_sysnew[j]]-sys_space->Num_hole_block[H_sys[j]]-1+2*m_n, -sysnew_space->Num_hole_block[H_sysnew[j]]+sys_space->Num_hole_block[H_sys[j]]+1-2*(m_n+m_1)+sys_space->Value_J_block[H_sys[j]][J_sys[j]])*gsl_sf_coupling_3j(env_space->Value_J_block[H_env[j]][J_env[j]], 1-envnew_space->Num_hole_block[H_envnew[j]]+env_space->Num_hole_block[H_env[j]], envnew_space->Value_J_block[H_envnew[j]][J_envnew[j]], -env_space->Value_J_block[H_env[j]][J_env[j]]+2*m_2, envnew_space->Num_hole_block[H_envnew[j]]-env_space->Num_hole_block[H_env[j]]-1+2*m_e, -envnew_space->Num_hole_block[H_envnew[j]]+env_space->Num_hole_block[H_env[j]]+1-2*(m_e+m_2)+env_space->Value_J_block[H_env[j]][J_env[j]])*gsl_sf_coupling_3j(sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]], envnew_space->Value_J_block[H_envnew[i]][J_envnew[i]], QuantumNumber_J, sysnew_space->Num_hole_block[H_sysnew[i]]-sys_space->Num_hole_block[H_sys[i]]-1+2*(m_e+m_1)-sys_space->Value_J_block[H_sys[i]][J_sys[i]], envnew_space->Num_hole_block[H_envnew[i]]-env_space->Num_hole_block[H_env[i]]-1+2*(m_n+m_2)-env_space->Value_J_block[H_env[i]][J_env[i]], -QuantumNumber_J);

			index++;

		}

	}

	//---- 6j coefficient
        six_j_J_config_3 = new double * [ BlockNumber_for_TargetBlock_config_3 ];

        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i ++ )
	if( sys_space -> Num_hole_block[ H_sys_config_3[ i ] ] != sys_space -> TotSiteNo  &&  env_space -> Num_hole_block[ H_env_config_3[ i ] ] != env_space -> TotSiteNo ) {

		index = 0;

		for(int j=0; j<BlockNumber_for_TargetBlock_config_3; j++)
		if( H_sys_config_3[i] == H_sys_config_3[j]  &&  H_env_config_3[i] == H_env_config_3[j]  &&  H_ns_config_3[j] == H_ns_config_3[i]  &&  H_ne_config_3[j] == H_ne_config_3[i]  &&  J_sysnew_config_3[i] == J_sysnew_config_3[j]  &&  J_envnew_config_3[i] == J_envnew_config_3[j]  &&  abs( sys_space -> Value_J_block[ H_sys_config_3[i] ][ J_sys_config_3[i] ] - sys_space -> Value_J_block[ H_sys_config_3[j] ][ J_sys_config_3[j] ] ) <= 2  &&  abs( env_space -> Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ] - env_space -> Value_J_block[ H_env_config_3[j] ][ J_env_config_3[j] ] ) <= 2 )

				index++;

                six_j_J_config_3[i] = new double [index];

                for(int n=0; n<index; n++) {

                        six_j_J_config_3[i][n]=(double) 0;

		}

	}

  //------Initialize six_j_J_config_3------------------------------------------
        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i ++) 
	if( sys_space -> Num_hole_block[ H_sys_config_3[ i ] ] != sys_space -> TotSiteNo  &&  env_space -> Num_hole_block[ H_env_config_3[ i ] ] != env_space -> TotSiteNo ) {

		index = 0;

		for(int j=0; j<BlockNumber_for_TargetBlock_config_3; j++)
		if( H_sys_config_3[ i ] == H_sys_config_3[ j ]  &&  H_env_config_3[ i ] == H_env_config_3[ j ]  &&  H_ns_config_3[ j ] == H_ns_config_3[ i ]  &&  H_ne_config_3[ j ] == H_ne_config_3[ i ]  &&  J_sysnew_config_3[ i ] == J_sysnew_config_3[ j ]  &&  J_envnew_config_3[ i ] == J_envnew_config_3[ j ]  &&  abs( sys_space -> Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] - sys_space -> Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] ) <= 2  &&  abs( env_space -> Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] - env_space -> Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ] ) <= 2 ) {

			six_j_J_config_3[ i ][ index++ ] = pow( -1.0, ( sys_space->Value_J_block[ H_sys_config_3[j] ][ J_sys_config_3[j] ] + env_space->Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ] + J_sysnew_config_3[i] ) / 2 ) * gsl_sf_coupling_6j( J_sysnew_config_3[i], env_space->Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ], sys_space->Value_J_block[ H_sys_config_3[i] ][ J_sys_config_3[i] ], 2, sys_space->Value_J_block[ H_sys_config_3[j] ][ J_sys_config_3[j] ], env_space->Value_J_block[ H_env_config_3[j] ][ J_env_config_3[j] ] );

		}

	}

	//---six_j_1_sys_n and six_j_1_e_env
	int i_g, i_f;
	int Length = BlockNumber_for_TargetBlock*BlockNumber_for_TargetBlock;

	six_j_1_sys_n = new double [Length];
	six_j_1_e_env = new double [Length];

	for(int i=0; i<Length; i++) {
		six_j_1_sys_n[i]=0.0;  six_j_1_e_env[i]=0.0;
	}

        for(int i=0; i<Length; i++) {

                i_g=i/BlockNumber_for_TargetBlock;//stands for the index of "g" function, "g" function is bra
                i_f=i%BlockNumber_for_TargetBlock;//stands for the index of "f" function  "f" function is ket

		if((sys_space->TotSiteNo - sys_space->Num_hole_block[H_sys[i_g]])>0 && sys_space->Num_hole_block[H_sys[i_g]] == sysnew_space->Num_hole_block[H_sysnew[i_g]] && H_sysnew[i_f]==H_sysnew[i_g] && H_envnew[i_f]==H_envnew[i_g] && H_env[i_f]==H_env[i_g] && H_sys[i_f]==H_sys[i_g] && J_sysnew[i_f]==J_sysnew[i_g] && J_envnew[i_f]==J_envnew[i_g] && J_env[i_f]==J_env[i_g] && (sys_space->Value_J_block[H_sys[i_f]][J_sys[i_f]] - sys_space->Value_J_block[H_sys[i_g]][J_sys[i_g]])<=2) {

                	six_j_1_sys_n[i]=gsl_sf_coupling_6j(sysnew_space->Value_J_block[H_sysnew[i_f]][J_sysnew[i_f]], 1, sys_space->Value_J_block[H_sys[i_g]][J_sys[i_g]], 2, sys_space->Value_J_block[H_sys[i_f]][J_sys[i_f]], 1) * pow(-1.0, (1+sysnew_space->Value_J_block[H_sysnew[i_f]][J_sysnew[i_f]]+sys_space->Value_J_block[H_sys[i_f]][J_sys[i_f]])/2)*constant; 

		}
	
		if((env_space->TotSiteNo-env_space->Num_hole_block[H_env[i_g]])>0 && env_space->Num_hole_block[H_env[i_g]]==envnew_space->Num_hole_block[H_envnew[i_g]] && H_sysnew[i_f]==H_sysnew[i_g] && H_envnew[i_f]==H_envnew[i_g] && H_env[i_f]==H_env[i_g] && H_sys[i_f]==H_sys[i_g] && J_sysnew[i_f]==J_sysnew[i_g] && J_envnew[i_f]==J_envnew[i_g] && J_sys[i_f]==J_sys[i_g] && (env_space->Value_J_block[H_env[i_f]][J_env[i_f]]-env_space->Value_J_block[H_env[i_g]][J_env[i_g]])<=2) {

			six_j_1_e_env[i]=gsl_sf_coupling_6j(envnew_space->Value_J_block[H_envnew[i_f]][J_envnew[i_f]], 1, env_space->Value_J_block[H_env[i_g]][J_env[i_g]], 2, env_space->Value_J_block[H_env[i_f]][J_env[i_f]], 1)*pow(-1.0, (1+envnew_space->Value_J_block[H_envnew[i_f]][J_envnew[i_f]]+env_space->Value_J_block[H_env[i_f]][J_env[i_f]])/2)*constant; 

		}

	}

	//------Create and Initialize six_j_2_config_3
        six_j_2_config_3 = new double [ BlockNumber_for_TargetBlock_config_3 ];

        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ ) 

		six_j_2_config_3[ i ] = 0.0;

        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ ) 
	if( H_ns_config_3[ i ] == 0  &&  H_ne_config_3[ i ] == 0 ) {

                if( J_envnew_config_3[ i ] == 0 )      

			six_j_2_config_3[ i ] = -0.75;

                else                            

			six_j_2_config_3[ i ] = 0.25;

        }

}

//=============Allocate block information of wavefunction for configuration 3===============
inline void Super::AllocateBlockNumber_config_3(Parameter &para) {
	
	//Find the total block number of the wavefunction matrix for a given targetspin
        BlockNumber_for_TargetBlock_config_3 = 0;

	for(int n_1=0; n_1<sys->Block_Num_hole; n_1++)
	for(int n_2=0; n_2<env->Block_Num_hole; n_2++)
	for(int n_n=0; n_n<2; n_n++)
	for(int n_e=0; n_e<2; n_e++)
	if( sys->Num_hole_block[n_1] + env->Num_hole_block[n_2] + n_n + n_e == QuantumNumber_hole ) {

	        for(int jsys=0; jsys<sys->Num_J_hole_block[n_1]; jsys++)
        	for(int jenv=0; jenv<env->Num_J_hole_block[n_2]; jenv++) {

                	int J_min=abs(sys->Value_J_block[n_1][jsys]-env->Value_J_block[n_2][jenv]);
                	int J_max=(sys->Value_J_block[n_1][jsys]+env->Value_J_block[n_2][jenv]);
                	int J_num=(J_max-J_min)/2+1;

                	for(int n=0; n<J_num; n++) {

                        	int J_mid=(J_min+2*n);

				int J_ne_min=abs(n_n-n_e);
				int J_ne_max=(2-n_n-n_e);
				int J_ne_num=(J_ne_max-J_ne_min)/2+1;
		
				for(int m=0; m<J_ne_num; m++) {

					int J_ne=(J_ne_min+2*m);

					int J_total_min=abs(J_mid-J_ne);
					int J_total_max=(J_mid+J_ne);
					int J_total_num=(J_total_max-J_total_min)/2+1;
				
					for(int l=0; l<J_total_num; l++) {
						int J_total=(J_total_min+2*l);
						if(J_total== QuantumNumber_J ) 
							BlockNumber_for_TargetBlock_config_3++;
					}

                        	}

                	}

        	}

	}
//	cout<<"\n BlockNumber_config_3="<<BlockNumber_for_TargetBlock_config_3<<endl;

  //------Create Space for sys_3, env_3, sysnew_3, envnew_3, Dim_block_config_3-------------------------------------------
        H_sys_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        H_env_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        H_ns_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        H_ne_config_3=new int [BlockNumber_for_TargetBlock_config_3];

        J_sys_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        J_env_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        J_sysnew_config_3=new int [BlockNumber_for_TargetBlock_config_3];
        J_envnew_config_3=new int [BlockNumber_for_TargetBlock_config_3];

        Dim_block_config_3=new int [BlockNumber_for_TargetBlock_config_3];

        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) {
		H_sys_config_3[i]=0;     H_env_config_3[i]=0;     H_ns_config_3[i]=0;  H_ne_config_3[i]=0;
                J_sys_config_3[i]=0;     J_env_config_3[i]=0;     J_sysnew_config_3[i]=0;  J_envnew_config_3[i]=0;  Dim_block_config_3[i]=0;
        }

  //------Initialize the quantum number arraies
	index=0;

	for(int n_1=0; n_1<sys->Block_Num_hole; n_1++)
	for(int n_2=0; n_2<env->Block_Num_hole; n_2++)
	for(int n_n=0; n_n<2; n_n++)
	for(int n_e=0; n_e<2; n_e++)
	if(sys->Num_hole_block[n_1]+env->Num_hole_block[n_2]+n_n+n_e== QuantumNumber_hole ) {

	        for(int jsys=0; jsys<sys->Num_J_hole_block[n_1]; jsys++)
        	for(int jenv=0; jenv<env->Num_J_hole_block[n_2]; jenv++) {

                	int J_min=abs(sys->Value_J_block[n_1][jsys]-env->Value_J_block[n_2][jenv]);
                	int J_max=(sys->Value_J_block[n_1][jsys]+env->Value_J_block[n_2][jenv]);
                	int J_num=(J_max-J_min)/2+1;

                	for(int n=0; n<J_num; n++) {

                        	int J_mid=(J_min+2*n);

				int J_ne_min=abs(n_n-n_e);
				int J_ne_max=(2-n_n-n_e);
				int J_ne_num=(J_ne_max-J_ne_min)/2+1;
		
				for(int m=0; m<J_ne_num; m++) {

					int J_ne=(J_ne_min+2*m);

					int J_total_min=abs(J_mid-J_ne);
					int J_total_max=(J_mid+J_ne);
					int J_total_num=(J_total_max-J_total_min)/2+1;
				
					for(int l=0; l<J_total_num; l++) {
						int J_total=(J_total_min+2*l);
						if(J_total== QuantumNumber_J ) {
							H_sys_config_3[index]=n_1;	H_env_config_3[index]=n_2;
							H_ns_config_3[index]=n_n;	H_ne_config_3[index]=n_e;
							J_sys_config_3[index]=jsys;	J_env_config_3[index]=jenv;
							J_sysnew_config_3[index]=J_mid;
							J_envnew_config_3[index]=J_ne;
							Dim_block_config_3[index++]=sys->Dim_J_block[n_1][jsys]*env->Dim_J_block[n_2][jenv];
						}
					}

                        	}

                	}

        	}

	}

//	for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++)
//		cout<<"\n sys_J["<<i<<"]="<<sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]]<<"\t env_J["<<i<<"]="<<env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]];

  //------Check the Dimensions of the different configurations!!!------------------------------------------------
        Dim_config_3=0;

        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++)
                Dim_config_3+=Dim_block_config_3[i];

        if(Dim_config_3!=Dim)   cout<<"\n Error: Dim_config_3 is not equal to Dim_config_1,2"<<endl;

  //------Create space for wave function f3 and g3---------------------------------------------------------------
        f3=new double * [BlockNumber_for_TargetBlock_config_3];
        g3=new double * [BlockNumber_for_TargetBlock_config_3];
        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) {
                f3[i]=new double [Dim_block_config_3[i]];        g3[i]=new double [Dim_block_config_3[i]];
                for(int j=0; j<Dim_block_config_3[i]; j++) {
                        f3[i][j]=(double) 0;            g3[i][j]=(double) 0;
                }
        }

}

//=============================Alocate the 6j coefficients for configurations 1 and 2=============================
inline void Super::Allocate_6j_config_1_2(Super &space) {

	int i_g, i_f;
	int Length=BlockNumber_for_TargetBlock*BlockNumber_for_TargetBlock;

	six_j_1_sys_n=new double [Length];
	six_j_1_e_env=new double [Length];

	for(int i=0; i<Length; i++) {
		six_j_1_sys_n[i]=0.0;  six_j_1_e_env[i]=0.0;
	}

        for(int i=0; i<Length; i++) {

                i_g=i/BlockNumber_for_TargetBlock;//stands for the index of "g" function, "g" function is bra
                i_f=i%BlockNumber_for_TargetBlock;//stands for the index of "f" function  "f" function is ket

		if((sys->TotSiteNo-sys->Num_hole_block[H_sys[i_g]])>0 && sys->Num_hole_block[H_sys[i_g]]==space.sysnew_space->Num_hole_block[H_sysnew[i_g]] && H_sysnew[i_f]==H_sysnew[i_g] && H_envnew[i_f]==H_envnew[i_g] && H_env[i_f]==H_env[i_g] && H_sys[i_f]==H_sys[i_g] && J_sysnew[i_f]==J_sysnew[i_g] && J_envnew[i_f]==J_envnew[i_g] && J_env[i_f]==J_env[i_g] && (sys->Value_J_block[H_sys[i_f]][J_sys[i_f]]-sys->Value_J_block[H_sys[i_g]][J_sys[i_g]])<=2) {

                	six_j_1_sys_n[i]=gsl_sf_coupling_6j(space.sysnew_space->Value_J_block[H_sysnew[i_f]][J_sysnew[i_f]], 1, sys->Value_J_block[H_sys[i_g]][J_sys[i_g]], 2, sys->Value_J_block[H_sys[i_f]][J_sys[i_f]], 1)*pow(-1.0, (1+space.sysnew_space->Value_J_block[H_sysnew[i_f]][J_sysnew[i_f]]+sys->Value_J_block[H_sys[i_f]][J_sys[i_f]])/2)*constant; 

		}
	
		if((env->TotSiteNo-env->Num_hole_block[H_env[i_g]])>0 && env->Num_hole_block[H_env[i_g]]==envnew->Num_hole_block[H_envnew[i_g]] && H_sysnew[i_f]==H_sysnew[i_g] && H_envnew[i_f]==H_envnew[i_g] && H_env[i_f]==H_env[i_g] && H_sys[i_f]==H_sys[i_g] && J_sysnew[i_f]==J_sysnew[i_g] && J_envnew[i_f]==J_envnew[i_g] && J_sys[i_f]==J_sys[i_g] && (env->Value_J_block[H_env[i_f]][J_env[i_f]]-env->Value_J_block[H_env[i_g]][J_env[i_g]])<=2) {

			six_j_1_e_env[i]=gsl_sf_coupling_6j(envnew->Value_J_block[H_envnew[i_f]][J_envnew[i_f]], 1, env->Value_J_block[H_env[i_g]][J_env[i_g]], 2, env->Value_J_block[H_env[i_f]][J_env[i_f]], 1)*pow(-1.0, (1+envnew->Value_J_block[H_envnew[i_f]][J_envnew[i_f]]+env->Value_J_block[H_env[i_f]][J_env[i_f]])/2)*constant; 

		}

	}

}

//==============================Allocate the 9j coefficients for configuration 2==================================
//Define the 9j coefficient for the basis transformation between configurations 1 and 2!!! From (sys+ns & env+ne)
//to (sys+ne & env+ns)!!! Obtained directly by the definition of 9j coefficient from 3j coefficients!!!
//================================================================================================================
inline void Super::Allocate_9j_config_2(Parameter &para, Super &space) {
  //------Create space--------------------------------------------------------------------------------------------
        nine_j_config_2=new double * [BlockNumber_for_TargetBlock];
        nine_j_config_2_inverse=new double * [BlockNumber_for_TargetBlock];

        for(int i=0; i<BlockNumber_for_TargetBlock; i++) {	//"i" denote the basis of configuration II

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock; j++)
		if(H_sys[i]==H_sys[j] && H_env[i]==H_env[j] && J_sys[i]==J_sys[j] && J_env[i]==J_env[j] && (space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]])==(envnew->Num_hole_block[H_envnew[j]]-env->Num_hole_block[H_env[j]]) && (space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]])==(envnew->Num_hole_block[H_envnew[i]]-env->Num_hole_block[H_env[i]]))
			index++;

                nine_j_config_2[i]=new double [index];
                nine_j_config_2_inverse[i]=new double [index];

                for(int n=0; n<index; n++) {
                        nine_j_config_2[i][n]=0.0;  nine_j_config_2_inverse[i][n]=0.0;

		}

        }

  //------Initialize 9-j coefficient for config_2----------------------------------------------------------------
        for(int i=0; i<BlockNumber_for_TargetBlock; i++) {

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock; j++)
		if(H_sys[i]==H_sys[j] && H_env[i]==H_env[j] && J_sys[i]==J_sys[j] && J_env[i]==J_env[j] && (space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]])==(envnew->Num_hole_block[H_envnew[j]]-env->Num_hole_block[H_env[j]]) && (space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]])==(envnew->Num_hole_block[H_envnew[i]]-env->Num_hole_block[H_env[i]])) {

			alpha=sqrt((space.sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]]+1.0)*(envnew->Value_J_block[H_envnew[i]][J_envnew[i]]+1.0)*(space.sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]]+1.0)*(envnew->Value_J_block[H_envnew[j]][J_envnew[j]]+1.0))*(QuantumNumber_J+1.0)*pow(-1.0,(-space.sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]]+envnew->Value_J_block[H_envnew[i]][J_envnew[i]]-space.sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]]+envnew->Value_J_block[H_envnew[j]][J_envnew[j]]-2*(sys->Value_J_block[H_sys[i]][J_sys[i]]+env->Value_J_block[H_env[i]][J_env[i]])+2*(2-space.sysnew_space->Num_hole_block[H_sysnew[i]]+sys->Num_hole_block[H_sys[i]]-space.sysnew_space->Num_hole_block[H_sysnew[j]]+sys->Num_hole_block[H_sys[j]])-4*QuantumNumber_J)/2);

			if(((space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]])==0 && (space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]])==0) || ((space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]])==0 && (space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]])==1 && (env->TotSiteNo-env->Num_hole_block[H_env[i]])%2==1) || ((space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]])==1 && (space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]])==0 && (env->TotSiteNo-env->Num_hole_block[H_env[i]])%2==1))
				alpha=-alpha;

                        for(int m_n=0; m_n<(2-space.sysnew_space->Num_hole_block[H_sysnew[j]]+sys->Num_hole_block[H_sys[j]]); m_n++)
                        for(int m_e=0; m_e<(2-space.sysnew_space->Num_hole_block[H_sysnew[i]]+sys->Num_hole_block[H_sys[i]]); m_e++)
                        for(int m_1=0; m_1<sys->Value_J_block[H_sys[i]][J_sys[i]]+1; m_1++)
                        for(int m_2=0; m_2<env->Value_J_block[H_env[i]][J_env[i]]+1; m_2++)
			if((space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]]+space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]]+2*(m_n+m_e+m_1+m_2-1)-sys->Value_J_block[H_sys[i]][J_sys[i]]-env->Value_J_block[H_env[i]][J_env[i]])==QuantumNumber_J)
				nine_j_config_2[i][index]+=alpha*gsl_sf_coupling_3j(sys->Value_J_block[H_sys[i]][J_sys[i]], 1-space.sysnew_space->Num_hole_block[H_sysnew[i]]+sys->Num_hole_block[H_sys[i]], space.sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]], -sys->Value_J_block[H_sys[i]][J_sys[i]]+2*m_1, space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]]-1+2*m_e, -space.sysnew_space->Num_hole_block[H_sysnew[i]]+sys->Num_hole_block[H_sys[i]]+1-2*(m_e+m_1)+sys->Value_J_block[H_sys[j]][J_sys[j]])*gsl_sf_coupling_3j(env->Value_J_block[H_env[j]][J_env[j]], 1-envnew->Num_hole_block[H_envnew[i]]+env->Num_hole_block[H_env[i]], envnew->Value_J_block[H_envnew[i]][J_envnew[i]], -env->Value_J_block[H_env[j]][J_env[j]]+2*m_2, envnew->Num_hole_block[H_envnew[i]]-env->Num_hole_block[H_env[i]]-1+2*m_n, -envnew->Num_hole_block[H_envnew[i]]+env->Num_hole_block[H_env[i]]+1-2*(m_n+m_2)+env->Value_J_block[H_env[j]][J_env[j]])*gsl_sf_coupling_3j(space.sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]], envnew->Value_J_block[H_envnew[j]][J_envnew[j]], QuantumNumber_J, space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]]-1+2*(m_n+m_1)-sys->Value_J_block[H_sys[j]][J_sys[j]], envnew->Num_hole_block[H_envnew[j]]-env->Num_hole_block[H_env[j]]-1+2*(m_e+m_2)-env->Value_J_block[H_env[j]][J_env[j]], -QuantumNumber_J)*gsl_sf_coupling_3j(sys->Value_J_block[H_sys[j]][J_sys[j]], 1-space.sysnew_space->Num_hole_block[H_sysnew[j]]+sys->Num_hole_block[H_sys[j]], space.sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]], -sys->Value_J_block[H_sys[j]][J_sys[j]]+2*m_1, space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]]-1+2*m_n, -space.sysnew_space->Num_hole_block[H_sysnew[j]]+sys->Num_hole_block[H_sys[j]]+1-2*(m_n+m_1)+sys->Value_J_block[H_sys[j]][J_sys[j]])*gsl_sf_coupling_3j(env->Value_J_block[H_env[j]][J_env[j]], 1-envnew->Num_hole_block[H_envnew[j]]+env->Num_hole_block[H_env[j]], envnew->Value_J_block[H_envnew[j]][J_envnew[j]], -env->Value_J_block[H_env[j]][J_env[j]]+2*m_2, envnew->Num_hole_block[H_envnew[j]]-env->Num_hole_block[H_env[j]]-1+2*m_e, -envnew->Num_hole_block[H_envnew[j]]+env->Num_hole_block[H_env[j]]+1-2*(m_e+m_2)+env->Value_J_block[H_env[j]][J_env[j]])*gsl_sf_coupling_3j(space.sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]], envnew->Value_J_block[H_envnew[i]][J_envnew[i]], QuantumNumber_J, space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]]-1+2*(m_e+m_1)-sys->Value_J_block[H_sys[i]][J_sys[i]], envnew->Num_hole_block[H_envnew[i]]-env->Num_hole_block[H_env[i]]-1+2*(m_n+m_2)-env->Value_J_block[H_env[i]][J_env[i]], -QuantumNumber_J);

                        for(int m_n=0; m_n<(2-space.sysnew_space->Num_hole_block[H_sysnew[i]]+sys->Num_hole_block[H_sys[i]]); m_n++)
                        for(int m_e=0; m_e<(2-space.sysnew_space->Num_hole_block[H_sysnew[j]]+sys->Num_hole_block[H_sys[j]]); m_e++)
                        for(int m_1=0; m_1<sys->Value_J_block[H_sys[i]][J_sys[i]]+1; m_1++)
                        for(int m_2=0; m_2<env->Value_J_block[H_env[i]][J_env[i]]+1; m_2++)
			if((space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]]+space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]]+2*(m_n+m_e+m_1+m_2-1)-sys->Value_J_block[H_sys[i]][J_sys[i]]-env->Value_J_block[H_env[i]][J_env[i]])==QuantumNumber_J)
				nine_j_config_2_inverse[i][index]+=alpha*gsl_sf_coupling_3j(sys->Value_J_block[H_sys[i]][J_sys[i]], 1-space.sysnew_space->Num_hole_block[H_sysnew[j]]+sys->Num_hole_block[H_sys[j]], space.sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]], -sys->Value_J_block[H_sys[i]][J_sys[i]]+2*m_1, space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]]-1+2*m_e, -space.sysnew_space->Num_hole_block[H_sysnew[j]]+sys->Num_hole_block[H_sys[j]]+1-2*(m_e+m_1)+sys->Value_J_block[H_sys[i]][J_sys[i]])*gsl_sf_coupling_3j(env->Value_J_block[H_env[i]][J_env[i]], 1-envnew->Num_hole_block[H_envnew[j]]+env->Num_hole_block[H_env[j]], envnew->Value_J_block[H_envnew[j]][J_envnew[j]], -env->Value_J_block[H_env[i]][J_env[i]]+2*m_2, envnew->Num_hole_block[H_envnew[j]]-env->Num_hole_block[H_env[j]]-1+2*m_n, -envnew->Num_hole_block[H_envnew[j]]+env->Num_hole_block[H_env[j]]+1-2*(m_n+m_2)+env->Value_J_block[H_env[i]][J_env[i]])*gsl_sf_coupling_3j(space.sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]], envnew->Value_J_block[H_envnew[i]][J_envnew[i]], QuantumNumber_J, space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]]-1+2*(m_n+m_1)-sys->Value_J_block[H_sys[i]][J_sys[i]], envnew->Num_hole_block[H_envnew[i]]-env->Num_hole_block[H_env[i]]-1+2*(m_e+m_2)-env->Value_J_block[H_env[i]][J_env[i]], -QuantumNumber_J)*gsl_sf_coupling_3j(sys->Value_J_block[H_sys[i]][J_sys[i]], 1-space.sysnew_space->Num_hole_block[H_sysnew[i]]+sys->Num_hole_block[H_sys[i]], space.sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]], -sys->Value_J_block[H_sys[i]][J_sys[i]]+2*m_1, space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]]-1+2*m_n, -space.sysnew_space->Num_hole_block[H_sysnew[i]]+sys->Num_hole_block[H_sys[i]]+1-2*(m_n+m_1)+sys->Value_J_block[H_sys[j]][J_sys[i]])*gsl_sf_coupling_3j(env->Value_J_block[H_env[i]][J_env[i]], 1-envnew->Num_hole_block[H_envnew[i]]+env->Num_hole_block[H_env[i]], envnew->Value_J_block[H_envnew[i]][J_envnew[i]], -env->Value_J_block[H_env[i]][J_env[i]]+2*m_2, envnew->Num_hole_block[H_envnew[i]]-env->Num_hole_block[H_env[i]]-1+2*m_e, -envnew->Num_hole_block[H_envnew[i]]+env->Num_hole_block[H_env[i]]+1-2*(m_e+m_2)+env->Value_J_block[H_env[i]][J_env[i]])*gsl_sf_coupling_3j(space.sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]], envnew->Value_J_block[H_envnew[j]][J_envnew[j]], QuantumNumber_J, space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]]-1+2*(m_e+m_1)-sys->Value_J_block[H_sys[j]][J_sys[j]], envnew->Num_hole_block[H_envnew[j]]-env->Num_hole_block[H_env[j]]-1+2*(m_n+m_2)-env->Value_J_block[H_env[j]][J_env[j]], -QuantumNumber_J);

			index++;

		}

	}

}

//====================================Allocate 6j coefficients for configuration 3===============================
inline void Super::Allocate_6j_config_3() {

//------Create space for six_j_T_config_3------------------------------------------------------------------------

        six_j_T_config_3 = new double * [ BlockNumber_for_TargetBlock_config_3 ];

        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i ++ )
	if( ( ( sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] ) > 0 )  &&  ( sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] ) < ( sys -> TotSiteNo + env -> TotSiteNo ) ) {

		index = 0;

		for( int j = 0; j < BlockNumber_for_TargetBlock_config_3; j ++ )
		if( ( sys -> Num_hole_block[ H_sys_config_3[ j ] ] + env -> Num_hole_block[ H_env_config_3[ j ] ]  ==  sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] )  &&  H_ns_config_3[ j ] == H_ns_config_3[ i ]  &&  H_ne_config_3[ j ] == H_ne_config_3[ i ]  &&  abs( sys -> Num_hole_block[ H_sys_config_3[ i ] ] - sys -> Num_hole_block[ H_sys_config_3[ j ] ] ) == 1  &&  abs( env -> Num_hole_block[ H_env_config_3[ i ] ] - env -> Num_hole_block[ H_env_config_3[ j ] ] ) == 1  &&  J_sysnew_config_3[ i ] == J_sysnew_config_3[ j ]  &&  J_envnew_config_3[ i ] == J_envnew_config_3[ j ]  &&  abs( sys -> Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] - sys -> Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] ) == 1  &&  abs( env -> Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] - env -> Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ] ) == 1 )

			index++;

                six_j_T_config_3[ i ] = new double [ index ];

                for( int n = 0; n < index; n ++ ) 
                        six_j_T_config_3[ i ][ n ] = 0.0;

	}

  //------Initialize six_j_T_config_3------------------------------------------
        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) 
	if( ( ( sys->Num_hole_block[ H_sys_config_3[i] ] + env->Num_hole_block[ H_env_config_3[i] ] ) > 0 )  &&  ( sys->Num_hole_block[ H_sys_config_3[i] ] + env->Num_hole_block[ H_env_config_3[i] ] ) < ( sys->TotSiteNo + env->TotSiteNo ) ) {

		index = 0;

		for(int j=0; j<BlockNumber_for_TargetBlock_config_3; j++)
		if( ( sys -> Num_hole_block[ H_sys_config_3[ j ] ] + env -> Num_hole_block[ H_env_config_3[ j ] ]  ==  sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] )  &&  H_ns_config_3[ j ] == H_ns_config_3[ i ]  &&  H_ne_config_3[ j ] == H_ne_config_3[ i ]  &&  abs( sys -> Num_hole_block[ H_sys_config_3[ i ] ] - sys -> Num_hole_block[ H_sys_config_3[ j ] ] ) == 1  &&  abs( env -> Num_hole_block[ H_env_config_3[ i ] ] - env -> Num_hole_block[ H_env_config_3[ j ] ] ) == 1  &&  J_sysnew_config_3[ i ] == J_sysnew_config_3[ j ]  &&  J_envnew_config_3[ i ] == J_envnew_config_3[ j ]  &&  abs( sys -> Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] - sys -> Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] ) == 1  &&  abs( env -> Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] - env -> Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ] ) == 1 )

			six_j_T_config_3[ i ][ index++ ] = pow( -1.0, ( sys -> Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] + env -> Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] + J_sysnew_config_3[ i ] ) / 2 ) * gsl_sf_coupling_6j( sys -> Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], sys -> Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ], 1, env -> Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ], env -> Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], J_sysnew_config_3[ i ] );

	}

//------Create space for six_j_J_config_3------------------------------------------------------------------------
        six_j_J_config_3 = new double * [BlockNumber_for_TargetBlock_config_3];

        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i ++ )
	if( sys -> Num_hole_block[ H_sys_config_3[ i ] ] != sys -> TotSiteNo  &&  env -> Num_hole_block[ H_env_config_3[ i ] ] != env -> TotSiteNo ) {

		index = 0;

		for(int j=0; j<BlockNumber_for_TargetBlock_config_3; j++)
		if( H_sys_config_3[i] == H_sys_config_3[j]  &&  H_env_config_3[i] == H_env_config_3[j]  &&  H_ns_config_3[j] == H_ns_config_3[i]  &&  H_ne_config_3[j] == H_ne_config_3[i]  &&  J_sysnew_config_3[i] == J_sysnew_config_3[j]  &&  J_envnew_config_3[i] == J_envnew_config_3[j]  &&  abs( sys -> Value_J_block[ H_sys_config_3[i] ][ J_sys_config_3[i] ] - sys -> Value_J_block[ H_sys_config_3[j] ][ J_sys_config_3[j] ] ) <= 2  &&  abs( env -> Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ] - env -> Value_J_block[ H_env_config_3[j] ][ J_env_config_3[j] ] ) <= 2 )

				index++;

                six_j_J_config_3[i] = new double [index];
                for(int n=0; n<index; n++)
                        six_j_J_config_3[i][n]=(double) 0;

	}

  //------Initialize six_j_J_config_3------------------------------------------
        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i ++) 
	if( sys -> Num_hole_block[ H_sys_config_3[ i ] ] != sys -> TotSiteNo  &&  env -> Num_hole_block[ H_env_config_3[ i ] ] != env -> TotSiteNo ) {

		index = 0;

		for(int j=0; j<BlockNumber_for_TargetBlock_config_3; j++)
		if( H_sys_config_3[ i ] == H_sys_config_3[ j ]  &&  H_env_config_3[ i ] == H_env_config_3[ j ]  &&  H_ns_config_3[ j ] == H_ns_config_3[ i ]  &&  H_ne_config_3[ j ] == H_ne_config_3[ i ]  &&  J_sysnew_config_3[ i ] == J_sysnew_config_3[ j ]  &&  J_envnew_config_3[ i ] == J_envnew_config_3[ j ]  &&  abs( sys -> Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] - sys -> Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] ) <= 2  &&  abs( env -> Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] - env -> Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ] ) <= 2 ) {

			six_j_J_config_3[ i ][ index++ ] = pow( -1.0, ( sys->Value_J_block[ H_sys_config_3[j] ][ J_sys_config_3[j] ] + env->Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ] + J_sysnew_config_3[i] ) / 2 ) * gsl_sf_coupling_6j( J_sysnew_config_3[i], env->Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ], sys->Value_J_block[ H_sys_config_3[i] ][ J_sys_config_3[i] ], 2, sys->Value_J_block[ H_sys_config_3[j] ][ J_sys_config_3[j] ], env->Value_J_block[ H_env_config_3[j] ][ J_env_config_3[j] ] );

		}

	}

//------Create space for six_j_J_config_3_diaele for the use in calculating the diagonal elements of H matrix----
        six_j_J_config_3_diaele = new double [ BlockNumber_for_TargetBlock_config_3 ];
        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ ) 

		six_j_J_config_3_diaele[i] = pow( -1.0, ( sys -> Value_J_block[ H_sys_config_3[i] ][ J_sys_config_3[i] ] + env -> Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ] + J_sysnew_config_3[i] ) / 2 ) * gsl_sf_coupling_6j( J_sysnew_config_3[i], env -> Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ], sys -> Value_J_block[ H_sys_config_3[i] ][ J_sys_config_3[i] ], 2, sys -> Value_J_block[ H_sys_config_3[i] ][ J_sys_config_3[i] ], env -> Value_J_block[ H_env_config_3[i] ][ J_env_config_3[i] ] );

//------Create and Initialize six_j_2_config_3
        six_j_2_config_3 = new double [ BlockNumber_for_TargetBlock_config_3 ];

        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ ) 

		six_j_2_config_3[ i ] = 0.0;

        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ ) 
	if( H_ns_config_3[ i ] == 0  &&  H_ne_config_3[ i ] == 0 ) {

                if( J_envnew_config_3[ i ] == 0 )      

			six_j_2_config_3[ i ] = -0.75;

                else                            

			six_j_2_config_3[ i ] = 0.25;

        }

}

//============================Allocate 9j coefficients for the configuration 3===================================
inline void Super::Allocate_9j_config_3(Parameter &para, Super &space) {
//-----------------------------------------9j coefficient for config 3--------------------------------------------
  //------Create space

        nine_j_config_3=new double * [BlockNumber_for_TargetBlock_config_3];

        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) {

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock; j++)
		if(H_sys[j]==H_sys_config_3[i] && H_env[j]==H_env_config_3[i] && J_sys[j]==J_sys_config_3[i] && J_env[j]==J_env_config_3[i] && H_ns_config_3[i]==(space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]]) && H_ne_config_3[i]==(envnew->Num_hole_block[H_envnew[j]]-env->Num_hole_block[H_env[j]]))
			index++;

                nine_j_config_3[i]=new double [index];

                for(int n=0; n<index; n++)
                        nine_j_config_3[i][n]=(double) 0;

        }

  //------Initialize 9-j coefficient
//-------------------------------------Obtained by the 9-j coefficient formula-----------------------------------
        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) {

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock; j++)
		if(H_sys[j]==H_sys_config_3[i] && H_env[j]==H_env_config_3[i] && J_sys[j]==J_sys_config_3[i] && J_env[j]==J_env_config_3[i] && H_ns_config_3[i]==(space.sysnew_space->Num_hole_block[H_sysnew[j]]-sys->Num_hole_block[H_sys[j]]) && H_ne_config_3[i]==(envnew->Num_hole_block[H_envnew[j]]-env->Num_hole_block[H_env[j]])) {

			alpha=1.0;

			if(H_ns_config_3[i]==0 && (env->TotSiteNo-env->Num_hole_block[H_env[j]])%2==1)	alpha=-1.0;
		
			nine_j_config_3[i][index++]=alpha*gsl_sf_coupling_9j(sys->Value_J_block[H_sys[j]][J_sys[j]], env->Value_J_block[H_env[j]][J_env[j]], J_sysnew_config_3[i], 1-H_ns_config_3[i], 1-H_ne_config_3[i], J_envnew_config_3[i], space.sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]], envnew->Value_J_block[H_envnew[j]][J_envnew[j]], QuantumNumber_J)*sqrt((J_sysnew_config_3[i]+1.0)*(J_envnew_config_3[i]+1.0)*(space.sysnew_space->Value_J_block[H_sysnew[j]][J_sysnew[j]]+1.0)*(envnew->Value_J_block[H_envnew[j]][J_envnew[j]]+1.0));
                        
		}

        }

//---------------------------------------------------------------------------------------------------------------
//----------------------------------------inverse 9-j coefficient for config 3------------------------------------
  //------Create space
        nine_j_config_3_inverse=new double * [BlockNumber_for_TargetBlock];

        for(int i=0; i<BlockNumber_for_TargetBlock; i++) {

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock_config_3; j++) 
		if(H_sys[i]==H_sys_config_3[j] && H_env[i]==H_env_config_3[j] && J_sys[i]==J_sys_config_3[j] && J_env[i]==J_env_config_3[j] && H_ns_config_3[j]==(space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]]) && H_ne_config_3[j]==(envnew->Num_hole_block[H_envnew[i]]-env->Num_hole_block[H_env[i]])) 
			index++;

                nine_j_config_3_inverse[i]=new double [index];

                for(int n=0; n<index; n++) {
                        nine_j_config_3_inverse[i][n]=(double) 0;

                }

        }

//-------------------------------------Obtained by the 9-j formula------------------------------------------------
  //------Initialize inverse 9-j coefficient
        for(int i=0; i<BlockNumber_for_TargetBlock; i++) {

                index=0;

                for(int j=0; j<BlockNumber_for_TargetBlock_config_3; j++) 
		if(H_sys[i]==H_sys_config_3[j] && H_env[i]==H_env_config_3[j] && J_sys[i]==J_sys_config_3[j] && J_env[i]==J_env_config_3[j] && H_ns_config_3[j]==(space.sysnew_space->Num_hole_block[H_sysnew[i]]-sys->Num_hole_block[H_sys[i]]) && H_ne_config_3[j]==(envnew->Num_hole_block[H_envnew[i]]-env->Num_hole_block[H_env[i]])) {

			alpha=1.0;
			if(H_ns_config_3[j]==0 && (env->TotSiteNo-env->Num_hole_block[H_env[i]])%2==1) alpha=-1.0;

			nine_j_config_3_inverse[i][index++]=alpha*gsl_sf_coupling_9j(sys->Value_J_block[H_sys[i]][J_sys[i]], 1-H_ns_config_3[j], space.sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]], env->Value_J_block[H_env[i]][J_env[i]], 1-H_ne_config_3[j], envnew->Value_J_block[H_envnew[i]][J_envnew[i]], J_sysnew_config_3[j], J_envnew_config_3[j], QuantumNumber_J)*sqrt((J_sysnew_config_3[j]+1.0)*(J_envnew_config_3[j]+1.0)*(space.sysnew_space->Value_J_block[H_sysnew[i]][J_sysnew[i]]+1.0)*(envnew->Value_J_block[H_envnew[i]][J_envnew[i]]+1.0));
 
		}

        }

}

//===========================getMatrixDiagonalElement================================
void Super::getMatrixDiagElement(double *f, const int &dim) {

//----------------------Realization by direct multiplication-------------------------------------
//      for(int i=0; i<dim; i++)
//              f[i]+=1.0;

	//entries of the diagonal elements of the Hamiltonian matrix
        int a_sys, a_env, Num_six_j, old_sys, old_env;

        //Contributions from H_sys, H_env, H_sys_ns, and H_env_ne
        for( int i = 0; i < BlockNumber_for_TargetBlock; i ++ )
        for( int j = 0; j < Dim_block[ i ]; j ++ ) {

		Num_six_j = i * BlockNumber_for_TargetBlock + i;

        	a_sys = j % sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ];
                a_env = j / sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ];

		old_sys = a_sys * sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ] + a_sys;
                old_env = a_env * env -> Dim_J_block[ H_env[ i ] ][ J_env[ i ] ] + a_env;

                f[ Table_2to1[ i ][ j ] ] += ( sys -> H[ H_sys[ i ] ][ J_sys[ i ] ][ old_sys ] + env -> H[ H_env[ i ] ][ J_env[ i ] ][ old_env ] );
                 
		//contributions from H_sys_ns_J, and H_sys_ns_N
       		SiteNum = N * sys -> TotSiteNo;

		if( ( sys -> TotSiteNo - sys -> Num_hole_block[ H_sys[i] ] ) > 0  &&  sys -> Num_hole_block[H_sys[i]] == sysnew -> Num_hole_block[ H_sysnew[i] ] ) {	//H_sys_ns_J

	                for( site_s = 0; site_s < operator_number_J_sys; site_s ++ )
        	        if( Table_J[ SiteNum + Table_J_sys[ site_s ] ] == 1 )

        			f[ Table_2to1[ i ][ j ] ] += Interaction_J[ SiteNum + Table_J_sys[ site_s ] ] * six_j_1_sys_n[ Num_six_j ] * sys -> S_Dia[ site_s ][ H_sys[ i ] ][ J_sys[ i ] ][ old_sys ];

		        for( site_s = 0; site_s < operator_number_N_sys; site_s ++ )
		        if( Table_N[ SiteNum + Table_N_sys[ site_s ] ] == 1 )

        			f[ Table_2to1[ i ][ j ] ] += Interaction_N[ SiteNum + Table_N_sys[ site_s ] ] * sys -> NN[ site_s ][ H_sys[ i ] ][ J_sys[ i ] ][ old_sys ] / sqrt( sys -> Value_J_block[ H_sys[i] ][ J_sys[i] ] + 1.0 );

		}

		//contributions from H_env_ne J couplings
        	SiteNum = N * ( StartSite - env -> TotSiteNo ) + StartSite;

        	if( (env -> TotSiteNo - env -> Num_hole_block[H_env[i]] ) > 0   &&  env -> Num_hole_block[H_env[i]] == envnew -> Num_hole_block[H_envnew[i]] ) {

                	for( site_e = 0; site_e < operator_number_J_env; site_e ++ )
	                if( Table_J[ SiteNum - Table_J_env[ site_e ] ] == 1 )

	                	f[ Table_2to1[ i ][ j ] ] += Interaction_J[ SiteNum - Table_J_env[ site_e ] ] * six_j_1_e_env[ Num_six_j ] * env -> S_Dia[ site_e ][ H_env[ i ] ][ J_env[ i ] ][ old_env ];

			for( site_e = 0; site_e < operator_number_N_env; site_e ++ )
	                if( Table_N[ SiteNum - Table_N_env[ site_e ] ] == 1 )

	                	f[ Table_2to1[ i ][ j ] ] += Interaction_N[ SiteNum - Table_N_env[ site_e ] ] * env -> NN[ site_e ][ H_env[ i ] ][ J_env[ i ] ][ old_env ] / sqrt( env -> Value_J_block[ H_env[i] ][ J_env[i] ] + 1.0 );

		}

	}

	//contributions from H_sys_ne and H_env_ns
        for( int i = 0; i < BlockNumber_for_TargetBlock; i ++ ) {

                index = 0;

                for( int i_p = 0; i_p < BlockNumber_for_TargetBlock; i_p ++ )
                if( H_sys[ i_p ] == H_sys[ i ]  &&  H_env[ i_p ] == H_env[ i ]  &&  J_sys[ i_p ] == J_sys[ i ]  &&  J_env[ i_p ] == J_env[ i ] ) {

                        Num_six_j = i_p * BlockNumber_for_TargetBlock + i_p;

                        for( int j = 0; j < Dim_block[ i ]; j ++ ) {

                                a_sys = j % sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ];
                                a_env = j / sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ];

                                old_sys = a_sys * sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ] + a_sys;
                                old_env = a_env * env -> Dim_J_block[ H_env[ i ] ][ J_env[ i ] ] + a_env;
 
                                SiteNum = N * ( StartSite - env -> TotSiteNo );

                                for( site_s = 0; site_s < operator_number_J_sys; site_s ++ )
                                if( Table_J[ SiteNum + Table_J_sys[ site_s ] ] == 1 )
                                        f[ Table_2to1[ i ][ j ] ] += Interaction_J[ SiteNum + Table_J_sys[ site_s ] ] * nine_j_config_2[ i ][ index ] * nine_j_config_2[ i ][ index ] * six_j_1_sys_n[ Num_six_j ] * sys -> S_Dia[ site_s ][ H_sys[ i ]][ J_sys[ i ] ][ old_sys ];
 
			        for( site_s = 0; site_s < operator_number_N_sys; site_s ++ )
			        if( Table_N[ SiteNum + Table_N_sys[ site_s ] ] == 1 )
        				f[ Table_2to1[ i ][ j ] ] += Interaction_N[ SiteNum + Table_N_sys[ site_s ] ] * nine_j_config_2[ i ][ index ] * nine_j_config_2[ i ][ index ] * sys -> NN[ site_s ][ H_sys[ i ] ][ J_sys[ i ] ][ old_sys ] / sqrt( sys -> Value_J_block[ H_sys[i] ][ J_sys[i] ] + 1.0 );
                              

                                SiteNum = N * sys -> TotSiteNo + StartSite;

                                for( site_e = 0; site_e < operator_number_J_env; site_e ++ )
                                if( Table_J[ SiteNum - Table_J_env[ site_e ] ] == 1 )
                                        f[ Table_2to1[ i ][ j ] ] += Interaction_J[ SiteNum - Table_J_env[ site_e ] ] * nine_j_config_2[ i ][ index ] * nine_j_config_2[ i ][ index ] * six_j_1_e_env[ Num_six_j ] * env -> S_Dia[ site_e ][ H_env[ i ] ][ J_env[ i ] ][ old_env ];

                                for( site_e = 0; site_e < operator_number_N_env; site_e ++ )
                                if( Table_N[ SiteNum - Table_N_env[ site_e ] ] == 1 )
                                        f[ Table_2to1[ i ][ j ] ] += Interaction_N[ SiteNum - Table_N_env[ site_e ] ] * nine_j_config_2[ i ][ index ] * nine_j_config_2[ i ][ index ] * env -> NN[ site_e ][ H_env[ i ] ][ J_env[ i ] ][ old_env ] / sqrt( env -> Value_J_block[ H_env[i] ][ J_env[i] ] + 1.0 );

                        }

                        index++;
                }

        }

	//Contributions from H_sys_env and H_n_e
        for( int i = 0; i < BlockNumber_for_TargetBlock; i ++ ) {

                index = 0;

                for( int i_p = 0; i_p < BlockNumber_for_TargetBlock_config_3; i_p ++ )
                if( H_sys_config_3[ i_p ] == H_sys[ i ]  &&  H_env_config_3[ i_p ] == H_env[ i ]  &&  J_sys_config_3[ i_p ] == J_sys[ i ]  &&  J_env_config_3[ i_p ] == J_env[ i ] ) {

                        for( int j = 0; j < Dim_block[ i ]; j ++ ) {

                                a_sys = j % sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ];
                                a_env = j / sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ];

                                old_sys = a_sys * sys -> Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ] + a_sys;
                                old_env = a_env * env -> Dim_J_block[ H_env[ i ] ][ J_env[ i ] ] + a_env;

                             //------contributions from Interactions between system and environment--------------
                                for( site_s = 0; site_s < operator_number_J_sys; site_s ++ )
                                for( site_e = 0; site_e < operator_number_J_env; site_e ++ )
                                if( Table_J[ N * Table_J_sys[ site_s ] + StartSite - Table_J_env[ site_e ] ] == 1 )
                                        f[ Table_2to1[ i ][ j ] ] += Interaction_J[ N * Table_J_sys[ site_s ] + StartSite - Table_J_env[ site_e ] ] * nine_j_config_3_inverse[ i ][ index ] * nine_j_config_3_inverse[ i ][ index ] * six_j_J_config_3_diaele[ i_p ] * sys -> S_Dia[ site_s ][ H_sys[ i ] ][ J_sys[ i ] ][ old_sys ] * env -> S_Dia[ site_e ][ H_env[ i ] ][ J_env[ i ] ][ old_env ];
 
                                for( site_s = 0; site_s < operator_number_N_sys; site_s ++ )
                                for( site_e = 0; site_e < operator_number_N_env; site_e ++ )
                                if( Table_N[ N * Table_N_sys[ site_s ] + StartSite - Table_N_env[ site_e ] ] == 1 )
                                        f[ Table_2to1[ i ][ j ] ] += Interaction_N[ N * Table_N_sys[ site_s ] + StartSite - Table_N_env[ site_e ] ] * nine_j_config_3_inverse[ i ][ index ] * nine_j_config_3_inverse[ i ][ index ] * sys -> NN[ site_s ][ H_sys[ i ] ][ J_sys[ i ] ][ old_sys ] * env -> NN[ site_e ][ H_env[ i ] ][ J_env[ i ] ][ old_env ] / sqrt( ( sys -> Value_J_block[ H_sys[i] ][ J_sys[i] ] + 1.0 ) * ( env -> Value_J_block[ H_env[i] ][ J_env[i] ] + 1.0 ) );
                            
                             //------contributions from Interaction between ns and ne----------------------------
                                if( Table_J[ N * sys -> TotSiteNo + StartSite - env -> TotSiteNo ] == 1 )
                                        f[ Table_2to1[ i ][ j ] ] += Interaction_J[ N * sys -> TotSiteNo + StartSite - env -> TotSiteNo ] * nine_j_config_3_inverse[ i ][ index ] * nine_j_config_3_inverse[ i ][ index ] * six_j_2_config_3[ i_p ];
 
                                if( Table_N[ N * sys -> TotSiteNo + StartSite - env -> TotSiteNo ] == 1 )
                                        f[ Table_2to1[ i ][ j ] ] += Interaction_N[ N * sys -> TotSiteNo + StartSite - env -> TotSiteNo ] * nine_j_config_3_inverse[ i ][ index ] * nine_j_config_3_inverse[ i ][ index ];
                            
			}

			index++;

		}	

	}

}

//=======================Matrix-Vector Multiplication================================
void Super::H_V(const double *f, double *g, const int &dim) {

	//Initialize the vectors
        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ )
        for( int j = 0; j < Dim_block[i]; j++ ) {

                f1[i][j] = 0.0;   
		g1[i][j] = 0.0;

		f2[i][j] = 0.0;	
		g2[i][j] = 0.0;

        }

        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ )
        for( int j = 0; j < Dim_block_config_3[i]; j++ ) {

                f3[i][j] = 0.0;    
		g3[i][j] = 0.0;

        }

	//Read from f to f1 and initialize g to zero
        for(int i=0; i<dim; i++) {

                g[i] = 0.0;       

		f1[ Table_1to2_Num[i] ][ Table_1to2_Site[i] ] = f[i];

        }

	//Read from f1 to f2
        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) {

                index = 0;

                for( int j = 0; j < BlockNumber_for_TargetBlock; j++ )
		if( H_sys[i] == H_sys[j]  &&  H_env[i] == H_env[j]  &&  J_sys[i] == J_sys[j]  &&  J_env[i] == J_env[j]  &&  ( sysnew->Num_hole_block[H_sysnew[i]] - sys->Num_hole_block[H_sys[i]] ) == ( envnew->Num_hole_block[H_envnew[j]] - env->Num_hole_block[H_env[j]] )  &&  ( sysnew->Num_hole_block[H_sysnew[j]] - sys->Num_hole_block[H_sys[j]] ) == ( envnew->Num_hole_block[H_envnew[i]] - env->Num_hole_block[H_env[i]] ) )  {

                        daxpy_( &Dim_block[i], &nine_j_config_2[i][index++], f1[j], &inc, f2[i], &inc );

		}

        }

	//Read from f1 to f3
        for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ ) {

                index = 0;

                for( int j = 0; j < BlockNumber_for_TargetBlock; j++ )
		if( H_sys[j] == H_sys_config_3[i]  &&  H_env[j] == H_env_config_3[i]  &&  J_sys[j] == J_sys_config_3[i]  &&  J_env[j] == J_env_config_3[i]  &&  H_ns_config_3[i] == ( sysnew->Num_hole_block[H_sysnew[j]] - sys->Num_hole_block[H_sys[j]] )  &&  H_ne_config_3[i] == ( envnew->Num_hole_block[H_envnew[j]] - env->Num_hole_block[H_env[j]]) ) {

			daxpy_( &Dim_block_config_3[i], &nine_j_config_3[i][index++], f1[j], &inc, f3[i], &inc );

		}

        }

	//The first configuration: sys+ns+ne+env
	H_Sys_and_Env(); 
	H_Sys_Ns_T(); 
	H_Env_Ne_T();	  
	H_Sys_Ns_J();
	H_Env_Ne_J();
	H_Sys_Ns_N();
	H_Env_Ne_N();

	//The second configuration: sys+ne+ns+env
	H_Sys_Ne_T(); 
	H_Env_Ns_T();	  
	H_Sys_Ne_J();
	H_Env_Ns_J();
	H_Sys_Ne_N();
	H_Env_Ns_N();

	//Read from g2 to g1
        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) {

                index = 0;

                for( int j = 0; j < BlockNumber_for_TargetBlock; j++ )
		if( H_sys[i] == H_sys[j]  &&  H_env[i] == H_env[j]  &&  J_sys[i] == J_sys[j]  &&  J_env[i] == J_env[j]  &&  ( sysnew->Num_hole_block[H_sysnew[i]] - sys->Num_hole_block[H_sys[i]] ) == ( envnew->Num_hole_block[H_envnew[j]] - env->Num_hole_block[H_env[j]] )  &&  ( sysnew->Num_hole_block[H_sysnew[j]] - sys->Num_hole_block[H_sys[j]] ) == ( envnew->Num_hole_block[H_envnew[i]] - env->Num_hole_block[H_env[i]] ) ) {

			daxpy_(&Dim_block[i], &nine_j_config_2_inverse[i][index++], g2[j], &inc, g1[i], &inc);

		}

        }

	//The third configuration: sys+env+ns+ne
	H_Ns_Ne_T();
	H_Ns_Ne_J();
	H_Ns_Ne_N();	
	H_Sys_Env_T();
	H_Sys_Env_J();	
	H_Sys_Env_N();

	//Read from g3 to g1
	for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) {

                index = 0;

                for( int j = 0; j < BlockNumber_for_TargetBlock_config_3; j++ )
		if( H_sys[i] == H_sys_config_3[j]  &&  H_env[i] == H_env_config_3[j]  &&  J_sys[i] == J_sys_config_3[j]  &&  J_env[i] == J_env_config_3[j]  &&  H_ns_config_3[j] == ( sysnew->Num_hole_block[H_sysnew[i]] - sys->Num_hole_block[H_sys[i]] )  &&  H_ne_config_3[j] == ( envnew->Num_hole_block[H_envnew[i]] - env->Num_hole_block[H_env[i]] ) ) {

			daxpy_(&Dim_block[i], &nine_j_config_3_inverse[i][index++], g3[j], &inc, g1[i], &inc);

		}

        }

	//Read from g1 to g
        for( int i = 0; i < BlockNumber_for_TargetBlock; i++ )
        for( int j = 0; j < Dim_block[i]; j++ )

                g[ Table_2to1[i][j] ] = g1[i][j];

}

//=============================H Hamiltonian acting on the vector f===============================
//         g[a1,a2]=H_sys[a1,a1']*f[a1',a2]                g[a1,a2]=f[a1,a2']*H_env[a2',a2]
//================================================================================================
inline void Super::H_Sys_and_Env() {

        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) {

                dsymm_( &side_L, &uplo, &sys->Dim_J_block[ H_sys[l] ][ J_sys[l] ], &env->Dim_J_block[ H_env[l] ][ J_env[l] ], &alpha_p, sys->H[ H_sys[l] ][ J_sys[l] ], &sys->Dim_J_block[ H_sys[l] ][ J_sys[l] ], f1[l], &sys->Dim_J_block[ H_sys[l] ][ J_sys[l] ], &beta, g1[l], &sys->Dim_J_block[ H_sys[l] ][ J_sys[l] ] );

                dsymm_( &side_R, &uplo, &sys->Dim_J_block[ H_sys[l] ][ J_sys[l] ], &env->Dim_J_block[ H_env[l] ][ J_env[l] ], &alpha_p, env->H[ H_env[l] ][ J_env[l] ], &env->Dim_J_block[ H_env[l] ][ J_env[l] ], f1[l], &sys->Dim_J_block[ H_sys[l] ][ J_sys[l] ], &beta, g1[l], &sys->Dim_J_block[ H_sys[l] ][ J_sys[l] ] );

        }

}

//=================H between systems sites and ns acting on the vector f============================
inline void Super::H_Sys_Ns_T() {

  //----T_sys operator
        SiteNum = N * sys->TotSiteNo;

        for( int i = 0; i < sys->Block_Num_hole - 1; i++ ) {

	        for( int j = 0; j < sys->Num_block_T[i]; j++ ) {

                	dimension = sys->Dim_J_block[i][ sys->J_block_T_bra[i][j] ] * sys->Dim_J_block[ i + 1 ][ sys->J_block_T_ket[i][j] ];

                        for( int k = 0; k < dimension; k++ )

                        	T_sys[i][j][k] = 0.0;

                }

	}

        for( site_s = 0; site_s < operator_number_T_sys; site_s ++ )
        if( Table_T[ SiteNum + Table_T_sys[ site_s ] ] == 1 ) {

		for( int i = 0; i < sys->Block_Num_hole - 1; i++ )
		for( int j = 0; j < sys->Num_block_T[i]; j++ ) {

			alpha = - Interaction_T[ SiteNum + Table_T_sys[ site_s ] ] * hoping;

			dimension = sys->Dim_J_block[i][ sys->J_block_T_bra[i][j] ] * sys->Dim_J_block[ i + 1 ][ sys->J_block_T_ket[i][j] ];

			daxpy_( &dimension, &alpha, sys->T[site_s][i][j], &inc, T_sys[i][j], &inc );

		}

	}

  //----Calculations
        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ )  // g1 label
	if( sysnew->Num_hole_block[ H_sysnew[ l ] ] != 0  &&  sysnew->Num_hole_block[ H_sysnew[ l ] ] != sysnew->TotSiteNo ) {  // not either half-filling or empty

                for( int i = 0; i < BlockNumber_for_TargetBlock; i++ )	 // f1 label
		if( H_sysnew[ i ] == H_sysnew[ l ]  &&  H_envnew[ i ] == H_envnew[ l ]  &&  J_envnew[ i ] == J_envnew[ l ]  &&  J_sysnew[ i ] == J_sysnew[ l ]  &&  abs( H_sys[ i ] - H_sys[ l ] ) == 1  &&  H_env[ i ] == H_env[ l ]  &&  sysnew->Num_hole_block[ H_sysnew[ l ] ] - sys->Num_hole_block[ H_sys[ l ] ] + sysnew->Num_hole_block[ H_sysnew[ i ] ] - sys->Num_hole_block[ H_sys[ i ] ] == 1  &&  J_env[ i ] == J_env[ l ]  &&  abs( sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] - sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] ) == 1 ) {

			double six_j = 1.0 / sqrt( 2.0 * ( sysnew->Value_J_block[ H_sysnew[ l ] ][ J_sysnew[ l ] ] + 1.0 ) );

			//---- T^{\dagger}_i T_j
			if(  sys->Num_hole_block[ H_sys[ i ] ] - sys->Num_hole_block[ H_sys[ l ] ]  == 1 ) {

				//----
				alpha = six_j;

				//----
				if( ( sys->TotSiteNo - sys->Num_hole_block[ H_sys[ i ] ] ) % 2 == 1 )  

					alpha = - alpha;  //phase

				//----
                        	for( int old_i = 0; old_i < sys->Num_block_T[ H_sys[ l ] ]; old_i++ )
           		        if( sys->J_block_T_bra[ H_sys[ l ] ][ old_i ] == J_sys[ l ]  &&  sys->J_block_T_ket[ H_sys[ l ] ][ old_i ] == J_sys[ i ] )

                               		old_J = old_i;

				//----
				dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &alpha, T_sys[ H_sys[ l ] ][ old_J ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], f1[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

			}

			//---- T_i T^{\dagger}_j
			else if( sys->Num_hole_block[ H_sys[ l ] ] - sys->Num_hole_block[ H_sys[ i ] ] == 1 ) {

				//----
				alpha = six_j;

				if( ( sys->TotSiteNo - sys->Num_hole_block[ H_sys[ l ] ] ) % 2 == 1 )  

					alpha = - alpha;  //phase

				//----
	                       	for( int old_i = 0; old_i < sys->Num_block_T[ H_sys[ i ] ]; old_i++ )
               	  	        if( sys->J_block_T_bra[ H_sys[ i ] ][ old_i ] == J_sys[ i ]  &&  sys->J_block_T_ket[ H_sys[ i ] ][ old_i ] == J_sys[ l ] )

                             		old_J = old_i;
	
				//---- complex conjugate if the model is complex !!!
				dgemm_( &trans_T, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &alpha, T_sys[ H_sys[ i ] ][ old_J ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], f1[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

			}

		}

	}

}

//==============================================================================
inline void Super::H_Env_Ne_T() {

 //----T_env operator
	SiteNum = N * ( StartSite - env->TotSiteNo ) + StartSite;

        for( int i = 0; i < env->Block_Num_hole - 1; i++ ) {

	        for( int j = 0; j < env->Num_block_T[i]; j++ ) {

                	dimension = env->Dim_J_block[i][ env->J_block_T_bra[i][j] ] * env->Dim_J_block[ i + 1 ][ env->J_block_T_ket[i][j] ];

                        for( int k = 0; k < dimension; k++ )

                        	T_env[i][j][k] = 0.0;

                }

	}

        for( site_e = 0; site_e < operator_number_T_env; site_e++ )
        if( Table_T[ SiteNum - Table_T_env[ site_e ] ] == 1 ) {

		for( int i = 0; i < env->Block_Num_hole - 1; i++ )
		for( int j = 0; j < env->Num_block_T[i]; j++ ) {

			alpha = - Interaction_T[ SiteNum - Table_T_env[ site_e ] ] * hoping;

			dimension = env->Dim_J_block[i][ env->J_block_T_bra[i][j] ] * env->Dim_J_block[ i + 1 ][ env->J_block_T_ket[i][j] ];

			daxpy_( &dimension, &alpha, env -> T[site_e][i][j], &inc, T_env[i][j], &inc );

		}

	}

  //----Calculation
        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ )
	if( envnew->Num_hole_block[ H_envnew[ l ] ] != 0  &&  envnew->Num_hole_block[ H_envnew[ l ] ] != envnew->TotSiteNo ) { //g1

                for( int i = 0; i < BlockNumber_for_TargetBlock; i++ )	//"f1"
                if( H_sysnew[ i ] == H_sysnew[ l ]  &&  H_envnew[ i ] == H_envnew[ l ]  &&  H_sys[ i ] == H_sys[ l ]  &&  abs( H_env[ i ] - H_env[ l ] ) == 1  &&  envnew->Num_hole_block[ H_envnew[ l ] ] - env->Num_hole_block[ H_env[ l ] ] + envnew->Num_hole_block[ H_envnew[ i ] ] - env->Num_hole_block[ H_env[ i ] ] == 1  &&  J_envnew[ i ] == J_envnew[ l ]  &&  J_sys[ i ] == J_sys[ l ]  &&  J_sysnew[ i ] == J_sysnew[ l ]  &&  abs( env->Value_J_block[ H_env[ i ] ][ J_env[ i ] ] - env->Value_J_block[ H_env[ l ] ][ J_env[ l ] ] ) == 1 ) {

			//----
			double six_j = 1.0 / sqrt( 2.0 * (envnew->Value_J_block[ H_envnew[ l ] ][ J_envnew[ l ] ] + 1.0 ) );

			//----
			if( env->Num_hole_block[ H_env[ i ] ] - env->Num_hole_block[ H_env[ l ] ] == 1 ) { 

				//----
				alpha = six_j;

				//----
				if( ( env->TotSiteNo - env->Num_hole_block[ H_env[ i ] ] ) % 2 == 1 )  

					alpha = - alpha;  //phase

	                       	for( int old_i = 0; old_i < env->Num_block_T[ H_env[ l ] ]; old_i++ )
               	        	if( env->J_block_T_bra[ H_env[ l ] ][ old_i ] == J_env[ l ]  &&  env->J_block_T_ket[ H_env[ l ] ][ old_i ] == J_env[ i ] )

                              		old_J = old_i;

				//----
				dgemm_( &trans_N, &trans_T, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &env->Dim_J_block[ H_env[ i ] ][ J_env[ i ] ], &alpha, f1[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], T_env[ H_env[ l ] ][ old_J ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

			}

			else if( env->Num_hole_block[ H_env[ l ] ] - env->Num_hole_block[ H_env[ i ] ] == 1 ) {

				//----
				alpha = six_j;

				//----
				if( ( env->TotSiteNo - env->Num_hole_block[ H_env[ l ] ] ) % 2 == 1 )

					alpha = - alpha;//phase

				//----
	                       	for( int old_i = 0; old_i < env->Num_block_T[ H_env[ i ] ]; old_i++ )
               		       	if( env->J_block_T_bra[ H_env[ i ] ][ old_i ] == J_env[ i ]  &&  env->J_block_T_ket[ H_env[ i ] ][ old_i ] == J_env[ l ] )

                               		old_J = old_i;

				//---- complex conjugate of T_env if the Hamiltonian is complex
				dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &env->Dim_J_block[ H_env[ i ] ][ J_env[ i ] ], &alpha, f1[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], T_env[ H_env[ i ] ][ old_J ], &env->Dim_J_block[ H_env[ i ] ][ J_env[ i ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

			}

		}

	}

}

//========================================================================
inline void Super::H_Sys_Ns_J() {

	SiteNum = N * sys->TotSiteNo;  //The added site n_s with the number N*sys->TotSiteNo in Table

	//----
        for( int i = 0; i < sys->Block_Num_hole; i++ ) {

                for( int j = 0; j < sys->Num_J_hole_block[ i ]; j++ ) {

                	dimension = sys->Dim_J_block[ i ][ j ] * sys->Dim_J_block[ i ][ j ];

                        for( int k = 0; k < dimension; k++ )

                        	S_Dia_sys[ i ][ j ][ k ] = (double) 0;
	
		}

	}

        for( int i = 0; i < sys->Block_Num_hole; i++ ) {

                for( int j = 0; j < sys->Num_J_hole_block[ i ] - 1; j++ ) {

                	dimension = sys->Dim_J_block[ i ][ j ] * sys->Dim_J_block[ i ][ j + 1 ];

                        for( int k =0; k < dimension; k++ )

                        	S_M_Dia_sys[ i ][ j ][ k ] = (double) 0;

                }

	}

	//----
        for( site_s = 0; site_s < operator_number_J_sys; site_s++ )
        if( Table_J[ SiteNum + Table_J_sys[ site_s ] ] == 1 ) {

                for( int i = 0; i < sys->Block_Num_hole; i++ ) 
		for( int j = 0; j < sys->Num_J_hole_block[ i ]; j++ ) {

                        alpha = Interaction_J[ SiteNum + Table_J_sys[ site_s ] ];

                        dimension = sys->Dim_J_block[ i ][ j ] * sys->Dim_J_block[ i ][ j ];

                        daxpy_( &dimension, &alpha, sys->S_Dia[ site_s ][ i ][ j ], &inc, S_Dia_sys[ i ][ j ], &inc );

                }

        }

        for( site_s = 0; site_s < operator_number_J_sys; site_s++ )
        if( Table_J[ SiteNum + Table_J_sys[ site_s ] ] == 1 ) {

                for( int i = 0; i < sys->Block_Num_hole; i++ )
		for( int j = 0; j < sys->Num_J_hole_block[ i ] - 1; j++ ) {

                        alpha = Interaction_J[ SiteNum + Table_J_sys[ site_s ] ];

                        dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[ i ][ j + 1 ];

                        daxpy_( &dimension, &alpha, sys->S_M_Dia[ site_s ][ i ][ j ], &inc, S_M_Dia_sys[ i ][ j ], &inc );

                }

        }

	//----
        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( sys->TotSiteNo != sys->Num_hole_block[ H_sys[ l ] ] )  //with electron
	if( sys->Num_hole_block[ H_sys[ l ] ] == sysnew->Num_hole_block[ H_sysnew[ l ] ] ) { //no hole

		for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) 
		if( H_sysnew[ i ] == H_sysnew[ l ]  &&  H_envnew[ i ] == H_envnew[ l ]  &&  H_env[ i ] == H_env[ l ]  &&  H_sys[ i ] == H_sys[ l ]  &&  J_sysnew[ i ] == J_sysnew[ l ]  &&  J_envnew[ i ] == J_envnew[ l ]  &&  J_env[ i ] == J_env[ l ]  &&  abs( sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] - sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] ) <= 2 ) { //difference of Value[J_sys] must be 0 or 1 !

			//----
			if( J_sys[ i ] == J_sys[ l ]  &&  sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] != 0 ) {

                               	alpha = six_j_1_sys_n[ l * BlockNumber_for_TargetBlock + i ];

                                dsymm_( &side_L, &uplo, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &alpha, S_Dia_sys[ H_sys[ l ] ][ J_sys[ l ] ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], f1[ i ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

                        }

			//----
                        else if( J_sys[ i ] - J_sys[ l ] == 1  &&  sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] - sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] == 2 ) {

	                        alpha = six_j_1_sys_n[ l * BlockNumber_for_TargetBlock + i ];

                                dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &alpha, S_M_Dia_sys[ H_sys[ l ] ][ J_sys[ l ] ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], f1[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

                        }

			//----
                        else if( J_sys[ l ] - J_sys[ i ] == 1  &&  sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] - sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] == 2 ) {

                                alpha = - six_j_1_sys_n[ l * BlockNumber_for_TargetBlock + i ]; //"-" from transpose "S_M_Dia"

                                dgemm_( &trans_T, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &alpha, S_M_Dia_sys[ H_sys[ i ] ][ J_sys[ i ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], f1[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

                        }

		}

	}

}

//===================================================================================
inline void Super::H_Env_Ne_J() {

        SiteNum = N * ( StartSite - env->TotSiteNo ) + StartSite;

	//----
        for( int i = 0; i < env->Block_Num_hole; i++ ) {

                for( int j = 0; j < env->Num_J_hole_block[ i ]; j++ ) {

                	dimension = env -> Dim_J_block[ i ][ j ] * env -> Dim_J_block[ i ][ j ];

                        for( int k = 0; k < dimension; k++ )

                        	S_Dia_env[ i ][ j ][ k ] = (double) 0;
	
	
		}

	}

        for( int i = 0; i < env->Block_Num_hole; i++ ) {

                for( int j = 0; j < env->Num_J_hole_block[ i ] - 1; j++ ) {

                	dimension = env->Dim_J_block[ i ][ j ] * env->Dim_J_block[ i ][ j + 1 ];

                        for( int k =0; k < dimension; k++ )

                        	S_M_Dia_env[ i ][ j ][ k ] = (double) 0;

                }

	}

	//----
        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
        if( Table_J[ SiteNum - Table_J_env[ site_e ] ] == 1 ) {

                for( int i = 0; i < env->Block_Num_hole; i++ ) 
		for( int j = 0; j < env->Num_J_hole_block[ i ]; j++ ) {

                        alpha = Interaction_J[ SiteNum - Table_J_env[ site_e ] ];

                        dimension = env->Dim_J_block[ i ][ j ] * env->Dim_J_block[ i ][ j ];

                        daxpy_( &dimension, &alpha, env->S_Dia[ site_e ][ i ][ j ], &inc, S_Dia_env[ i ][ j ], &inc );

                }

        }

	//----
        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
        if( Table_J[ SiteNum - Table_J_env[ site_e ] ] == 1 ) {

                for( int i = 0; i < env->Block_Num_hole; i++ )
		for( int j = 0; j < env->Num_J_hole_block[ i ] - 1; j++ ) {

                        alpha = Interaction_J[ SiteNum - Table_J_env[ site_e ] ];

                        dimension = env->Dim_J_block[ i ][ j ] * env->Dim_J_block[ i ][ j + 1 ];

                        daxpy_( &dimension, &alpha, env->S_M_Dia[ site_e ][ i ][ j ], &inc, S_M_Dia_env[ i ][ j ], &inc );

                }

        }

	//----
	for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( env->TotSiteNo != env->Num_hole_block[ H_env[ l ] ] )   
	if( env->Num_hole_block[ H_env[ l ] ] == envnew->Num_hole_block[ H_envnew[ l ] ] ) {

		for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) 
		if( H_sysnew[ i ] == H_sysnew[ l ]  &&  H_envnew[ i ] == H_envnew[ l ]  &&  H_env[ i ] == H_env[ l ]  &&  H_sys[ i ] == H_sys[ l ]  &&  J_sysnew[ i ] == J_sysnew[ l ]  &&  J_envnew[ i ] == J_envnew[ l ]  &&  J_sys[ i ] == J_sys[ l ]  &&  abs( env->Value_J_block[ H_env[ i ] ][ J_env[ i ] ] - env->Value_J_block[ H_env[ l ] ][ J_env[ l ] ] ) <= 2 ) {

			//----
			if( J_env[ i ] == J_env[ l ]  &&  env->Value_J_block[ H_env[ i ] ][ J_env[ i ] ] != 0 )  {

	                        alpha = six_j_1_e_env[ l * BlockNumber_for_TargetBlock + i ];

                                dsymm_( &side_R, &uplo, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &alpha, S_Dia_env[ H_env[ l ] ][ J_env[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], f1[ i ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

                        }   

			//----
                        else if( J_env[i] - J_env[l] == 1  &&  env->Value_J_block[H_env[i]][J_env[i]] - env->Value_J_block[H_env[l]][J_env[l]] == 2 ) {

	                        alpha = six_j_1_e_env[l*BlockNumber_for_TargetBlock+i];

                                dgemm_(&trans_N, &trans_T, &sys->Dim_J_block[H_sys[l]][J_sys[l]], &env->Dim_J_block[H_env[l]][J_env[l]], &env->Dim_J_block[H_env[i]][J_env[i]], &alpha, f1[i], &sys->Dim_J_block[H_sys[l]][J_sys[l]], S_M_Dia_env[H_env[l]][J_env[l]], &env->Dim_J_block[H_env[l]][J_env[l]], &beta, g1[l], &sys->Dim_J_block[H_sys[l]][J_sys[l]]);

                        }

			//----
                        else if( J_env[l] - J_env[i] == 1  &&  env->Value_J_block[H_env[l]][J_env[l]] - env->Value_J_block[H_env[i]][J_env[i]] == 2 ) {

	                        alpha = -six_j_1_e_env[l*BlockNumber_for_TargetBlock+i];

                                dgemm_(&trans_N, &trans_N, &sys->Dim_J_block[H_sys[l]][J_sys[l]], &env->Dim_J_block[H_env[l]][J_env[l]], &env->Dim_J_block[H_env[i]][J_env[i]], &alpha, f1[i], &sys->Dim_J_block[H_sys[l]][J_sys[l]], S_M_Dia_env[H_env[i]][J_env[i]], &env->Dim_J_block[H_env[i]][J_env[i]], &beta, g1[l], &sys->Dim_J_block[H_sys[l]][J_sys[l]]);

                        }

		}

	}

}

//==========================================================================
inline void Super::H_Sys_Ns_N() {

        SiteNum = N * sys->TotSiteNo;       //The added site n_s with the number N*sys->TotSiteNo in Table

	//----
	for( int i = 0; i < sys->Block_Num_hole; i++ ) {

                for( int j = 0; j < sys->Num_J_hole_block[i]; j++ ) {

                	dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[i][j];

                        for( int k = 0; k < dimension; k++ )

                        	NN_sys[i][j][k] = (double) 0;
	
		}

	}

	//----
        for( site_s = 0; site_s < operator_number_J_sys; site_s++ )
        if( Table_N[ SiteNum + Table_N_sys[ site_s ] ] == 1 ) {

                for( int i = 0; i < sys->Block_Num_hole; i++ ) 
		for( int j = 0; j < sys->Num_J_hole_block[i]; j++ ) {

                        alpha = Interaction_N[ SiteNum + Table_N_sys[ site_s ] ];

                        dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[i][j];

                        daxpy_( &dimension, &alpha, sys->NN[site_s][i][j], &inc, NN_sys[i][j], &inc );

                }

        }

	//----
        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( sys->Num_hole_block[ H_sys[ l ] ] != sys->TotSiteNo  &&  sys->Num_hole_block[ H_sys[ l ] ] == sysnew->Num_hole_block[ H_sysnew[ l ] ] ) {

		alpha = 1.0 / sqrt( sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] + 1.0 );

		dsymm_( &side_L, &uplo, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &alpha, NN_sys[ H_sys[ l ] ][ J_sys[ l ] ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], f1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

	}

}

//=============================================================================================
inline void Super::H_Env_Ne_N() {

        SiteNum = N * ( StartSite - env->TotSiteNo ) + StartSite;

	//----
	for( int i = 0; i < env->Block_Num_hole; i++ ) {

                for( int j = 0; j < env->Num_J_hole_block[i]; j++ ) {

                	dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j];

                        for( int k = 0; k < dimension; k++ )

                        	NN_env[i][j][k] = (double) 0;
	
		}

	}

	//----
        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
        if( Table_N[ SiteNum - Table_N_env[ site_e ] ] == 1 ) {

                for( int i = 0; i < env->Block_Num_hole; i++ ) 
		for( int j = 0; j < env->Num_J_hole_block[i]; j++ ) {

                        alpha = Interaction_N[ SiteNum - Table_N_env[ site_e ] ];

                        dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j];

                        daxpy_( &dimension, &alpha, env->NN[site_e][i][j], &inc, NN_env[i][j], &inc );

                }

        }

	//----
	for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( env->TotSiteNo != env->Num_hole_block[ H_env[ l ] ]  &&  env->Num_hole_block[ H_env[ l ] ] == envnew->Num_hole_block[ H_envnew[ l ] ] ) {

		alpha = 1.0 / sqrt( env->Value_J_block[ H_env[ l ] ][ J_env[ l ] ] + 1.0 );

		dsymm_( &side_R, &uplo, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &alpha, NN_env[ H_env[ l ] ][ J_env[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], f1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &beta, g1[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

	}

}

//=========================================================================================
inline void Super::H_Sys_Ne_T() {

  //----T_sys operator
	SiteNum = N * ( StartSite - env->TotSiteNo );

        for( int i = 0; i < sys->Block_Num_hole - 1; i++ ) {

	        for( int j = 0; j < sys->Num_block_T[i]; j++ ) {

                	dimension = sys->Dim_J_block[i][ sys->J_block_T_bra[i][j] ] * sys->Dim_J_block[ i + 1 ][ sys->J_block_T_ket[i][j] ];

                        for( int k = 0; k < dimension; k++ )

                        	T_sys[i][j][k] = 0.0;

                }

	}

        for( site_s = 0; site_s < operator_number_T_sys; site_s ++ )
        if( Table_T[ SiteNum + Table_T_sys[ site_s ] ] == 1 ) {

		for( int i = 0; i < sys->Block_Num_hole - 1; i++ )
		for( int j = 0; j < sys->Num_block_T[i]; j++ ) {

			alpha = - Interaction_T[ SiteNum + Table_T_sys[ site_s ] ] * hoping;

			dimension = sys->Dim_J_block[i][ sys->J_block_T_bra[i][j] ] * sys->Dim_J_block[ i + 1 ][ sys->J_block_T_ket[i][j] ];

			daxpy_( &dimension, &alpha, sys->T[site_s][i][j], &inc, T_sys[i][j], &inc );

		}

	}

  //----calculation
        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( sysnew->Num_hole_block[ H_sysnew[ l ] ] != 0  &&  sysnew->Num_hole_block[ H_sysnew[ l ] ] != sysnew->TotSiteNo ) {//g2

                for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) //f2
                if( H_sysnew[ i ] == H_sysnew[ l ]  &&  H_envnew[ i ] == H_envnew[ l ]  &&  J_envnew[ i ] == J_envnew[ l ]  &&  J_sysnew[ i ] == J_sysnew[ l ]  &&  H_env[ i ] == H_env [ l ]  &&  abs( H_sys[ i ] - H_sys[ l ] ) == 1  && sysnew->Num_hole_block[ H_sysnew[ l ] ] - sys->Num_hole_block[ H_sys[ l ] ] + sysnew->Num_hole_block[ H_sysnew[ i ] ] - sys->Num_hole_block[ H_sys[ i ] ] == 1  &&  J_env[ i ] == J_env[ l ]  &&  abs( sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] - sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] ) == 1 ) {

			double six_j = 1.0 / sqrt( 2.0 * ( sysnew->Value_J_block[ H_sysnew[ l ] ][ J_sysnew[ l ] ] + 1.0 ) );
			
                        //---- T^{\dagger}_i T_j
			if( sys->Num_hole_block[ H_sys[ i ] ] - sys->Num_hole_block[ H_sys[ l ] ] == 1 ) { 

				//----
				alpha = six_j;

				//----
				if( ( sys->TotSiteNo - sys->Num_hole_block[ H_sys[ i ] ] ) % 2 == 1 )  

					alpha = - alpha;  //phase

				//----
	                       	for( int old_i = 0; old_i < sys->Num_block_T[ H_sys[ l ] ]; old_i++ )
               	        	if( sys->J_block_T_bra[ H_sys[ l ] ][ old_i ] == J_sys[ l ]  &&  sys->J_block_T_ket[ H_sys[ l ] ][ old_i ] == J_sys[ i ] )

                            		old_J = old_i;

				//----
				dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &alpha, T_sys[ H_sys[ l ] ][ old_J ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], f2[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

			}

                        //---- T_i T^{\dagger}_j
			else if( sys->Num_hole_block[ H_sys[ l ] ] - sys->Num_hole_block[ H_sys[ i ] ] == 1 ) {

				//----
				alpha = six_j;

				if( ( sys->TotSiteNo - sys->Num_hole_block[ H_sys[ l ] ] ) % 2 == 1 )  

					alpha = - alpha;	//phase

				//----
	                        for( int old_i = 0; old_i < sys->Num_block_T[ H_sys[ i ] ]; old_i++ )
               		  	if( sys->J_block_T_bra[ H_sys[ i ] ][ old_i ] == J_sys[ i ]  &&  sys->J_block_T_ket[ H_sys[ i ] ][ old_i ] == J_sys[ l ] )

                               		old_J = old_i;

                                //---- complex conjugate if the model is complex !!!
				dgemm_( &trans_T, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &alpha, T_sys[ H_sys[ i ] ][ old_J ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], f2[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

			}

		}

	}

}

//=====================================================================================
inline void Super::H_Env_Ns_T() {

 //----T_env operator
        SiteNum = N * sys->TotSiteNo + StartSite;

        for( int i = 0; i < env->Block_Num_hole - 1; i++ ) {

	        for( int j = 0; j < env->Num_block_T[i]; j++ ) {

                	dimension = env->Dim_J_block[i][ env->J_block_T_bra[i][j] ] * env->Dim_J_block[ i + 1 ][ env->J_block_T_ket[i][j] ];

                        for( int k = 0; k < dimension; k++ )

                        	T_env[i][j][k] = 0.0;

                }

	}

        for( site_e = 0; site_e < operator_number_T_env; site_e++ )
        if( Table_T[ SiteNum - Table_T_env[ site_e ] ] == 1 ) {

		for( int i = 0; i < env->Block_Num_hole - 1; i++ )
		for( int j = 0; j < env->Num_block_T[i]; j++ ) {

			alpha = - Interaction_T[ SiteNum - Table_T_env[ site_e ] ] * hoping;

			dimension = env->Dim_J_block[i][ env->J_block_T_bra[i][j] ] * env->Dim_J_block[ i + 1 ][ env->J_block_T_ket[i][j] ];

			daxpy_( &dimension, &alpha, env -> T[site_e][i][j], &inc, T_env[i][j], &inc );

		}

	}

  //----calculation
        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( envnew->Num_hole_block[ H_envnew[ l ] ] != 0  &&  envnew->Num_hole_block[ H_envnew[ l ] ] != envnew->TotSiteNo ) { //g2

                for( int i = 0; i < BlockNumber_for_TargetBlock; i++ )  //"f2"
		if( H_sysnew[ i ] == H_sysnew[ l ]  &&  H_envnew[ i ] == H_envnew[ l ]  &&  H_sys[ i ] == H_sys[ l ]  &&  abs( H_env[ i ] - H_env[ l ] ) == 1  &&  envnew->Num_hole_block[ H_envnew[ l ] ] - env->Num_hole_block[ H_env[ l ] ] + envnew->Num_hole_block[ H_envnew[ i ] ] - env->Num_hole_block[ H_env[ i ] ] == 1  &&  J_envnew[ i ] == J_envnew[ l ]  &&  J_sys[ i ] == J_sys[ l ]  &&  J_sysnew[ i ] == J_sysnew[ l ]  &&  abs( env->Value_J_block[ H_env[ i ] ][ J_env[ i ] ] - env->Value_J_block[ H_env[ l ] ][ J_env[ l ] ] ) == 1 ) {

			//----
			double six_j = 1.0 / sqrt( 2.0 * ( envnew->Value_J_block[ H_envnew[ l ] ][ J_envnew[ l ] ] + 1.0 ) );

			//----
			if( env->Num_hole_block[ H_env[ i ] ] - env->Num_hole_block[ H_env[ l ] ] == 1 ) {

				//----
				alpha = six_j;

				//----
				if( ( env->TotSiteNo - env->Num_hole_block[ H_env[ i ] ] ) % 2 == 1 )  

					alpha = - alpha;  //phase

                        	for( int old_i = 0; old_i < env->Num_block_T[ H_env[ l ] ]; old_i++ )
            		        if( env->J_block_T_bra[ H_env[ l ] ][ old_i ] == J_env[ l ]  &&  env->J_block_T_ket[ H_env[ l ] ][ old_i ] == J_env[ i ] )

                              		old_J = old_i;

				//----
				dgemm_( &trans_N, &trans_T, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &env->Dim_J_block[ H_env[ i ] ][ J_env[ i ] ], &alpha, f2[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], T_env[ H_env[ l ] ][ old_J ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

			}

			//----
			else if( env->Num_hole_block[ H_env[ l ] ] - env->Num_hole_block[ H_env[ i ] ] == 1 ) {

				//----
				alpha = six_j;

				//----
				if( ( env->TotSiteNo - env->Num_hole_block[ H_env[ l ] ] ) % 2 == 1 )  

					alpha = - alpha;  //phase

				//----
	                       	for( int old_i = 0; old_i < env->Num_block_T[ H_env[ i ] ]; old_i++ )
               		       	if( env->J_block_T_bra[ H_env[ i ] ][ old_i ] == J_env[ i ]  &&  env->J_block_T_ket[ H_env[ i ] ][ old_i ] == J_env[ l ] )

                               		old_J = old_i;

                                //---- complex conjugate of T_env if the Hamiltonian is complex
				dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &env->Dim_J_block[ H_env[ i ] ][ J_env[ i ] ], &alpha, f2[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], T_env[ H_env[ i ] ][ old_J ], &env->Dim_J_block[ H_env[ i ] ][ J_env[ i ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );  //sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ]

			}

		}

	}

}

//=======================================================================================
inline void Super::H_Sys_Ne_J() {

        SiteNum = N * ( StartSite - env->TotSiteNo );

	//----
        for( int i = 0; i < sys->Block_Num_hole; i++ ) {

                for( int j = 0; j < sys->Num_J_hole_block[i]; j++ ) {

                	dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[i][j];

                        for( int k = 0; k < dimension; k++ )

                        	S_Dia_sys[i][j][k] = (double) 0;
	
		}

	}

        for( int i = 0; i < sys->Block_Num_hole; i++ ) {

                for( int j = 0; j < sys->Num_J_hole_block[i] - 1; j++ ) {

                	dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[i][j+1];

                        for( int k =0; k < dimension; k++ )

                        	S_M_Dia_sys[i][j][k] = (double) 0;

                }

	}

	//----
        for( site_s = 0; site_s < operator_number_J_sys; site_s++ )
        if( Table_J[ SiteNum + Table_J_sys[ site_s ] ] == 1 ) {

                for( int i = 0; i < sys->Block_Num_hole; i++ ) 
		for( int j = 0; j < sys->Num_J_hole_block[i]; j++ ) {

                        alpha = Interaction_J[ SiteNum + Table_J_sys[ site_s ] ];

                        dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[i][j];

                        daxpy_( &dimension, &alpha, sys->S_Dia[site_s][i][j], &inc, S_Dia_sys[i][j], &inc );

                }

        }

        for( site_s = 0; site_s < operator_number_J_sys; site_s++ )
        if( Table_J[ SiteNum + Table_J_sys[site_s] ] == 1 ) {

                for( int i = 0; i < sys->Block_Num_hole; i++ )
		for( int j = 0; j < sys->Num_J_hole_block[i] - 1; j++ ) {

                        alpha = Interaction_J[ SiteNum + Table_J_sys[site_s] ];

                        dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[i][j+1];

                        daxpy_( &dimension, &alpha, sys->S_M_Dia[site_s][i][j], &inc, S_M_Dia_sys[i][j], &inc );

                }

        }

	//----
        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( sys->TotSiteNo != sys->Num_hole_block[ H_sys[ l ] ] )
	if( sys->Num_hole_block[ H_sys[ l ] ] == sysnew->Num_hole_block[ H_sysnew[ l ] ] ) {

		for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) 
		if( H_sysnew[ i ] == H_sysnew[ l ]  &&  H_envnew[ i ] == H_envnew[ l ]  &&  H_env[ i ] == H_env[ l ]  &&  H_sys[ i ] == H_sys[ l ]  &&  J_sysnew[ i ] == J_sysnew[ l ]  &&  J_envnew[ i ] == J_envnew[ l ]  &&  J_env[ i ] == J_env[ l ]  &&  abs( sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] - sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] ) <= 2 ) {

			//----
			if( J_sys[ i ] == J_sys[ l ]  &&  sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] != 0 ) {

                               	alpha = six_j_1_sys_n[ l * BlockNumber_for_TargetBlock + i ];

                                dsymm_( &side_L, &uplo, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &alpha, S_Dia_sys[ H_sys[ l ] ][ J_sys[ l ] ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], f2[ i ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

                        }

			//----
                        else if( J_sys[ i ] - J_sys[ l ] == 1  &&  sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] - sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] == 2 ) {

                                alpha = six_j_1_sys_n[ l * BlockNumber_for_TargetBlock + i ];

                                dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &alpha, S_M_Dia_sys[ H_sys[ l ] ][ J_sys[ l ] ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], f2[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

                        }

			//----
                        else if( J_sys[ l ] - J_sys[ i ] == 1  &&  sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] - sys->Value_J_block[ H_sys[ i ] ][ J_sys[ i ] ] == 2 ) {

	                        alpha = -six_j_1_sys_n[ l * BlockNumber_for_TargetBlock + i ];//"-" from transpose of "S_M_Dia"

	                        dgemm_( &trans_T, &trans_N, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &alpha, S_M_Dia_sys[ H_sys[ i ] ][ J_sys[ i ] ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], f2[ i ], &sys->Dim_J_block[ H_sys[ i ] ][ J_sys[ i ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

			}

		}

	}

}

//======================================================================================
inline void Super::H_Env_Ns_J() {

        SiteNum = N * sys->TotSiteNo + StartSite;

	//----
        for( int i = 0; i < env->Block_Num_hole; i++ ) {

                for( int j = 0; j < env->Num_J_hole_block[i]; j++ ) {

                	dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j];

                        for( int k = 0; k < dimension; k++ )

                        	S_Dia_env[i][j][k] = (double) 0;
	
		}

	}

        for( int i = 0; i < env->Block_Num_hole; i++ ) {

                for( int j = 0; j < env->Num_J_hole_block[i] - 1; j++ ) {

                	dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j+1];

                        for( int k =0; k < dimension; k++ )

                        	S_M_Dia_env[i][j][k] = (double) 0;

                }

	}
	
	//----
        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
        if( Table_J[ SiteNum - Table_J_env[ site_e ] ] == 1 ) {

                for( int i = 0; i < env->Block_Num_hole; i++ ) 
		for( int j = 0; j < env->Num_J_hole_block[i]; j++ ) {

                        alpha = Interaction_J[ SiteNum - Table_J_env[site_e] ];

                        dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j];

                        daxpy_( &dimension, &alpha, env->S_Dia[site_e][i][j], &inc, S_Dia_env[i][j], &inc );

                }

        }

	//----
        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
        if( Table_J[ SiteNum - Table_J_env[site_e] ] == 1 ) {

                for( int i = 0; i < env->Block_Num_hole; i++ )
		for( int j = 0; j < env->Num_J_hole_block[i] - 1; j++ ) {

                        alpha = Interaction_J[ SiteNum - Table_J_env[site_e] ];

                        dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j+1];

                        daxpy_( &dimension, &alpha, env->S_M_Dia[site_e][i][j], &inc, S_M_Dia_env[i][j], &inc );

                }

        }

	//----
	for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( env->TotSiteNo != env->Num_hole_block[ H_env[ l ] ] )
	if( env->Num_hole_block[ H_env[ l ] ] == envnew->Num_hole_block[ H_envnew[ l ] ] ) {

		for( int i = 0; i < BlockNumber_for_TargetBlock; i++ ) 
		if( H_sysnew[ i ] == H_sysnew[ l ]  &&  H_envnew[ i ] == H_envnew[ l ]  &&  H_env[ i ] == H_env[ l ]  &&  H_sys[ i ] == H_sys[ l ]  &&  J_sysnew[ i ] == J_sysnew[ l ]  &&  J_envnew[ i ] == J_envnew[ l ]  &&  J_sys[ i ] == J_sys[ l ]  &&  abs( env->Value_J_block[ H_env[ i ] ][ J_env[ i ] ] - env->Value_J_block[ H_env[ l ] ][ J_env[ l ] ] ) <= 2 ) {

			//----
			if( J_env[ i ] == J_env[ l ]  &&  env->Value_J_block[ H_env[ i ] ][ J_env[ i ] ] != 0 )  {

	                        alpha = six_j_1_e_env[ l * BlockNumber_for_TargetBlock + i ];

                                dsymm_( &side_R, &uplo, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &alpha, S_Dia_env[ H_env[ l ] ][ J_env[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], f2[ i ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

                        }   
	
			//----
                        else if( J_env[ i ] - J_env[ l ] == 1  &&  env->Value_J_block[ H_env[ i ] ][ J_env[ i ] ] - env->Value_J_block[ H_env[ l ] ][ J_env[ l ] ] == 2 ) {

	                        alpha = six_j_1_e_env[l*BlockNumber_for_TargetBlock+i];

                                dgemm_(&trans_N, &trans_T, &sys->Dim_J_block[H_sys[l]][J_sys[l]], &env->Dim_J_block[H_env[l]][J_env[l]], &env->Dim_J_block[H_env[i]][J_env[i]], &alpha, f2[i], &sys->Dim_J_block[H_sys[l]][J_sys[l]], S_M_Dia_env[H_env[l]][J_env[l]], &env->Dim_J_block[H_env[l]][J_env[l]], &beta, g2[l], &sys->Dim_J_block[H_sys[l]][J_sys[l]]);

                        }

			//----
                        else if( J_env[l] - J_env[i] == 1  &&  env->Value_J_block[H_env[l]][J_env[l]] - env->Value_J_block[H_env[i]][J_env[i]] == 2 ) {

	                        alpha = -six_j_1_e_env[l*BlockNumber_for_TargetBlock+i];

                                dgemm_(&trans_N, &trans_N, &sys->Dim_J_block[H_sys[l]][J_sys[l]], &env->Dim_J_block[H_env[l]][J_env[l]], &env->Dim_J_block[H_env[i]][J_env[i]], &alpha, f2[i], &sys->Dim_J_block[H_sys[l]][J_sys[l]], S_M_Dia_env[H_env[i]][J_env[i]], &env->Dim_J_block[H_env[i]][J_env[i]], &beta, g2[l], &sys->Dim_J_block[H_sys[l]][J_sys[l]]);

			}

		}

	}

}

//======================================================================================
inline void Super::H_Sys_Ne_N() {

        SiteNum = N * ( StartSite - env->TotSiteNo );

	//----
	for( int i = 0; i < sys->Block_Num_hole; i++ ) {

                for( int j = 0; j < sys->Num_J_hole_block[i]; j++ ) {

                	dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[i][j];

                        for( int k = 0; k < dimension; k++ )

                        	NN_sys[i][j][k] = (double) 0;
	
		}

	}

	//----
        for( site_s = 0; site_s < operator_number_J_sys; site_s++ )
        if( Table_N[ SiteNum + Table_N_sys[site_s] ] == 1 ) {

                for( int i = 0; i < sys->Block_Num_hole; i++ ) 
		for( int j = 0; j < sys->Num_J_hole_block[i]; j++ ) {

                        alpha = Interaction_N[ SiteNum + Table_N_sys[site_s] ];

                        dimension = sys->Dim_J_block[i][j] * sys->Dim_J_block[i][j];

                        daxpy_( &dimension, &alpha, sys->NN[site_s][i][j], &inc, NN_sys[i][j], &inc );

                }

        }

	//----
        for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( sys->TotSiteNo != sys->Num_hole_block[ H_sys[ l ] ]  &&  sys->Num_hole_block[ H_sys[ l ] ] == sysnew->Num_hole_block[ H_sysnew[ l ] ] ) {

		alpha = 1.0 / sqrt( sys->Value_J_block[ H_sys[ l ] ][ J_sys[ l ] ] + 1.0 );

		dsymm_( &side_L, &uplo, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &alpha, NN_sys[ H_sys[ l ] ][ J_sys[ l ] ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], f2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

	}

}

//=============================================================================================
inline void Super::H_Env_Ns_N() {

        SiteNum = N * sys->TotSiteNo + StartSite;

	//----
	for( int i = 0; i < env->Block_Num_hole; i++ ) {

                for( int j = 0; j < env->Num_J_hole_block[i]; j++ ) {

                	dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j];

                        for( int k = 0; k < dimension; k++ )

                        	NN_env[i][j][k] = (double) 0;
	
		}

	}

	//----
        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
        if( Table_N[ SiteNum - Table_N_env[site_e] ] == 1 ) {

                for( int i = 0; i < env->Block_Num_hole; i++ ) 
		for( int j = 0; j < env->Num_J_hole_block[i]; j++ ) {

                        alpha = Interaction_N[ SiteNum - Table_N_env[site_e] ];

                        dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j];

                        daxpy_( &dimension, &alpha, env->NN[site_e][i][j], &inc, NN_env[i][j], &inc );

                }

        }

	//----
	for( int l = 0; l < BlockNumber_for_TargetBlock; l++ ) 
	if( env->TotSiteNo != env->Num_hole_block[ H_env[ l ] ]  &&  env->Num_hole_block[ H_env[ l ] ] == envnew->Num_hole_block[ H_envnew[ l ] ] ) {

		alpha = 1.0 / sqrt( env->Value_J_block[ H_env[ l ] ][ J_env[ l ] ] + 1.0 );

		dsymm_( &side_R, &uplo, &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], &alpha, NN_env[ H_env[ l ] ][ J_env[ l ] ], &env->Dim_J_block[ H_env[ l ] ][ J_env[ l ] ], f2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ], &beta, g2[ l ], &sys->Dim_J_block[ H_sys[ l ] ][ J_sys[ l ] ] );

	}

}

//===============================================================================
inline void Super::H_Ns_Ne_T() {

	alpha = - Interaction_T[ N * sys->TotSiteNo + StartSite - env->TotSiteNo ];
	
        if( Table_T[ N * sys->TotSiteNo + StartSite - env->TotSiteNo ] == 1 ) {

	        for( int l = 0; l < BlockNumber_for_TargetBlock_config_3; l++ ) 
		if( H_ns_config_3[ l ] + H_ne_config_3[ l ]  == 1 ) { //g3

	                for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ ) //f3
                        if( H_sys_config_3[ i ] == H_sys_config_3[ l ]  &&  H_env_config_3[ i ] == H_env_config_3[ l ]  &&   H_ns_config_3[ i ] + H_ns_config_3[ l ] == 1  &&  H_ne_config_3[ i ] + H_ne_config_3[ l ] == 1  &&  J_sys_config_3[ i ] == J_sys_config_3[ l ]  &&  J_env_config_3[ i ] == J_env_config_3[ l ]  &&  J_sysnew_config_3[ i ] == J_sysnew_config_3[ l ] ) {	//constraints

 				daxpy_( &Dim_block_config_3[ l ], &alpha, f3[ i ], &inc, g3[ l ], &inc );

			}

		}

	}

}

//===============================================================================
inline void Super::H_Ns_Ne_J() {

        if( Table_J[ N * sys->TotSiteNo + StartSite - env->TotSiteNo ] == 1 ) {

	        for( int l = 0; l < BlockNumber_for_TargetBlock_config_3; l++ ) 
		if( H_ns_config_3[ l ] == 0  &&  H_ne_config_3[ l ] == 0 ) {

			alpha = six_j_2_config_3[ l ] * Interaction_J[ N * sys->TotSiteNo + StartSite - env->TotSiteNo ];

			daxpy_( &Dim_block_config_3[ l ], &alpha, f3[ l ], &inc, g3[ l ], &inc);

		}

	}

}

//=======================================================================================
inline void Super::H_Ns_Ne_N() {

        if( Table_N[ N * sys->TotSiteNo + StartSite - env->TotSiteNo ] == 1 ) {

	        for( int l = 0; l < BlockNumber_for_TargetBlock_config_3; l++ ) 
		if( H_ns_config_3[ l ] == 0  &&  H_ne_config_3[ l ] == 0 ) {

			alpha = Interaction_N[ N * sys->TotSiteNo + StartSite - env->TotSiteNo ];

			daxpy_( &Dim_block_config_3[ l ], &alpha, f3[ l ], &inc, g3[ l ], &inc );

		}

	}

}

//=========================================================================================
inline void Super::H_Sys_Env_T() {

	int mul;

        for( site_s = 0; site_s < operator_number_T_sys; site_s ++ ) {

		//----
	        for( int i = 0; i < env -> Block_Num_hole - 1; i++ ) {

		        for( int j = 0; j < env -> Num_block_T[ i ]; j++ ) {

               			dimension = env -> Dim_J_block[ i ][ env -> J_block_T_bra[ i ][ j ] ] * env -> Dim_J_block[ i + 1 ][ env -> J_block_T_ket[ i ][ j ] ];

                       		for( int k = 0; k < dimension; k++ ) {

               	        		T_env[ i ][ j ][ k ] = 0.0;

				}

                	}

		}

	        for( site_e = 0; site_e < operator_number_T_env; site_e ++ )
		if( Table_T[ N * Table_T_sys[ site_s ] + StartSite - Table_T_env[ site_e ] ] == 1 ) {

			for( int i = 0; i < env -> Block_Num_hole - 1; i++ )
			for( int j = 0; j < env -> Num_block_T[i]; j++ ) {

				alpha = Interaction_T[ N * Table_T_sys[ site_s ] + StartSite - Table_T_env[ site_e ] ];

				dimension = env -> Dim_J_block[ i ][ env -> J_block_T_bra[ i ][ j ] ] * env -> Dim_J_block[ i + 1 ][ env -> J_block_T_ket[ i ][ j ] ];

				daxpy_( &dimension, &alpha, env -> T[ site_e ][ i ][ j ], &inc, T_env[ i ][ j ], &inc );

			}

		}

		//----
                for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ )	// g vector
		if( ( ( sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] ) != 0 )  &&  ( ( sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] ) != ( sys -> TotSiteNo + env -> TotSiteNo ) ) ) {

			index = 0;

			for( int j = 0; j < BlockNumber_for_TargetBlock_config_3; j ++ )
			if( ( sys -> Num_hole_block[ H_sys_config_3[ j ] ] + env -> Num_hole_block[ H_env_config_3[ j ] ]  ==  sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] )  &&  H_ns_config_3[ j ] == H_ns_config_3[ i ]  &&  H_ne_config_3[ j ] == H_ne_config_3[ i ]  &&  abs( sys -> Num_hole_block[ H_sys_config_3[ i ] ] - sys -> Num_hole_block[ H_sys_config_3[ j ] ] ) == 1  &&  abs( env -> Num_hole_block[ H_env_config_3[ i ] ] - env -> Num_hole_block[ H_env_config_3[ j ] ] ) == 1  &&  J_sysnew_config_3[ i ] == J_sysnew_config_3[ j ]  &&  J_envnew_config_3[ i ] == J_envnew_config_3[ j ]  &&  abs( sys -> Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] - sys -> Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] ) == 1  &&  abs( env -> Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] - env -> Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ] ) == 1 ) {

				if( ( sys -> Num_hole_block[ H_sys_config_3[ j ] ] - sys -> Num_hole_block[ H_sys_config_3[ i ] ] ) == 1  &&  ( env -> Num_hole_block[ H_env_config_3[ i ] ] - env -> Num_hole_block[ H_env_config_3[ j ] ] ) == 1 ) {

					alpha = -six_j_T_config_3[ i ][ index ];	
	
					if( ( sys -> TotSiteNo - sys -> Num_hole_block[ H_sys_config_3[ i ] ] ) % 2 == 1 )  
						alpha = -alpha;
					
					if( ( env -> Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] - env -> Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ] ) == 1 )  
						alpha = -alpha;

					mul = sys -> Dim_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] * env -> Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ];

					double * f_sub = new double [ mul ];
					for( int s = 0; s < mul; s ++ )  f_sub[ s ] = 0.0;

	                        	for( int old_i = 0; old_i < env -> Num_block_T[ H_env_config_3[ j ] ]; old_i ++ )
               		        	if( env -> J_block_T_bra[ H_env_config_3[ j ] ][ old_i ] == J_env_config_3[ j ]  &&  env -> J_block_T_ket[ H_env_config_3[ j ] ][ old_i ] == J_env_config_3[ i ] )
                               			old_J = old_i;

					dgemm_( &trans_N, &trans_N, &sys -> Dim_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ], &env -> Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &env -> Dim_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ], &alpha_p, f3[ j ], &sys -> Dim_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ], T_env[ H_env_config_3[ j ] ][ old_J ], &env -> Dim_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ], &beta_p, f_sub, &sys -> Dim_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] );

	                        	for( int old_i = 0; old_i < sys->Num_block_T[ H_sys_config_3[ i ] ]; old_i ++ )
               		        	if( sys -> J_block_T_bra[ H_sys_config_3[ i ] ][ old_i ] == J_sys_config_3[ i ]  &&  sys -> J_block_T_ket[ H_sys_config_3[ i ] ][ old_i ] == J_sys_config_3[ j ] )
                               			old_J = old_i;

					dgemm_( &trans_N, &trans_N, &sys -> Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &env -> Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &sys -> Dim_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ], &alpha, sys -> T[ site_s ][ H_sys_config_3[ i ] ][ old_J ], &sys -> Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], f_sub, &sys -> Dim_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ], &beta, g3[ i ], &sys -> Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] );

					delete [] f_sub;

				}

				else if( ( sys -> Num_hole_block[ H_sys_config_3[ i ] ] - sys -> Num_hole_block[ H_sys_config_3[ j ] ] ) == 1  &&  ( env -> Num_hole_block[ H_env_config_3[ j ] ] - env -> Num_hole_block[ H_env_config_3[ i ] ] ) == 1 ) {
	
					alpha = six_j_T_config_3[ i ][ index ];

					if( ( sys -> TotSiteNo - sys -> Num_hole_block[ H_sys_config_3[ j ] ] ) % 2 == 1 )  
						alpha = -alpha;
					if( ( sys -> Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] - sys -> Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] ) == 1 )  
						alpha = -alpha;

					mul = sys -> Dim_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] * env -> Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ];
					double *f_sub=new double [mul];
					for(int s=0; s<mul; s++)	f_sub[s]=0.0;

	                        	for( int old_i = 0; old_i < env -> Num_block_T[ H_env_config_3[ i ] ]; old_i ++ )
               		        	if( env -> J_block_T_bra[ H_env_config_3[ i ] ][ old_i ] == J_env_config_3[ i ]  &&  env -> J_block_T_ket[ H_env_config_3[ i ] ][ old_i ] == J_env_config_3[ j ] ) 
                               			old_J=old_i;

					dgemm_(&trans_N, &trans_T, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &alpha_p, f3[j], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], T_env[H_env_config_3[i]][old_J], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]]);

	                        	for( int old_i = 0; old_i < sys -> Num_block_T[ H_sys_config_3[ j ] ]; old_i ++ )
               		        	if( sys -> J_block_T_bra[ H_sys_config_3[ j ] ][ old_i ] == J_sys_config_3[ j ]  &&  sys -> J_block_T_ket[ H_sys_config_3[ j ] ][ old_i ] == J_sys_config_3[ i ] )
                               			old_J=old_i;

					dgemm_(&trans_T, &trans_N, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &alpha, sys->T[site_s][H_sys_config_3[j]][old_J], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]]);

					delete [] f_sub;

				}
	
				index++;

			}

		}

	}

}

//=======================================================================================
inline void Super::H_Sys_Env_J() {

        int mul;

        for( site_s = 0; site_s < operator_number_J_sys; site_s++ ) { //we add env operators together

		//----S_Dia operators
		for( int i = 0; i < env->Block_Num_hole; i++ ) {

	                for( int j = 0; j < env->Num_J_hole_block[i]; j++ ) {

        	        	dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j];

                	        for( int k = 0; k < dimension; k++ )

                        		S_Dia_env[i][j][k] = (double) 0;
	
			}

		}

		//----
	        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
	        if( Table_J[ N * Table_J_sys[site_s] + StartSite - Table_J_env[site_e]] == 1 ) {

       	                alpha = Interaction_J[ N * Table_J_sys[site_s] + StartSite - Table_J_env[site_e] ];

                	for( int i = 0; i < env->Block_Num_hole; i++ ) 
			for( int j = 0; j < env->Num_J_hole_block[i]; j++ ) {

                	        dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j];

                        	daxpy_(&dimension, &alpha, env->S_Dia[site_e][i][j], &inc, S_Dia_env[i][j], &inc);

                	}

        	}
	
		//----S_M_Dia operators
	        for( int i = 0; i < env->Block_Num_hole; i++ ) {

        	        for( int j = 0; j < env->Num_J_hole_block[i] - 1; j++ ) {

                		dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j+1];

                        	for( int k =0; k < dimension; k++ )

                        		S_M_Dia_env[i][j][k] = (double) 0;

	                }

		}

		//----
	        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
	        if( Table_J[ N * Table_J_sys[site_s] + StartSite - Table_J_env[site_e]] == 1 ) {

			alpha = Interaction_J[ N * Table_J_sys[site_s] + StartSite - Table_J_env[site_e] ];

                	for( int i = 0; i < env->Block_Num_hole; i++ )
			for( int j = 0; j < env->Num_J_hole_block[i] - 1; j++ ) {

                	        dimension = env->Dim_J_block[i][j] * env->Dim_J_block[i][j+1];

                        	daxpy_( &dimension, &alpha, env->S_M_Dia[site_e][i][j], &inc, S_M_Dia_env[i][j], &inc );

                	}

        	}

		//----calculation
                for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i++ )
		if( sys->Num_hole_block[ H_sys_config_3[ i ] ] != sys->TotSiteNo )
		if( env->Num_hole_block[ H_env_config_3[ i ] ] != env->TotSiteNo ) {

			index = 0;

			for( int j = 0; j < BlockNumber_for_TargetBlock_config_3; j++ )
			if( H_sys_config_3[ i ] == H_sys_config_3[ j ]  &&  H_env_config_3[ i ] == H_env_config_3[ j ] )
			if( H_ns_config_3[ i ] == H_ns_config_3[ j ]  &&  H_ne_config_3[ i ] == H_ne_config_3[ j ] )
			if( J_sysnew_config_3[ i ] == J_sysnew_config_3[ j ]  &&  J_envnew_config_3[ i ] == J_envnew_config_3[ j ] )
			if( abs( sys->Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] - sys->Value_J_block[ H_sys_config_3[ j ] ][ J_sys_config_3[ j ] ] ) <= 2  &&  abs( env->Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] - env->Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ] ) <= 2 ) {

                        //------1
                                if( i == j  &&  sys->Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] != 0  &&  env->Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] != 0 ) {

                                        alpha = six_j_J_config_3[ i ][ index ];

                                        mul = sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] * env->Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ];

                                        double * f_sub = new double [mul];

                                        for( int s = 0; s < mul; s++ )  f_sub[ s ] = (double) 0;

                                        dsymm_( &side_R, &uplo, &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &env->Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &alpha_p, S_Dia_env[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &env->Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], f3[ j ], &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &beta_p, f_sub, &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] );

                                        dsymm_( &side_L, &uplo, &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &env->Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &alpha, sys->S_Dia[ site_s ][ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], f_sub, &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &beta, g3[ i ], &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] );

                                        delete [] f_sub;

                                }

                        //------2
                                else if( J_sys_config_3[ i ] == J_sys_config_3[ j ]  &&  sys->Value_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] != 0  &&  J_env_config_3[ j ] - J_env_config_3[ i ] == 1  &&  env->Value_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ] - env->Value_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ] == 2 ) {

                                        alpha = six_j_J_config_3[ i ][ index ];

                                        mul = sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] * env->Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ];

                                        double * f_sub = new double [ mul ];

                                        for( int s = 0; s < mul; s++ )  f_sub[ s ] = (double) 0;

                                        dgemm_( &trans_N, &trans_T, &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &env->Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &env->Dim_J_block[ H_env_config_3[ j ] ][ J_env_config_3[ j ] ], &alpha_p, f3[ j ], &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], S_M_Dia_env[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &env->Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &beta_p, f_sub, &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] );

                                        dsymm_( &side_L, &uplo, &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &env->Dim_J_block[ H_env_config_3[ i ] ][ J_env_config_3[ i ] ], &alpha, sys->S_Dia[ site_s ][ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], f_sub, &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ], &beta, g3[ i ], &sys->Dim_J_block[ H_sys_config_3[ i ] ][ J_sys_config_3[ i ] ] );

                                        delete [] f_sub;

                                }

                        //------3
                                else if( J_sys_config_3[i] == J_sys_config_3[j]  &&  sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]] != 0  &&  J_env_config_3[i] - J_env_config_3[j] == 1  &&  env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]] - env->Value_J_block[H_env_config_3[j]][J_env_config_3[j]] == 2 ) {

                                        alpha = - six_j_J_config_3[i][index];

                                        mul = sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]] * env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]];

                                        double * f_sub=new double [mul];

                                        for(int s=0; s<mul; s++)    f_sub[s]=(double) 0;

                                        dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &alpha_p, f3[j], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], S_M_Dia_env[H_env_config_3[j]][J_env_config_3[j]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]] );

                                        dsymm_( &side_L, &uplo, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &alpha, sys->S_Dia[site_s][H_sys_config_3[i]][J_sys_config_3[i]], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], f_sub, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]] );

                                        delete [] f_sub;

                                }

                        //------4
                                else if( J_sys_config_3[j] - J_sys_config_3[i] == 1  &&  sys->Value_J_block[H_sys_config_3[j]][J_sys_config_3[j]] - sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]] == 2  &&  J_env_config_3[i] == J_env_config_3[j]  &&  env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]] != 0 ) {

                                        alpha = six_j_J_config_3[i][index];

                                        mul = sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] * env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]];

                                        double * f_sub=new double [mul];

                                        for(int s=0; s<mul; s++)    f_sub[s]=(double) 0;

                                        dsymm_( &side_R, &uplo, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &alpha_p, S_Dia_env[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], f3[j], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] );

                                        dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &alpha, sys->S_M_Dia[site_s][H_sys_config_3[i]][J_sys_config_3[i]], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]] );

                                        delete [] f_sub;

                                }

                        //------5
                                else if( J_sys_config_3[j] - J_sys_config_3[i] == 1  &&  sys->Value_J_block[H_sys_config_3[j]][J_sys_config_3[j]] - sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]] == 2  &&  J_env_config_3[j] - J_env_config_3[i] == 1  &&  env->Value_J_block[H_env_config_3[j]][J_env_config_3[j]] - env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]] == 2 ) {

                                        alpha = six_j_J_config_3[i][index];

                                        mul = sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] * env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]];

                                        double * f_sub=new double [mul];

                                        for(int s=0; s<mul; s++)    f_sub[s]=(double) 0;

                                        dgemm_( &trans_N, &trans_T, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &alpha_p, f3[j], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], S_M_Dia_env[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] );

                                        dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &alpha, sys->S_M_Dia[site_s][H_sys_config_3[i]][J_sys_config_3[i]], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]] );

                                        delete [] f_sub;

                                }

                        //------6
                                else if( J_sys_config_3[j] - J_sys_config_3[i]==1  &&  sys->Value_J_block[H_sys_config_3[j]][J_sys_config_3[j]] - sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]] == 2  &&  J_env_config_3[i] - J_env_config_3[j] == 1  &&  env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]] - env->Value_J_block[H_env_config_3[j]][J_env_config_3[j]] == 2 ) {

                                        alpha = -six_j_J_config_3[i][index];

                                        mul = sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] * env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]];

                                        double * f_sub=new double [mul];

                                        for(int s=0; s<mul; s++)    f_sub[s]=(double) 0;

                                        dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &alpha_p, f3[j], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], S_M_Dia_env[H_env_config_3[j]][J_env_config_3[j]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] );

                                        dgemm_(&trans_N, &trans_N, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &alpha, sys->S_M_Dia[site_s][H_sys_config_3[i]][J_sys_config_3[i]], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]]);

                                        delete [] f_sub;

                                }

                        //------7
                                else if( J_sys_config_3[i] - J_sys_config_3[j] == 1  &&  sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]] - sys->Value_J_block[H_sys_config_3[j]][J_sys_config_3[j]] == 2  &&  J_env_config_3[i] == J_env_config_3[j] && env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]] != 0 ) {

                                        alpha = -six_j_J_config_3[i][index];

                                        mul = sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] * env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]];

                                        double * f_sub=new double [mul];

                                        for(int s=0; s<mul; s++)    f_sub[s]=(double) 0;

                                        dsymm_( &side_R, &uplo, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &alpha_p, S_Dia_env[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], f3[j], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] );

                                        dgemm_( &trans_T, &trans_N, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &alpha, sys->S_M_Dia[site_s][H_sys_config_3[j]][J_sys_config_3[j]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]] );

                                        delete [] f_sub;

                                }

                        //------8
                                else if( J_sys_config_3[i] - J_sys_config_3[j] == 1  &&  sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]] - sys->Value_J_block[H_sys_config_3[j]][J_sys_config_3[j]] == 2  &&  J_env_config_3[j] - J_env_config_3[i] == 1  &&  env->Value_J_block[H_env_config_3[j]][J_env_config_3[j]] - env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]] == 2 ) {

                                        alpha = -six_j_J_config_3[i][index];

                                        mul = sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] * env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]];

                                        double * f_sub=new double [mul];

                                        for(int s=0; s<mul; s++)    f_sub[s]=(double) 0;

                                        dgemm_( &trans_N, &trans_T, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &alpha_p, f3[j], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], S_M_Dia_env[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] );

                                        dgemm_( &trans_T, &trans_N, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &alpha, sys->S_M_Dia[site_s][H_sys_config_3[j]][J_sys_config_3[j]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]] );

                                        delete [] f_sub;

                                }

                        //------9
                                else if( J_sys_config_3[i] - J_sys_config_3[j] == 1  &&  sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]] - sys->Value_J_block[H_sys_config_3[j]][J_sys_config_3[j]] == 2  &&  J_env_config_3[i] - J_env_config_3[j] == 1  &&  env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]] - env->Value_J_block[H_env_config_3[j]][J_env_config_3[j]] == 2 ) {

                                        alpha = six_j_J_config_3[i][index];

                                        mul = sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] * env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]];

                                        double * f_sub=new double [mul];

                                        for(int s=0; s<mul; s++)    f_sub[s]=(double) 0;

                                        dgemm_( &trans_N, &trans_N, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &alpha_p, f3[j], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], S_M_Dia_env[H_env_config_3[j]][J_env_config_3[j]], &env->Dim_J_block[H_env_config_3[j]][J_env_config_3[j]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]] );

                                        dgemm_( &trans_T, &trans_N, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &alpha, sys->S_M_Dia[site_s][H_sys_config_3[j]][J_sys_config_3[j]], &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], f_sub, &sys->Dim_J_block[H_sys_config_3[j]][J_sys_config_3[j]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]] );

                                        delete [] f_sub;

                                }

				index++;

			}

		}

	}

}

//=======================================================================================
inline void Super::H_Sys_Env_N() {

        for(site_s=0; site_s<operator_number_J_sys; site_s++) {

		for( int i = 0; i < env -> Block_Num_hole; i++ ) {

        	        for( int j = 0; j < env -> Num_J_hole_block[i]; j++ ) {

                		dimension = env -> Dim_J_block[i][j] * env -> Dim_J_block[i][j];

                        	for( int k = 0; k < dimension; k++ )

                        		NN_env[i][j][k] = (double) 0;
	
			}

		}

	        for( site_e = 0; site_e < operator_number_J_env; site_e++ )
        	if(Table_N[ N * Table_N_sys[site_s] + StartSite - Table_N_env[site_e] ] == 1 ) {

                	for( int i = 0; i < env -> Block_Num_hole; i++ ) 
			for( int j = 0; j < env -> Num_J_hole_block[i]; j++ ) {

        	                alpha = Interaction_N[ N * Table_N_sys[site_s] + StartSite - Table_N_env[site_e]];
                	        dimension = env -> Dim_J_block[i][j] * env -> Dim_J_block[i][j];
                        	daxpy_( &dimension, &alpha, env -> NN[site_e][i][j], &inc, NN_env[i][j], &inc );

                	}

        	}


               	for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++)
		if(sys->Num_hole_block[H_sys_config_3[i]]!=sys->TotSiteNo && env->Num_hole_block[H_env_config_3[i]]!=env->TotSiteNo) {

			alpha = 1.0 / sqrt((sys->Value_J_block[H_sys_config_3[i]][J_sys_config_3[i]]+1.0)*(env->Value_J_block[H_env_config_3[i]][J_env_config_3[i]]+1.0));

			double * f_sub=new double [Dim_block_config_3[i]];
			for(int s=0; s<Dim_block_config_3[i]; s++)  f_sub[s]=(double) 0;
		
			dsymm_(&side_R, &uplo, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &alpha_p, NN_env[H_env_config_3[i]][J_env_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], f3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &beta_p, f_sub, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]]);

			dsymm_(&side_L, &uplo, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &env->Dim_J_block[H_env_config_3[i]][J_env_config_3[i]], &alpha, sys->NN[site_s][H_sys_config_3[i]][J_sys_config_3[i]], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], f_sub, &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]], &beta, g3[i], &sys->Dim_J_block[H_sys_config_3[i]][J_sys_config_3[i]]);

			delete [] f_sub;

		}

	}

}

//====================================Free the space of constructing the Hamiltonian==============================
Super::~Super() {

        if(destruct == 's') {

		delete [] H_sys;        delete [] H_env;        delete [] H_sysnew;     delete [] H_envnew;
                delete [] J_sys;        delete [] J_env;        delete [] J_sysnew;     delete [] J_envnew;
                delete [] Dim_block;

                delete [] WaveFunction;                 //delete [] WaveFunction_excited;

                for(int i=0; i<BlockNumber_for_TargetBlock; i++) {
                        delete [] WaveFunction_block[i];        //delete [] WaveFunction_excited_block[i];
                }
                delete [] WaveFunction_block;                   //delete [] WaveFunction_excited_block;

        }

        else if(destruct == 'd') { 

          //------Delete the tables
                delete [] Table_T;  delete [] Table_J;  delete [] Table_N;
		delete [] Interaction_T;  delete [] Interaction_J;  delete [] Interaction_N;
                delete [] Table_T_sys;  delete [] Table_T_env;  delete [] Table_J_sys;  delete [] Table_J_env;
		delete [] Table_N_sys;  delete [] Table_N_env;

                delete [] Table_1to2_Num;       delete [] Table_1to2_Site;

                for(int i=0; i<BlockNumber_for_TargetBlock; i++) 
                        delete [] Table_2to1[i];
                delete [] Table_2to1;

          //------Delete the variables of configurations 1 and 2
             //------block numbers
		delete [] H_sys;  delete [] H_env;  delete [] H_sysnew;  delete [] H_envnew;
                delete [] J_sys;  delete [] J_env;  delete [] J_sysnew;  delete [] J_envnew;
                delete [] Dim_block;

             //------wave functions
                for(int i=0; i<BlockNumber_for_TargetBlock; i++) {
                        delete [] f1[i];        delete [] g1[i];        delete [] f2[i];        delete [] g2[i];
                }
                delete [] f1;   delete [] g1;   delete [] f2;   delete [] g2;

             //------6j coefficients
                delete [] six_j_1_sys_n;        delete [] six_j_1_e_env;

             //------9j coefficient
                for(int i=0; i<BlockNumber_for_TargetBlock; i++) {
                        delete [] nine_j_config_2[i];  delete [] nine_j_config_2_inverse[i];
		}
                delete [] nine_j_config_2;  delete [] nine_j_config_2_inverse;

          //------Delete the variables of configuration 3
              //------9j coefficients
                for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++)
                        delete [] nine_j_config_3[i];
                delete [] nine_j_config_3;

                for(int i=0; i<BlockNumber_for_TargetBlock; i++)
                        delete [] nine_j_config_3_inverse[i];
                delete [] nine_j_config_3_inverse;

              //------6j coefficients
                delete [] six_j_J_config_3_diaele;
                delete [] six_j_2_config_3;
             
                for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++)
		if( sys -> Num_hole_block[ H_sys_config_3[ i ] ] < sys -> TotSiteNo  &&  env -> Num_hole_block[ H_env_config_3[ i ] ] < env -> TotSiteNo )
                        delete [] six_j_J_config_3[i];

                delete [] six_j_J_config_3;

        	for( int i = 0; i < BlockNumber_for_TargetBlock_config_3; i ++ )
		if( ( ( sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] ) > 0 )  &&  ( sys -> Num_hole_block[ H_sys_config_3[ i ] ] + env -> Num_hole_block[ H_env_config_3[ i ] ] ) < ( sys -> TotSiteNo + env -> TotSiteNo ) )
			delete [] six_j_T_config_3[i];

		delete [] six_j_T_config_3;

     //------wave functions
                for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++) {
                        delete [] f3[i];        delete [] g3[i];
                }
                delete [] f3;                   delete [] g3;
         
             //------block numbers	
		delete [] H_sys_config_3;  delete [] H_env_config_3;  delete [] H_ns_config_3;  delete [] H_ne_config_3;
                delete [] J_sys_config_3;  delete [] J_env_config_3;       
		delete [] J_sysnew_config_3;  delete [] J_envnew_config_3;    delete [] Dim_block_config_3;

    //------delete operators
		for( int i = 0; i < sys -> Block_Num_hole - 1; i++ ) {
		
			for( int j = 0; j < sys -> Num_block_T[i]; j++ ) {

				delete [] T_sys[i][j];
		
			}

			delete [] T_sys[i];

		}

		delete [] T_sys;


		for( int i = 0; i < env -> Block_Num_hole - 1; i++ ) {
			
			for( int j = 0; j < env -> Num_block_T[i]; j++ ) {

				delete [] T_env[i][j];
		
			}

			delete [] T_env[i];

		}

		delete [] T_env;


                for( int i = 0; i < sys -> Block_Num_hole; i++ ) {

                        for( int j = 0; j < sys -> Num_J_hole_block[i]; j++ ) {
        
                                delete [] S_Dia_sys[i][j];
                
                        }

                        delete [] S_Dia_sys[i];

                }
        
                delete [] S_Dia_sys;


                for( int i = 0; i < sys -> Block_Num_hole; i++ ) {

                        for( int j = 0; j < sys -> Num_J_hole_block[i] - 1; j++ ) {
        
                                delete [] S_M_Dia_sys[i][j];
                
                        }

                        delete [] S_M_Dia_sys[i];

                }
        
                delete [] S_M_Dia_sys;


		for( int i = 0; i < env -> Block_Num_hole; i++ ) {

                        for( int j = 0; j < env -> Num_J_hole_block[i]; j++ ) {
        
                                delete [] S_Dia_env[i][j];
                
                        }

                        delete [] S_Dia_env[i];

                }
        
                delete [] S_Dia_env;


                for( int i = 0; i < env -> Block_Num_hole; i++ ) {

                        for( int j = 0; j < env -> Num_J_hole_block[i] - 1; j++ ) {
        
                                delete [] S_M_Dia_env[i][j];
                
                        }

                        delete [] S_M_Dia_env[i];

                }
        
                delete [] S_M_Dia_env;


                for( int i = 0; i < sys -> Block_Num_hole; i++ ) {

                        for( int j = 0; j < sys -> Num_J_hole_block[i]; j++ ) {
        
                                delete [] NN_sys[i][j];
                
                        }

                        delete [] NN_sys[i];

                }
        
                delete [] NN_sys;


		for( int i = 0; i < env -> Block_Num_hole; i++ ) {

                        for( int j = 0; j < env -> Num_J_hole_block[i]; j++ ) {
        
                                delete [] NN_env[i][j];
                
                        }

                        delete [] NN_env[i];

                }
        
                delete [] NN_env;

	}

	else if( destruct == 'm' ) {

		delete [] six_j_1_sys_n;  delete [] six_j_1_e_env;

		delete [] six_j_2_config_3;

		//----
                for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++)

                        delete [] nine_j_config_3[i];

                delete [] nine_j_config_3;

                for(int i=0; i<BlockNumber_for_TargetBlock; i++)

                        delete [] nine_j_config_2[i];

                delete [] nine_j_config_2;

		//----
	        for(int i=0; i<BlockNumber_for_TargetBlock_config_3; i++)
		if( sys_space -> Num_hole_block[ H_sys_config_3[ i ] ] < sys_space -> TotSiteNo  &&  env_space -> Num_hole_block[ H_env_config_3[ i ] ] < env_space -> TotSiteNo ) {

                        delete [] six_j_J_config_3[i];

		}

                delete [] six_j_J_config_3;

		//----
		delete [] H_sys_config_3;  delete [] H_env_config_3;  
		delete [] H_ns_config_3;  delete [] H_ne_config_3;
                delete [] J_sys_config_3;  delete [] J_env_config_3;       
		delete [] J_sysnew_config_3;  delete [] J_envnew_config_3;    
		delete [] Dim_block_config_3;

		//----
		delete [] H_sys;        delete [] H_env;        
		delete [] H_sysnew;     delete [] H_envnew;
                delete [] J_sys;        delete [] J_env;        
		delete [] J_sysnew;     delete [] J_envnew;
                delete [] Dim_block;

	}

}

//=====================================Normalize the WaveFunction================================================
void Super::NormalizedCopy(const double *f, double *g) {
        double res=(double) 0 ;
        for(int i=0; i<Dim; i++) res+=f[i]*f[i];
                res=sqrt(res);
        for(int j=0; j<Dim; j++)
                g[j]=f[j]/res;
}
//=====================================================END=======================================================
