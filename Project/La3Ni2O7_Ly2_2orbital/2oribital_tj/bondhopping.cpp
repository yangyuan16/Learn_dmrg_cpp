#include<iostream>
using namespace std;
#include<math.h>
#include<time.h>
#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_sf_coupling.h>

#include"common.h"
#include"sub.h"
#include"super.h"
#include"bondhopping.h"

#define constant sqrt(6.0)*0.5
#define hoping -sqrt(2.0)
#define ref_site 0	// You may change the reference site here
//===================================================================================
//				BLAS ROUTINES
//===================================================================================
extern "C" {

void daxpy_( const int *n, const double *alpha, double *x, const int *incx, double *y, const int *incy );

void dsymm_( char *side, char *uplo, const int *m, const int *n, const double *alpha, double *a, const int *lda, double *b, const int *ldb, const double *beta, double *c, const int *ldc );

void dgemm_( char *transa, char *transb, const int *m, const int *n, const int *k, const double *alpha, double *a, const int *lda, double *b, const int *ldb, const double *beta, double *c, const int *ldc );

}
//====================================================================================

//====================================================================================
//			Calculating the electron density
//====================================================================================
BondHopping::BondHopping( const int &lsys, Parameter &para, Super &sup ) {

	//----Define variables
	trans_N = 'N';	trans_T = 'T';
	side_R = 'R';	side_L = 'L';
	uplo = 'U';

	CreateSpace( lsys, para );  //matrices A, B and angular coupling coefficients
	CreateFunction( sup );	    //The function in which the physical quantities 

	//----
	correlation = new double * [ para.Total_N ];

	for( int i = 0; i < para.Total_N; i++ ) {

		correlation[ i ] = new double [ para.Total_N ];

		for( int j = 0; j < para.Total_N; j++ ) {

			correlation[ i ][ j ] = 0.0;

		}

	}

	//----Create filename
        stringstream ss;
        ss<<"measurement_bond_hopping.dat";
        filename=ss.str();
        fout.precision(12);
        fout.open(filename.c_str()); //start writing file

  //----Calculate the spin correlations-----------------------------
	Correlation_sys( lsys, para, sup );
	Correlation_env( para.Total_N - lsys - 2, para, sup );
	Correlation_ns_ne( lsys, para, sup );
	Correlation_sys_env( lsys, para, sup );

        fout.close(); //finish writing file

  //----Delete matrices A and B-------------------------------------
	for( int i = 0; i < para.Total_N; i++ ) {

		delete [] correlation[ i ];

	}

	delete [] correlation;

	//----
	DeleteFunction( sup );

	DeleteSpace( lsys, para );

}

//============================================================================
//		     Create Space for variables
//============================================================================
inline void BondHopping::CreateSpace( const int &lsys, Parameter &para ) {

	int A_site = lsys;
	int B_site = para.Total_N - lsys - 2;

//------Create space for A_N-matrices----------------------------------------------------------------------------
	A_N_Block_Num_hole = new int [ A_site ];

	for( int i = 0; i < A_site; i++ ) {

		A_N_Block_Num_hole[i] = 0;

	}

	A_N_Num_hole_block = new int * [ A_site ];
	A_N_Num_J_hole_block = new int * [ A_site ];

	A_N_Value_J_block = new int ** [ A_site ];
	A_N_Dim_J_block = new int ** [ A_site ];

	A_N_Hole_blockOld = new int *** [ A_site ];
	A_N_J_blockOld = new int *** [ A_site ];
	A_N_Start = new int *** [ A_site ];

  //----The first site-----------------------------------------------
	A_N_Block_Num_hole[0] = 2;

	//----
	A_N_Num_hole_block[0] = new int [ A_N_Block_Num_hole[0] ];
	A_N_Num_J_hole_block[0] = new int [ A_N_Block_Num_hole[0] ];

	for( int i = 0; i < A_N_Block_Num_hole[0]; i++ ) {

		A_N_Num_hole_block[0][i] = i;  
		A_N_Num_J_hole_block[0][i] = 1;

	}
 
	//----
	A_N_Value_J_block[0] = new int * [ A_N_Block_Num_hole[0] ];
	A_N_Dim_J_block[0] = new int * [ A_N_Block_Num_hole[0] ];

	for( int i = 0; i < A_N_Block_Num_hole[0]; i++ ) {

		A_N_Value_J_block[0][i] =  new int [ A_N_Num_J_hole_block[0][i] ];
		A_N_Dim_J_block[0][i] =  new int [ A_N_Num_J_hole_block[0][i] ];

		for( int j = 0; j < A_N_Num_J_hole_block[0][i]; j++ ) {

			A_N_Value_J_block[0][i][j] = 1 - i;
			A_N_Dim_J_block[0][i][j] = 1;

		}

	}

	//----
	A_N_Hole_blockOld[0] = new int ** [ A_N_Block_Num_hole[0] ];
	A_N_J_blockOld[0] = new int ** [ A_N_Block_Num_hole[0] ];
	A_N_Start[0] = new int ** [ A_N_Block_Num_hole[0] ];

	for( int i = 0; i < A_N_Block_Num_hole[0]; i++ ) {

		A_N_Hole_blockOld[0][i] = new int * [3];
		A_N_J_blockOld[0][i] = new int * [3];
		A_N_Start[0][i] = new int * [3];

		for( int j = 0; j < 3; j++ ) {

			A_N_Hole_blockOld[0][i][j] = new int [ A_N_Num_J_hole_block[0][i] ];
			A_N_J_blockOld[0][i][j] = new int [ A_N_Num_J_hole_block[0][i] ];
			A_N_Start[0][i][j] = new int [ A_N_Num_J_hole_block[0][i] ];

			for( int k = 0; k < A_N_Num_J_hole_block[0][i]; k++ ) {

				A_N_Hole_blockOld[0][i][j][k] = 0;
				A_N_J_blockOld[0][i][j][k] = 0;
				A_N_Start[0][i][j][k] = 0;

			}

		}

	}

  //----Read other sites from disk------------------------------------
	for( int i = 1; i < A_site; i++ ) {

		FILE *fr = fopen( Combine ( Combine ( "new_block/", 1 ), i+1 ), "rb" );	//"i"=1 stand for sys

		//----
		fread( &A_N_Block_Num_hole[i], sizeof(int), 1, fr );

		//----
		A_N_Num_hole_block[i] = new int [ A_N_Block_Num_hole[i] ];
		A_N_Num_J_hole_block[i] = new int [ A_N_Block_Num_hole[i] ];

		fread( A_N_Num_hole_block[i], sizeof(int), A_N_Block_Num_hole[i], fr );
		fread( A_N_Num_J_hole_block[i], sizeof(int), A_N_Block_Num_hole[i], fr );

		//----
		A_N_Value_J_block[i] = new int * [ A_N_Block_Num_hole[i] ];
		A_N_Dim_J_block[i] = new int * [ A_N_Block_Num_hole[i] ];

		for( int j = 0; j < A_N_Block_Num_hole[i]; j++ ) {

			A_N_Value_J_block[i][j] = new int [ A_N_Num_J_hole_block[i][j] ];
			A_N_Dim_J_block[i][j] = new int [ A_N_Num_J_hole_block[i][j] ];

		}

		for( int j = 0; j < A_N_Block_Num_hole[i]; j++ ) {

			fread( A_N_Value_J_block[i][j], sizeof(int), A_N_Num_J_hole_block[i][j], fr );
			fread( A_N_Dim_J_block[i][j], sizeof(int), A_N_Num_J_hole_block[i][j], fr );

		}

		//----
		A_N_Hole_blockOld[i] = new int ** [ A_N_Block_Num_hole[i] ];
		A_N_J_blockOld[i] = new int ** [ A_N_Block_Num_hole[i] ];
		A_N_Start[i] = new int ** [ A_N_Block_Num_hole[i] ];

		for( int j = 0; j < A_N_Block_Num_hole[i]; j++ ) {

			A_N_Hole_blockOld[i][j] = new int * [3];
			A_N_J_blockOld[i][j] = new int * [3];
			A_N_Start[i][j] = new int * [3];

			for( int k = 0; k < 3; k++ ) {

				A_N_Hole_blockOld[i][j][k] = new int [ A_N_Num_J_hole_block[i][j] ];
				A_N_J_blockOld[i][j][k] = new int [ A_N_Num_J_hole_block[i][j] ];
				A_N_Start[i][j][k] = new int [ A_N_Num_J_hole_block[i][j] ];

			}

			for( int k = 0; k < 3; k++ ) {

				fread( A_N_Hole_blockOld[i][j][k], sizeof(int), A_N_Num_J_hole_block[i][j], fr );
				fread( A_N_J_blockOld[i][j][k], sizeof(int), A_N_Num_J_hole_block[i][j], fr );
				fread( A_N_Start[i][j][k], sizeof(int), A_N_Num_J_hole_block[i][j], fr );

			}

		}

		fclose(fr);

	}

  //----T operator block
	A_N_Num_block_T = new int * [ A_site ];

	A_N_J_block_T_bra = new int ** [ A_site ];
	A_N_J_block_T_ket = new int ** [ A_site ];

	//----The first site

	//----
	A_N_Num_block_T[ 0 ] = new int [ 1 ];

	for( int i = 0; i < 1; i++ ) {

		A_N_Num_block_T[ 0 ][ i ] = 1;

	}

	//----
	A_N_J_block_T_bra[ 0 ] = new int * [ 1 ];
	A_N_J_block_T_ket[ 0 ] = new int * [ 1 ];

	for( int i = 0; i < 1; i++ ) {

		A_N_J_block_T_bra[ 0 ][ i ] =  new int [ 1 ];
		A_N_J_block_T_ket[ 0 ][ i ] =  new int [ 1 ];

		for( int j = 0; j < 1; j++ ) {

			A_N_J_block_T_bra[ 0 ][ i ][ j ] =  0;
			A_N_J_block_T_ket[ 0 ][ i ][ j ] =  0;

		}

	}

	//----Read from the disk
	for( int i = 1; i < A_site; i++ ) {

		FILE *fr = fopen( Combine ( Combine ( "new_block_T/", 1 ), i+1 ), "rb" );

		//----
		A_N_Num_block_T[ i ] = new int [ A_N_Block_Num_hole[ i ] - 1 ];

		fread( A_N_Num_block_T[ i ], sizeof(int), A_N_Block_Num_hole[ i ] - 1, fr );

		//----
		A_N_J_block_T_bra[ i ] = new int * [ A_N_Block_Num_hole[ i ] - 1 ];
		A_N_J_block_T_ket[ i ] = new int * [ A_N_Block_Num_hole[ i ] - 1 ];

		for( int j = 0; j < A_N_Block_Num_hole[ i ] - 1; j++ ) {

			A_N_J_block_T_bra[ i ][ j ] = new int [ A_N_Num_block_T[ i ][ j ] ];
			A_N_J_block_T_ket[ i ][ j ] = new int [ A_N_Num_block_T[ i ][ j ] ];

		}

		for( int j = 0; j < A_N_Block_Num_hole[ i ] - 1; j++ ) {

			fread( A_N_J_block_T_bra[ i ][ j ], sizeof(int), A_N_Num_block_T[ i ][ j ], fr );
			fread( A_N_J_block_T_ket[ i ][ j ], sizeof(int), A_N_Num_block_T[ i ][ j ], fr );

		}

		fclose(fr);

	}

//------Create space for B_N-matrices----------------------------------------------------------------------------
	B_N_Block_Num_hole = new int [ B_site ];

       for( int i = 0; i < B_site; i++ ) {

                B_N_Block_Num_hole[i] = 0;

        }

	B_N_Num_hole_block = new int * [ B_site ];
	B_N_Num_J_hole_block = new int * [ B_site ];

	B_N_Value_J_block = new int ** [ B_site ];
	B_N_Dim_J_block = new int ** [ B_site ];

	B_N_Hole_blockOld = new int *** [ B_site ];
	B_N_J_blockOld = new int *** [ B_site ];
	B_N_Start = new int *** [ B_site ];

  //----The first site--------------------------------------------
	B_N_Block_Num_hole[0] = 2;

	//----
	B_N_Num_hole_block[0] = new int [ B_N_Block_Num_hole[0] ];
	B_N_Num_J_hole_block[0] = new int [ B_N_Block_Num_hole[0] ];

	for( int i = 0; i < B_N_Block_Num_hole[0]; i++ ) {

		B_N_Num_hole_block[0][i] = i;  
		B_N_Num_J_hole_block[0][i] = 1;

	}
 
	//----
	B_N_Value_J_block[0] = new int * [ B_N_Block_Num_hole[0] ];
	B_N_Dim_J_block[0] = new int * [ B_N_Block_Num_hole[0] ];

	for( int i = 0; i < B_N_Block_Num_hole[0]; i++ ) {

		B_N_Value_J_block[0][i] =  new int [ B_N_Num_J_hole_block[0][i] ];
		B_N_Dim_J_block[0][i] =  new int [ B_N_Num_J_hole_block[0][i] ];

		for( int j = 0; j < B_N_Num_J_hole_block[0][i]; j++ ) {

			B_N_Value_J_block[0][i][j] = 1 - i;
			B_N_Dim_J_block[0][i][j] = 1;

		}

	}

	//----
	B_N_Hole_blockOld[0] = new int ** [ B_N_Block_Num_hole[0] ];
	B_N_J_blockOld[0] = new int ** [ B_N_Block_Num_hole[0] ];
	B_N_Start[0] = new int ** [ B_N_Block_Num_hole[0] ];

	for( int i = 0; i < B_N_Block_Num_hole[0]; i++ ) {

		B_N_Hole_blockOld[0][i] = new int * [ 3 ];
		B_N_J_blockOld[0][i] = new int * [ 3 ];
		B_N_Start[0][i] = new int * [ 3 ];

		for( int j = 0; j < 3; j++ ) {

			B_N_Hole_blockOld[0][i][j] = new int [ B_N_Num_J_hole_block[0][i] ];
			B_N_J_blockOld[0][i][j] = new int [ B_N_Num_J_hole_block[0][i] ];
			B_N_Start[0][i][j] = new int [ B_N_Num_J_hole_block[0][i] ];

			for( int k = 0; k < B_N_Num_J_hole_block[0][i]; k++ ) {

				B_N_Hole_blockOld[0][i][j][k] = 0;
				B_N_J_blockOld[0][i][j][k] = 0;
				B_N_Start[0][i][j][k] = 0;

			}

		}

	}

  //----Read other sites from disk---------------------------------
	for( int i = 1; i < B_site; i++ ) {

		FILE *fr = fopen( Combine ( Combine ( "new_block/", 2 ), i+1 ), "rb" );	//"i"=2 stand for sys

		//----
		fread( &B_N_Block_Num_hole[i], sizeof(int), 1, fr );

		//----
		B_N_Num_hole_block[i] = new int [ B_N_Block_Num_hole[i] ];
		B_N_Num_J_hole_block[i] = new int [ B_N_Block_Num_hole[i] ];

		fread( B_N_Num_hole_block[i], sizeof(int), B_N_Block_Num_hole[i], fr );
		fread( B_N_Num_J_hole_block[i], sizeof(int), B_N_Block_Num_hole[i], fr );

		//----
		B_N_Value_J_block[i] = new int * [ B_N_Block_Num_hole[i] ];
		B_N_Dim_J_block[i] = new int * [ B_N_Block_Num_hole[i] ];

		for( int j = 0; j < B_N_Block_Num_hole[i]; j++ ) {

			B_N_Value_J_block[i][j] = new int [ B_N_Num_J_hole_block[i][j] ];
			B_N_Dim_J_block[i][j] = new int [ B_N_Num_J_hole_block[i][j] ];

		}

		for( int j = 0; j < B_N_Block_Num_hole[i]; j++ ) {

			fread( B_N_Value_J_block[i][j], sizeof(int), B_N_Num_J_hole_block[i][j], fr );
			fread( B_N_Dim_J_block[i][j], sizeof(int), B_N_Num_J_hole_block[i][j], fr );

		}

		//----
		B_N_Hole_blockOld[i] = new int ** [ B_N_Block_Num_hole[i] ];
		B_N_J_blockOld[i] = new int ** [ B_N_Block_Num_hole[i] ];
		B_N_Start[i] = new int ** [ B_N_Block_Num_hole[i] ];

		for( int j = 0; j < B_N_Block_Num_hole[i]; j++ ) {

			B_N_Hole_blockOld[i][j] = new int * [ 3 ];
			B_N_J_blockOld[i][j] = new int * [ 3 ];
			B_N_Start[i][j] = new int * [ 3 ];

			for( int k = 0; k < 3; k++ ) {

				B_N_Hole_blockOld[i][j][k] = new int [ B_N_Num_J_hole_block[i][j] ];
				B_N_J_blockOld[i][j][k] = new int [ B_N_Num_J_hole_block[i][j] ];
				B_N_Start[i][j][k] = new int [ B_N_Num_J_hole_block[i][j] ];

			}

			for( int k = 0; k < 3; k++ ) {

				fread( B_N_Hole_blockOld[i][j][k], sizeof(int), B_N_Num_J_hole_block[i][j], fr );
				fread( B_N_J_blockOld[i][j][k], sizeof(int), B_N_Num_J_hole_block[i][j], fr );
				fread( B_N_Start[i][j][k], sizeof(int), B_N_Num_J_hole_block[i][j], fr );

			}

		}

		fclose(fr);

	}

  //----T operator block
	B_N_Num_block_T = new int * [ B_site ];

	B_N_J_block_T_bra = new int ** [ B_site ];
	B_N_J_block_T_ket = new int ** [ B_site ];

	//----The first site

	//----
	B_N_Num_block_T[ 0 ] = new int [ 1 ];

	for( int i = 0; i < 1; i++ ) {

		B_N_Num_block_T[ 0 ][ i ] = 1;

	}

	//----
	B_N_J_block_T_bra[ 0 ] = new int * [ 1 ];
	B_N_J_block_T_ket[ 0 ] = new int * [ 1 ];

	for( int i = 0; i < 1; i++ ) {

		B_N_J_block_T_bra[ 0 ][ i ] =  new int [ 1 ];
		B_N_J_block_T_ket[ 0 ][ i ] =  new int [ 1 ];

		for( int j = 0; j < 1; j++ ) {

			B_N_J_block_T_bra[ 0 ][ i ][ j ] =  0;
			B_N_J_block_T_ket[ 0 ][ i ][ j ] =  0;

		}

	}

	//----Read from the disk
	for( int i = 1; i < B_site; i++ ) {

		FILE *fr = fopen( Combine ( Combine ( "new_block_T/", 2 ), i+1 ), "rb" );

		//----
		B_N_Num_block_T[ i ] = new int [ B_N_Block_Num_hole[ i ] - 1 ];

		fread( B_N_Num_block_T[ i ], sizeof(int), B_N_Block_Num_hole[ i ] - 1, fr );

		//----
		B_N_J_block_T_bra[ i ] = new int * [ B_N_Block_Num_hole[ i ] - 1 ];
		B_N_J_block_T_ket[ i ] = new int * [ B_N_Block_Num_hole[ i ] - 1 ];

		for( int j = 0; j < B_N_Block_Num_hole[ i ] - 1; j++ ) {

			B_N_J_block_T_bra[ i ][ j ] = new int [ B_N_Num_block_T[ i ][ j ] ];
			B_N_J_block_T_ket[ i ][ j ] = new int [ B_N_Num_block_T[ i ][ j ] ];

		}

		for( int j = 0; j < B_N_Block_Num_hole[ i ] - 1; j++ ) {

			fread( B_N_J_block_T_bra[ i ][ j ], sizeof(int), B_N_Num_block_T[ i ][ j ], fr );
			fread( B_N_J_block_T_ket[ i ][ j ], sizeof(int), B_N_Num_block_T[ i ][ j ], fr );

		}

		fclose(fr);

	}

//------Create space for A-matrices-----------------------------------------------------------------------------
	A_Block_Num_hole = new int [ A_site ];

	A_Num_hole_block =  new int * [ A_site ];
	A_Num_J_hole_block = new int * [ A_site ];

	A_Value_J_block =  new int ** [ A_site ];
	A_Dim_J_block = new int ** [ A_site ];

	A_density_dim = new int ** [ A_site ];
	A = new double *** [ A_site ];

	A_Old_hole = new int * [ A_site ];
	A_Old_J = new int ** [ A_site ];

  //----The first site----------------------------------------------
  	A_Block_Num_hole[0] = A_N_Block_Num_hole[0];

	//----
	A_Num_hole_block[0] = new int [ A_Block_Num_hole[0] ];
	A_Num_J_hole_block[0] = new int [ A_Block_Num_hole[0] ];

	for( int i = 0; i < A_Block_Num_hole[0]; i++ ) {

		A_Num_hole_block[0][i] = i;  A_Num_J_hole_block[0][i] = 1;

	}

	//----
	A_Value_J_block[0] = new int * [ A_Block_Num_hole[0] ];
	A_Dim_J_block[0] = new int * [ A_Block_Num_hole[0] ];
	A_density_dim[0] = new int * [ A_Block_Num_hole[0] ];

	for( int i = 0; i < A_Block_Num_hole[0]; i++ ) {

		A_Value_J_block[0][i] =  new int [ A_Num_J_hole_block[0][i] ];
		A_Dim_J_block[0][i] =  new int [ A_Num_J_hole_block[0][i] ];
		A_density_dim[0][i] =  new int [ A_Num_J_hole_block[0][i] ];

		for( int j = 0; j < A_Num_J_hole_block[0][i]; j++ ) {

			A_Value_J_block[0][i][j] = 1 - i;
			A_Dim_J_block[0][i][j] = 1;
			A_density_dim[0][i][j] = 1;

		}

	}

	//----
	A[0] = new double ** [ A_Block_Num_hole[0] ];

	for( int i = 0; i < A_Block_Num_hole[0]; i++ ) {

		A[0][i] = new double * [ A_Num_J_hole_block[0][i] ];

		for ( int j = 0; j < A_Num_J_hole_block[0][i]; j++ ) {

			A[0][i][j] =  new double [ A_density_dim[0][i][j] ];

			for( int k = 0; k < A_density_dim[0][i][j]; k++ ) {

				A[0][i][j][k] = 1.0;

			}

		}

	}

	//----
	A_Old_hole[0] = new int [ A_Block_Num_hole[0] ];

	for( int i = 0; i < A_Block_Num_hole[0]; i++ ) {

		A_Old_hole[0][i] = i;

	}

	//----
	A_Old_J[0] = new int * [ A_Block_Num_hole[0] ];

	for( int i = 0; i < A_Block_Num_hole[0]; i++ ) {

		A_Old_J[0][i] = new int [ A_Num_J_hole_block[0][i] ];

		for( int j = 0; j < A_Num_J_hole_block[0][i]; j++ ) {

			A_Old_J[0][i][j] = j;

		}

	}

  //----Read other sites from disk----------------------------------
	for( int i = 1; i < A_site; i++ ) {

		FILE *fp = fopen( Combine ( Combine ( "truncated_density_eigenvector/", 1 ), i+1 ), "rb" );

		//----
		fread( &A_Block_Num_hole[i], sizeof(int), 1, fp );
	
		//----
		A_Num_hole_block[i] = new int [ A_Block_Num_hole[i] ];
		A_Num_J_hole_block[i] = new int [ A_Block_Num_hole[i] ];

		fread( A_Num_hole_block[i], sizeof(int), A_Block_Num_hole[i], fp );
		fread( A_Num_J_hole_block[i], sizeof(int), A_Block_Num_hole[i], fp );

		//----
		A_Value_J_block[i] = new int * [ A_Block_Num_hole[i] ];
		A_Dim_J_block[i] = new int * [ A_Block_Num_hole[i] ];
		A_density_dim[i] = new int * [ A_Block_Num_hole[i] ];

		for( int j = 0; j < A_Block_Num_hole[i]; j++ ) {

			A_Value_J_block[i][j] = new int [ A_Num_J_hole_block[i][j] ];
			A_Dim_J_block[i][j] = new int [ A_Num_J_hole_block[i][j] ];
			A_density_dim[i][j] = new int [ A_Num_J_hole_block[i][j] ];

			for( int k = 0; k < A_Num_J_hole_block[i][j]; k++ ) {

				A_Value_J_block[i][j][k] = 0;
				A_Dim_J_block[i][j][k] = 0;
				A_density_dim[i][j][k] = 0;

			}

		}

		for( int j = 0; j < A_Block_Num_hole[i]; j++ ) {

			fread( A_Value_J_block[i][j], sizeof(int), A_Num_J_hole_block[i][j], fp );
			fread( A_Dim_J_block[i][j], sizeof(int), A_Num_J_hole_block[i][j], fp );
			fread( A_density_dim[i][j], sizeof(int), A_Num_J_hole_block[i][j], fp );

		}

		//----
		A[i] = new double ** [ A_Block_Num_hole[i] ];

		for( int j = 0; j < A_Block_Num_hole[i]; j++ ) {

			A[i][j] = new double * [ A_Num_J_hole_block[i][j] ];

			for( int k = 0; k < A_Num_J_hole_block[i][j]; k++ ) {

				A[i][j][k] =  new double [ A_density_dim[i][j][k] ];

				for( int l = 0; l < A_density_dim[i][j][k]; l++ ) {

					A[i][j][k][l] = 0.0;

				}

			}

		}

		for( int j = 0; j < A_Block_Num_hole[i]; j++ )
		for( int k = 0; k < A_Num_J_hole_block[i][j]; k++ ) {

	                fread( A[i][j][k], sizeof(double), A_density_dim[i][j][k], fp );

		}

		//----
		A_Old_hole[i] = new int [ A_Block_Num_hole[i] ];

		fread( A_Old_hole[i], sizeof(int), A_Block_Num_hole[i], fp );


		//----
		A_Old_J[i] = new int * [ A_Block_Num_hole[i] ];

		for( int j = 0; j < A_Block_Num_hole[i]; j++ ) {

			A_Old_J[i][j] = new int [ A_Num_J_hole_block[i][j] ];

			for( int k = 0; k < A_Num_J_hole_block[i][j]; k++ ) {

				A_Old_J[i][j][k] = 0.0;

			}

		}

		for( int j = 0; j < A_Block_Num_hole[i]; j++ ) {

			fread( A_Old_J[i][j], sizeof(int), A_Num_J_hole_block[i][j], fp );

		}

		fclose(fp);

	}

	//----
	A_Num_block_T = new int * [ A_site ];

	for( int n = 0; n < A_site; n++ ) {

		A_Num_block_T[ n ] = new int [ A_Block_Num_hole[ n ] - 1 ];

		for( int i = 0; i < A_Block_Num_hole[ n ] - 1; i++ ) {
	
			A_Num_block_T[ n ][ i ] = 0;	

		}

	}

	int counter = 0;

	for( int n = 0; n < A_site; n++ ) {

		for( int i = 0; i < A_Block_Num_hole[ n ] - 1; i++ )
        	if( A_Num_hole_block[ n ][ i + 1 ] - A_Num_hole_block[ n ][ i ] == 1 ) {

			counter = 0;

	                for( int j = 0; j < A_Num_J_hole_block[ n ][ i ]; j++ )  // J number of the bra vector
        	        for( int k = 0; k < A_Num_J_hole_block[ n ][ i + 1 ]; k++ )  // J number of the ket vector
			if( abs( A_Value_J_block[ n ][ i ][ j ] - A_Value_J_block[ n ][ i + 1 ][ k ] ) == 1 ) {

                        	counter++;

			}

                	A_Num_block_T[ n ][ i ] = counter;  // counter denotes the number of the subblocks

		}

	}

	//----
	A_J_block_T_bra = new int ** [ A_site ];
	A_J_block_T_ket = new int ** [ A_site ];

	for( int n = 0; n < A_site; n++ ) {

		A_J_block_T_bra[ n ] = new int * [ A_Block_Num_hole[ n ] - 1 ];
		A_J_block_T_ket[ n ] = new int * [ A_Block_Num_hole[ n ] - 1 ];

		for( int i = 0; i < A_Block_Num_hole[ n ] - 1; i++ ) {

			A_J_block_T_bra[ n ][ i ] = new int [ A_Num_block_T[ n ][ i ] ];
			A_J_block_T_ket[ n ][ i ] = new int [ A_Num_block_T[ n ][ i ] ];

	                for( int j = 0; j < A_Num_block_T[ n ][ i ]; j++ ) {

	                        A_J_block_T_bra[ n ][ i ][ j ] = 0;  

				A_J_block_T_ket[ n ][ i ][ j ] = 0;

        	        }

		}

	}

	for( int n = 0; n < A_site; n++ ) {

		for( int i = 0; i < A_Block_Num_hole[ n ] - 1; i++ ) 
		if( A_Num_hole_block[ n ][ i + 1 ] - A_Num_hole_block[ n ][ i ] == 1 ) {

			counter = 0;

	                for( int j = 0; j < A_Num_J_hole_block[ n ][ i ]; j++ ) 
	                for( int k = 0; k < A_Num_J_hole_block[ n ][ i + 1 ]; k++ ) 
			if( abs( A_Value_J_block[ n ][ i ][ j ] - A_Value_J_block[ n ][ i + 1 ][ k ] ) == 1 ) {

	                        A_J_block_T_bra[ n ][ i ][ counter ] = j;  

				A_J_block_T_ket[ n ][ i ][ counter ] = k;

				counter++;

        	        }

		}

	}

//------Create space for B-matrices------------------------------------------------------------------------------
	B_Block_Num_hole = new int [ B_site ];

	B_Num_hole_block =  new int * [ B_site ];
	B_Num_J_hole_block = new int * [ B_site ];

	B_Value_J_block =  new int ** [ B_site ];
	B_Dim_J_block = new int ** [ B_site ];

	B_density_dim = new int ** [ B_site ];
	B = new double *** [ B_site ];

	B_Old_hole = new int * [ B_site ];
	B_Old_J = new int ** [ B_site ];

  //----The first site------------------------------------------------
	B_Block_Num_hole[0] = 2;

	//----
	B_Num_hole_block[0] = new int [ B_Block_Num_hole[0] ];
	B_Num_J_hole_block[0] = new int [ B_Block_Num_hole[0] ];

	for( int i = 0; i < B_Block_Num_hole[0]; i++ ) {

		B_Num_hole_block[0][i] = i;  B_Num_J_hole_block[0][i] = 1;

	}

	//----
	B_Value_J_block[0] = new int * [ B_Block_Num_hole[0] ];
	B_Dim_J_block[0] = new int * [ B_Block_Num_hole[0] ];
	B_density_dim[0] = new int * [ B_Block_Num_hole[0] ];

	for( int i = 0; i < B_Block_Num_hole[0]; i++ ) {

		B_Value_J_block[0][i] =  new int [ B_Num_J_hole_block[0][i] ];
		B_Dim_J_block[0][i] =  new int [ B_Num_J_hole_block[0][i] ];
		B_density_dim[0][i] =  new int [ B_Num_J_hole_block[0][i] ];

		for( int j = 0; j < B_Num_J_hole_block[0][i]; j++ ) {

			B_Value_J_block[0][i][j] = 1 - i;
			B_Dim_J_block[0][i][j] = 1;
			B_density_dim[0][i][j] = 1;

		}

	}

	//----
	B[0] = new double ** [ B_Block_Num_hole[0] ];

	for( int i = 0; i < B_Block_Num_hole[0]; i++ ) {

		B[0][i] = new double * [ B_Num_J_hole_block[0][i] ];

		for ( int j = 0; j < B_Num_J_hole_block[0][i]; j++ ) {

			B[0][i][j] =  new double [ B_density_dim[0][i][j] ];

			for( int k = 0; k < B_density_dim[0][i][j]; k++ ) {

				B[0][i][j][k] = 1.0;

			}

		}

	}

	//----
	B_Old_hole[0] = new int [ B_Block_Num_hole[0] ];

	for( int i = 0; i < B_Block_Num_hole[0]; i++ ) {

		B_Old_hole[0][i] = i;

	}

	B_Old_J[0] = new int * [ B_Block_Num_hole[0] ];

	for( int i = 0; i < B_Block_Num_hole[0]; i++ ) {

		B_Old_J[0][i] = new int [ B_Num_J_hole_block[0][i] ];

		for( int j = 0; j < B_Num_J_hole_block[0][i]; j++ ) {

			B_Old_J[0][i][j] = j;

		}

	}

  //----Read other sites from disk-------------------------------------
	for( int i = 1; i < B_site; i++ ) {

		FILE *fp = fopen( Combine ( Combine ( "truncated_density_eigenvector/", 2 ), i+1 ), "rb" );

		//----
		fread( &B_Block_Num_hole[i], sizeof(int), 1, fp );
	
		//----
		B_Num_hole_block[i] = new int [ B_Block_Num_hole[i] ];
		B_Num_J_hole_block[i] = new int [ B_Block_Num_hole[i] ];

		fread( B_Num_hole_block[i], sizeof(int), B_Block_Num_hole[i], fp );
		fread( B_Num_J_hole_block[i], sizeof(int), B_Block_Num_hole[i], fp );

		//----
		B_Value_J_block[i] = new int * [ B_Block_Num_hole[i] ];
		B_Dim_J_block[i] = new int * [ B_Block_Num_hole[i] ];
		B_density_dim[i] = new int * [ B_Block_Num_hole[i] ];

		for( int j = 0; j < B_Block_Num_hole[i]; j++ ) {

			B_Value_J_block[i][j] = new int [ B_Num_J_hole_block[i][j] ];
			B_Dim_J_block[i][j] = new int [ B_Num_J_hole_block[i][j] ];
			B_density_dim[i][j] = new int [ B_Num_J_hole_block[i][j] ];

			for( int k = 0; k < B_Num_J_hole_block[i][j]; k++ ) {

				B_Value_J_block[i][j][k] = 0.0;
				B_Dim_J_block[i][j][k] = 0.0;
				B_density_dim[i][j][k] = 0.0;

			}

		}

		for( int j = 0; j < B_Block_Num_hole[i]; j++ ) {

			fread( B_Value_J_block[i][j], sizeof(int), B_Num_J_hole_block[i][j], fp );
			fread( B_Dim_J_block[i][j], sizeof(int), B_Num_J_hole_block[i][j], fp );
			fread( B_density_dim[i][j], sizeof(int), B_Num_J_hole_block[i][j], fp );

		}

		//----
		B[i] = new double ** [ B_Block_Num_hole[i] ];

		for( int j = 0; j < B_Block_Num_hole[i]; j++ ) {

			B[i][j] = new double * [ B_Num_J_hole_block[i][j] ];

			for( int k = 0; k < B_Num_J_hole_block[i][j]; k++ ) {

				B[i][j][k] =  new double [ B_density_dim[i][j][k] ];

				for( int l = 0; l < B_density_dim[i][j][k]; l++ ) {

					B[i][j][k][l] = 0.0;

				}

			}

		}

		for( int j = 0; j < B_Block_Num_hole[i]; j++ )
		for( int k = 0; k < B_Num_J_hole_block[i][j]; k++ ) {

	                fread( B[i][j][k], sizeof(double), B_density_dim[i][j][k], fp );

		}

		//----
		B_Old_hole[i] = new int [ B_Block_Num_hole[i] ];

		fread( B_Old_hole[i], sizeof(int), B_Block_Num_hole[i], fp );

		//----
		B_Old_J[i] = new int * [ B_Block_Num_hole[i] ];

		for( int j = 0; j < B_Block_Num_hole[i]; j++ ) {

			B_Old_J[i][j] = new int [ B_Num_J_hole_block[i][j] ];

			for( int k = 0; k < B_Num_J_hole_block[i][j]; k++ ) {

				B_Old_J[i][j][k] = 0.0;

			}

		}

		for( int j = 0; j < B_Block_Num_hole[i]; j++ ) {

			fread( B_Old_J[i][j], sizeof(int), B_Num_J_hole_block[i][j], fp );

		}

                fclose(fp);

        }

//----
	B_Num_block_T = new int * [ B_site ];

	for( int n = 0; n < B_site; n++ ) {

		B_Num_block_T[ n ] = new int [ B_Block_Num_hole[ n ] - 1 ];

		for( int i = 0; i < B_Block_Num_hole[ n ] - 1; i++ ) {
	
			B_Num_block_T[ n ][ i ] = 0;	

		}

	}

	for( int n = 0; n < B_site; n++ ) {

		for( int i = 0; i < B_Block_Num_hole[ n ] - 1; i++ )
        	if( B_Num_hole_block[ n ][ i + 1 ] - B_Num_hole_block[ n ][ i ] == 1 ) {

			counter = 0;

	                for( int j = 0; j < B_Num_J_hole_block[ n ][ i ]; j++ )  // J number of the bra vector
        	        for( int k = 0; k < B_Num_J_hole_block[ n ][ i + 1 ]; k++ )  // J number of the ket vector
			if( abs( B_Value_J_block[ n ][ i ][ j ] - B_Value_J_block[ n ][ i + 1 ][ k ] ) == 1 ) {

                        	counter++;

			}

                	B_Num_block_T[ n ][ i ] = counter;  // counter denotes the number of the subblocks

		}

	}

	//----
	B_J_block_T_bra = new int ** [ B_site ];
	B_J_block_T_ket = new int ** [ B_site ];

	for( int n = 0; n < B_site; n++ ) {

		B_J_block_T_bra[ n ] = new int * [ B_Block_Num_hole[ n ] - 1 ];
		B_J_block_T_ket[ n ] = new int * [ B_Block_Num_hole[ n ] - 1 ];

		for( int i = 0; i < B_Block_Num_hole[ n ] - 1; i++ ) {

			B_J_block_T_bra[ n ][ i ] = new int [ B_Num_block_T[ n ][ i ] ];
			B_J_block_T_ket[ n ][ i ] = new int [ B_Num_block_T[ n ][ i ] ];

	                for( int j = 0; j < B_Num_block_T[ n ][ i ]; j++ ) {

	                        B_J_block_T_bra[ n ][ i ][ j ] = 0;  

				B_J_block_T_ket[ n ][ i ][ j ] = 0;

        	        }

		}

	}

	for( int n = 0; n < B_site; n++ ) {

		for( int i = 0; i < B_Block_Num_hole[ n ] - 1; i++ ) 
		if( B_Num_hole_block[ n ][ i + 1 ] - B_Num_hole_block[ n ][ i ] == 1 ) {

			counter = 0;

	                for( int j = 0; j < B_Num_J_hole_block[ n ][ i ]; j++ ) 
	                for( int k = 0; k < B_Num_J_hole_block[ n ][ i + 1 ]; k++ ) 
			if( abs( B_Value_J_block[ n ][ i ][ j ] - B_Value_J_block[ n ][ i + 1 ][ k ] ) == 1 ) {

	                        B_J_block_T_bra[ n ][ i ][ counter ] = j;  

				B_J_block_T_ket[ n ][ i ][ counter ] = k;

				counter++;

        	        }

		}

	}

//------6j coefficients for A-blocks----------------------------------------------------------------------------- 
   //---A_six_j_S_Dia_old----------------------------------------------------------------------------------------
	A_six_j_S_Dia_old = new double **** [ lsys - 1 ];

        for( int i = 0; i < lsys - 1; i++ ) {

                A_six_j_S_Dia_old[i] = new double *** [ A_N_Block_Num_hole[ i + 1 ] ];

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ ) {

			A_six_j_S_Dia_old[i][j] = new double ** [ A_N_Num_J_hole_block[ i + 1 ][ j ] ];

			for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

				A_six_j_S_Dia_old[i][j][k] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					A_six_j_S_Dia_old[i][j][k][l] = new double [ 3 ];

					for( int m = 0; m < 3; m++ ) {

						A_six_j_S_Dia_old[i][j][k][l][m] = 0.0;

					}

				}

			}

		}

		//----
		FILE * fs = fopen( Combine ( Combine ( "6j_factor/S_Dia_old/", 1 ),  i+2 ), "rb" );

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ]; k++ )
		for( int l = 0; l < 3; l++ ) {

			fread( A_six_j_S_Dia_old[i][j][k][l], sizeof(double), 3, fs );

		}
 
		fclose(fs);

        }

   //---A_six_j_S_Dia_n------------------------------------------------------------------------------------------
	A_six_j_S_Dia_n = new double *** [ lsys - 1 ];

	for( int i = 0; i < lsys - 1; i++ ) {

		A_six_j_S_Dia_n[i] = new double ** [ A_N_Block_Num_hole[ i + 1 ] ];

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ ) {

			A_six_j_S_Dia_n[i][j] = new double * [ A_N_Num_J_hole_block[ i + 1 ][ j ] ];

			for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

				A_six_j_S_Dia_n[i][j][k] = new double [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					A_six_j_S_Dia_n[i][j][k][l] = 0.0;

				}

			}

		}

		//----
		FILE *fa = fopen( Combine ( Combine ( "6j_factor/S_Dia_n/", 1 ), i+2 ), "rb" );

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

			fread( A_six_j_S_Dia_n[i][j][k], sizeof(double), 3, fa );

		}

		fclose(fa);

	}

   //---A_six_j_S_M_Dia_old--------------------------------------------------------------------------------------
	A_six_j_S_M_Dia_old = new double **** [ lsys - 1 ];

        for( int i = 0; i < lsys - 1; i++ ) {

                A_six_j_S_M_Dia_old[i] = new double *** [ A_N_Block_Num_hole[ i + 1 ] ];

                for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ ) {

                        A_six_j_S_M_Dia_old[i][j] = new double ** [ A_N_Num_J_hole_block[ i + 1][ j ] - 1 ];

                        for( int k = 0; k < A_N_Num_J_hole_block[ i +  1 ][ j ] - 1; k++ ) {

                                A_six_j_S_M_Dia_old[i][j][k] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					A_six_j_S_M_Dia_old[i][j][k][l] = new double [ 3 ]; 

					for( int m = 0; m < 3; m++ ) {

						A_six_j_S_M_Dia_old[i][j][k][l][m] = 0.0;

					}

				}

			}

                }

		//----
                FILE *fp = fopen( Combine ( Combine ("6j_factor/S_M_Dia_old/", 1 ), i+2 ), "rb" );

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ] - 1; k++ )
		for( int l = 0; l < 3; l++ ) {

                        fread( A_six_j_S_M_Dia_old[i][j][k][l], sizeof(double), 3, fp);

		}

                fclose(fp);

        }

   //---A_six_j_S_M_Dia_n----------------------------------------------------------------------------------------
	A_six_j_S_M_Dia_n = new double **** [ lsys - 1 ];

        for( int i = 0; i < lsys - 1; i++ ) {

	        A_six_j_S_M_Dia_n[i] = new double *** [ A_N_Block_Num_hole[ i + 1 ] ];

                for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ ) {

                        A_six_j_S_M_Dia_n[i][j] = new double ** [ A_N_Num_J_hole_block[ i + 1][ j ] - 1 ];

                        for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ] - 1; k++ ) {

                                A_six_j_S_M_Dia_n[i][j][k] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					A_six_j_S_M_Dia_n[i][j][k][l] = new double [ 3 ]; 

					for( int m = 0; m < 3; m++ ) {

						A_six_j_S_M_Dia_n[i][j][k][l][m] = 0.0;

					}

				}

			}

                }

		//----
                FILE *fp = fopen( Combine ( Combine ("6j_factor/S_M_Dia_n/", 1 ), i+2 ), "rb" );

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ] - 1; k++ )
		for( int l = 0; l < 3; l++ ) {

                        fread( A_six_j_S_M_Dia_n[i][j][k][l], sizeof(double), 3, fp);

		}

                fclose(fp);

        }

   //---A_six_j_T----------------------------------------------------------------------------------------
	A_six_j_T = new double **** [ lsys - 1 ];

        for( int i = 0; i < lsys - 1; i++ ) {

                A_six_j_T[ i ] = new double *** [ A_N_Block_Num_hole[ i + 1 ] - 1 ];

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ] - 1; j++ ) {

			A_six_j_T[ i ][ j ] = new double ** [ A_N_Num_block_T[ i + 1 ][ j ] ];

			for( int k = 0; k < A_N_Num_block_T[ i + 1 ][ j ]; k++ ) {

				A_six_j_T[ i ][ j ][ k ] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					A_six_j_T[ i ][ j ][ k ][ l ] = new double [ 3 ];

					for( int m = 0; m < 3; m++ ) {

						A_six_j_T[ i ][ j ][ k ][ l ][ m ] = 0.0;

					}

				}

			}

		}

		//----
		FILE * fs = fopen( Combine ( Combine ( "6j_factor/T_old/", 1 ),  i+2 ), "rb" );

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ] - 1; j++ )
		for( int k = 0; k < A_N_Num_block_T[ i + 1 ][ j ]; k++ )
		for( int l = 0; l < 3; l++ ) {

			fread( A_six_j_T[ i ][ j ][ k ][ l ], sizeof(double), 3, fs );

		}
 
		fclose(fs);

        }

   //---A_six_j_H------------------------------------------------------------------------------------------------
   	A_six_j_H = new double **** [ lsys - 1 ];

	for( int i = 0; i < lsys - 1; i++ ) {

                A_six_j_H[i] = new double *** [ A_N_Block_Num_hole[ i + 1 ] ];

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ ) {

			A_six_j_H[i][j] = new double ** [ A_N_Num_J_hole_block[ i + 1 ][ j ] ];

			for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

				A_six_j_H[i][j][k] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					A_six_j_H[i][j][k][l] = new double [ 3 ];

					for( int m = 0; m < 3; m++ ) {

						A_six_j_H[i][j][k][l][m] = 0.0;

					}

				}

			}

		}

		//----
		FILE * fs = fopen( Combine ( Combine ( "6j_factor/H/", 1 ),  i+2 ), "rb" );

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ]; k++ )
		for( int l = 0; l < 3; l++ ) {

			fread( A_six_j_H[i][j][k][l], sizeof(double), 3, fs );

		}
 
		fclose(fs);

	}

//------6j coefficients for B-block------------------------------------------------------------------------------
   //---B_six_j_S_Dia_old----------------------------------------------------------------------------------------
	B_six_j_S_Dia_old = new double **** [ B_site - 1 ];

        for( int i = 0; i < B_site - 1; i++ ) {

                B_six_j_S_Dia_old[i] = new double *** [ B_N_Block_Num_hole[ i + 1 ] ];

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ ) {

			B_six_j_S_Dia_old[i][j] = new double ** [ B_N_Num_J_hole_block[ i + 1 ][ j ] ];

			for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

				B_six_j_S_Dia_old[i][j][k] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					B_six_j_S_Dia_old[i][j][k][l] = new double [ 3 ];

					for( int m = 0; m < 3; m++ ) {

                                                B_six_j_S_Dia_old[i][j][k][l][m] = 0.0;

					}

				}

			}

		}

		//----
		FILE * fs = fopen( Combine ( Combine ( "6j_factor/S_Dia_old/", 2 ),  i+2 ), "rb" );

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ]; k++ )
		for( int l = 0; l < 3; l++ ) {

			fread( B_six_j_S_Dia_old[i][j][k][l], sizeof(double), 3, fs );

		}
 
		fclose(fs);

        }

   //---B_six_j_S_Dia_n------------------------------------------------------------------------------------------
	B_six_j_S_Dia_n = new double *** [ B_site - 1 ];

	for( int i = 0; i < B_site - 1; i++ ) {

		B_six_j_S_Dia_n[i] = new double ** [ B_N_Block_Num_hole[ i + 1 ] ];

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ ) {

			B_six_j_S_Dia_n[i][j] = new double * [ B_N_Num_J_hole_block[ i + 1 ][ j ] ];

			for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

				B_six_j_S_Dia_n[i][j][k] = new double [ 3 ];

				for( int l = 0; l < 3; l++ ) {

	                                B_six_j_S_Dia_n[i][j][k][l] = 0.0;

				}

			}

		}

		//----
		FILE *fa = fopen( Combine ( Combine ( "6j_factor/S_Dia_n/", 2 ), i+2 ), "rb" );

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

			fread( B_six_j_S_Dia_n[i][j][k], sizeof(double), 3, fa );

		}

		fclose(fa);

	}

   //---B_six_j_S_M_Dia_old--------------------------------------------------------------------------------------
	B_six_j_S_M_Dia_old = new double **** [ B_site - 1 ];

        for( int i = 0; i < B_site - 1; i++ ) {

                B_six_j_S_M_Dia_old[i] = new double *** [ B_N_Block_Num_hole[ i + 1 ] ];

                for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ ) {

                        B_six_j_S_M_Dia_old[i][j] = new double ** [ B_N_Num_J_hole_block[ i + 1][ j ] - 1 ];

                        for( int k = 0; k < B_N_Num_J_hole_block[ i +  1 ][ j ] - 1; k++ ) {

                                B_six_j_S_M_Dia_old[i][j][k] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					B_six_j_S_M_Dia_old[i][j][k][l] = new double [ 3 ]; 

					for( int m = 0; m < 3; m++ ) {

                                                B_six_j_S_M_Dia_old[i][j][k][l][m] = 0.0;

					}

				}

			}

                }

		//----
                FILE *fp = fopen( Combine ( Combine ("6j_factor/S_M_Dia_old/", 2 ), i+2 ), "rb" );

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < B_N_Num_J_hole_block[ i +  1 ][ j ] - 1; k++ )
		for( int l = 0; l < 3; l++ ) {

                        fread( B_six_j_S_M_Dia_old[i][j][k][l], sizeof(double), 3, fp);

		}

                fclose(fp);

        }

   //---B_six_j_S_M_Dia_n----------------------------------------------------------------------------------------
	B_six_j_S_M_Dia_n = new double **** [ B_site - 1 ];

        for( int i = 0; i < B_site - 1; i++ ) {

	        B_six_j_S_M_Dia_n[i] = new double *** [ B_N_Block_Num_hole[ i + 1 ] ];

                for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ ) {

                        B_six_j_S_M_Dia_n[i][j] = new double ** [ B_N_Num_J_hole_block[ i + 1][ j ] - 1 ];

                        for( int k = 0; k < B_N_Num_J_hole_block[ i +  1 ][ j ] - 1; k++ ) {

                                B_six_j_S_M_Dia_n[i][j][k] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					B_six_j_S_M_Dia_n[i][j][k][l] = new double [ 3 ]; 

                                        for( int m = 0; m < 3; m++ ) {

                                                B_six_j_S_M_Dia_n[i][j][k][l][m] = 0.0;

					}

				}

			}

                }

		//----
                FILE *fp = fopen( Combine ( Combine ("6j_factor/S_M_Dia_n/", 2 ), i+2 ), "rb" );

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ] - 1; k++ )
		for( int l = 0; l < 3; l++ ) {

                        fread( B_six_j_S_M_Dia_n[i][j][k][l], sizeof(double), 3, fp);

		}

                fclose(fp);

        }

   //---B_six_j_T----------------------------------------------------------------------------------------
	B_six_j_T = new double **** [ B_site - 1 ];

        for( int i = 0; i < B_site - 1; i++ ) {

                B_six_j_T[ i ] = new double *** [ B_N_Block_Num_hole[ i + 1 ] - 1 ];

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ] - 1; j++ ) {

			B_six_j_T[ i ][ j ] = new double ** [ B_N_Num_block_T[ i + 1 ][ j ] ];

			for( int k = 0; k < B_N_Num_block_T[ i + 1 ][ j ]; k++ ) {

				B_six_j_T[ i ][ j ][ k ] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					B_six_j_T[ i ][ j ][ k ][ l ] = new double [ 3 ];

					for( int m = 0; m < 3; m++ ) {

						B_six_j_T[ i ][ j ][ k ][ l ][ m ] = 0.0;

					}

				}

			}

		}

		//----
		FILE * fs = fopen( Combine ( Combine ( "6j_factor/T_old/", 2 ),  i+2 ), "rb" );

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ] - 1; j++ )
		for( int k = 0; k < B_N_Num_block_T[ i + 1 ][ j ]; k++ )
		for( int l = 0; l < 3; l++ ) {

			fread( B_six_j_T[ i ][ j ][ k ][ l ], sizeof(double), 3, fs );

		}
 
		fclose(fs);

        }

   //---B_six_j_H------------------------------------------------------------------------------------------------
 	B_six_j_H = new double **** [ B_site - 1 ];

	for( int i = 0; i < B_site - 1; i++ ) {

                B_six_j_H[i] = new double *** [ B_N_Block_Num_hole[ i + 1 ] ];

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ ) {

			B_six_j_H[i][j] = new double ** [ B_N_Num_J_hole_block[ i + 1 ][ j ] ];

			for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

				B_six_j_H[i][j][k] = new double * [ 3 ];

				for( int l = 0; l < 3; l++ ) {

					B_six_j_H[i][j][k][l] = new double [ 3 ];

                                        for( int m = 0; m < 3; m++ ) {

                                                B_six_j_H[i][j][k][l][m] = 0.0;

					}

				}

			}

		}

		//----
		FILE * fs = fopen( Combine ( Combine ( "6j_factor/H/", 2 ),  i+2 ), "rb" );

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ )
		for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ]; k++ )
		for( int l = 0; l < 3; l++ ) {

			fread( B_six_j_H[i][j][k][l], sizeof(double), 3, fs );

		}
 
		fclose(fs);

	}
   
}

//=============================================================================
//			Create space for wavefunctions
//=============================================================================
inline void BondHopping::CreateFunction( Super &sup ) {

	int inc = 1;

	//----
        WaveFunction_block = new double * [ sup.BlockNumber_for_TargetBlock ];

        for( int i = 0; i < sup.BlockNumber_for_TargetBlock; i++ ) {

                WaveFunction_block[i] = new double [ sup.Dim_block[i] ];

                for( int j = 0; j < sup.Dim_block[i]; j++ ) {

                        WaveFunction_block[i][j] = 0.0;

		}

        }

	//----
        FILE *f = fopen( "wavefunction/wavefunction_ground", "r+" );

        for( int i = 0; i < sup.BlockNumber_for_TargetBlock; i++ )
        for( int j = 0; j < sup.Dim_block[i]; j++ ) {

                fscanf( f, "%lf\n", &WaveFunction_block[i][j] );

	}

        fclose(f);

	//----WaveFunction_config_2
        WaveFunction_config_2 = new double * [ sup.BlockNumber_for_TargetBlock ];

        for( int i = 0; i < sup.BlockNumber_for_TargetBlock; i++ ) {

                WaveFunction_config_2[i] = new double [ sup.Dim_block[i] ];

                for( int j = 0; j < sup.Dim_block[i]; j++ ) {

                        WaveFunction_config_2[i][j] = 0.0;

		}

        }

        for( int i = 0; i < sup.BlockNumber_for_TargetBlock; i++ ) {

         	int index = 0;

                for( int j = 0; j < sup.BlockNumber_for_TargetBlock; j++ )
                if( sup.H_sys[i] == sup.H_sys[j]  &&  sup.H_env[i] == sup.H_env[j]  &&  sup.J_sys[i] == sup.J_sys[j]  &&  sup.J_env[i] == sup.J_env[j]  &&  ( sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[i] ] - sup.sys_space -> Num_hole_block[ sup.H_sys[i] ] ) == ( sup.envnew_space -> Num_hole_block[ sup.H_envnew[j] ] - sup.env_space -> Num_hole_block[ sup.H_env[j] ] )  &&  ( sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[j] ] - sup.sys_space -> Num_hole_block[ sup.H_sys[j] ] ) == ( sup.envnew_space -> Num_hole_block[ sup.H_envnew[i] ] - sup.env_space -> Num_hole_block[ sup.H_env[i] ] ) )  {

                        daxpy_( &sup.Dim_block[i], &sup.nine_j_config_2[i][index++], WaveFunction_block[j], &inc, WaveFunction_config_2[i], &inc );

                }

        }

	//----WaveFunction_config_3
        WaveFunction_config_3 = new double * [ sup.BlockNumber_for_TargetBlock_config_3 ];

        for( int i = 0; i < sup.BlockNumber_for_TargetBlock_config_3; i++ ) {

        	WaveFunction_config_3[ i ] = new double [ sup.Dim_block_config_3[ i ] ];

                for( int j = 0; j < sup.Dim_block_config_3[ i ]; j++ ) {

                	WaveFunction_config_3[ i ][ j ] = 0.0; 

		}

	}

        for( int i = 0; i < sup.BlockNumber_for_TargetBlock_config_3; i++ ) {
        
		int index = 0;

                for( int j = 0; j < sup.BlockNumber_for_TargetBlock; j++ )
                if( sup.H_sys[j] == sup.H_sys_config_3[i]  &&  sup.H_env[j] == sup.H_env_config_3[i]  &&  sup.J_sys[j] == sup.J_sys_config_3[i]  &&  sup.J_env[j] == sup.J_env_config_3[i]  &&  sup.H_ns_config_3[i] == ( sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j ] ] - sup.sys_space -> Num_hole_block[ sup.H_sys[ j ] ] )  &&  sup.H_ne_config_3[i] == ( sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j ] ] - sup.env_space -> Num_hole_block[ sup.H_env[ j ] ] ) ) {

			daxpy_( &sup.Dim_block_config_3[i], &sup.nine_j_config_3[i][ index++ ], WaveFunction_block[j], &inc, WaveFunction_config_3[i], &inc);

                }
       }

}

//========================================================================
//				system 
//========================================================================
inline void BondHopping::Correlation_sys( const int &lsys, Parameter &para, Super &sup ) {

	//----Calculate the inner product of two wavefunctions 
        double ***tmp_wavefunction;

        tmp_wavefunction = new double **[ sup.sys_space -> Block_Num_hole ];

        for( int i = 0; i < sup.sys_space -> Block_Num_hole; i++ ) {

                tmp_wavefunction[i] = new double *[ sup.sys_space -> Num_J_hole_block[i] ];

                for( int j = 0; j < sup.sys_space -> Num_J_hole_block[i]; j++ ) {

                        int dim = sup.sys_space -> Dim_J_block[i][j] * sup.sys_space -> Dim_J_block[i][j];

                        tmp_wavefunction[i][j] = new double [dim];

                        for( int k = 0; k < dim; k++ ) {

                                tmp_wavefunction[i][j][k] = (double) 0;

                        }

                }

        }

	//----Initialize tmp_wavefunction
        alpha = 1.0;      beta = 1.0;

        //----
        for( int n = 0; n < sup.sys_space -> Block_Num_hole; n++ )
        for( int i = 0; i < sup.sys_space -> Num_J_hole_block[ n ]; i++ ) {

                for( int j = 0; j < sup.BlockNumber_for_TargetBlock; j++ )
                if( sup.H_sys[ j ] == n  &&  sup.J_sys[ j ] == i ) {

                        int Dim_sys = sup.sys_space -> Dim_J_block[ n ][ i ];

                        int Dim_env = sup.env_space -> Dim_J_block[ sup.H_env[ j ] ][ sup.J_env[ j ] ];

                        dgemm_( &trans_N, &trans_T, &Dim_sys, &Dim_sys, &Dim_env, &alpha, WaveFunction_block[ j ], &Dim_sys, WaveFunction_block[ j ], &Dim_sys, &beta, tmp_wavefunction[ n ][ i ], &Dim_sys );

                }

        }

	//----calculate correlation
	for( int i = 0; i < lsys; i++ )
	if( i >= ref_site  ) {

		//----Initialize T^{ \dagger (1/2) }_i
		if( i == 0 )	New_A_Si_initial( 0 );

		else {

			New_A_Si_initial( i );

			Truncate_A_Si( i );

		}

		//----Enlarge T^{ \dagger (1/2) }_i
		for( int j = i + 1; j < lsys; j++ ) {

        if( para.Table_T[i*para.Total_N + j] == 1 ) {

			//----Initialize T^{ \dagger (1/2) }_i T^{ (1/2) }_j
			Initial_A_SiSj( j );

			Truncate_A_SiSj( j );

                        //----Transform SiSj operator
                        for( int k = j + 1; k < lsys; k++ ) {

	                        New_A_SiSj( k );

                                Truncate_A_SiSj( k );

                        }

			//----correlation value
                        beta = 0.0;

                        for( int n = 0; n < sup.sys_space -> Block_Num_hole; n++ )
                        for( int j_1 = 0; j_1 < sup.sys_space -> Num_J_hole_block[ n ]; j_1++ ) {

				alpha = sqrt( 2.0 / ( sup.sys_space -> Value_J_block[ n ][ j_1 ] + 1.0 ) );

                                int dim = sup.sys_space -> Dim_J_block[ n ][ j_1 ] * sup.sys_space -> Dim_J_block[ n ][ j_1 ];

                                double *matrix = new double [ dim ];

                                for( int p = 0; p < dim; p++ )  matrix[ p ] = (double) 0;

				dgemm_( &trans_N, &trans_N, &sup.sys_space -> Dim_J_block[ n ][ j_1 ], &sup.sys_space -> Dim_J_block[ n ][ j_1 ], &sup.sys_space -> Dim_J_block[ n ][ j_1 ], &alpha, tmp_wavefunction[ n ][ j_1 ], &sup.sys_space -> Dim_J_block[ n ][ j_1 ], SiSj_old_A[ n ][ j_1 ], &sup.sys_space -> Dim_J_block[ n ][ j_1 ], &beta, matrix, &sup.sys_space -> Dim_J_block[ n ][ j_1 ] );

                                for( int p = 0; p < sup.sys_space -> Dim_J_block[ n ][ j_1 ]; p++ ) {

                                        correlation[ i ][ j ] += matrix[ p * ( 1 + sup.sys_space -> Dim_J_block[ n ][ j_1 ] ) ];

                                        correlation[ j ][ i ] += matrix[ p * ( 1 + sup.sys_space -> Dim_J_block[ n ][ j_1 ] ) ];

                                 }

                                 delete [] matrix;

                        }

                        fout<<"\n"<<i<<"\t"<<j<<"\t"<<correlation[ i ][ j ];

                        //----Delete space of SiSj_old
                        for( int k = 0; k < A_Block_Num_hole[ lsys - 1 ]; k++ ) {

                                for( int l = 0; l < A_Num_J_hole_block[ lsys - 1 ][ k ]; l++ ) {

                                        delete [] SiSj_old_A[ k ][ l ];

                                }

                                delete [] SiSj_old_A[ k ];

                        }

                        delete [] SiSj_old_A;

		    }  //!!!!!!

                	//----New A_Si after shift a site to right
                        New_A_Si_new( j );

                        Truncate_A_Si( j );

                }

		//----Compute the correlation functions between sites in system and ns
        	for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock; j_i++ )  //bra
		if( sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_i ] ] != 0  &&  sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_i ] ] != sup.sysnew_space -> TotSiteNo  &&  sup.sys_space -> Num_hole_block[ sup.H_sys[ j_i ] ] != sup.sys_space -> TotSiteNo  &&  sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_i ] ] - sup.sys_space -> Num_hole_block[ sup.H_sys[ j_i ] ] == 1 ) {

	        	for( int j_j = 0; j_j < sup.BlockNumber_for_TargetBlock; j_j++ ) 
			if( sup.H_sysnew[ j_i ] == sup.H_sysnew[ j_j ]  && sup.H_envnew[ j_i ] == sup.H_envnew[ j_j ]  &&  sup.J_envnew[ j_i ] == sup.J_envnew[ j_j ]  &&  sup.J_sysnew[ j_i ] == sup.J_sysnew[ j_j ]  &&  sup.H_env[ j_i ] == sup.H_env[ j_j ]  && sup.J_env[ j_i ] == sup.J_env[ j_j ]  &&  sup.sysnew_space -> Value_J_block[ sup.H_sysnew[ j_i ] ][ sup.J_sysnew[ j_i ] ] == sup.sysnew_space -> Value_J_block[ sup.H_sysnew[ j_j ] ][ sup.J_sysnew[ j_j ] ]  &&  sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_j ] ] == sup.sys_space -> Num_hole_block[ sup.H_sys[ j_j ] ]  &&  abs( sup.sys_space -> Value_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ] - sup.sys_space -> Value_J_block[ sup.H_sys[ j_j ] ][ sup.J_sys[ j_j ] ] ) == 1 ) {

				int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ];

				int Dim_sys_j = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_j ] ][ sup.J_sys[ j_j ] ];

				int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ];

				int Dim_env_j = sup.env_space -> Dim_J_block[ sup.H_env[ j_j ] ][ sup.J_env[ j_j ] ];

				//----
				double *f_1 = new double [ Dim_sys_i * Dim_env_j ];
				for( int l = 0; l < Dim_sys_i * Dim_env_j; l++ )  f_1[ l ] = 0.0;

				double *f_2 = new double [ Dim_sys_i * Dim_sys_i ];
				for( int l = 0; l < Dim_sys_i * Dim_sys_i; l++ )  f_2[ l ] = 0.0;

				//----
				int a_J = 0;

				for( int l = 0; l < A_Num_block_T[ lsys - 1 ][ sup.H_sys[ j_i ] ]; l++ )
				if( A_J_block_T_bra[ lsys - 1 ][ sup.H_sys[ j_i ] ][ l ] ==  sup.J_sys[ j_i ]  &&  A_J_block_T_ket[ lsys - 1 ][ sup.H_sys[ j_i ] ][ l ] == sup.J_sys[ j_j ] )

					a_J = l;

				//----
				alpha = 1.0;  beta = 0.0;

				dgemm_( &trans_N, &trans_N, &Dim_sys_i, &Dim_env_j, &Dim_sys_j, &alpha, T_old_A[ sup.H_sys[ j_i ] ][ a_J ], &Dim_sys_i, WaveFunction_block[ j_j ], &Dim_sys_j, &beta, f_1, &Dim_sys_i );

				alpha =  - 1.0 / sqrt( sup.sysnew_space -> Value_J_block[ sup.H_sysnew[ j_j ] ][ sup.J_sysnew[ j_j ] ] + 1.0 );

	                        //---- fermion phase factor
                                if( ( sup.sys_space -> TotSiteNo - sup.sys_space -> Num_hole_block[ sup.H_sys[ j_j ] ] ) % 2 == 1 )

					alpha = -alpha;

				dgemm_( &trans_N, &trans_T, &Dim_sys_i, &Dim_sys_i, &Dim_env_i, &alpha, f_1, &Dim_sys_i, WaveFunction_block[ j_i ], &Dim_sys_i, &beta, f_2, &Dim_sys_i );

				//----
				for( int l = 0; l < Dim_sys_i; l++ ) {

	                        	correlation[ i ][ lsys ] += f_2[ l * ( 1 + Dim_sys_i ) ];
	
               	                        correlation[ lsys ][ i ] += f_2[ l * ( 1 + Dim_sys_i ) ];

                        	}

				delete [] f_1;  delete [] f_2; 

			}

		}

                fout<<"\n"<<i<<"\t"<<lsys<<"\t"<<correlation[ i ][ lsys ];

		//----Compute the correlation functions between sites in system and ne
        	for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock; j_i++ )  //bra
		if( sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_i ] ] != 0  &&  sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_i ] ] != sup.sysnew_space -> TotSiteNo  &&  sup.sys_space -> Num_hole_block[ sup.H_sys[ j_i ] ] != sup.sys_space -> TotSiteNo  &&  sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_i ] ] - sup.sys_space -> Num_hole_block[ sup.H_sys[ j_i ] ] == 1 ) {

	        	for( int j_j = 0; j_j < sup.BlockNumber_for_TargetBlock; j_j++ ) 
			if( sup.H_sysnew[ j_i ] == sup.H_sysnew[ j_j ]  && sup.H_envnew[ j_i ] == sup.H_envnew[ j_j ]  &&  sup.J_envnew[ j_i ] == sup.J_envnew[ j_j ]  &&  sup.J_sysnew[ j_i ] == sup.J_sysnew[ j_j ]  &&  sup.H_env[ j_i ] == sup.H_env[ j_j ]  && sup.J_env[ j_i ] == sup.J_env[ j_j ]  &&  sup.sysnew_space -> Value_J_block[ sup.H_sysnew[ j_i ] ][ sup.J_sysnew[ j_i ] ] == sup.sysnew_space -> Value_J_block[ sup.H_sysnew[ j_j ] ][ sup.J_sysnew[ j_j ] ]  &&  sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_j ] ] == sup.sys_space -> Num_hole_block[ sup.H_sys[ j_j ] ]  &&  abs( sup.sys_space -> Value_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ] - sup.sys_space -> Value_J_block[ sup.H_sys[ j_j ] ][ sup.J_sys[ j_j ] ] ) == 1 ) {

				int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ];

				int Dim_sys_j = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_j ] ][ sup.J_sys[ j_j ] ];

				int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ];

				int Dim_env_j = sup.env_space -> Dim_J_block[ sup.H_env[ j_j ] ][ sup.J_env[ j_j ] ];

				//----
				double *f_1 = new double [ Dim_sys_i * Dim_env_j ];
				for( int l = 0; l < Dim_sys_i * Dim_env_j; l++ )  f_1[ l ] = 0.0;

				double *f_2 = new double [ Dim_sys_i * Dim_sys_i ];
				for( int l = 0; l < Dim_sys_i * Dim_sys_i; l++ )  f_2[ l ] = 0.0;

				//----
				int a_J = 0;

				for( int l = 0; l < A_Num_block_T[ lsys - 1 ][ sup.H_sys[ j_i ] ]; l++ )
				if( A_J_block_T_bra[ lsys - 1 ][ sup.H_sys[ j_i ] ][ l ] ==  sup.J_sys[ j_i ]  &&  A_J_block_T_ket[ lsys - 1 ][ sup.H_sys[ j_i ] ][ l ] == sup.J_sys[ j_j ] )

					a_J = l;

				//----
				alpha = 1.0;  beta = 0.0;

				dgemm_( &trans_N, &trans_N, &Dim_sys_i, &Dim_env_j, &Dim_sys_j, &alpha, T_old_A[ sup.H_sys[ j_i ] ][ a_J ], &Dim_sys_i, WaveFunction_config_2[ j_j ], &Dim_sys_j, &beta, f_1, &Dim_sys_i );

				alpha =  - 1.0 / sqrt( sup.sysnew_space -> Value_J_block[ sup.H_sysnew[ j_j ] ][ sup.J_sysnew[ j_j ] ] + 1.0 );

	                        //---- fermion phase factor
                                if( ( sup.sys_space -> TotSiteNo - sup.sys_space -> Num_hole_block[ sup.H_sys[ j_j ] ] ) % 2 == 1 )

					alpha = -alpha;

				dgemm_( &trans_N, &trans_T, &Dim_sys_i, &Dim_sys_i, &Dim_env_i, &alpha, f_1, &Dim_sys_i, WaveFunction_config_2[ j_i ], &Dim_sys_i, &beta, f_2, &Dim_sys_i );

				//----
				for( int l = 0; l < Dim_sys_i; l++ ) {

	                        	correlation[ i ][ lsys + 1 ] += f_2[ l * ( 1 + Dim_sys_i ) ];
	
               	                        correlation[ lsys + 1 ][ i ] += f_2[ l * ( 1 + Dim_sys_i ) ];

                        	}

				delete [] f_1;  delete [] f_2; 

			}

		}

                fout<<"\n"<<i<<"\t"<<lsys+1<<"\t"<<correlation[ i ][ lsys + 1 ];

		//----Free Ti
		for( int n = 0; n < A_Block_Num_hole[ lsys - 1 ] - 1; n++ ) 
		if( A_Num_hole_block[ lsys - 1 ][ n + 1 ] - A_Num_hole_block[ lsys - 1 ][ n ] == 1 ) {

       			for( int m = 0; m < A_Num_block_T[ lsys - 1 ][ n ]; m++ ) {

			        delete [] T_old_A[ n ][ m ];

			}

			delete [] T_old_A[ n ];

		}

	        delete [] T_old_A;

	}
	
	//----free space
        for( int n = 0; n < sup.sys_space -> Block_Num_hole; n++ ) {

                for( int i = 0; i < sup.sys_space -> Num_J_hole_block[ n ]; i++ ) {

                        delete [] tmp_wavefunction[ n ][ i ];

                }

                delete [] tmp_wavefunction[ n ];

        }

        delete [] tmp_wavefunction;

}

//========================================================================
//				environment
//========================================================================
inline void BondHopping::Correlation_env( const int &lsys, Parameter &para, Super &sup ) {

	//----Calculate the inner product of two wavefunctions 
        double ***tmp_wavefunction;

        tmp_wavefunction = new double **[ sup.env_space -> Block_Num_hole ];

        for( int i = 0; i < sup.env_space -> Block_Num_hole; i++ ) {

                tmp_wavefunction[i] = new double *[ sup.env_space -> Num_J_hole_block[i] ];

                for( int j = 0; j < sup.env_space -> Num_J_hole_block[i]; j++ ) {

                        int dim = sup.env_space -> Dim_J_block[i][j] * sup.env_space -> Dim_J_block[i][j];

                        tmp_wavefunction[i][j] = new double [dim];

                        for( int k = 0; k < dim; k++ ) {

                                tmp_wavefunction[i][j][k] = (double) 0;

                        }

                }

        }

	//----Initialize tmp_wavefunction
        alpha = 1.0;      beta = 1.0;

        //----
        for( int n = 0; n < sup.env_space -> Block_Num_hole; n++ )
        for( int i = 0; i < sup.env_space -> Num_J_hole_block[ n ]; i++ ) {

                for( int j = 0; j < sup.BlockNumber_for_TargetBlock; j++ )
                if( sup.H_env[ j ] == n  &&  sup.J_env[ j ] == i ) {

                        int Dim_env = sup.env_space -> Dim_J_block[ n ][ i ];

                        int Dim_sys = sup.sys_space -> Dim_J_block[ sup.H_sys[ j ] ][ sup.J_sys[ j ] ];

                        dgemm_( &trans_T, &trans_N, &Dim_env, &Dim_env, &Dim_sys, &alpha, WaveFunction_block[ j ], &Dim_sys, WaveFunction_block[ j ], &Dim_sys, &beta, tmp_wavefunction[ n ][ i ], &Dim_env );

                }

        }

	//----calculate correlation
	for( int i = 0; i < lsys; i++ )
	if( i >= ref_site ) {

		//----Initialize T^{ \dagger (1/2) }_i
		if( i == 0 )	New_B_Si_initial( 0 );

		else {

			New_B_Si_initial( i );

			Truncate_B_Si( i );

		}

		//----Enlarge T^{ \dagger (1/2) }_i
		for( int j = i + 1; j < lsys; j++ ) {

        if( para.Table_T[ para.Total_N * (para.Total_N-1-i) + para.Total_N-1-j] == 1 ) {

			//----Initialize T^{ \dagger (1/2) }_i T^{ (1/2) }_j
			Initial_B_SiSj( j );

			Truncate_B_SiSj( j );

                        //----Transform SiSj operator
                        for( int k = j + 1; k < lsys; k++ ) {

	                        New_B_SiSj( k );

                                Truncate_B_SiSj( k );

                        }

			//----correlation value
                        beta = 0.0;

                        for( int n = 0; n < sup.env_space -> Block_Num_hole; n++ )
                        for( int j_1 = 0; j_1 < sup.env_space -> Num_J_hole_block[ n ]; j_1++ ) {

				alpha = sqrt( 2.0 / ( sup.env_space -> Value_J_block[ n ][ j_1 ] + 1.0 ) );

                                int dim = sup.env_space -> Dim_J_block[ n ][ j_1 ] * sup.env_space -> Dim_J_block[ n ][ j_1 ];

                                double *matrix = new double [ dim ];

                                for( int p = 0; p < dim; p++ )  matrix[ p ] = (double) 0;

				dgemm_( &trans_T, &trans_N, &sup.env_space -> Dim_J_block[ n ][ j_1 ], &sup.env_space -> Dim_J_block[ n ][ j_1 ], &sup.env_space -> Dim_J_block[ n ][ j_1 ], &alpha, tmp_wavefunction[ n ][ j_1 ], &sup.env_space -> Dim_J_block[ n ][ j_1 ], SiSj_old_B[ n ][ j_1 ], &sup.env_space -> Dim_J_block[ n ][ j_1 ], &beta, matrix, &sup.env_space -> Dim_J_block[ n ][ j_1 ] );

                                for( int p = 0; p < sup.env_space -> Dim_J_block[ n ][ j_1 ]; p++ ) {

                                        correlation[ para.Total_N - 1 - i ][ para.Total_N - 1 - j ] += matrix[ p * ( 1 + sup.env_space -> Dim_J_block[ n ][ j_1 ] ) ];

                                        correlation[ para.Total_N - 1 - j ][ para.Total_N - 1 - i ] += matrix[ p * ( 1 + sup.env_space -> Dim_J_block[ n ][ j_1 ] ) ];

                                 }

                                 delete [] matrix;

                        }

                        fout<<"\n"<<para.Total_N - 1 - j<<"\t"<<para.Total_N - 1 - i<<"\t"<<correlation[ para.Total_N - 1 - j ][ para.Total_N - 1 - i ];

                        //----Delete space of SiSj_old
                        for( int k = 0; k < B_Block_Num_hole[ lsys - 1 ]; k++ ) {

                                for( int l = 0; l < B_Num_J_hole_block[ lsys - 1 ][ k ]; l++ ) {

                                        delete [] SiSj_old_B[ k ][ l ];

                                }

                                delete [] SiSj_old_B[ k ];

                        }

                        delete [] SiSj_old_B;

		    }  //!!!!!!

                	//----New B_Si after shift a site to right
                        New_B_Si_new( j );

                        Truncate_B_Si( j );

                }

		//----Compute the correlation functions between sites in env  and ne
        	for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock; j_i++ )  //bra
		if( sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j_i ] ] != 0  &&  sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j_i ] ] != sup.envnew_space -> TotSiteNo  &&  sup.env_space -> Num_hole_block[ sup.H_env[ j_i ] ] != sup.env_space -> TotSiteNo  &&  sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j_i ] ] - sup.env_space -> Num_hole_block[ sup.H_env[ j_i ] ] == 1 ) {

	        	for( int j_j = 0; j_j < sup.BlockNumber_for_TargetBlock; j_j++ ) 
			if( sup.H_sysnew[ j_i ] == sup.H_sysnew[ j_j ]  && sup.H_envnew[ j_i ] == sup.H_envnew[ j_j ]  &&  sup.J_envnew[ j_i ] == sup.J_envnew[ j_j ]  &&  sup.J_sysnew[ j_i ] == sup.J_sysnew[ j_j ]  &&  sup.H_sys[ j_i ] == sup.H_sys[ j_j ]  && sup.J_sys[ j_i ] == sup.J_sys[ j_j ]  &&  abs( sup.env_space -> Value_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ] - sup.env_space -> Value_J_block[ sup.H_env[ j_j ] ][ sup.J_env[ j_j ] ] ) == 1 ) {

				int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ];

				int Dim_sys_j = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_j ] ][ sup.J_sys[ j_j ] ];

				int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ];

				int Dim_env_j = sup.env_space -> Dim_J_block[ sup.H_env[ j_j ] ][ sup.J_env[ j_j ] ];

				//----
				double *f_1 = new double [ Dim_env_i * Dim_sys_j ];
				for( int l = 0; l < Dim_sys_j * Dim_env_i; l++ )  f_1[ l ] = 0.0;

				double *f_2 = new double [ Dim_env_i * Dim_env_i ];
				for( int l = 0; l < Dim_env_i * Dim_env_i; l++ )  f_2[ l ] = 0.0;

				//----
				int a_J = 0;

				for( int l = 0; l < B_Num_block_T[ lsys - 1 ][ sup.H_env[ j_i ] ]; l++ )
				if( B_J_block_T_bra[ lsys - 1 ][ sup.H_env[ j_i ] ][ l ] ==  sup.J_env[ j_i ]  &&  B_J_block_T_ket[ lsys - 1 ][ sup.H_env[ j_i ] ][ l ] == sup.J_env[ j_j ] )

					a_J = l;

				//----
				alpha = 1.0;  beta = 0.0;

				dgemm_( &trans_N, &trans_T, &Dim_env_i, &Dim_sys_j, &Dim_env_j, &alpha, T_old_B[ sup.H_env[ j_i ] ][ a_J ], &Dim_env_i, WaveFunction_block[ j_j ], &Dim_sys_j, &beta, f_1, &Dim_env_i );

				alpha =  - 1.0 / sqrt( sup.envnew_space -> Value_J_block[ sup.H_envnew[ j_j ] ][ sup.J_envnew[ j_j ] ] + 1.0 );

	                        //---- fermion phase factor
                                if( ( sup.env_space -> TotSiteNo - sup.env_space -> Num_hole_block[ sup.H_env[ j_j ] ] ) % 2 == 1 )

					alpha = -alpha;

				dgemm_( &trans_N, &trans_N, &Dim_env_i, &Dim_env_i, &Dim_sys_i, &alpha, f_1, &Dim_env_i, WaveFunction_block[ j_i ], &Dim_sys_i, &beta, f_2, &Dim_env_i );

				//----
				for( int l = 0; l < Dim_env_i; l++ ) {

	                        	correlation[ para.Total_N - 1 - i ][ para.Total_N - lsys - 1 ] += f_2[ l * ( 1 + Dim_env_i ) ];
	
               	                        correlation[ para.Total_N - lsys - 1 ][ para.Total_N - 1 - i ] += f_2[ l * ( 1 + Dim_env_i ) ];

                        	}

				delete [] f_1;  delete [] f_2; 

			}

		}

                fout<<"\n"<<para.Total_N - lsys -  1<<"\t"<<para.Total_N - 1 - i<<"\t"<<correlation[ para.Total_N - i - 1 ][ para.Total_N - lsys - 1 ];

		//----Compute the correlation functions between sites in env and ns
        	for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock; j_i++ )  //bra
		if( sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j_i ] ] != 0  &&  sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j_i ] ] != sup.envnew_space -> TotSiteNo  &&  sup.env_space -> Num_hole_block[ sup.H_env[ j_i ] ] != sup.env_space -> TotSiteNo  &&  sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j_i ] ] - sup.env_space -> Num_hole_block[ sup.H_env[ j_i ] ] == 1 ) {

	        	for( int j_j = 0; j_j < sup.BlockNumber_for_TargetBlock; j_j++ ) 
			if( sup.H_sysnew[ j_i ] == sup.H_sysnew[ j_j ]  && sup.H_envnew[ j_i ] == sup.H_envnew[ j_j ]  &&  sup.J_envnew[ j_i ] == sup.J_envnew[ j_j ]  &&  sup.J_sysnew[ j_i ] == sup.J_sysnew[ j_j ]  &&  sup.H_sys[ j_i ] == sup.H_sys[ j_j ]  && sup.J_sys[ j_i ] == sup.J_sys[ j_j ]  &&  abs( sup.env_space -> Value_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ] - sup.env_space -> Value_J_block[ sup.H_env[ j_j ] ][ sup.J_env[ j_j ] ] ) == 1 ) {

				int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ];

				int Dim_sys_j = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_j ] ][ sup.J_sys[ j_j ] ];

				int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ];

				int Dim_env_j = sup.env_space -> Dim_J_block[ sup.H_env[ j_j ] ][ sup.J_env[ j_j ] ];

				//----
				double *f_1 = new double [ Dim_env_i * Dim_sys_j ];
				for( int l = 0; l < Dim_sys_j * Dim_env_i; l++ )  f_1[ l ] = 0.0;

				double *f_2 = new double [ Dim_env_i * Dim_env_i ];
				for( int l = 0; l < Dim_env_i * Dim_env_i; l++ )  f_2[ l ] = 0.0;

				//----
				int a_J = 0;

				for( int l = 0; l < B_Num_block_T[ lsys - 1 ][ sup.H_env[ j_i ] ]; l++ )
				if( B_J_block_T_bra[ lsys - 1 ][ sup.H_env[ j_i ] ][ l ] ==  sup.J_env[ j_i ]  &&  B_J_block_T_ket[ lsys - 1 ][ sup.H_env[ j_i ] ][ l ] == sup.J_env[ j_j ] )

					a_J = l;

				//----
				alpha = 1.0;  beta = 0.0;

				dgemm_( &trans_N, &trans_T, &Dim_env_i, &Dim_sys_j, &Dim_env_j, &alpha, T_old_B[ sup.H_env[ j_i ] ][ a_J ], &Dim_env_i, WaveFunction_config_2[ j_j ], &Dim_sys_j, &beta, f_1, &Dim_env_i );

				alpha =  - 1.0 / sqrt( sup.envnew_space -> Value_J_block[ sup.H_envnew[ j_j ] ][ sup.J_envnew[ j_j ] ] + 1.0 );

	                        //---- fermion phase factor
                                if( ( sup.env_space -> TotSiteNo - sup.env_space -> Num_hole_block[ sup.H_env[ j_j ] ] ) % 2 == 1 )

					alpha = -alpha;

				dgemm_( &trans_N, &trans_N, &Dim_env_i, &Dim_env_i, &Dim_sys_i, &alpha, f_1, &Dim_env_i, WaveFunction_config_2[ j_i ], &Dim_sys_i, &beta, f_2, &Dim_env_i );

				//----
				for( int l = 0; l < Dim_env_i; l++ ) {

	                        	correlation[ para.Total_N - 1 - i ][ para.Total_N - lsys - 2 ] += f_2[ l * ( 1 + Dim_env_i ) ];
	
               	                        correlation[ para.Total_N - lsys - 2 ][ para.Total_N - 1 - i ] += f_2[ l * ( 1 + Dim_env_i ) ];

                        	}

				delete [] f_1;  delete [] f_2; 

			}

		}

                fout<<"\n"<<para.Total_N - 2 - lsys<<"\t"<<para.Total_N - 1 - i<<"\t"<<correlation[ para.Total_N - i - 1 ][ para.Total_N - 2 - lsys ];

		//----Free Ti
		for( int n = 0; n < B_Block_Num_hole[ lsys - 1 ] - 1; n++ ) 
		if( B_Num_hole_block[ lsys - 1 ][ n + 1 ] - B_Num_hole_block[ lsys - 1 ][ n ] == 1 ) {

       			for( int m = 0; m < B_Num_block_T[ lsys - 1 ][ n ]; m++ ) {

			        delete [] T_old_B[ n ][ m ];

			}

			delete [] T_old_B[ n ];

		}

	        delete [] T_old_B;

	}
	
	//----free space
        for( int n = 0; n < sup.env_space -> Block_Num_hole; n++ ) {

                for( int i = 0; i < sup.env_space -> Num_J_hole_block[ n ]; i++ ) {

                        delete [] tmp_wavefunction[ n ][ i ];

                }

                delete [] tmp_wavefunction[ n ];

        }

        delete [] tmp_wavefunction;

}

//========================================================================
//				ns and ne
//========================================================================
inline void BondHopping::Correlation_ns_ne( const int &lsys, Parameter &para, Super &sup ) {

	//----
	alpha = 1.0;  beta = 0.0;

	//----
	for( int l = 0; l < sup.BlockNumber_for_TargetBlock_config_3; l++ ) //bra
        if( sup.H_ns_config_3[ l ] == 0  &&  sup.H_ne_config_3[ l ] == 1 ) { 

	        for( int i = 0; i < sup.BlockNumber_for_TargetBlock_config_3; i++ )
		if( sup.H_sys_config_3[ i ] == sup.H_sys_config_3[ l ]  &&  sup.H_env_config_3[ i ] == sup.H_env_config_3[ l ]  &&  sup.H_ns_config_3[ i ]  == 1  &&  sup.H_ne_config_3[ i ]  == 0  &&  sup.J_sys_config_3[ i ] == sup.J_sys_config_3[ l ]  &&  sup.J_env_config_3[ i ] == sup.J_env_config_3[ l ]  &&  sup.J_sysnew_config_3[ i ] == sup.J_sysnew_config_3[ l ] ) {   

			int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys_config_3[ i ] ][ sup.J_sys_config_3[ i ] ];

			int Dim_sys_l = sup.sys_space -> Dim_J_block[ sup.H_sys_config_3[ l ] ][ sup.J_sys_config_3[ l ] ];

			int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env_config_3[ i ] ][ sup.J_env_config_3[ i ] ];

			int Dim_env_l = sup.env_space -> Dim_J_block[ sup.H_env_config_3[ l ] ][ sup.J_env_config_3[ l ] ];

			//----
			double *f = new double [ Dim_sys_i * Dim_sys_l ];
			for( int n = 0; n < Dim_sys_i * Dim_sys_l; n++ )  f[ n ] = 0.0;

			dgemm_( &trans_N, &trans_T, &Dim_sys_i, &Dim_sys_l, &Dim_env_i, &alpha, WaveFunction_config_3[ i ], &Dim_sys_i, WaveFunction_config_3[ l ], &Dim_sys_l, &beta, f, &Dim_sys_i );

			//----
			for( int n = 0; n < Dim_sys_i; n++ ) {

                        	correlation[ lsys ][ lsys + 1 ] += f[ n * ( 1 + Dim_sys_i ) ];
	
    	                        correlation[ lsys + 1 ][ lsys ] += f[ n * ( 1 + Dim_sys_i ) ];

                       	}

			delete [] f;

		}

	}

	fout<<"\n"<<lsys<<"\t"<<lsys+1<<"\t"<<correlation[ lsys ][ lsys + 1 ];

}

//========================================================================
//		system and environment pairing correlation
//========================================================================
inline void BondHopping::Correlation_sys_env( const int &lsys, Parameter &para, Super &sup ) {

	for( int i = 0; i < lsys ; i++ )
	if( i >= ref_site ) {

		//----Initialize T^{ \dagger (1/2) }_i
		if( i == 0 )	New_A_Si_initial( 0 );

		else {

			New_A_Si_initial( i );
			Truncate_A_Si( i );

		}

		//----Enlarge T^{ \dagger (1/2) }_i
		for( int n = i + 1; n < lsys; n++ ) {

			New_A_Si_new( n );
			Truncate_A_Si( n );

		}

		for( int k = 0; k < para.Total_N - 2 - lsys; k++ )
            if( k >= ref_site && para.Table_T[para.Total_N * i + para.Total_N-1-k] == 1 ) {

			//----Initialize T^{ (1/2) }_k
			if( k == 0 )	New_B_Si_initial( 0 );

			else {

				New_B_Si_initial( k );
				Truncate_B_Si( k );

			}

			//----Enlarge T^{ (1/2) }_k
			for( int n = k + 1; n < para.Total_N - 2 - lsys; n++ ) {

				New_B_Si_new( n );
				Truncate_B_Si( n );

			}

		//----computing the correlation
        	        for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock_config_3; j_i++ ) 
			if( ( sup.sys_space -> Num_hole_block[ sup.H_sys_config_3[ j_i ] ] + sup.env_space -> Num_hole_block[ sup.H_env_config_3[ j_i ] ] ) != 0  &&  ( sup.sys_space -> Num_hole_block[ sup.H_sys_config_3[ j_i ] ] + sup.env_space -> Num_hole_block[ sup.H_env_config_3[ j_i ] ] ) != ( sup.sys_space -> TotSiteNo + sup.env_space -> TotSiteNo ) ) {  //contribution to the operator

		                for( int j_j = 0; j_j < sup.BlockNumber_for_TargetBlock_config_3; j_j++ )
				if( ( sup.sys_space -> Num_hole_block[ sup.H_sys_config_3[ j_j ] ] + sup.env_space -> Num_hole_block[ sup.H_env_config_3[ j_j ] ] == sup.sys_space -> Num_hole_block[ sup.H_sys_config_3[ j_i ] ] + sup.env_space -> Num_hole_block[ sup.H_env_config_3[ j_i ] ] )  &&  sup.H_ns_config_3[ j_j ] == sup.H_ns_config_3[ j_i ]  &&  sup.H_ne_config_3[ j_j ] == sup.H_ne_config_3[ j_i ]  &&  sup.sys_space -> Num_hole_block[ sup.H_sys_config_3[ j_j ] ] - sup.sys_space -> Num_hole_block[ sup.H_sys_config_3[ j_i ] ] == 1  &&  sup.env_space -> Num_hole_block[ sup.H_env_config_3[ j_i ] ] - sup.env_space -> Num_hole_block[ sup.H_env_config_3[ j_j ] ] == 1  &&  sup.J_sysnew_config_3[ j_i ] == sup.J_sysnew_config_3[ j_j ]  &&  sup.J_envnew_config_3[ j_i ] == sup.J_envnew_config_3[ j_j ]  &&  abs( sup.sys_space -> Value_J_block[ sup.H_sys_config_3[ j_i ] ][ sup.J_sys_config_3[ j_i ] ] - sup.sys_space -> Value_J_block[ sup.H_sys_config_3[ j_j ] ][ sup.J_sys_config_3[ j_j ] ] ) == 1  &&  abs( sup.env_space -> Value_J_block[ sup.H_env_config_3[ j_i ] ][ sup.J_env_config_3[ j_i ] ] - sup.env_space -> Value_J_block[ sup.H_env_config_3[ j_j ] ][ sup.J_env_config_3[ j_j ] ] ) == 1 ) {

					int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys_config_3[ j_i ] ][ sup.J_sys_config_3[ j_i ] ];

					int Dim_sys_j = sup.sys_space -> Dim_J_block[ sup.H_sys_config_3[ j_j ] ][ sup.J_sys_config_3[ j_j ] ];

					int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env_config_3[ j_i ] ][ sup.J_env_config_3[ j_i ] ];

					int Dim_env_j = sup.env_space -> Dim_J_block[ sup.H_env_config_3[ j_j ] ][ sup.J_env_config_3[ j_j ] ];

					//----
					double six_j = - pow( -1.0, ( sup.sys_space -> Value_J_block[ sup.H_sys_config_3[ j_j ] ][ sup.J_sys_config_3[ j_j ] ] + sup.env_space -> Value_J_block[ sup.H_env_config_3[ j_i ] ][ sup.J_env_config_3[ j_i ] ] + sup.J_sysnew_config_3[ j_i ] + 1 ) / 2 ) * gsl_sf_coupling_6j( sup.sys_space -> Value_J_block[ sup.H_sys_config_3[ j_i ] ][ sup.J_sys_config_3[ j_i ] ], sup.sys_space -> Value_J_block[ sup.H_sys_config_3[ j_j ] ][ sup.J_sys_config_3[ j_j ] ], 1, sup.env_space -> Value_J_block[ sup.H_env_config_3[ j_j ] ][ sup.J_env_config_3[ j_j ] ], sup.env_space -> Value_J_block[ sup.H_env_config_3[ j_i ] ][ sup.J_env_config_3[ j_i ] ], sup.J_sysnew_config_3[ j_i ] );

                        	        if( ( sup.sys_space -> TotSiteNo - sup.sys_space -> Num_hole_block[ sup.H_sys_config_3[ j_i ] ] ) % 2 == 1 )

	                        	        six_j = -six_j;

	                                if( ( sup.env_space -> Value_J_block[ sup.H_env_config_3[ j_i ] ][ sup.J_env_config_3[ j_i ] ] - sup.env_space -> Value_J_block[ sup.H_env_config_3[ j_j ] ][ sup.J_env_config_3[ j_j ] ] ) == 1 )

		                                six_j = - six_j;

					//----
					double *f_1 = new double [ Dim_sys_i * Dim_env_j ];
					for( int l = 0; l < Dim_sys_i * Dim_env_j; l++ )  f_1[ l ] = 0.0;

					double *f_2 = new double [ Dim_sys_i * Dim_env_i ];
					for( int l = 0; l < Dim_sys_i * Dim_env_i; l++ )  f_2[ l ] = 0.0;

					double *f_3 = new double [ Dim_sys_i * Dim_sys_i ];
					for( int l = 0; l < Dim_sys_i * Dim_sys_i; l++ )  f_3[ l ] = 0.0;

					//----
					int a_J = 0;

					for( int l = 0; l < A_Num_block_T[ lsys - 1 ][ sup.H_sys_config_3[ j_i ] ]; l++ )
					if( A_J_block_T_bra[ lsys - 1 ][ sup.H_sys_config_3[ j_i ] ][ l ] ==  sup.J_sys_config_3[ j_i ]  &&  A_J_block_T_ket[ lsys - 1 ][ sup.H_sys_config_3[ j_i ] ][ l ] == sup.J_sys_config_3[ j_j ] )

						a_J = l;

					int b_J = 0;

					for( int l = 0; l < B_Num_block_T[ para.Total_N - lsys - 3 ][ sup.H_env_config_3[ j_j ] ]; l++ )
					if( B_J_block_T_bra[ para.Total_N - lsys - 3 ][ sup.H_env_config_3[ j_j ] ][ l ] ==  sup.J_env_config_3[ j_j ]  &&  B_J_block_T_ket[ para.Total_N - lsys - 3 ][ sup.H_env_config_3[ j_j ] ][ l ] == sup.J_env_config_3[ j_i ] )

						b_J = l;

					//----
					alpha = 1.0;  beta = 0.0;

					dgemm_( &trans_N, &trans_N, &Dim_sys_i, &Dim_env_j, &Dim_sys_j, &alpha, T_old_A[ sup.H_sys_config_3[ j_i ] ][ a_J ], &Dim_sys_i, WaveFunction_config_3[ j_j ], &Dim_sys_j, &beta, f_1, &Dim_sys_i );

					dgemm_( &trans_N, &trans_N, &Dim_sys_i, &Dim_env_i, &Dim_env_j, &alpha, f_1, &Dim_sys_i, T_old_B[ sup.H_env_config_3[ j_j ] ][ b_J ], &Dim_env_j, &beta, f_2, &Dim_sys_i );

					dgemm_( &trans_N, &trans_T, &Dim_sys_i, &Dim_sys_i, &Dim_env_i, &six_j, f_2, &Dim_sys_i, WaveFunction_config_3[ j_i ], &Dim_sys_i, &beta, f_3, &Dim_sys_i );

					//----
					for( int l = 0; l < Dim_sys_i; l++ ) {

		                        	correlation[ i ][ para.Total_N - k - 1 ] += f_3[ l * ( 1 + Dim_sys_i ) ];
	
                	                        correlation[ para.Total_N - k - 1 ][ i ] += f_3[ l * ( 1 + Dim_sys_i ) ];

                        	        }

					delete [] f_1;  delete [] f_2;  delete [] f_3;

				}

			}		

			fout<<"\n"<<i<<"\t"<<para.Total_N - 1 - k<<"\t"<<correlation[ i ][ para.Total_N - 1 - k ];

			//----Free Tj
			for( int n = 0; n < B_Block_Num_hole[ para.Total_N - lsys - 3 ] - 1; n++ ) 
			if( B_Num_hole_block[ para.Total_N - lsys - 3 ][ n + 1 ] - B_Num_hole_block[ para.Total_N - lsys - 3 ][ n ] == 1 ) {

	       			for( int m = 0; m < B_Num_block_T[ para.Total_N - lsys - 3 ][ n ]; m++ ) {

			        	delete [] T_old_B[ n ][ m ];

				}

				delete [] T_old_B[ n ];

			}

		       	delete [] T_old_B;

		}

		//----Free Ti
		for( int n = 0; n < A_Block_Num_hole[ lsys - 1 ] - 1; n++ ) 
		if( A_Num_hole_block[ lsys - 1 ][ n + 1 ] - A_Num_hole_block[ lsys - 1 ][ n ] == 1 ) {

       			for( int m = 0; m < A_Num_block_T[ lsys - 1 ][ n ]; m++ ) {

			        delete [] T_old_A[ n ][ m ];

			}

			delete [] T_old_A[ n ];

		}

	        delete [] T_old_A;

	}
	
}

//=========================================================================================
//                                     New_A_Si_initial
//=========================================================================================
inline void BondHopping::New_A_Si_initial( const int &num ) {

	if( num == 0 ) {

		//----
		T_old_A = new double ** [ 1 ];

		for( int i = 0; i < 1; i++ ) {

			T_old_A[ i ] = new double * [ 1 ];

			for( int j = 0; j < 1; j++ ) {

				T_old_A[ i ][ j ] = new double [ 1 ];

				for( int k = 0; k < 1; k++ ) {

					T_old_A[ i ][ j ][ k ] = hoping;
			
				}

			}

		}

	}

	//----
	else if( num > 0 ) {

		//----
		T_old_A = new double ** [ A_N_Block_Num_hole[ num ] - 1 ];

		for( int i = 0; i < A_N_Block_Num_hole[ num ] - 1; i++ ) 
		if( A_N_Num_hole_block[ num ][ i + 1 ] - A_N_Num_hole_block[ num ][ i ] == 1 ) {

			T_old_A[ i ] = new double * [ A_N_Num_block_T[ num ][ i ] ];

                        for( int j = 0; j < A_N_Num_block_T[ num ][ i ]; j++ ) {

	                        int Dim = A_N_Dim_J_block[ num ][ i ][ A_N_J_block_T_bra[ num ][ i ][ j ] ] * A_N_Dim_J_block[ num ][ i + 1 ][ A_N_J_block_T_ket[ num ][ i ][ j ] ];

		        	T_old_A[ i ][ j ] = new double [ Dim ];

			        for( int k = 0; k < Dim; k++ ) {

			                T_old_A[ i ][ j ][ k ] = (double) 0;

				}

			}
	
		}

		//----
		int oldhm, oldhl, oldJm, oldJl, position, a_new;
		double factor;

		//----
		for( int n = 0; n < A_N_Block_Num_hole[ num ] - 1; n++ )
		for( int i = 0; i < A_N_Num_block_T[ num ][ n ]; i++ ) {

	                for( int m = 0; m < 3; m ++ )  // bra label
                        for( int l = 0; l < 3; l ++ )  // ket label
                        if( m != 0  &&  l == 0  &&  ( oldhm = A_N_Hole_blockOld[ num ][ n ][ m ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] ) != -1  &&  ( oldhl = A_N_Hole_blockOld[ num ][ n + 1 ][ l ][ A_N_J_block_T_ket[ num ][ n ][ i ] ] ) != -1  &&  ( oldJm = A_N_J_blockOld[ num ][ n ][ m ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] ) != -1  &&  ( oldJl = A_N_J_blockOld[ num ][ n + 1 ][ l ][ A_N_J_block_T_ket[ num ][ n ][ i ] ] ) != -1  &&  oldhm == oldhl  &&  oldJm == oldJl ) {

				//----
				position = A_N_Dim_J_block[ num ][ n ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] * A_N_Start[ num ][ n + 1 ][ l ][ A_N_J_block_T_ket[ num ][ n ][ i ] ] + A_N_Start[ num ][ n ][ m ][ A_N_J_block_T_bra[ num ][ n ][ i ] ];

                                factor = - sqrt( A_N_Value_J_block[ num ][ n ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] + 1.0 );

				//---phase factor
                                if( ( num - A_Num_hole_block[ num - 1 ][ oldhl ] ) % 2 == 1 ) {

	                                factor = - factor;

				}

				//----
                                for( int a_old = 0; a_old < A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ]; a_old++ ) {

	                                a_new = position + ( A_N_Dim_J_block[ num ][ n ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] + 1 ) * a_old;

                                        T_old_A[ n ][ i ][ a_new ] += factor;

                                }

			}

		}

	}

}

//=========================================================================================
//                                     New_B_Si_initial
//=========================================================================================
inline void BondHopping::New_B_Si_initial( const int &num ) {

	if( num == 0 ) {

		//----
		T_old_B = new double ** [ 1 ];

		for( int i = 0; i < 1; i++ ) {

			T_old_B[ i ] = new double * [ 1 ];

			for( int j = 0; j < 1; j++ ) {

				T_old_B[ i ][ j ] = new double [ 1 ];

				for( int k = 0; k < 1; k++ ) {

					T_old_B[ i ][ j ][ k ] = hoping;
			
				}

			}

		}

	}

	//----
	else if( num > 0 ) {

		//----
		T_old_B = new double ** [ B_N_Block_Num_hole[ num ] - 1 ];

		for( int i = 0; i < B_N_Block_Num_hole[ num ] - 1; i++ ) 
		if( B_N_Num_hole_block[ num ][ i + 1 ] - B_N_Num_hole_block[ num ][ i ] == 1 ) {

			T_old_B[ i ] = new double * [ B_N_Num_block_T[ num ][ i ] ];

                        for( int j = 0; j < B_N_Num_block_T[ num ][ i ]; j++ ) {

	                        int Dim = B_N_Dim_J_block[ num ][ i ][ B_N_J_block_T_bra[ num ][ i ][ j ] ] * B_N_Dim_J_block[ num ][ i + 1 ][ B_N_J_block_T_ket[ num ][ i ][ j ] ];

		        	T_old_B[ i ][ j ] = new double [ Dim ];

			        for( int k = 0; k < Dim; k++ ) {

			                T_old_B[ i ][ j ][ k ] = (double) 0;

				}

			}
	
		}

		//----
		int oldhm, oldhl, oldJm, oldJl, position, a_new;
		double factor;

		//----
		for( int n = 0; n < B_N_Block_Num_hole[ num ] - 1; n++ )
		for( int i = 0; i < B_N_Num_block_T[ num ][ n ]; i++ ) {

	                for( int m = 0; m < 3; m ++ )  // bra label
                        for( int l = 0; l < 3; l ++ )  // ket label
                        if( m != 0  &&  l == 0  &&  ( oldhm = B_N_Hole_blockOld[ num ][ n ][ m ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] ) != -1  &&  ( oldhl = B_N_Hole_blockOld[ num ][ n + 1 ][ l ][ B_N_J_block_T_ket[ num ][ n ][ i ] ] ) != -1  &&  ( oldJm = B_N_J_blockOld[ num ][ n ][ m ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] ) != -1  &&  ( oldJl = B_N_J_blockOld[ num ][ n + 1 ][ l ][ B_N_J_block_T_ket[ num ][ n ][ i ] ] ) != -1  &&  oldhm == oldhl  &&  oldJm == oldJl ) {

				//----
				position = B_N_Dim_J_block[ num ][ n ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] * B_N_Start[ num ][ n + 1 ][ l ][ B_N_J_block_T_ket[ num ][ n ][ i ] ] + B_N_Start[ num ][ n ][ m ][ B_N_J_block_T_bra[ num ][ n ][ i ] ];

                                factor = - sqrt( B_N_Value_J_block[ num ][ n ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] + 1.0 );

				//---phase factor
                                if( ( num - B_Num_hole_block[ num - 1 ][ oldhl ] ) % 2 == 1 ) {

	                                factor = - factor;

				}

				//----
                                for( int a_old = 0; a_old < B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ]; a_old++ ) {

	                                a_new = position + ( B_N_Dim_J_block[ num ][ n ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] + 1 ) * a_old;

                                        T_old_B[ n ][ i ][ a_new ] += factor;

                                }

			}

		}

	}

}

//=====================================================================================
//                                      New_A_Ti
//=====================================================================================
inline void BondHopping::New_A_Si_new( const int &num ) {

//------Create space and define the operators of Ti
	T_new_A = new double ** [ A_N_Block_Num_hole[ num ] - 1 ];

	for( int i = 0; i < A_N_Block_Num_hole[ num ] - 1; i++ ) 
	if( A_N_Num_hole_block[ num ][ i + 1 ] - A_N_Num_hole_block[ num ][ i ] == 1 ) {

		T_new_A[ i ] = new double * [ A_N_Num_block_T[ num ][ i ] ];

                for( int j = 0; j < A_N_Num_block_T[ num ][ i ]; j++ ) {

                	int Dim = A_N_Dim_J_block[ num ][ i ][ A_N_J_block_T_bra[ num ][ i ][ j ] ] * A_N_Dim_J_block[ num ][ i + 1 ][ A_N_J_block_T_ket[ num ][ i ][ j ] ];

	        	T_new_A[ i ][ j ] = new double [ Dim ];

		        for( int k = 0; k < Dim; k++ ) {

		                T_new_A[ i ][ j ][ k ] = (double) 0;

			}

		}
	
	}

	//----
	int oldhm, oldhl, oldJm, oldJl, oldJ, SubDim, position, a_new, sign;
	double factor;

        for( int n = 0; n < A_N_Block_Num_hole[ num ] - 1; n++ )
	for( int i = 0; i < A_N_Num_block_T[ num ][ n ]; i++ ) {

                for( int m = 0; m < 3; m++ )
	        for( int l = 0; l < 3; l++ )
                if( ( oldhm = A_N_Hole_blockOld[ num ][ n ][ m ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] ) != -1  &&  ( oldhl = A_N_Hole_blockOld[ num ][ n + 1 ][ l ][ A_N_J_block_T_ket[ num ][ n ][ i ] ] ) != -1  &&  ( oldJm = A_N_J_blockOld[ num ][ n ][ m ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] ) != -1  &&  ( oldJl = A_N_J_blockOld[ num ][ n + 1 ][ l ][ A_N_J_block_T_ket[ num ][ n ][ i ] ] ) != -1  &&  ( A_Num_hole_block[ num - 1 ][ oldhl ] - A_Num_hole_block[ num - 1 ][ oldhm ] ) == 1  &&  abs( A_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] - A_Value_J_block[ num - 1 ][ oldhl ][ oldJl ] ) == 1 ) {

			//---- J numbers of the old T-block
			oldJ = 0;

                        for( int old_i = 0; old_i < A_Num_block_T[ num - 1 ][ oldhm ]; old_i++ )
                        if( A_J_block_T_bra[ num - 1 ][ oldhm ][ old_i ] == oldJm  &&  A_J_block_T_ket[ num - 1 ][ oldhm ][ old_i ] == oldJl ) {

	                        oldJ = old_i;

			}

                        //---- sign is the contribution of "1/2 + s" in the "factor"
                        sign = 2;

                        if( m == 0  &&  l == 0 )  sign = 1;

                        //----
                        position = A_N_Dim_J_block[ num ][ n ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] * A_N_Start[ num ][ n + 1 ][ l ][ A_N_J_block_T_ket[ num ][ n ][ i ] ] + A_N_Start[ num ][ n ][ m ][ A_N_J_block_T_bra[ num ][ n ][ i ] ];

                        factor = pow( -1.0, ( A_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] + A_N_Value_J_block[ num ][ n + 1 ][ A_N_J_block_T_ket[ num ][ n ][ i ] ] + sign ) / 2 ) * sqrt( ( A_N_Value_J_block[ num ][ n ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] + 1.0 ) * ( A_N_Value_J_block[ num ][ n + 1 ][ A_N_J_block_T_ket[ num ][ n ][ i ] ] + 1.0 ) ) * A_six_j_T[ num - 1 ][ n ][ i ][ m ][ l ];

                        SubDim = A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * A_Dim_J_block[ num - 1 ][ oldhl ][ oldJl ];

                        //----
                        for( int a_old = 0; a_old < SubDim; a_old++ ) {

	                        a_new = position + A_N_Dim_J_block[ num ][ n ][ A_N_J_block_T_bra[ num ][ n ][ i ] ] * ( a_old / A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_old % A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                                T_new_A[ n ][ i ][ a_new ] += factor * T_old_A[ oldhm ][ oldJ ][ a_old ];

			}

		}

	}

//------Delete T_old
	for( int n = 0; n < A_Block_Num_hole[ num - 1 ] - 1; n++ ) 
	if( A_Num_hole_block[ num - 1 ][ n + 1 ] - A_Num_hole_block[ num - 1 ][ n ] == 1 ) {

        	for( int i = 0; i < A_Num_block_T[ num - 1 ][ n ]; i++ ) {

		        delete [] T_old_A[ n ][ i ];

		}

		delete [] T_old_A[ n ];

	}

        delete [] T_old_A;

//------Create space and define the operators of T
	T_old_A = new double ** [ A_N_Block_Num_hole[ num ] - 1 ];

	for( int i = 0; i < A_N_Block_Num_hole[ num ] - 1; i++ ) 
	if( A_N_Num_hole_block[ num ][ i + 1 ] - A_N_Num_hole_block[ num ][ i ] == 1 ) {

		T_old_A[ i ] = new double * [ A_N_Num_block_T[ num ][ i ] ];

                for( int j = 0; j < A_N_Num_block_T[ num ][ i ]; j++ ) {

                	int Dim = A_N_Dim_J_block[ num ][ i ][ A_N_J_block_T_bra[ num ][ i ][ j ] ] * A_N_Dim_J_block[ num ][ i + 1 ][ A_N_J_block_T_ket[ num ][ i ][ j ] ];

	        	T_old_A[ i ][ j ] = new double [ Dim ];

		        for( int k = 0; k < Dim; k++ ) {

		                T_old_A[ i ][ j ][ k ] = T_new_A[ i ][ j ][ k ];

			}

		}
	
	}

//------Delete T_new
	for( int n = 0; n < A_N_Block_Num_hole[ num ] - 1; n++ ) 
	if( A_N_Num_hole_block[ num ][ n + 1 ] - A_N_Num_hole_block[ num ][ n ] == 1 ) {

        	for( int i = 0; i < A_N_Num_block_T[ num ][ n ]; i++ ) {

		        delete [] T_new_A[ n ][ i ];

		}

		delete [] T_new_A[ n ];

	}

        delete [] T_new_A;

}

//=====================================================================================
//                                      New_B_Ti
//=====================================================================================
inline void BondHopping::New_B_Si_new( const int &num ) {

//------Create space and define the operators of Ti
	T_new_B = new double ** [ B_N_Block_Num_hole[ num ] - 1 ];

	for( int i = 0; i < B_N_Block_Num_hole[ num ] - 1; i++ ) 
	if( B_N_Num_hole_block[ num ][ i + 1 ] - B_N_Num_hole_block[ num ][ i ] == 1 ) {

		T_new_B[ i ] = new double * [ B_N_Num_block_T[ num ][ i ] ];

                for( int j = 0; j < B_N_Num_block_T[ num ][ i ]; j++ ) {

                	int Dim = B_N_Dim_J_block[ num ][ i ][ B_N_J_block_T_bra[ num ][ i ][ j ] ] * B_N_Dim_J_block[ num ][ i + 1 ][ B_N_J_block_T_ket[ num ][ i ][ j ] ];

	        	T_new_B[ i ][ j ] = new double [ Dim ];

		        for( int k = 0; k < Dim; k++ ) {

		                T_new_B[ i ][ j ][ k ] = (double) 0;

			}

		}
	
	}

	//----
	int oldhm, oldhl, oldJm, oldJl, oldJ, SubDim, position, a_new, sign;
	double factor;

        for( int n = 0; n < B_N_Block_Num_hole[ num ] - 1; n++ )
	for( int i = 0; i < B_N_Num_block_T[ num ][ n ]; i++ ) {

                for( int m = 0; m < 3; m++ )
	        for( int l = 0; l < 3; l++ )
                if( ( oldhm = B_N_Hole_blockOld[ num ][ n ][ m ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] ) != -1  &&  ( oldhl = B_N_Hole_blockOld[ num ][ n + 1 ][ l ][ B_N_J_block_T_ket[ num ][ n ][ i ] ] ) != -1  &&  ( oldJm = B_N_J_blockOld[ num ][ n ][ m ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] ) != -1  &&  ( oldJl = B_N_J_blockOld[ num ][ n + 1 ][ l ][ B_N_J_block_T_ket[ num ][ n ][ i ] ] ) != -1  &&  ( B_Num_hole_block[ num - 1 ][ oldhl ] - B_Num_hole_block[ num - 1 ][ oldhm ] ) == 1  &&  abs( B_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] - B_Value_J_block[ num - 1 ][ oldhl ][ oldJl ] ) == 1 ) {

			//---- J numbers of the old T-block
			oldJ = 0;

                        for( int old_i = 0; old_i < B_Num_block_T[ num - 1 ][ oldhm ]; old_i++ )
                        if( B_J_block_T_bra[ num - 1 ][ oldhm ][ old_i ] == oldJm  &&  B_J_block_T_ket[ num - 1 ][ oldhm ][ old_i ] == oldJl ) {

	                        oldJ = old_i;

			}

                        //---- sign is the contribution of "1/2 + s" in the "factor"
                        sign = 2;

                        if( m == 0  &&  l == 0 )  sign = 1;

                        //----
                        position = B_N_Dim_J_block[ num ][ n ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] * B_N_Start[ num ][ n + 1 ][ l ][ B_N_J_block_T_ket[ num ][ n ][ i ] ] + B_N_Start[ num ][ n ][ m ][ B_N_J_block_T_bra[ num ][ n ][ i ] ];

                        factor = pow( -1.0, ( B_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] + B_N_Value_J_block[ num ][ n + 1 ][ B_N_J_block_T_ket[ num ][ n ][ i ] ] + sign ) / 2 ) * sqrt( ( B_N_Value_J_block[ num ][ n ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] + 1.0 ) * ( B_N_Value_J_block[ num ][ n + 1 ][ B_N_J_block_T_ket[ num ][ n ][ i ] ] + 1.0 ) ) * B_six_j_T[ num - 1 ][ n ][ i ][ m ][ l ];

                        SubDim = B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * B_Dim_J_block[ num - 1 ][ oldhl ][ oldJl ];

                        //----
                        for( int a_old = 0; a_old < SubDim; a_old++ ) {

	                        a_new = position + B_N_Dim_J_block[ num ][ n ][ B_N_J_block_T_bra[ num ][ n ][ i ] ] * ( a_old / B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_old % B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                                T_new_B[ n ][ i ][ a_new ] += factor * T_old_B[ oldhm ][ oldJ ][ a_old ];

			}

		}

	}

//------Delete T_old
	for( int n = 0; n < B_Block_Num_hole[ num - 1 ] - 1; n++ ) 
	if( B_Num_hole_block[ num - 1 ][ n + 1 ] - B_Num_hole_block[ num - 1 ][ n ] == 1 ) {

        	for( int i = 0; i < B_Num_block_T[ num - 1 ][ n ]; i++ ) {

		        delete [] T_old_B[ n ][ i ];

		}

		delete [] T_old_B[ n ];

	}

        delete [] T_old_B;

//------Create space and define the operators of T
	T_old_B = new double ** [ B_N_Block_Num_hole[ num ] - 1 ];

	for( int i = 0; i < B_N_Block_Num_hole[ num ] - 1; i++ ) 
	if( B_N_Num_hole_block[ num ][ i + 1 ] - B_N_Num_hole_block[ num ][ i ] == 1 ) {

		T_old_B[ i ] = new double * [ B_N_Num_block_T[ num ][ i ] ];

                for( int j = 0; j < B_N_Num_block_T[ num ][ i ]; j++ ) {

                	int Dim = B_N_Dim_J_block[ num ][ i ][ B_N_J_block_T_bra[ num ][ i ][ j ] ] * B_N_Dim_J_block[ num ][ i + 1 ][ B_N_J_block_T_ket[ num ][ i ][ j ] ];

	        	T_old_B[ i ][ j ] = new double [ Dim ];

		        for( int k = 0; k < Dim; k++ ) {

		                T_old_B[ i ][ j ][ k ] = T_new_B[ i ][ j ][ k ];

			}

		}
	
	}

//------Delete T_new
	for( int n = 0; n < B_N_Block_Num_hole[ num ] - 1; n++ ) 
	if( B_N_Num_hole_block[ num ][ n + 1 ] - B_N_Num_hole_block[ num ][ n ] == 1 ) {

        	for( int i = 0; i < B_N_Num_block_T[ num ][ n ]; i++ ) {

		        delete [] T_new_B[ n ][ i ];

		}

		delete [] T_new_B[ n ];

	}

        delete [] T_new_B;

}

//========================================================================
//                       Truncate A_Si operators
//========================================================================
inline void BondHopping::Truncate_A_Si( const int &num ) {

	int Dim, Dim_old, Dim_old_p, Dim_new, Dim_new_p;
	char trans_N = 'N', trans_T = 'T';
        double alpha = 1.0, beta = 0.0;

//------Create space and define the operators of Ti

	//----
	T_new_A = new double ** [ A_Block_Num_hole[ num ] - 1 ];

	for( int i = 0; i < A_Block_Num_hole[ num ] - 1; i++ ) 
	if( A_Num_hole_block[ num ][ i + 1 ] - A_Num_hole_block[ num ][ i ] == 1 ) {

		T_new_A[ i ] = new double * [ A_Num_block_T[ num ][ i ] ];

                for( int j = 0; j < A_Num_block_T[ num ][ i ]; j++ ) {

                        Dim = A_Dim_J_block[ num ][ i ][ A_J_block_T_bra[ num ][ i ][ j ] ] * A_Dim_J_block[ num ][ i + 1 ][ A_J_block_T_ket[ num ][ i ][ j ] ];

	        	T_new_A[ i ][ j ] = new double [ Dim ];

		        for( int k = 0; k < Dim; k++ ) {

		                T_new_A[ i ][ j ][ k ] = (double) 0;

			}

		}
	
	}

//------Truncate operator
	int old_T, old_N;

        for( int n = 0; n < A_Block_Num_hole[ num ] - 1; n++ )
        for( int i = 0; i < A_Num_block_T[ num ][ n ]; i++ ) {

		old_T = 0;  old_N = 0;

	        for( int n_old = 0; n_old < A_N_Block_Num_hole[ num ] - 1; n_old++ )
                for( int i_old = 0; i_old < A_N_Num_block_T[ num ][ n_old ]; i_old++ )
                if( A_Num_hole_block[ num ][ n ] == A_N_Num_hole_block[ num ][ n_old ]  &&  A_Value_J_block[ num ][ n ][ A_J_block_T_bra[ num ][ n ][ i ] ] == A_N_Value_J_block[ num ][ n_old ][ A_N_J_block_T_bra[ num ][ n_old ][ i_old ] ]  &&  A_Value_J_block[ num ][ n + 1 ][ A_J_block_T_ket[ num ][ n ][ i ] ] == A_N_Value_J_block[ num ][ n_old + 1 ][ A_N_J_block_T_ket[ num ][ n_old ][ i_old ] ] ) {

	                old_T = n_old;    old_N = i_old;

		}

		Dim_old = A_N_Dim_J_block[ num ][ old_T ][ A_N_J_block_T_bra[ num ][ old_T ][ old_N ] ];
		Dim_old_p = A_N_Dim_J_block[ num ][ old_T + 1 ][ A_N_J_block_T_ket[ num ][ old_T ][ old_N ] ];

		Dim_new = A_Dim_J_block[ num ][ n + 1 ][ A_J_block_T_ket[ num ][ n ][ i ] ];
		Dim_new_p = A_Dim_J_block[ num ][ n ][ A_J_block_T_bra[ num ][ n ][ i ] ];

	        Dim = Dim_old * Dim_new;

                double * f_sub = new double [ Dim ];

                for( int j = 0;  j < Dim; j++ )  f_sub[ j ] = (double) 0;

                dgemm_( &trans_N, &trans_N, &Dim_old, &Dim_new, &Dim_old_p, &alpha, T_old_A[ old_T ][ old_N ], &Dim_old, A[ num ][ n + 1 ][ A_J_block_T_ket[ num ][ n ][ i ] ], &Dim_old_p, &beta, f_sub, &Dim_old );

                dgemm_( &trans_T, &trans_N, &Dim_new_p, &Dim_new, &Dim_old, &alpha, A[ num ][ n ][ A_J_block_T_bra[ num ][ n ][ i ] ], &Dim_old, f_sub, &Dim_old, &beta, T_new_A[ n ][ i ], &Dim_new_p );

                delete [] f_sub;

        }

//------Delete T_old
	for( int n = 0; n < A_N_Block_Num_hole[ num ] - 1; n++ ) 
	if( A_N_Num_hole_block[ num ][ n + 1 ] - A_N_Num_hole_block[ num ][ n ] == 1 ) {

        	for( int i = 0; i < A_N_Num_block_T[ num ][ n ]; i++ ) {

		        delete [] T_old_A[ n ][ i ];

		}

		delete [] T_old_A[ n ];

	}

        delete [] T_old_A;

//------Create space and define the operators of Si
        T_old_A = new double ** [ A_Block_Num_hole[ num ] - 1 ];

	for( int n = 0; n < A_Block_Num_hole[ num ] - 1; n++ )
	if( A_Num_hole_block[ num ][ n + 1 ] - A_Num_hole_block[ num ][ n ] == 1 ) {

	        T_old_A[ n ] = new double * [ A_Num_block_T[ num ][ n ] ];

		for( int i = 0; i < A_Num_block_T[ num ][ n ]; i++ ) {

	                Dim = A_Dim_J_block[ num ][ n ][ A_J_block_T_bra[ num ][ n ][ i ] ] * A_Dim_J_block[ num ][ n + 1 ][ A_J_block_T_ket[ num ][ n ][ i ] ];

	                T_old_A[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	T_old_A[ n ][ i ][ j ] = T_new_A[ n ][ i ][ j ];

			}

		}

        }

//------Delete T_new
	for( int n = 0; n < A_Block_Num_hole[ num ] - 1; n++ )
 	if( A_Num_hole_block[ num ][ n + 1 ] - A_Num_hole_block[ num ][ n ] == 1 ) {

        	for( int i = 0; i < A_Num_block_T[ num ][ n ]; i++ ) {

		        delete [] T_new_A[ n ][ i ];

		}

		delete [] T_new_A[ n ];

	}

        delete [] T_new_A;

}

//==========================================================================================
//                                Truncate B_Si operators
//==========================================================================================
inline void BondHopping::Truncate_B_Si( const int &num ) {

	int Dim, Dim_old, Dim_old_p, Dim_new, Dim_new_p;
	char trans_N = 'N', trans_T = 'T';
        double alpha = 1.0, beta = 0.0;

//------Create space and define the operators of Ti

	//----
	T_new_B = new double ** [ B_Block_Num_hole[ num ] - 1 ];

	for( int i = 0; i < B_Block_Num_hole[ num ] - 1; i++ ) 
	if( B_Num_hole_block[ num ][ i + 1 ] - B_Num_hole_block[ num ][ i ] == 1 ) {

		T_new_B[ i ] = new double * [ B_Num_block_T[ num ][ i ] ];

                for( int j = 0; j < B_Num_block_T[ num ][ i ]; j++ ) {

                        Dim = B_Dim_J_block[ num ][ i ][ B_J_block_T_bra[ num ][ i ][ j ] ] * B_Dim_J_block[ num ][ i + 1 ][ B_J_block_T_ket[ num ][ i ][ j ] ];

	        	T_new_B[ i ][ j ] = new double [ Dim ];

		        for( int k = 0; k < Dim; k++ ) {

		                T_new_B[ i ][ j ][ k ] = (double) 0;

			}

		}
	
	}

//------Truncate operator
	int old_T, old_N;

        for( int n = 0; n < B_Block_Num_hole[ num ] - 1; n++ )
        for( int i = 0; i < B_Num_block_T[ num ][ n ]; i++ ) {

		old_T = 0;  old_N = 0;

	        for( int n_old = 0; n_old < B_N_Block_Num_hole[ num ] - 1; n_old++ )
                for( int i_old = 0; i_old < B_N_Num_block_T[ num ][ n_old ]; i_old++ )
                if( B_Num_hole_block[ num ][ n ] == B_N_Num_hole_block[ num ][ n_old ]  &&  B_Value_J_block[ num ][ n ][ B_J_block_T_bra[ num ][ n ][ i ] ] == B_N_Value_J_block[ num ][ n_old ][ B_N_J_block_T_bra[ num ][ n_old ][ i_old ] ]  &&  B_Value_J_block[ num ][ n + 1 ][ B_J_block_T_ket[ num ][ n ][ i ] ] == B_N_Value_J_block[ num ][ n_old + 1 ][ B_N_J_block_T_ket[ num ][ n_old ][ i_old ] ] ) {

	                old_T = n_old;    old_N = i_old;

		}

		Dim_old = B_N_Dim_J_block[ num ][ old_T ][ B_N_J_block_T_bra[ num ][ old_T ][ old_N ] ];
		Dim_old_p = B_N_Dim_J_block[ num ][ old_T + 1 ][ B_N_J_block_T_ket[ num ][ old_T ][ old_N ] ];

		Dim_new = B_Dim_J_block[ num ][ n + 1 ][ B_J_block_T_ket[ num ][ n ][ i ] ];
		Dim_new_p = B_Dim_J_block[ num ][ n ][ B_J_block_T_bra[ num ][ n ][ i ] ];

	        Dim = Dim_old * Dim_new;

                double * f_sub = new double [ Dim ];

                for( int j = 0;  j < Dim; j++ )  f_sub[ j ] = (double) 0;

                dgemm_( &trans_N, &trans_N, &Dim_old, &Dim_new, &Dim_old_p, &alpha, T_old_B[ old_T ][ old_N ], &Dim_old, B[ num ][ n + 1 ][ B_J_block_T_ket[ num ][ n ][ i ] ], &Dim_old_p, &beta, f_sub, &Dim_old );

                dgemm_( &trans_T, &trans_N, &Dim_new_p, &Dim_new, &Dim_old, &alpha, B[ num ][ n ][ B_J_block_T_bra[ num ][ n ][ i ] ], &Dim_old, f_sub, &Dim_old, &beta, T_new_B[ n ][ i ], &Dim_new_p );

                delete [] f_sub;

        }

//------Delete T_old
	for( int n = 0; n < B_N_Block_Num_hole[ num ] - 1; n++ ) 
	if( B_N_Num_hole_block[ num ][ n + 1 ] - B_N_Num_hole_block[ num ][ n ] == 1 ) {

        	for( int i = 0; i < B_N_Num_block_T[ num ][ n ]; i++ ) {

		        delete [] T_old_B[ n ][ i ];

		}

		delete [] T_old_B[ n ];

	}

        delete [] T_old_B;

//------Create space and define the operators of Si
        T_old_B = new double ** [ B_Block_Num_hole[ num ] - 1 ];

	for( int n = 0; n < B_Block_Num_hole[ num ] - 1; n++ )
	if( B_Num_hole_block[ num ][ n + 1 ] - B_Num_hole_block[ num ][ n ] == 1 ) {

	        T_old_B[ n ] = new double * [ B_Num_block_T[ num ][ n ] ];

		for( int i = 0; i < B_Num_block_T[ num ][ n ]; i++ ) {

	                Dim = B_Dim_J_block[ num ][ n ][ B_J_block_T_bra[ num ][ n ][ i ] ] * B_Dim_J_block[ num ][ n + 1 ][ B_J_block_T_ket[ num ][ n ][ i ] ];

	                T_old_B[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	T_old_B[ n ][ i ][ j ] = T_new_B[ n ][ i ][ j ];

			}

		}

        }

//------Delete T_new
	for( int n = 0; n < B_Block_Num_hole[ num ] - 1; n++ )
 	if( B_Num_hole_block[ num ][ n + 1 ] - B_Num_hole_block[ num ][ n ] == 1 ) {

        	for( int i = 0; i < B_Num_block_T[ num ][ n ]; i++ ) {

		        delete [] T_new_B[ n ][ i ];

		}

		delete [] T_new_B[ n ];

	}

        delete [] T_new_B;

}

//========================================================================
//                           Initial_A_SiSj
//========================================================================
inline void BondHopping::Initial_A_SiSj( const int &num ) {

	int oldhm, oldhl, oldJl, oldJm, old_J, Dim, a_new, position;
        double factor;

//------Create space for SiSj_old
        SiSj_old_A = new double ** [ A_N_Block_Num_hole[ num ] ];

	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

	        SiSj_old_A[ n ] = new double * [ A_N_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = A_N_Dim_J_block[ num ][ n ][ i ] * A_N_Dim_J_block[ num ][ n ][ i ];

	                SiSj_old_A[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_old_A[ n ][ i ][ j ] = (double) 0;

			}

		}

        }

//------Initialize SiSj
        for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ )
	for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ )
	for( int m = 0; m < 3; m++ )
	for( int l = 0; l < 3; l++ )
	if( m == 0  &&  l != 0  &&  (oldhm = A_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldhl = A_N_Hole_blockOld[ num ][ n ][ l ][ i ] ) != -1  &&  ( A_Num_hole_block[ num - 1 ][ oldhm ] + 1 ) == A_Num_hole_block[ num - 1 ][ oldhl ]  &&  ( oldJm = A_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJl = A_N_J_blockOld[ num ][ n ][ l ][ i ] ) != -1  &&  abs( A_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] - A_Value_J_block[ num - 1 ][ oldhl ][ oldJl ] ) == 1 ) {

		position = A_N_Dim_J_block[ num ][ n ][ i ] * A_N_Start[ num ][ n ][ l ][ i ] + A_N_Start[ num ][ n ][ m ][ i ];

                factor = - 1.0 / sqrt( 2.0 );

                Dim = A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * A_Dim_J_block[ num - 1 ][ oldhl ][ oldJl ];

                //---- fermion phase factor
                if( ( num - A_Num_hole_block[ num - 1 ][ oldhl ] ) % 2 == 1 )

	                factor = - factor; 

                //----
                for( int old_i = 0; old_i < A_Num_block_T[ num - 1 ][ oldhm ]; old_i++ )
                if( A_J_block_T_bra[ num - 1 ][ oldhm ][ old_i ] == oldJm  &&  A_J_block_T_ket[ num - 1 ][ oldhm ][ old_i ] == oldJl )

                                old_J = old_i;

		//----
                for( int a_sys = 0; a_sys < Dim; a_sys++ ) {

	                a_new = position + A_N_Dim_J_block[ num ][ n ][ i ] * ( a_sys / A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_sys % A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                        SiSj_old_A[ n ][ i ][ a_new ] += factor * T_old_A[ oldhm ][ old_J ][ a_sys ];

                }

	}

}

//========================================================================
//                           Initial_B_SiSj
//========================================================================
inline void BondHopping::Initial_B_SiSj( const int &num ) {

	int oldhm, oldhl, oldJl, oldJm, old_J, Dim, a_new, position;
        double factor;

//------Create space for SiSj_old
        SiSj_old_B = new double ** [ B_N_Block_Num_hole[ num ] ];

	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

	        SiSj_old_B[ n ] = new double * [ B_N_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = B_N_Dim_J_block[ num ][ n ][ i ] * B_N_Dim_J_block[ num ][ n ][ i ];

	                SiSj_old_B[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_old_B[ n ][ i ][ j ] = (double) 0;

			}

		}

        }

//------Initialize SiSj
        for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ )
	for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ )
	for( int m = 0; m < 3; m++ )
	for( int l = 0; l < 3; l++ )
	if( m == 0  &&  l != 0  &&  (oldhm = B_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldhl = B_N_Hole_blockOld[ num ][ n ][ l ][ i ] ) != -1  &&  ( B_Num_hole_block[ num - 1 ][ oldhm ] + 1 ) == B_Num_hole_block[ num - 1 ][ oldhl ]  &&  ( oldJm = B_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJl = B_N_J_blockOld[ num ][ n ][ l ][ i ] ) != -1  &&  abs( B_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] - B_Value_J_block[ num - 1 ][ oldhl ][ oldJl ] ) == 1 ) {

		position = B_N_Dim_J_block[ num ][ n ][ i ] * B_N_Start[ num ][ n ][ l ][ i ] + B_N_Start[ num ][ n ][ m ][ i ];

                factor = - 1.0 / sqrt( 2.0 );

                Dim = B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * B_Dim_J_block[ num - 1 ][ oldhl ][ oldJl ];

                //---- fermion phase factor
                if( ( num - B_Num_hole_block[ num - 1 ][ oldhl ] ) % 2 == 1 )

	                factor = - factor; 

                //----
                for( int old_i = 0; old_i < B_Num_block_T[ num - 1 ][ oldhm ]; old_i++ )
                if( B_J_block_T_bra[ num - 1 ][ oldhm ][ old_i ] == oldJm  &&  B_J_block_T_ket[ num - 1 ][ oldhm ][ old_i ] == oldJl )

                                old_J = old_i;

		//----
                for( int a_sys = 0; a_sys < Dim; a_sys++ ) {

	                a_new = position + B_N_Dim_J_block[ num ][ n ][ i ] * ( a_sys / B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_sys % B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                        SiSj_old_B[ n ][ i ][ a_new ] += factor * T_old_B[ oldhm ][ old_J ][ a_sys ];

                }

	}

}

//===============================================================================
//                             Truncate A_SiSj operator
//===============================================================================
inline void BondHopping::Truncate_A_SiSj( const int &num ) {

	int Dim;
        char trans_N = 'N', trans_T = 'T';
        double alpha = 1.0, beta = 0.0;

//------Create space for A_SiSj_new
        SiSj_new_A = new double ** [ A_Block_Num_hole[ num ] ];

	for( int n = 0; n < A_Block_Num_hole[ num ]; n++ ) {

	        SiSj_new_A[ n ] = new double * [ A_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < A_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = A_Dim_J_block[ num ][ n ][ i ] * A_Dim_J_block[ num ][ n ][ i ];

	                SiSj_new_A[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_new_A[ n ][ i ][ j ] = (double) 0;

			}

		}

        }

//------Truncate SiSj_old
        for( int n = 0; n < A_Block_Num_hole[ num ]; n++ )
        for( int i = 0; i < A_Num_J_hole_block[ num ][ n ]; i++ ) {

	        Dim = A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ] * A_Dim_J_block[ num ][ n ][ i ];

                double * f_sub = new double [ Dim ];

                for( int j = 0;  j < Dim; j++ )  f_sub[ j ] = (double) 0;

                dgemm_( &trans_N, &trans_N, &A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], &A_Dim_J_block[ num ][ n ][ i ], &A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], &alpha, SiSj_old_A[ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], &A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], A[ num ][ n ][ i ], &A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], &beta, f_sub, &A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ] );

                dgemm_( &trans_T, &trans_N, &A_Dim_J_block[ num ][ n ][ i ], &A_Dim_J_block[ num ][ n ][ i ], &A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], &alpha, A[ num ][ n ][ i ], &A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], f_sub, &A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], &beta, SiSj_new_A[ n ][ i ], &A_Dim_J_block[ num ][ n ][ i ] );

                delete [] f_sub;

        }

//------Delete SiSj_old
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

        	for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ ) {

		        delete [] SiSj_old_A[ n ][ i ];

		}

		delete [] SiSj_old_A[ n ];

	}

        delete [] SiSj_old_A;

//------Create new SiSj
        SiSj_old_A = new double ** [ A_Block_Num_hole[ num ] ];

	for( int n = 0; n < A_Block_Num_hole[ num ]; n++ ) {

	        SiSj_old_A[ n ] = new double * [ A_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < A_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = A_Dim_J_block[ num ][ n ][ i ] * A_Dim_J_block[ num ][ n ][ i ];

	                SiSj_old_A[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_old_A[ n ][ i ][ j ] = SiSj_new_A[ n ][ i ][ j ];

			}

		}

        }

//------Delete SiSj_new
	for( int n = 0; n < A_Block_Num_hole[ num ]; n++ ) {

        	for( int i = 0; i < A_Num_J_hole_block[ num ][ n ]; i++ ) {

		        delete [] SiSj_new_A[ n ][ i ];

		}

		delete [] SiSj_new_A[ n ];

	}

        delete [] SiSj_new_A;

}

//===============================================================================
//                             Truncate B_SiSj operator
//===============================================================================
inline void BondHopping::Truncate_B_SiSj( const int &num ) {

	int Dim;
        char trans_N = 'N', trans_T = 'T';
        double alpha = 1.0, beta = 0.0;

//------Create space for B_SiSj_new
        SiSj_new_B = new double ** [ B_Block_Num_hole[ num ] ];

	for( int n = 0; n < B_Block_Num_hole[ num ]; n++ ) {

	        SiSj_new_B[ n ] = new double * [ B_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < B_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = B_Dim_J_block[ num ][ n ][ i ] * B_Dim_J_block[ num ][ n ][ i ];

	                SiSj_new_B[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_new_B[ n ][ i ][ j ] = (double) 0;

			}

		}

        }

//------Truncate SiSj_old
        for( int n = 0; n < B_Block_Num_hole[ num ]; n++ )
        for( int i = 0; i < B_Num_J_hole_block[ num ][ n ]; i++ ) {

	        Dim = B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ] * B_Dim_J_block[ num ][ n ][ i ];

                double * f_sub = new double [ Dim ];

                for( int j = 0;  j < Dim; j++ )  f_sub[ j ] = (double) 0;

                dgemm_( &trans_N, &trans_N, &B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], &B_Dim_J_block[ num ][ n ][ i ], &B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], &alpha, SiSj_old_B[ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], &B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], B[ num ][ n ][ i ], &B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], &beta, f_sub, &B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ] );

                dgemm_( &trans_T, &trans_N, &B_Dim_J_block[ num ][ n ][ i ], &B_Dim_J_block[ num ][ n ][ i ], &B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], &alpha, B[ num ][ n ][ i ], &B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], f_sub, &B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], &beta, SiSj_new_B[ n ][ i ], &B_Dim_J_block[ num ][ n ][ i ] );

                delete [] f_sub;

        }

//------Delete SiSj_old
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

        	for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ ) {

		        delete [] SiSj_old_B[ n ][ i ];

		}

		delete [] SiSj_old_B[ n ];

	}

        delete [] SiSj_old_B;

//------Create new SiSj
        SiSj_old_B = new double ** [ B_Block_Num_hole[ num ] ];

	for( int n = 0; n < B_Block_Num_hole[ num ]; n++ ) {

	        SiSj_old_B[ n ] = new double * [ B_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < B_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = B_Dim_J_block[ num ][ n ][ i ] * B_Dim_J_block[ num ][ n ][ i ];

	                SiSj_old_B[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_old_B[ n ][ i ][ j ] = SiSj_new_B[ n ][ i ][ j ];

			}

		}

        }

//------Delete SiSj_new
	for( int n = 0; n < B_Block_Num_hole[ num ]; n++ ) {

        	for( int i = 0; i < B_Num_J_hole_block[ num ][ n ]; i++ ) {

		        delete [] SiSj_new_B[ n ][ i ];

		}

		delete [] SiSj_new_B[ n ];

	}

        delete [] SiSj_new_B;

}

//===================================================================================
//                                    New A_SiSj
//===================================================================================
inline void BondHopping::New_A_SiSj( const int &num ) {

	int oldhm, oldJm, position, Dim, a_new;

	double factor;

//------Create space for new SiSj
        SiSj_new_A = new double ** [ A_N_Block_Num_hole[ num ] ];

	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

	        SiSj_new_A[ n ] = new double * [ A_N_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = A_N_Dim_J_block[ num ][ n ][ i ] * A_N_Dim_J_block[ num ][ n ][ i ];

	                SiSj_new_A[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_new_A[ n ][ i ][ j ] = (double) 0;

			}

		}

        }

//------New A_SiSj
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ )
	if( A_N_Num_hole_block[ num ][ n ] != num + 1 ) {

	        for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ )
        	for( int m = 0; m < 3; m++ )
	        if( ( oldhm = A_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = A_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1 ) {

			position = ( A_N_Dim_J_block[ num ][ n ][ i ] + 1 ) * A_N_Start[ num ][ n ][ m ][ i ];

        	        Dim = A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                        factor = sqrt( ( A_N_Value_J_block[ num ][ n ][ i ] + 1.0 ) / ( A_Value_J_block[ num -1 ][ oldhm ][ oldJm ] + 1.0 ) );

                	for( int a_old = 0; a_old < Dim; a_old++ ) {

	                	a_new = position + A_N_Dim_J_block[ num ][ n ][ i ] * ( a_old / A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_old % A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                        	SiSj_new_A[ n ][ i ][ a_new ] += factor * SiSj_old_A[ oldhm ][ oldJm ][ a_old ];

                	}

		}

	}

//------Delete SiSj_old
	for( int n = 0; n < A_Block_Num_hole[ num - 1 ]; n++ ) {

        	for( int i = 0; i < A_Num_J_hole_block[ num - 1 ][ n ]; i++ ) {

		        delete [] SiSj_old_A[ n ][ i ];

		}

		delete [] SiSj_old_A[ n ];

	}

        delete [] SiSj_old_A;

//------Create new SiSj
        SiSj_old_A = new double ** [ A_N_Block_Num_hole[ num ] ];

	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

	        SiSj_old_A[ n ] = new double * [ A_N_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = A_N_Dim_J_block[ num ][ n ][ i ] * A_N_Dim_J_block[ num ][ n ][ i ];

	                SiSj_old_A[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_old_A[ n ][ i ][ j ] = SiSj_new_A[ n ][ i ][ j ];

			}

		}

        }

//------Delete SiSj_new
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

        	for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ ) {

		        delete [] SiSj_new_A[ n ][ i ];

		}

		delete [] SiSj_new_A[ n ];

	}

        delete [] SiSj_new_A;

}

//===================================================================================
//                                    New B_SiSj
//===================================================================================
inline void BondHopping::New_B_SiSj( const int &num ) {

	int oldhm, oldJm, position, Dim, a_new;

	double factor;

//------Create space for new SiSj
        SiSj_new_B = new double ** [ B_N_Block_Num_hole[ num ] ];

	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

	        SiSj_new_B[ n ] = new double * [ B_N_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = B_N_Dim_J_block[ num ][ n ][ i ] * B_N_Dim_J_block[ num ][ n ][ i ];

	                SiSj_new_B[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_new_B[ n ][ i ][ j ] = (double) 0;

			}

		}

        }

//------New B_SiSj
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ )
	if( B_N_Num_hole_block[ num ][ n ] != num + 1 ) {

	        for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ )
        	for( int m = 0; m < 3; m++ )
	        if( ( oldhm = B_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = B_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1 ) {

			position = ( B_N_Dim_J_block[ num ][ n ][ i ] + 1 ) * B_N_Start[ num ][ n ][ m ][ i ];

        	        Dim = B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                        factor = sqrt( ( B_N_Value_J_block[ num ][ n ][ i ] + 1.0 ) / ( B_Value_J_block[ num -1 ][ oldhm ][ oldJm ] + 1.0 ) );

                	for( int a_old = 0; a_old < Dim; a_old++ ) {

	                	a_new = position + B_N_Dim_J_block[ num ][ n ][ i ] * ( a_old / B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_old % B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                        	SiSj_new_B[ n ][ i ][ a_new ] += factor * SiSj_old_B[ oldhm ][ oldJm ][ a_old ];

                	}

		}

	}

//------Delete SiSj_old
	for( int n = 0; n < B_Block_Num_hole[ num - 1 ]; n++ ) {

        	for( int i = 0; i < B_Num_J_hole_block[ num - 1 ][ n ]; i++ ) {

		        delete [] SiSj_old_B[ n ][ i ];

		}

		delete [] SiSj_old_B[ n ];

	}

        delete [] SiSj_old_B;

//------Create new SiSj
        SiSj_old_B = new double ** [ B_N_Block_Num_hole[ num ] ];

	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

	        SiSj_old_B[ n ] = new double * [ B_N_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                Dim = B_N_Dim_J_block[ num ][ n ][ i ] * B_N_Dim_J_block[ num ][ n ][ i ];

	                SiSj_old_B[ n ][ i ] = new double [ Dim ];

             		for( int j = 0; j < Dim; j++ ) {

                        	SiSj_old_B[ n ][ i ][ j ] = SiSj_new_B[ n ][ i ][ j ];

			}

		}

        }

//------Delete SiSj_new
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

        	for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ ) {

		        delete [] SiSj_new_B[ n ][ i ];

		}

		delete [] SiSj_new_B[ n ];

	}

        delete [] SiSj_new_B;

}

//=================================================================
//			Delete WaveFunctions
//=================================================================
inline void BondHopping::DeleteFunction( Super &sup ) {

        for( int i = 0; i < sup.BlockNumber_for_TargetBlock; i++ )

                delete [] WaveFunction_block[i];

        delete [] WaveFunction_block;

}

//==================================================================
//			Delete Space
//==================================================================
inline void BondHopping::DeleteSpace( const int &lsys, Parameter &para ) {

	int A_site = lsys;
	int B_site = para.Total_N - lsys - 2;

//------Delete 6j coefficient
   //---A_six_j_S_Dia_old----------------------------------------------------------------------------------------
        for( int i = 0; i < A_site - 1; i++ ) {

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ ) {

			for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

				for( int l = 0; l < 3; l++ ) {

					delete [] A_six_j_S_Dia_old[i][j][k][l];  delete [] A_six_j_H[i][j][k][l];

				}

				delete [] A_six_j_S_Dia_old[i][j][k];  delete [] A_six_j_H[i][j][k];

			}

			delete [] A_six_j_S_Dia_old[i][j];  delete [] A_six_j_H[i][j];

		}

		delete [] A_six_j_S_Dia_old[i];  delete [] A_six_j_H[i];

	}

	delete [] A_six_j_S_Dia_old;  delete [] A_six_j_H;

  //----delete A_six_j_S_Dia_n
	for( int i = 0; i < A_site - 1; i++ ) {

		for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ ) {

			for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ]; k++ )

				delete [] A_six_j_S_Dia_n[i][j][k];

			delete [] A_six_j_S_Dia_n[i][j];

		}

		delete [] A_six_j_S_Dia_n[i];

	}

	delete [] A_six_j_S_Dia_n;

  //----delete A_six_j_S_M_Dia_old(n)
	for( int i = 0; i < A_site - 1; i++ ) {

                for( int j = 0; j < A_N_Block_Num_hole[ i + 1 ]; j++ ) {

                        for( int k = 0; k < A_N_Num_J_hole_block[ i + 1 ][ j ] - 1; k++ ) {

				for( int l = 0; l < 3; l++ ) {

					delete [] A_six_j_S_M_Dia_old[i][j][k][l];  delete [] A_six_j_S_M_Dia_n[i][j][k][l];

				}

				delete [] A_six_j_S_M_Dia_old[i][j][k];  delete [] A_six_j_S_M_Dia_n[i][j][k];

			}

			delete [] A_six_j_S_M_Dia_old[i][j];  delete [] A_six_j_S_M_Dia_n[i][j];

                }

		delete [] A_six_j_S_M_Dia_old[i];  delete [] A_six_j_S_M_Dia_n[i];

	}

	delete [] A_six_j_S_M_Dia_old;  delete [] A_six_j_S_M_Dia_n;

   //---B_six_j_S_Dia_old----------------------------------------------------------------------------------------
        for( int i = 0; i < B_site - 1; i++ ) {

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ ) {

			for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ]; k++ ) {

				for( int l = 0; l < 3; l++ ) {

					delete [] B_six_j_S_Dia_old[i][j][k][l];  delete [] B_six_j_H[i][j][k][l];

				}

				delete [] B_six_j_S_Dia_old[i][j][k];  delete [] B_six_j_H[i][j][k];

			}

			delete [] B_six_j_S_Dia_old[i][j];  delete [] B_six_j_H[i][j];

		}

		delete [] B_six_j_S_Dia_old[i];  delete [] B_six_j_H[i];

	}

	delete [] B_six_j_S_Dia_old;  delete [] B_six_j_H;

  //----delete B_six_j_S_Dia_n
	for( int i = 0; i < B_site - 1; i++ ) {

		for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ ) {

			for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ]; k++ )

				delete [] B_six_j_S_Dia_n[i][j][k];

			delete [] B_six_j_S_Dia_n[i][j];

		}

		delete [] B_six_j_S_Dia_n[i];

	}

	delete [] B_six_j_S_Dia_n;

  //----delete B_six_j_S_M_Dia_old(n)
	for( int i = 0; i < B_site - 1; i++ ) {

                for( int j = 0; j < B_N_Block_Num_hole[ i + 1 ]; j++ ) {

                        for( int k = 0; k < B_N_Num_J_hole_block[ i + 1 ][ j ] - 1; k++ ) {

				for( int l = 0; l < 3; l++ ) {

					delete [] B_six_j_S_M_Dia_old[i][j][k][l];  delete [] B_six_j_S_M_Dia_n[i][j][k][l];

				}

				delete [] B_six_j_S_M_Dia_old[i][j][k];  delete [] B_six_j_S_M_Dia_n[i][j][k];

			}

			delete [] B_six_j_S_M_Dia_old[i][j];  delete [] B_six_j_S_M_Dia_n[i][j];

                }

		delete [] B_six_j_S_M_Dia_old[i];  delete [] B_six_j_S_M_Dia_n[i];

	}

	delete [] B_six_j_S_M_Dia_old;  delete [] B_six_j_S_M_Dia_n;


//------Delete A
  //----delete A_Old_hole, A_Old_J
	for( int i = 0; i < A_site; i++ )

		delete [] A_Old_hole[i];

	delete [] A_Old_hole;

	for( int i = 0; i < A_site; i++ ) {

		for( int j = 0; j < A_Block_Num_hole[i]; j++ )

			delete [] A_Old_J[i][j];

		delete [] A_Old_J[i];

	}

	delete [] A_Old_J;

  //----delete A
	for( int i = 0; i < A_site; i++ ) {

		for( int j = 0; j < A_Block_Num_hole[i]; j++ ) {

			for( int k = 0; k < A_Num_J_hole_block[i][j]; k++ )

				delete [] A[i][j][k];

			delete [] A[i][j];

		}

		delete [] A[i];

	}

	delete [] A;

  //----delete A_Value (Dim)
 	for( int i = 0; i < A_site; i++ ) {

		for( int j = 0; j < A_Block_Num_hole[i]; j++ ) {

			delete [] A_Value_J_block[i][j];  delete [] A_Dim_J_block[i][j];  delete [] A_density_dim[i][j];

		}

		delete [] A_Value_J_block[i];  delete [] A_Dim_J_block[i];  delete [] A_density_dim[i];

	}

	delete [] A_Value_J_block;  delete [] A_Dim_J_block;  delete [] A_density_dim;
 
  //----delete A_Num_hole_block
 	for( int i = 0; i < A_site; i++ ) {

		delete [] A_Num_hole_block[i];  delete [] A_Num_J_hole_block[i];

	}

	delete [] A_Num_hole_block;  delete [] A_Num_J_hole_block;

	delete [] A_Block_Num_hole;

//------Delete B
  //----delete B_Old_hole, B_Old_J
	for( int i = 0; i < B_site; i++ )

		delete [] B_Old_hole[i];

	delete [] B_Old_hole;

	for( int i = 0; i < B_site; i++ ) {

		for( int j = 0; j < B_Block_Num_hole[i]; j++ )

			delete [] B_Old_J[i][j];

		delete [] B_Old_J[i];

	}

	delete [] B_Old_J;

  //----delete B
	for( int i = 0; i < B_site; i++ ) {

		for( int j = 0; j < B_Block_Num_hole[i]; j++ ) {

			for( int k = 0; k < B_Num_J_hole_block[i][j]; k++ )

				delete [] B[i][j][k];

			delete [] B[i][j];

		}

		delete [] B[i];

	}

	delete [] B;

  //----delete B_Value (Dim)
 	for( int i = 0; i < B_site; i++ ) {

		for( int j = 0; j < B_Block_Num_hole[i]; j++ ) {

			delete [] B_Value_J_block[i][j];  delete [] B_Dim_J_block[i][j];  delete [] B_density_dim[i][j];

		}

		delete [] B_Value_J_block[i];  delete [] B_Dim_J_block[i];  delete [] B_density_dim[i];

	}

	delete [] B_Value_J_block;  delete [] B_Dim_J_block;  delete [] B_density_dim;
 
  //----delete B_Num_hole_block
 	for( int i = 0; i < B_site; i++ ) {

		delete [] B_Num_hole_block[i];  delete [] B_Num_J_hole_block[i];

	}

	delete [] B_Num_hole_block;  delete [] B_Num_J_hole_block;

	delete [] B_Block_Num_hole;

//------Delete A_N
  //----delete A_N_holeblockOld,......
	for( int i = 0; i < A_site; i++ ) {

		for( int j = 0; j < A_N_Block_Num_hole[i]; j++ ) {

			for( int k = 0; k < 3; k++ ) {
			
				delete [] A_N_Hole_blockOld[i][j][k];  delete [] A_N_J_blockOld[i][j][k];  delete [] A_N_Start[i][j][k];

			}
			
			delete [] A_N_Hole_blockOld[i][j];  delete [] A_N_J_blockOld[i][j];  delete [] A_N_Start[i][j];

		}

		delete [] A_N_Hole_blockOld[i];  delete [] A_N_J_blockOld[i];  delete [] A_N_Start[i];

	}

	delete [] A_N_Hole_blockOld;  delete [] A_N_J_blockOld;  delete [] A_N_Start;

  //----delete A_N_Value(Dim)
	for( int i = 0; i < A_site; i++ ) {

		for( int j = 0; j < A_N_Block_Num_hole[i]; j++ ) {

			delete [] A_N_Value_J_block[i][j];  delete [] A_N_Dim_J_block[i][j];

		}

		delete [] A_N_Value_J_block[i];  delete [] A_N_Dim_J_block[i];

	}

	delete [] A_N_Value_J_block;  delete [] A_N_Dim_J_block;

  //----delete A_N_Num_hole_block
	for( int i = 0; i < A_site; i++ ) {

		delete [] A_N_Num_hole_block[i];  delete [] A_N_Num_J_hole_block[i];

	}

	delete [] A_N_Num_hole_block;  delete [] A_N_Num_J_hole_block;

  //----delete A_N_Block_Num_hole
	delete [] A_N_Block_Num_hole;

//------Delete B_N
  //----delete B_N_holeblockOld,......
	for( int i = 0; i < B_site; i++ ) {

		for( int j = 0; j < B_N_Block_Num_hole[i]; j++ ) {

			for( int k = 0; k < 3; k++ ) {
			
				delete [] B_N_Hole_blockOld[i][j][k];  delete [] B_N_J_blockOld[i][j][k];  delete [] B_N_Start[i][j][k];

			}
			
			delete [] B_N_Hole_blockOld[i][j];  delete [] B_N_J_blockOld[i][j];  delete [] B_N_Start[i][j];

		}

		delete [] B_N_Hole_blockOld[i];  delete [] B_N_J_blockOld[i];  delete [] B_N_Start[i];

	}

	delete [] B_N_Hole_blockOld;  delete [] B_N_J_blockOld;  delete [] B_N_Start;

  //----delete B_N_Value(Dim)
	for( int i = 0; i < B_site; i++ ) {

		for( int j = 0; j < B_N_Block_Num_hole[i]; j++ ) {

			delete [] B_N_Value_J_block[i][j];  delete [] B_N_Dim_J_block[i][j];

		}

		delete [] B_N_Value_J_block[i];  delete [] B_N_Dim_J_block[i];

	}

	delete [] B_N_Value_J_block;  delete [] B_N_Dim_J_block;

  //----delete B_N_Num_hole_block
	for( int i = 0; i < B_site; i++ ) {

		delete [] B_N_Num_hole_block[i];  delete [] B_N_Num_J_hole_block[i];

	}

	delete [] B_N_Num_hole_block;  delete [] B_N_Num_J_hole_block;

  //----delete B_N_Block_Num_hole
	delete [] B_N_Block_Num_hole;

}

//===============================================================================================================
BondHopping::~BondHopping() {}
//==================================================END==========================================================
