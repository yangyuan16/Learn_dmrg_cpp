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
#include"bonddensity.h"

#define constant sqrt(6.0)*0.5
#define hoping -sqrt(2.0)
#define ref_site 0    // You may change the reference site

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
BondDensity::BondDensity( const int &lsys, Parameter &para, Super &sup ) {

  //----Define variables
	trans_N = 'N';	trans_T = 'T';
	side_R = 'R';	side_L = 'L';
	uplo = 'U';

  //----Create space for transformation matrices A and B, as well as other variables
	CreateSpace( lsys, para ); 
	CreateFunction( sup );	   

	electron_density = new double * [ para.Total_N ];

	for( int n = 0; n < para.Total_N; n++ ) {

		electron_density[ n ] = new double [ para.Total_N ];

		for( int m = 0; m < para.Total_N; m++ ) {

			electron_density[ n ][ m ] = 0.0;

		}

	}

  //----Create filename
        stringstream ss;
        ss<<"measurement_bond_density.dat";
        filename=ss.str();
        fout.precision(12);
        fout.open(filename.c_str()); //start writing file

  //----Calculate the electron density
	Density_sys( lsys, para, sup );
	Density_ns_ne( lsys, para, sup );
	Density_env( para.Total_N - lsys - 2, para, sup );
	Density_sys_env( lsys, para, sup );

        fout.close(); //finish writing file

  //----Delete matrices A and B
	for( int n = 0; n < para.Total_N; n++ ) {

		delete [] electron_density[ n ];

	}

	delete [] electron_density;

	DeleteFunction( sup );

	DeleteSpace( lsys, para );

}

//==============================================================
//		     Create Space for variables
//===============================================================
inline void BondDensity::CreateSpace( const int &lsys, Parameter &para ) {

	int A_site = lsys;
	int B_site = para.Total_N - lsys - 2;

//------Create space for A_N-matrices
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

   //---A_six_j_S_M_Dia_old
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

   //---A_six_j_S_M_Dia_n
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

   //---A_six_j_H
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

//------6j coefficients for B-block
   //---B_six_j_S_Dia_old
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

   //---B_six_j_S_Dia_n
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

   //---B_six_j_S_M_Dia_old
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

   //---B_six_j_S_M_Dia_n
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

   //---B_six_j_H
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

//===========================================================================
//			Create space for wavefunctions
//===========================================================================
inline void BondDensity::CreateFunction( Super &sup ) {

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

        //WaveFunction_config_2
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

//======================================================================
//			Electron density:system 
//======================================================================
inline void BondDensity::Density_sys( const int &lsys, Parameter &para, Super &sup ) {

	//----Calculate the inner product of two wavefunctions for the calculation
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

        //------Calculate spin correlation
	for( int i = 0; i < lsys; i++ )
	if( i >= ref_site ) {

		//----Initialize Ni operators
		if( i == 0 )  A_Initial_Ni(0);

                else  {

                        A_Initial_Ni(i);

                        A_Truncate_Ni(i);

                }

		//----Transform Si operators
                for( int j = i + 1; j < lsys; j++ ) {
                    
                    if( para.Table_N[i*para.Total_N + j] == 1 ) {

	                //----Initialize SiSj operator
                        Initial_A_SiSj( j );

                        Truncate_A_SiSj( j );

                        //----Transform SiSj operator
                        for( int k = j + 1; k < lsys; k++ ) {

                                New_A_SiSj( k );

                                Truncate_A_SiSj( k );

                        }

                        //----Calculate the expectation value of correlation function
                        beta = 0.0;

                        for( int n = 0; n < sup.sys_space -> Block_Num_hole; n++ )
                        for( int j_1 = 0; j_1 < sup.sys_space -> Num_J_hole_block[ n ]; j_1++ ) {

	                        alpha = 1.0 / sqrt( sup.sys_space -> Value_J_block[ n ][ j_1 ] + 1.0 );

                                int dim = sup.sys_space -> Dim_J_block[ n ][ j_1 ] * sup.sys_space -> Dim_J_block[ n ][ j_1 ];

                                double *matrix = new double [ dim ];

                                for( int p = 0; p < dim; p++ )  matrix[ p ] = (double) 0;

                                dsymm_( &side_R, &uplo, &sup.sys_space -> Dim_J_block[ n ][ j_1 ], &sup.sys_space -> Dim_J_block[ n ][ j_1 ], &alpha, SiSj_old_A[ n ][ j_1 ], &sup.sys_space -> Dim_J_block[ n ][ j_1 ], tmp_wavefunction[ n ][ j_1 ], &sup.sys_space -> Dim_J_block[ n ][ j_1 ], &beta, matrix, &sup.sys_space -> Dim_J_block[ n ][ j_1 ] );

                               for( int p = 0; p < sup.sys_space -> Dim_J_block[ n ][ j_1 ]; p++ ) {

                                        electron_density[ i ][ j ] += matrix[ p * ( 1 + sup.sys_space -> Dim_J_block[ n ][ j_1 ] ) ];

                                        electron_density[ j ][ i ] += matrix[ p * ( 1 + sup.sys_space -> Dim_J_block[ n ][ j_1 ] ) ];

                                 }

                                 delete [] matrix;

                        }

                        fout<<"\n"<<i<<"\t"<<j<<"\t"<<electron_density[ i ][ j ];

	                //----Delete space of SiSj_old
                        for( int k = 0; k < A_Block_Num_hole[ lsys - 1 ]; k++ ) {

                                for( int l = 0; l < A_Num_J_hole_block[ lsys - 1 ][ k ]; l++ ) {

                                        delete [] SiSj_old_A[ k ][ l ];

                                }

                                delete [] SiSj_old_A[ k ];

                        }

                        delete [] SiSj_old_A;

	         } //!!!

                    //----New A_Si after shift a site to right
			A_New_Ni( j );

			A_Truncate_Ni( j );

		}

                //----Calculate the expectation value of correlation function
                for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock; j_i++ )  
                if( sup.sys_space -> Num_hole_block[ sup.H_sys[ j_i ] ] != sup.sys_space -> TotSiteNo  &&  sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_i ] ] == sup.sys_space -> Num_hole_block[ sup.H_sys[ j_i ] ] ) {

                	int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ];

                        int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ];

                        //----
                        double *f_1 = new double [ Dim_sys_i * Dim_sys_i ];
                        for( int l = 0; l < Dim_sys_i * Dim_sys_i; l++ )  f_1[ l ] = 0.0;

                        double *f_2 = new double [ Dim_sys_i * Dim_sys_i ];
                        for( int l = 0; l < Dim_sys_i * Dim_sys_i; l++ )  f_2[ l ] = 0.0;

                        //----
                        alpha = 1.0;  beta = 0.0;

                        dgemm_( &trans_N, &trans_T, &Dim_sys_i, &Dim_sys_i, &Dim_env_i, &alpha, WaveFunction_block[ j_i ], &Dim_sys_i, WaveFunction_block[ j_i ], &Dim_sys_i, &beta, f_1, &Dim_sys_i );

                        alpha = 1.0 / sqrt( sup.sys_space -> Value_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ] + 1.0 );

                        dgemm_( &trans_N, &trans_N, &Dim_sys_i, &Dim_sys_i, &Dim_sys_i, &alpha, Ni_old_A[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ], &Dim_sys_i, f_1, &Dim_sys_i, &beta, f_2, &Dim_sys_i );

                        //----
                        for( int l = 0; l < Dim_sys_i; l++ ) {

	                        electron_density[ i ][ lsys ] += f_2[ l * ( 1 + Dim_sys_i ) ];

                                electron_density[ lsys ][ i ] += f_2[ l * ( 1 + Dim_sys_i ) ];

                        }

                        delete [] f_1;  delete [] f_2;

		}

                fout<<"\n"<<i<<"\t"<<lsys<<"\t"<<electron_density[ i ][ lsys ];

                //----Calculate the expectation value of correlation function
                for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock; j_i++ )  
                if( sup.sys_space -> Num_hole_block[ sup.H_sys[ j_i ] ] != sup.sys_space -> TotSiteNo  &&  sup.sysnew_space -> Num_hole_block[ sup.H_sysnew[ j_i ] ] == sup.sys_space -> Num_hole_block[ sup.H_sys[ j_i ] ] ) {

                	int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ];

                        int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ];

                        //----
                        double *f_1 = new double [ Dim_sys_i * Dim_sys_i ];
                        for( int l = 0; l < Dim_sys_i * Dim_sys_i; l++ )  f_1[ l ] = 0.0;

                        double *f_2 = new double [ Dim_sys_i * Dim_sys_i ];
                        for( int l = 0; l < Dim_sys_i * Dim_sys_i; l++ )  f_2[ l ] = 0.0;

                        //----
                        alpha = 1.0;  beta = 0.0;

                        dgemm_( &trans_N, &trans_T, &Dim_sys_i, &Dim_sys_i, &Dim_env_i, &alpha, WaveFunction_config_2[ j_i ], &Dim_sys_i, WaveFunction_config_2[ j_i ], &Dim_sys_i, &beta, f_1, &Dim_sys_i );

                        alpha = 1.0 / sqrt( sup.sys_space -> Value_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ] + 1.0 );

                        dgemm_( &trans_N, &trans_N, &Dim_sys_i, &Dim_sys_i, &Dim_sys_i, &alpha, Ni_old_A[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ], &Dim_sys_i, f_1, &Dim_sys_i, &beta, f_2, &Dim_sys_i );

                        //----
                        for( int l = 0; l < Dim_sys_i; l++ ) {

	                        electron_density[ i ][ lsys + 1 ] += f_2[ l * ( 1 + Dim_sys_i ) ];

                                electron_density[ lsys + 1 ][ i ] += f_2[ l * ( 1 + Dim_sys_i ) ];

                        }

                        delete [] f_1;  delete [] f_2;

		}

                fout<<"\n"<<i<<"\t"<<lsys+1<<"\t"<<electron_density[ i ][ lsys + 1 ];

                //--Delete space of Ni operators
		for( int n = 0; n < A_Block_Num_hole[ lsys - 1 ]; n++ ) {

	        	for( int i = 0; i < A_Num_J_hole_block[ lsys - 1 ][ n ]; i++ )

	                	delete [] Ni_old_A[ n ][ i ];
	
			delete [] Ni_old_A[ n ];

		}

		delete [] Ni_old_A;

	}

  //----delete tmp_wavefunction
        for( int i = 0; i < sup.sys_space -> Block_Num_hole; i++ ) {

		for( int j = 0; j < sup.sys_space -> Num_J_hole_block[i]; j++ ) {

			delete [] tmp_wavefunction[i][j];
		}

		delete [] tmp_wavefunction[i];
        }

	delete [] tmp_wavefunction;

}

//======================================================================
//		Electron density:system and environment
//======================================================================
inline void BondDensity::Density_sys_env( const int &lsys, Parameter &para, Super &sup ) {

	for( int i = 0; i < lsys; i++ )
	if( i >= ref_site ) {

          //----Initialize Ni operators---------------------------------------
		if( i == 0 )  A_Initial_Ni(0);

                else  {

                        A_Initial_Ni(i);
                        A_Truncate_Ni(i);

                }
          
         //-----Transform Ni operators--------------------------------------
         	for( int j = i + 1; j < lsys; j++ ) {

			A_New_Ni(j);
			A_Truncate_Ni(j);

		}

		for( int k = 0; k < para.Total_N - 2 - lsys; k++ )
            if( k >= ref_site && para.Table_N[para.Total_N * i + para.Total_N-1-k] == 1 ) {

		        //----Initialize Ni operators
			if( k == 0 )  B_Initial_Ni(0);

	                else  {

        	                B_Initial_Ni(k);
                	        B_Truncate_Ni(k);

                	}
          
		        //-----Transform Ni operators
	         	for( int l = k + 1; l < para.Total_N - lsys - 2; l++ ) {

				B_New_Ni(l);
				B_Truncate_Ni(l);

			}

			//----compute the average value
			for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock_config_3; j_i++ ) {

				//----
	                        int Dim_sys = sup.sys_space -> Dim_J_block[ sup.H_sys_config_3[ j_i ] ][ sup.J_sys_config_3[ j_i ] ];

                                int Dim_env = sup.env_space -> Dim_J_block[ sup.H_env_config_3[ j_i ] ][ sup.J_env_config_3[ j_i ] ];

                               //----
                                double *f_1 = new double [ Dim_sys * Dim_env ];
                                for( int p = 0; p < Dim_sys * Dim_env; p++ )  f_1[ p ] = 0.0;

                                double *f_2 = new double [ Dim_sys * Dim_sys ];
                                for( int p = 0; p < Dim_sys * Dim_sys; p++ )  f_2[ p ] = 0.0;

                                double *f_3 = new double [ Dim_sys * Dim_sys ];
                                for( int p = 0; p < Dim_sys * Dim_sys; p++ )  f_3[ p ] = 0.0;

                                //----
                                alpha = 1.0;  beta = 0.0;

				double factor =  1.0 / sqrt( ( sup.sys_space -> Value_J_block[ sup.H_sys_config_3[ j_i ] ][ sup.J_sys_config_3[ j_i ] ] + 1.0 ) * ( sup.env_space -> Value_J_block[ sup.H_env_config_3[ j_i ] ][ sup.J_env_config_3[ j_i ] ] + 1.0 ) );

				//
                                dgemm_( &trans_N, &trans_N, &Dim_sys, &Dim_env, &Dim_env, &alpha, WaveFunction_config_3[ j_i ], &Dim_sys, Ni_old_B[ sup.H_env_config_3[ j_i ] ][ sup.J_env_config_3[ j_i ] ], &Dim_env, &beta, f_1, &Dim_sys );

                                dgemm_( &trans_N, &trans_T, &Dim_sys, &Dim_sys, &Dim_env, &alpha, WaveFunction_config_3[ j_i ], &Dim_sys, f_1, &Dim_sys, &beta, f_2, &Dim_sys );

                                dgemm_( &trans_N, &trans_N, &Dim_sys, &Dim_sys, &Dim_sys, &factor, Ni_old_A[ sup.H_sys_config_3[ j_i ] ][ sup.J_sys_config_3[ j_i ] ], &Dim_sys, f_2, &Dim_sys, &beta, f_3, &Dim_sys );

                                //----
                                for( int p = 0; p < Dim_sys; p++ ) {

	                                electron_density[ i ][ para.Total_N - k - 1 ] += f_3[ p * ( 1 + Dim_sys ) ];

                                        electron_density[ para.Total_N - k - 1 ][ i ] += f_3[ p * ( 1 + Dim_sys ) ];

                                }

                                delete [] f_1;  delete [] f_2;  delete [] f_3;

			}

                        fout<<"\n"<<i<<"\t"<<para.Total_N - 1 - k<<"\t"<<electron_density[ i ][ para.Total_N - 1 - k ];

			//----delete Ni_old
			for( int n = 0; n < B_Block_Num_hole[ para.Total_N - lsys - 3 ]; n++ ) {

		        	for( int m = 0; m < B_Num_J_hole_block[ para.Total_N - lsys - 3 ][ n ]; m++ )

	        	        	delete [] Ni_old_B[ n ][ m ];
	
				delete [] Ni_old_B[ n ];

			}

			delete [] Ni_old_B;

		}

	//----delete Ni_old
		for( int n = 0; n < A_Block_Num_hole[ lsys - 1 ]; n++ ) {

	        	for( int m = 0; m < A_Num_J_hole_block[ lsys - 1 ][ n ]; m++ )

	                	delete [] Ni_old_A[ n ][ m ];
	
			delete [] Ni_old_A[ n ];

		}

		delete [] Ni_old_A;

	}

}

//===================================================================
//	Correlation function of sites in environment and ns, ne
//===================================================================
inline void BondDensity::Density_env( const int &lsys, Parameter &para, Super &sup ) {

  //----Calculate the inner product of two wavefunctions for the calculation
        double *** tmp_wavefunction;

        tmp_wavefunction = new double ** [ sup.env_space -> Block_Num_hole ];

        for( int i = 0; i < sup.env_space -> Block_Num_hole; i++ ) {

		tmp_wavefunction[i] = new double * [ sup.env_space -> Num_J_hole_block[i] ];

		for( int j = 0; j < sup.env_space -> Num_J_hole_block[i]; j++ ) {

                	int dim = sup.env_space -> Dim_J_block[i][j] * sup.env_space -> Dim_J_block[i][j];

                	tmp_wavefunction[i][j] = new double [dim];

                	for( int k = 0; k < dim; k++ ) {

                        	tmp_wavefunction[i][j][k] = 0.0;

			}

		}

        }

  //----Initialize tmp_wavefunction
        alpha = 1.0;      beta = 1.0;

	//----
        for( int n = 0; n < sup.env_space -> Block_Num_hole; n++ )
	for( int i = 0; i < sup.env_space -> Num_J_hole_block[n]; i++ ) {

	        for( int j = 0; j < sup.BlockNumber_for_TargetBlock; j++ )
        	if( sup.H_env[j] == n  &&  sup.J_env[j] == i ) {

			int Dim_sys = sup.sys_space -> Dim_J_block[ sup.H_sys[j] ][ sup.J_sys[j] ];

			int Dim_env = sup.env_space -> Dim_J_block[n][i];

        	        dgemm_( &trans_T, &trans_N, &Dim_env, &Dim_env, &Dim_sys, &alpha, WaveFunction_block[j], &Dim_sys, WaveFunction_block[j], &Dim_sys, &beta, tmp_wavefunction[n][i], &Dim_env );

		}

	}

//------Calculate the electron density
	for( int i = 0; i < lsys; i++ )
	if( i >= ref_site ) { 	//Loop of the site number

          //----Initialize Ni operators
		if( i == 0 )  B_Initial_Ni(0);

                else  {

                        B_Initial_Ni(i);
                        B_Truncate_Ni(i);

                }
          
         //-----Transform Ni operators
         	for( int j = i + 1; j < lsys; j++ ) {
                
                if( para.Table_N[ para.Total_N * (para.Total_N-1-i) + para.Total_N-1-j] == 1 ) {

                        //----Initialize SiSj operator
                        Initial_B_SiSj( j );

                        Truncate_B_SiSj( j );

                        //----Transform SiSj operator
                        for( int k = j + 1; k < lsys; k++ ) {

                                New_B_SiSj( k );

                                Truncate_B_SiSj( k );

                        }

			//------Calculate the expectation value of electron density 
			beta = 0.0;

			for( int n = 0; n < sup.env_space -> Block_Num_hole; n++ )
        	        for( int j_1 = 0; j_1 < sup.env_space -> Num_J_hole_block[n]; j_1++ ) {

				alpha = 1.0 / sqrt( sup.env_space -> Value_J_block[ n ][ j_1 ] + 1.0 );

	                        int dim = sup.env_space -> Dim_J_block[n][j_1] * sup.env_space -> Dim_J_block[n][j_1];

        	                double * matrix = new double [ dim ];

                	        for( int p = 0; p < dim; p++ ) 

					matrix[p] = 0.0;

	                        dsymm_( &side_R, &uplo, &sup.env_space -> Dim_J_block[n][j_1], &sup.env_space -> Dim_J_block[n][j_1], &alpha, SiSj_old_B[n][j_1], &sup.env_space -> Dim_J_block[n][j_1], tmp_wavefunction[n][j_1], &sup.env_space -> Dim_J_block[n][j_1], &beta, matrix, &sup.env_space -> Dim_J_block[n][j_1] );

	                        for( int p = 0; p < sup.env_space -> Dim_J_block[n][j_1]; p++ ) {

        	                	electron_density[ para.Total_N - 1 - j ][ para.Total_N - 1 - i ] += matrix[ p * ( 1 + sup.env_space -> Dim_J_block[n][j_1] ) ];

                	        	electron_density[ para.Total_N - 1 - i ][ para.Total_N - 1 - j ] += matrix[ p * ( 1 + sup.env_space -> Dim_J_block[n][j_1] ) ];

				}

	                        delete [] matrix;

       	        	}

			fout<<"\n"<<para.Total_N - j - 1<<"\t"<<para.Total_N - i - 1<<"\t"<<electron_density[ para.Total_N - i - 1 ][ para.Total_N - j - 1 ];

                        //----Delete space of SiSj_old
                        for( int k = 0; k < B_Block_Num_hole[ lsys - 1 ]; k++ ) {

                                for( int l = 0; l < B_Num_J_hole_block[ lsys - 1 ][ k ]; l++ ) {

                                        delete [] SiSj_old_B[ k ][ l ];

                                }

                                delete [] SiSj_old_B[ k ];

                        }

                        delete [] SiSj_old_B;

		    } //!!!

                    //----New B_Si after shift a site to right
                        B_New_Ni( j );

                        B_Truncate_Ni( j );

		}

                //----Calculate the expectation value of correlation function
                for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock; j_i++ )  
                if( sup.env_space -> Num_hole_block[ sup.H_env[ j_i ] ] != sup.env_space -> TotSiteNo  &&  sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j_i ] ] == sup.env_space -> Num_hole_block[ sup.H_env[ j_i ] ] ) {

                	int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ];

                        int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ];

                        //----
                        double *f_1 = new double [ Dim_env_i * Dim_env_i ];
                        for( int l = 0; l < Dim_env_i * Dim_env_i; l++ )  f_1[ l ] = 0.0;

                        double *f_2 = new double [ Dim_env_i * Dim_env_i ];
                        for( int l = 0; l < Dim_env_i * Dim_env_i; l++ )  f_2[ l ] = 0.0;

                        //----
                        alpha = 1.0;  beta = 0.0;

                        dgemm_( &trans_T, &trans_N, &Dim_env_i, &Dim_env_i, &Dim_sys_i, &alpha, WaveFunction_block[ j_i ], &Dim_sys_i, WaveFunction_block[ j_i ], &Dim_sys_i, &beta, f_1, &Dim_env_i );

                        alpha = 1.0 / sqrt( sup.env_space -> Value_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ] + 1.0 );

                        dgemm_( &trans_N, &trans_T, &Dim_env_i, &Dim_env_i, &Dim_env_i, &alpha, Ni_old_B[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ], &Dim_env_i, f_1, &Dim_env_i, &beta, f_2, &Dim_env_i );

                        //----
                        for( int l = 0; l < Dim_env_i; l++ ) {

	                        electron_density[ para.Total_N - 1 - i ][ para.Total_N - 1 - lsys ] += f_2[ l * ( 1 + Dim_env_i ) ];

                                electron_density[ para.Total_N - 1 - lsys ][ para.Total_N - 1 - i ] += f_2[ l * ( 1 + Dim_env_i ) ];

                        }

                        delete [] f_1;  delete [] f_2;

		}

                fout<<"\n"<<para.Total_N - 1 - lsys<<"\t"<<para.Total_N - 1 - i<<"\t"<<electron_density[ para.Total_N - 1 - i ][ para.Total_N -1 - lsys ];

                //----Calculate the expectation value of correlation function
                for( int j_i = 0; j_i < sup.BlockNumber_for_TargetBlock; j_i++ )  
                if( sup.env_space -> Num_hole_block[ sup.H_env[ j_i ] ] != sup.env_space -> TotSiteNo  &&  sup.envnew_space -> Num_hole_block[ sup.H_envnew[ j_i ] ] == sup.env_space -> Num_hole_block[ sup.H_env[ j_i ] ] ) {

                	int Dim_sys_i = sup.sys_space -> Dim_J_block[ sup.H_sys[ j_i ] ][ sup.J_sys[ j_i ] ];

                        int Dim_env_i = sup.env_space -> Dim_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ];

                        //----
                        double *f_1 = new double [ Dim_env_i * Dim_env_i ];
                        for( int l = 0; l < Dim_env_i * Dim_env_i; l++ )  f_1[ l ] = 0.0;

                        double *f_2 = new double [ Dim_env_i * Dim_env_i ];
                        for( int l = 0; l < Dim_env_i * Dim_env_i; l++ )  f_2[ l ] = 0.0;

                        //----
                        alpha = 1.0;  beta = 0.0;

                        dgemm_( &trans_T, &trans_N, &Dim_env_i, &Dim_env_i, &Dim_sys_i, &alpha, WaveFunction_config_2[ j_i ], &Dim_sys_i, WaveFunction_config_2[ j_i ], &Dim_sys_i, &beta, f_1, &Dim_env_i );

                        alpha = 1.0 / sqrt( sup.env_space -> Value_J_block[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ] + 1.0 );

                        dgemm_( &trans_N, &trans_T, &Dim_env_i, &Dim_env_i, &Dim_env_i, &alpha, Ni_old_B[ sup.H_env[ j_i ] ][ sup.J_env[ j_i ] ], &Dim_env_i, f_1, &Dim_env_i, &beta, f_2, &Dim_env_i );

                        //----
                        for( int l = 0; l < Dim_env_i; l++ ) {

	                        electron_density[ para.Total_N - 1 - i ][ para.Total_N - 2 - lsys ] += f_2[ l * ( 1 + Dim_env_i ) ];

                                electron_density[ para.Total_N - 2 - lsys ][ para.Total_N - 1 - i ] += f_2[ l * ( 1 + Dim_env_i ) ];

                        }

                        delete [] f_1;  delete [] f_2;

		}

                fout<<"\n"<<para.Total_N - 2 - lsys<<"\t"<<para.Total_N - 1 - i<<"\t"<<electron_density[ para.Total_N - 1 - i ][ para.Total_N - 2 - lsys ];

		//----delete Ni_old
		for( int n = 0; n < B_Block_Num_hole[ lsys - 1 ]; n++ ) {

	        	for( int i = 0; i < B_Num_J_hole_block[ lsys - 1 ][ n ]; i++ )

	                	delete [] Ni_old_B[ n ][ i ];
	
			delete [] Ni_old_B[ n ];

		}

		delete [] Ni_old_B;

	}

  //----delete tmp_wavefunction
        for( int i = 0; i < sup.env_space -> Block_Num_hole; i++ ) {

		for( int j = 0; j < sup.env_space -> Num_J_hole_block[i]; j++ ) {

			delete [] tmp_wavefunction[i][j];
		}

		delete [] tmp_wavefunction[i];
        }

	delete [] tmp_wavefunction;

}

//========================================================================
//		Correlation function between ns and ne
//========================================================================
inline void BondDensity::Density_ns_ne(const int &lsys, Parameter &para, Super &sup) {

	alpha = 1.0;
	beta = 0.0;

	//----
	for( int i = 0; i < sup.BlockNumber_for_TargetBlock_config_3; i++ ) 
	if( sup.H_ns_config_3[ i ] == 0  &&  sup.H_ne_config_3[ i ] == 0 ) {

		//----
		int Dim_sys = sup.sys_space -> Dim_J_block[ sup.H_sys_config_3[ i ] ][ sup.J_sys_config_3[ i ] ];

                int Dim_env = sup.env_space -> Dim_J_block[ sup.H_env_config_3[ i ] ][ sup.J_env_config_3[ i ] ];

		//----
                int dim = Dim_sys * Dim_sys;

                double *f  = new double [ dim ];

                for( int j = 0; j < dim; j++ )  f[ j ] = 0.0;

		//----
                dgemm_( &trans_N, &trans_T, &Dim_sys, &Dim_sys, &Dim_env, &alpha, WaveFunction_config_3[ i ], &Dim_sys, WaveFunction_config_3[ i ], &Dim_sys, &beta, f, &Dim_sys );

		//----
                for( int j = 0; j < Dim_sys; j++) {

                	electron_density[ lsys ][ lsys + 1 ] += f[ j * ( 1 + Dim_sys ) ];
                	electron_density[ lsys + 1 ][ lsys ] += f[ j * ( 1 + Dim_sys ) ];

		}

                delete [] f;

        }

        fout<<"\n"<<lsys<<"\t"<<lsys+1<<"\t"<<electron_density[ lsys ][ lsys + 1 ];

}

//==================================================================
//				     A_Initial_Ni
//===================================================================
inline void BondDensity::A_Initial_Ni( const int &num ) {

	if( num == 0 ) {

		//----
		Ni_old_A = new double ** [ A_N_Block_Num_hole[ num ] ];// A_N_Block_Num_hole[0] == 2

		for( int i = 0; i < A_N_Block_Num_hole[ num ]; i++ ) {
 
			Ni_old_A[i] = new double * [ A_N_Num_J_hole_block[ num ][ i ] ];

			for( int j = 0; j < A_N_Num_J_hole_block[ num ][ i ]; j++ ) {
		
				int Dim = A_N_Dim_J_block[ num ][ i ][ j ] * A_N_Dim_J_block[ num ][ i ][ j ];

				Ni_old_A[i][j] = new double [ Dim ];

				for( int k = 0; k < Dim; k++ ) {

					Ni_old_A[i][j][k] = (double) 0;

				}

			}

		}

		//----
		Ni_old_A[0][0][0] = -hoping;
		Ni_old_A[1][0][0] = 0.0;

	}

	else if( num > 0 ) {

		//----
		Ni_old_A = new double ** [ A_N_Block_Num_hole[ num ] ];

		for( int i = 0; i < A_N_Block_Num_hole[ num ]; i++ ) {
 
			Ni_old_A[i] = new double * [ A_N_Num_J_hole_block[ num ][ i ] ];

			for( int j = 0; j < A_N_Num_J_hole_block[ num ][ i ]; j++ ) {
		
				int Dim = A_N_Dim_J_block[ num ][ i ][ j ] * A_N_Dim_J_block[ num ][ i ][ j ];

				Ni_old_A[i][j] = new double [ Dim ];

				for( int k = 0; k < Dim; k++ ) {

					Ni_old_A[i][j][k] = (double) 0;

				}

			}

		}

		//----
		int oldhm, oldJm, a_new, position;
		double factor;

		for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ )
		if( A_N_Num_hole_block[ num ][ n ] != ( num + 1 ) ) {

			for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ )
			for( int m = 1; m < 3; m++ )
			if( ( oldhm = A_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = A_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1 ) {

				position = A_N_Dim_J_block[ num ][ n ][ i ] + 1;

				factor = sqrt( A_N_Value_J_block[ num ][ n ][ i ] + 1.0 );
	
				for( int a_old = 0; a_old < A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ]; a_old++ ) {

					a_new = position * ( A_N_Start[ num ][ n ][ m ][ i ] + a_old );

					Ni_old_A[ n ][ i ][ a_new ] += factor;

				}

			}

		}

	}

}

//======================================================================
//			     B_Initial_Ni
//=======================================================================
inline void BondDensity::B_Initial_Ni( const int &num ) {

	if( num == 0 ) {

		//----
		Ni_old_B = new double ** [ B_N_Block_Num_hole[ num ] ];// B_N_Block_Num_hole[0] == 2

		for( int i = 0; i < B_N_Block_Num_hole[ num ]; i++ ) {
 
			Ni_old_B[i] = new double * [ B_N_Num_J_hole_block[ num ][ i ] ];

			for( int j = 0; j < B_N_Num_J_hole_block[ num ][ i ]; j++ ) {
		
				int Dim = B_N_Dim_J_block[ num ][ i ][ j ] * B_N_Dim_J_block[ num ][ i ][ j ];

				Ni_old_B[i][j] = new double [ Dim ];

				for( int k = 0; k < Dim; k++ ) {

					Ni_old_B[i][j][k] = (double) 0;

				}

			}

		}

		//----
		Ni_old_B[0][0][0] = -hoping;
		Ni_old_B[1][0][0] = 0.0;

	}

	else if( num > 0 ) {

		//----
		Ni_old_B = new double ** [ B_N_Block_Num_hole[ num ] ];

		for( int i = 0; i < B_N_Block_Num_hole[ num ]; i++ ) {
 
			Ni_old_B[i] = new double * [ B_N_Num_J_hole_block[ num ][ i ] ];

			for( int j = 0; j < B_N_Num_J_hole_block[ num ][ i ]; j++ ) {
		
				int Dim = B_N_Dim_J_block[ num ][ i ][ j ] * B_N_Dim_J_block[ num ][ i ][ j ];

				Ni_old_B[i][j] = new double [ Dim ];

				for( int k = 0; k < Dim; k++ ) {

					Ni_old_B[i][j][k] = (double) 0;

				}

			}

		}

		//----
		int oldhm, oldJm, a_new, position;
		double factor;

		for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ )
		if( B_N_Num_hole_block[ num ][ n ] != ( num + 1 ) ) {

			for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ )
			for( int m = 1; m < 3; m++ )
			if( ( oldhm = B_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = B_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1 ) {

				position = B_N_Dim_J_block[ num ][ n ][ i ] + 1;

				factor = sqrt( B_N_Value_J_block[ num ][ n ][ i ] + 1.0 );
	
				for( int a_old = 0; a_old < B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ]; a_old++ ) {

					a_new = position * ( B_N_Start[ num ][ n ][ m ][ i ] + a_old );

					Ni_old_B[ n ][ i ][ a_new ] += factor;

				}

			}

		}

	}

}

//===========================================================================
//				     A_Truncate
//===========================================================================
inline void BondDensity::A_Truncate_Ni( const int &num ) {

	int Dim, Dim_old, Dim_new, SubDim;

        double alpha = 1.0;
	double beta = 0.0;

//------Create space
        Ni_new_A = new double ** [ A_Block_Num_hole[ num ] ];

	for( int i = 0; i < A_Block_Num_hole[ num ]; i++ ) {

		Ni_new_A[i] = new double * [ A_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < A_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = A_Dim_J_block[ num ][ i ][ j ] * A_Dim_J_block[ num ][ i ][ j ];

			Ni_new_A[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				Ni_new_A[i][j][k] = 0.0;

			}

		}

	}

//------Truncate Ni_old operator
        for( int n = 0; n < A_Block_Num_hole[ num ]; n++ )
        for( int i = 0; i < A_Num_J_hole_block[ num ][ n ]; i++ ) {

		Dim_old = A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ];
		Dim_new = A_Dim_J_block[ num ][ n ][ i ];	
		SubDim = Dim_old * Dim_new;

                double * f_sub = new double [ SubDim ];

                for( int j = 0; j < SubDim; j++ )     f_sub[j] = (double) 0;

                dgemm_( &trans_N, &trans_N, &Dim_old, &Dim_new, &Dim_old, &alpha, Ni_old_A[ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], &Dim_old, A[ num ][ n ][ i ], &Dim_old, &beta, f_sub, &Dim_old );

                dgemm_( &trans_T, &trans_N, &Dim_new, &Dim_new, &Dim_old, &alpha, A[ num ][ n ][ i ], &Dim_old, f_sub, &Dim_old, &beta, Ni_new_A[ n ][ i ], &Dim_new );

                delete [] f_sub;

        }

//------Delete Ni_old
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ )

	                delete [] Ni_old_A[ n ][ i ];
	
		delete [] Ni_old_A[ n ];

	}

	delete [] Ni_old_A;
	
//------Create Ni_old
	Ni_old_A = new double ** [ A_Block_Num_hole[ num ] ];

	for( int n = 0; n < A_Block_Num_hole[ num ]; n++ ) {

		Ni_old_A[ n ] = new double * [ A_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < A_Num_J_hole_block[ num ][ n ]; i++ ) {

			Dim = A_Dim_J_block[ num ][ n ][ i ] * A_Dim_J_block[ num ][ n ][ i ];

			Ni_old_A[ n ][ i ] =  new double [ Dim ];
	
			for( int j = 0; j < Dim; j++ ) {
	
				Ni_old_A[ n ][ i ][ j ] = Ni_new_A[ n ][ i ][ j ];

			}

		} 

	}

//------Delete Ni_new
	for( int n = 0; n < A_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < A_Num_J_hole_block[ num ][ n ]; i++ )

	                delete [] Ni_new_A[ n ][ i ];

		delete [] Ni_new_A[ n ];

	}

	delete [] Ni_new_A;

}

//===============================================================
//			     B_Truncate
//===============================================================
inline void BondDensity::B_Truncate_Ni( const int &num ) {

	int Dim, Dim_old, Dim_new, SubDim;

        double alpha = 1.0;
	double beta = 0.0;

//------Create space
        Ni_new_B = new double ** [ B_Block_Num_hole[ num ] ];

	for( int i = 0; i < B_Block_Num_hole[ num ]; i++ ) {

		Ni_new_B[i] = new double * [ B_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < B_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = B_Dim_J_block[ num ][ i ][ j ] * B_Dim_J_block[ num ][ i ][ j ];

			Ni_new_B[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				Ni_new_B[i][j][k] = 0.0;

			}

		}

	}

//------Truncate Ni_old operator
        for( int n = 0; n < B_Block_Num_hole[ num ]; n++ )
        for( int i = 0; i < B_Num_J_hole_block[ num ][ n ]; i++ ) {

		Dim_old = B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ];
		Dim_new = B_Dim_J_block[ num ][ n ][ i ];	
		SubDim = Dim_old * Dim_new;

                double * f_sub = new double [ SubDim ];

                for( int j = 0; j < SubDim; j++ )     f_sub[j] = (double) 0;

                dgemm_( &trans_N, &trans_N, &Dim_old, &Dim_new, &Dim_old, &alpha, Ni_old_B[ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], &Dim_old, B[ num ][ n ][ i ], &Dim_old, &beta, f_sub, &Dim_old );

                dgemm_( &trans_T, &trans_N, &Dim_new, &Dim_new, &Dim_old, &alpha, B[ num ][ n ][ i ], &Dim_old, f_sub, &Dim_old, &beta, Ni_new_B[ n ][ i ], &Dim_new );

                delete [] f_sub;

        }

//------Delete Ni_old
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ )

	                delete [] Ni_old_B[ n ][ i ];
	
		delete [] Ni_old_B[ n ];

	}

	delete [] Ni_old_B;
	
//------Create Ni_old
	Ni_old_B = new double ** [ B_Block_Num_hole[ num ] ];

	for( int n = 0; n < B_Block_Num_hole[ num ]; n++ ) {

		Ni_old_B[ n ] = new double * [ B_Num_J_hole_block[ num ][ n ] ];

		for( int i = 0; i < B_Num_J_hole_block[ num ][ n ]; i++ ) {

			Dim = B_Dim_J_block[ num ][ n ][ i ] * B_Dim_J_block[ num ][ n ][ i ];

			Ni_old_B[ n ][ i ] =  new double [ Dim ];
	
			for( int j = 0; j < Dim; j++ ) {
	
				Ni_old_B[ n ][ i ][ j ] = Ni_new_B[ n ][ i ][ j ];

			}

		} 

	}

//------Delete Ni_new
	for( int n = 0; n < B_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < B_Num_J_hole_block[ num ][ n ]; i++ )

	                delete [] Ni_new_B[ n ][ i ];

		delete [] Ni_new_B[ n ];

	}

	delete [] Ni_new_B;

}

//=============================================================
//		     Truncate_A_SiSj
//=============================================================
inline void BondDensity::Truncate_A_SiSj( const int &num ) {

	int Dim, Dim_old, Dim_new, SubDim;

        double alpha = 1.0;
	double beta = 0.0;

//------Create space
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

//------Truncate 
        for( int n = 0; n < A_Block_Num_hole[ num ]; n++ )
        for( int i = 0; i < A_Num_J_hole_block[ num ][ n ]; i++ ) {

		Dim_old = A_N_Dim_J_block[ num ][ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ];
		Dim_new = A_Dim_J_block[ num ][ n ][ i ];	
		SubDim = Dim_old * Dim_new;

                double * f_sub = new double [ SubDim ];

                for( int j = 0; j < SubDim; j++ )     f_sub[j] = (double) 0;

                dgemm_( &trans_N, &trans_N, &Dim_old, &Dim_new, &Dim_old, &alpha, SiSj_old_A[ A_Old_hole[ num ][ n ] ][ A_Old_J[ num ][ n ][ i ] ], &Dim_old, A[ num ][ n ][ i ], &Dim_old, &beta, f_sub, &Dim_old );

                dgemm_( &trans_T, &trans_N, &Dim_new, &Dim_new, &Dim_old, &alpha, A[ num ][ n ][ i ], &Dim_old, f_sub, &Dim_old, &beta, SiSj_new_A[ n ][ i ], &Dim_new );

                delete [] f_sub;

        }

//------Delete Ni_old
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                delete [] SiSj_old_A[ n ][ i ];

		}
	
		delete [] SiSj_old_A[ n ];

	}

	delete [] SiSj_old_A;
	
//------Create SiSj_old
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

//=============================================================
//		     Truncate_B_SiSj
//=============================================================
inline void BondDensity::Truncate_B_SiSj( const int &num ) {

	int Dim, Dim_old, Dim_new, SubDim;

        double alpha = 1.0;
	double beta = 0.0;

//------Create space
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

//------Truncate 
        for( int n = 0; n < B_Block_Num_hole[ num ]; n++ )
        for( int i = 0; i < B_Num_J_hole_block[ num ][ n ]; i++ ) {

		Dim_old = B_N_Dim_J_block[ num ][ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ];
		Dim_new = B_Dim_J_block[ num ][ n ][ i ];	
		SubDim = Dim_old * Dim_new;

                double * f_sub = new double [ SubDim ];

                for( int j = 0; j < SubDim; j++ )     f_sub[j] = (double) 0;

                dgemm_( &trans_N, &trans_N, &Dim_old, &Dim_new, &Dim_old, &alpha, SiSj_old_B[ B_Old_hole[ num ][ n ] ][ B_Old_J[ num ][ n ][ i ] ], &Dim_old, B[ num ][ n ][ i ], &Dim_old, &beta, f_sub, &Dim_old );

                dgemm_( &trans_T, &trans_N, &Dim_new, &Dim_new, &Dim_old, &alpha, B[ num ][ n ][ i ], &Dim_old, f_sub, &Dim_old, &beta, SiSj_new_B[ n ][ i ], &Dim_new );

                delete [] f_sub;

        }

//------Delete Ni_old
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                delete [] SiSj_old_B[ n ][ i ];

		}
	
		delete [] SiSj_old_B[ n ];

	}

	delete [] SiSj_old_B;
	
//------Create SiSj_old
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

//==================================================================
//			     A_New
//==================================================================
inline void BondDensity::A_New_Ni( const int &num ) {

	int Dim, oldhm, oldJm, position, SubDim, a_new;
	double factor;

//------Create space
        Ni_new_A = new double ** [ A_N_Block_Num_hole[ num ] ];

	for( int i = 0; i < A_N_Block_Num_hole[ num ]; i++ ) {

		Ni_new_A[i] = new double * [ A_N_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < A_N_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = A_N_Dim_J_block[ num ][ i ][ j ] * A_N_Dim_J_block[ num ][ i ][ j ];

			Ni_new_A[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				Ni_new_A[i][j][k] = (double) 0;

			}

		}

	}

	//----
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ )
	if( A_N_Num_hole_block[ num ][ n ] != ( num + 1 ) ) {

		for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ )
                for( int m = 0; m < 3; m++ )
                if( ( oldhm = A_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = A_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  A_Num_hole_block[ num - 1 ][ oldhm ] != num ) {

                	position = ( A_N_Dim_J_block[ num ][ n ][ i ] + 1 ) * A_N_Start[ num ][ n ][ m ][ i ];

                        SubDim = A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                        factor = sqrt( ( A_N_Value_J_block[ num ][ n ][ i ] + 1.0 ) / ( A_Value_J_block[ num -1 ][ oldhm ][ oldJm ] + 1.0 ) );

                        for( int a_old = 0; a_old < SubDim; a_old++ ) {

                        	a_new = position + A_N_Dim_J_block[ num ][ n ][ i ] * ( a_old / A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_old % A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                                Ni_new_A[ n ][ i ][ a_new ] += factor * Ni_old_A[ oldhm ][ oldJm ][ a_old ];

                        }

		}

	}

//------Delete Ni_old
	for( int n = 0; n < A_Block_Num_hole[ num - 1 ]; n++ ) {

	        for( int i = 0; i < A_Num_J_hole_block[ num - 1 ][ n ]; i++ )

	                delete [] Ni_old_A[ n ][ i ];
	
		delete [] Ni_old_A[ n ];

	}

	delete [] Ni_old_A;

//------Create Ni_old
        Ni_old_A = new double ** [ A_N_Block_Num_hole[ num ] ];

	for( int i = 0; i < A_N_Block_Num_hole[ num ]; i++ ) {

		Ni_old_A[i] = new double * [ A_N_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < A_N_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = A_N_Dim_J_block[ num ][ i ][ j ] * A_N_Dim_J_block[ num ][ i ][ j ];

			Ni_old_A[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				Ni_old_A[i][j][k] = Ni_new_A[i][j][k];

			}

		}

	}

//------Delete Ni_new
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ )

	                delete [] Ni_new_A[ n ][ i ];
	
		delete [] Ni_new_A[ n ];

	}

	delete [] Ni_new_A;

}

//==================================================================
//			     New_A_SiSj
//==================================================================
inline void BondDensity::New_A_SiSj( const int &num ) {

	int Dim, oldhm, oldJm, position, a_new;
	double factor;

//------Create space
        SiSj_new_A = new double ** [ A_N_Block_Num_hole[ num ] ];

	for( int i = 0; i < A_N_Block_Num_hole[ num ]; i++ ) {

		SiSj_new_A[i] = new double * [ A_N_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < A_N_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = A_N_Dim_J_block[ num ][ i ][ j ] * A_N_Dim_J_block[ num ][ i ][ j ];

			SiSj_new_A[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				SiSj_new_A[i][j][k] = (double) 0;

			}

		}

	}

	//----
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ )
	if( A_N_Num_hole_block[ num ][ n ] != ( num + 1 ) ) {

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

//------Delete Ni_old
	for( int n = 0; n < A_Block_Num_hole[ num - 1 ]; n++ ) {

	        for( int i = 0; i < A_Num_J_hole_block[ num - 1 ][ n ]; i++ ) {

	                delete [] SiSj_old_A[ n ][ i ];

		}
	
		delete [] SiSj_old_A[ n ];

	}

	delete [] SiSj_old_A;

//------Create SiSj_old
        SiSj_old_A = new double ** [ A_N_Block_Num_hole[ num ] ];

	for( int i = 0; i < A_N_Block_Num_hole[ num ]; i++ ) {

		SiSj_old_A[i] = new double * [ A_N_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < A_N_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = A_N_Dim_J_block[ num ][ i ][ j ] * A_N_Dim_J_block[ num ][ i ][ j ];

			SiSj_old_A[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				SiSj_old_A[i][j][k] = SiSj_new_A[i][j][k];

			}

		}

	}

//------Delete Ni_new
	for( int n = 0; n < A_N_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                delete [] SiSj_new_A[ n ][ i ];

		}
	
		delete [] SiSj_new_A[ n ];

	}

	delete [] SiSj_new_A;

}

//==================================================================
//			     New_B_SiSj
//==================================================================
inline void BondDensity::New_B_SiSj( const int &num ) {

	int Dim, oldhm, oldJm, position, a_new;
	double factor;

//------Create space
        SiSj_new_B = new double ** [ B_N_Block_Num_hole[ num ] ];

	for( int i = 0; i < B_N_Block_Num_hole[ num ]; i++ ) {

		SiSj_new_B[i] = new double * [ B_N_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < B_N_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = B_N_Dim_J_block[ num ][ i ][ j ] * B_N_Dim_J_block[ num ][ i ][ j ];

			SiSj_new_B[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				SiSj_new_B[i][j][k] = (double) 0;

			}

		}

	}

	//----
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ )
	if( B_N_Num_hole_block[ num ][ n ] != ( num + 1 ) ) {

		for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ )
                for( int m = 0; m < 3; m++ )
                if( ( oldhm = B_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = B_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  B_Num_hole_block[ num - 1 ][ oldhm ] != num ) {

                	position = ( B_N_Dim_J_block[ num ][ n ][ i ] + 1 ) * B_N_Start[ num ][ n ][ m ][ i ];

                        Dim = B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

			factor = sqrt( ( B_N_Value_J_block[ num ][ n ][ i ] + 1.0 ) / ( B_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] + 1.0 ) );

                        for( int a_old = 0; a_old < Dim; a_old++ ) {

                        	a_new = position + B_N_Dim_J_block[ num ][ n ][ i ] * ( a_old / B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_old % B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                                SiSj_new_B[ n ][ i ][ a_new ] += factor * SiSj_old_B[ oldhm ][ oldJm ][ a_old ];

                        }

		}

	}

//------Delete Ni_old
	for( int n = 0; n < B_Block_Num_hole[ num - 1 ]; n++ ) {

	        for( int i = 0; i < B_Num_J_hole_block[ num - 1 ][ n ]; i++ ) {

	                delete [] SiSj_old_B[ n ][ i ];

		}
	
		delete [] SiSj_old_B[ n ];

	}

	delete [] SiSj_old_B;

//------Create SiSj_old
        SiSj_old_B = new double ** [ B_N_Block_Num_hole[ num ] ];

	for( int i = 0; i < B_N_Block_Num_hole[ num ]; i++ ) {

		SiSj_old_B[i] = new double * [ B_N_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < B_N_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = B_N_Dim_J_block[ num ][ i ][ j ] * B_N_Dim_J_block[ num ][ i ][ j ];

			SiSj_old_B[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				SiSj_old_B[i][j][k] = SiSj_new_B[i][j][k];

			}

		}

	}

//------Delete Ni_new
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ ) {

	                delete [] SiSj_new_B[ n ][ i ];

		}
	
		delete [] SiSj_new_B[ n ];

	}

	delete [] SiSj_new_B;

}

//====================================================================
//			     B_New
//====================================================================
inline void BondDensity::B_New_Ni( const int &num ) {

	int Dim, oldhm, oldJm, position, SubDim, a_new;
	double factor;

//------Create space
        Ni_new_B = new double ** [ B_N_Block_Num_hole[ num ] ];

	for( int i = 0; i < B_N_Block_Num_hole[ num ]; i++ ) {

		Ni_new_B[i] = new double * [ B_N_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < B_N_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = B_N_Dim_J_block[ num ][ i ][ j ] * B_N_Dim_J_block[ num ][ i ][ j ];

			Ni_new_B[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				Ni_new_B[i][j][k] = (double) 0;

			}

		}

	}

	//----
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ )
	if( B_N_Num_hole_block[ num ][ n ] != ( num + 1 ) ) {

		for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ )
                for( int m = 0; m < 3; m++ )
                if( ( oldhm = B_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = B_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  B_Num_hole_block[ num - 1 ][ oldhm ] != num ) {

                	position = ( B_N_Dim_J_block[ num ][ n ][ i ] + 1 ) * B_N_Start[ num ][ n ][ m ][ i ];

                        SubDim = B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                        factor = sqrt( ( B_N_Value_J_block[ num ][ n ][ i ] + 1.0 ) / ( B_Value_J_block[ num -1 ][ oldhm ][ oldJm ] + 1.0 ) );

                        for( int a_old = 0; a_old < SubDim; a_old++ ) {

                        	a_new = position + B_N_Dim_J_block[ num ][ n ][ i ] * ( a_old / B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_old % B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                                Ni_new_B[ n ][ i ][ a_new ] += factor * Ni_old_B[ oldhm ][ oldJm ][ a_old ];

                        }

		}

	}

//------Delete Ni_old
	for( int n = 0; n < B_Block_Num_hole[ num - 1 ]; n++ ) {

	        for( int i = 0; i < B_Num_J_hole_block[ num - 1 ][ n ]; i++ )

	                delete [] Ni_old_B[ n ][ i ];
	
		delete [] Ni_old_B[ n ];

	}

	delete [] Ni_old_B;

//------Create Ni_old
        Ni_old_B = new double ** [ B_N_Block_Num_hole[ num ] ];

	for( int i = 0; i < B_N_Block_Num_hole[ num ]; i++ ) {

		Ni_old_B[i] = new double * [ B_N_Num_J_hole_block[ num ][ i ] ];

		for( int j = 0; j < B_N_Num_J_hole_block[ num ][ i ]; j++ ) {

			Dim = B_N_Dim_J_block[ num ][ i ][ j ] * B_N_Dim_J_block[ num ][ i ][ j ];

			Ni_old_B[i][j] = new double [ Dim ];

			for( int k = 0; k < Dim; k++ ) {

				Ni_old_B[i][j][k] = Ni_new_B[i][j][k];

			}

		}

	}

//------Delete Ni_new
	for( int n = 0; n < B_N_Block_Num_hole[ num ]; n++ ) {

	        for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ )

	                delete [] Ni_new_B[ n ][ i ];
	
		delete [] Ni_new_B[ n ];

	}

	delete [] Ni_new_B;

}

//==================================================================
//			Initial_A_NiNj
//===================================================================
inline void BondDensity::Initial_A_SiSj( const int &num ) {

        int oldhm, oldhl, oldJl, oldJm, Dim, a_new, a_new_t, position;
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
        if( A_N_Num_hole_block[ num ][ n ] != ( num + 1 ) ) {

	        for( int i = 0; i < A_N_Num_J_hole_block[ num ][ n ]; i++ )
                for( int m = 1; m < 3; m++ )
                if( ( oldhm = A_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = A_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1 ) 
		if( A_N_Num_hole_block[ num ][ n ] == A_Num_hole_block[ num - 1 ][ oldhm ]  &&  A_Num_hole_block[ num - 1 ][ oldhm ] != num ) {

	                //----
                        position = A_N_Dim_J_block[ num ][ n ][ i ] * A_N_Start[ num ][ n ][ m ][ i ] + A_N_Start[ num ][ n ][ m ][ i ];

                        factor = sqrt( A_N_Value_J_block[ num ][ n ][ i ] + 1.0 ) / sqrt( A_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] + 1.0 );

                        Dim = A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

			//----
                        for( int a_sys = 0; a_sys < Dim; a_sys++ ) {

	                        a_new = position + A_N_Dim_J_block[ num ][ n ][ i ] * ( a_sys / A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_sys % A_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                                SiSj_old_A[ n ][ i ][ a_new ] += factor * Ni_old_A[ oldhm ][ oldJm ][ a_sys ];

                         }

		}

	}

}

//==================================================================
//			Initial_B_NiNj
//===================================================================
inline void BondDensity::Initial_B_SiSj( const int &num ) {

        int oldhm, oldhl, oldJl, oldJm, Dim, a_new, a_new_t, position;
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
        if( B_N_Num_hole_block[ num ][ n ] != ( num + 1 ) ) {

	        for( int i = 0; i < B_N_Num_J_hole_block[ num ][ n ]; i++ )
                for( int m = 1; m < 3; m++ )
                if( ( oldhm = B_N_Hole_blockOld[ num ][ n ][ m ][ i ] ) != -1  &&  ( oldJm = B_N_J_blockOld[ num ][ n ][ m ][ i ] ) != -1 ) {

	                //----
                        position = B_N_Dim_J_block[ num ][ n ][ i ] * B_N_Start[ num ][ n ][ m ][ i ] + B_N_Start[ num ][ n ][ m ][ i ];

                        factor = sqrt( B_N_Value_J_block[ num ][ n ][ i ] + 1.0 ) / sqrt( B_Value_J_block[ num - 1 ][ oldhm ][ oldJm ] + 1.0 );

                        Dim = B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] * B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

			//----
                        for( int a_sys = 0; a_sys < Dim; a_sys++ ) {

	                        a_new = position + B_N_Dim_J_block[ num ][ n ][ i ] * ( a_sys / B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ] ) + a_sys % B_Dim_J_block[ num - 1 ][ oldhm ][ oldJm ];

                                SiSj_old_B[ n ][ i ][ a_new ] += factor * Ni_old_B[ oldhm ][ oldJm ][ a_sys ];

                         }

		}

	}

}

//=====================================================================
//			Delete WaveFunctions
//=====================================================================
inline void BondDensity::DeleteFunction( Super &sup ) {

        for( int i = 0; i < sup.BlockNumber_for_TargetBlock; i++ )

                delete [] WaveFunction_block[i];

        delete [] WaveFunction_block;

}

//======================================================================================
//						Delete Space
//======================================================================================
inline void BondDensity::DeleteSpace( const int &lsys, Parameter &para ) {

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
BondDensity::~BondDensity() {}
//==================================================END==========================================================
