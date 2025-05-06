#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<stdio.h>
using namespace std;

#include"subcommon.h"
#include"common.h"

//=====================Initialize the system and environment blocks==================
SubCommon::SubCommon() {}  //See the subroutine in Sub Class

//=====================CreateSpace_hole_block========================================
void SubCommon::CreateSpace_hole_block() {
	
  //----
	Num_hole_block = new int [ Block_Num_hole ];

	for( int i = 0; i < Block_Num_hole; i++ )

		Num_hole_block[i] = 0;

  //----
	Num_J_hole_block = new int [ Block_Num_hole ];

	for( int i = 0; i < Block_Num_hole; i++ )

		Num_J_hole_block[i] = 0;

}

//====================CreateSpace_J_block===============================================
void SubCommon::CreateSpace_J_block() {
	
	Value_J_block = new int * [ Block_Num_hole ];
	Dim_J_block = new int * [ Block_Num_hole ];
	
	for( int i = 0; i < Block_Num_hole; i++ ) {

		Value_J_block[i] = new int [ Num_J_hole_block[i] ];
		Dim_J_block[i] = new int [ Num_J_hole_block[i] ];

		for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

			Value_J_block[i][j] = 0;  Dim_J_block[i][j] = 0;

		}

	}

	if( IndexNo == 1 ) {

		Hole_blockOld = new int ** [ Block_Num_hole ];
		J_blockOld = new int ** [ Block_Num_hole ];
		Start = new int ** [ Block_Num_hole ];

		for( int i = 0; i < Block_Num_hole; i++ ) {

			Hole_blockOld[i] = new int * [3];
			J_blockOld[i] = new int * [3];	   
			Start[i] = new int * [3];

			for( int j = 0; j < 3; j++ ) {

				Hole_blockOld[i][j] = new int [ Num_J_hole_block[i] ];		
				J_blockOld[i][j] = new int [ Num_J_hole_block[i] ];
				Start[i][j] = new int [ Num_J_hole_block[i] ];

				for( int k = 0; k < Num_J_hole_block[i]; k++ ) {
		
					Hole_blockOld[i][j][k] = -1;
					J_blockOld[i][j][k] = -1;
					Start[i][j][k] = 0;

				}

			}

		}

	}

}

//=========================Create Space S_Dia=================================
void SubCommon::Create_S_Dia() {

	int Dim = 0;

	S_Dia = new double *** [ operator_number_J ];

        for( int n = 0; n < operator_number_J; n++ ) {

                S_Dia[n] = new double ** [ Block_Num_hole ];

                for( int i = 0; i < Block_Num_hole; i++ ) {

			S_Dia[n][i] = new double * [ Num_J_hole_block[i] ];

			for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

                        	Dim = Dim_J_block[i][j] * Dim_J_block[i][j];

                        	S_Dia[n][i][j] = new double [ Dim ];

                        	for( int k =0; k < Dim; k++ )

                                	S_Dia[n][i][j][k] = (double) 0;

			}

                }

        }

}

//=============================FreeSpace S_Dia===================================
inline void SubCommon::FreeSpace_S_Dia() {

	for( int n = 0; n < operator_number_J; n++ ) {

                for( int i = 0; i < Block_Num_hole; i++ ) {

			for( int j = 0; j < Num_J_hole_block[i]; j++ ) {
	
        	                delete [] S_Dia[n][i][j];
                
			}

        		delete [] S_Dia[n][i];

		}
	
		delete [] S_Dia[n];

        }

        delete [] S_Dia;

}

//===============================Create Space S_M_Dia===============================
void SubCommon::Create_S_M_Dia() {

	int Dim = 0;

        S_M_Dia = new double *** [ operator_number_J ];

        for( int n = 0; n < operator_number_J; n++ ) {

                S_M_Dia[n] = new double ** [ Block_Num_hole ];

                for( int i = 0; i < Block_Num_hole; i++ ) {

		 	S_M_Dia[n][i] = new double * [ Num_J_hole_block[i] - 1 ];

			for( int j = 0; j < Num_J_hole_block[i] - 1; j++ ) {
	
                        	Dim = Dim_J_block[i][j] * Dim_J_block[i][ j + 1 ];

                        	S_M_Dia[n][i][j] = new double [ Dim ];

                        	for( int k = 0; k < Dim; k++ )

                                	S_M_Dia[n][i][j][k] = (double) 0;

			}
	
                }

        }

}

//================================Free Space S_M_Dia===================================
inline void SubCommon::FreeSpace_S_M_Dia() {

	for( int n = 0; n < operator_number_J; n++ ) {

		for( int i = 0; i < Block_Num_hole; i++ ) {

			for( int j = 0; j < Num_J_hole_block[i] - 1; j++ ) {

				delete [] S_M_Dia[n][i][j];

			}

			delete [] S_M_Dia[n][i];

		}

		delete [] S_M_Dia[n];

	}

	delete [] S_M_Dia;

}

//=======================================Create space for N operators=======================================
void SubCommon::Create_NN() {	//store the irreducible tensor operator for NN operator

	int Dim = 0;

	NN = new double *** [ operator_number_J ];	//change to operator_number_J

        for( int n = 0; n < operator_number_J; n++ ) {

                NN[n] = new double ** [ Block_Num_hole ];

                for( int i = 0; i < Block_Num_hole; i++ ) {

			NN[n][i] = new double * [ Num_J_hole_block[i] ];

			for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

                        	Dim = Dim_J_block[i][j] * Dim_J_block[i][j];

                        	NN[n][i][j] = new double [ Dim ];

                        	for( int k = 0; k < Dim; k++ )

                                	NN[n][i][j][k] = (double) 0;

			}

                }

	}

}

//=====================================Delete space for N oeprator================================
inline void SubCommon::FreeSpace_NN() {

	for( int n = 0; n < operator_number_J; n++ ) {

                for( int i = 0; i < Block_Num_hole; i++ ) {

			for( int j = 0; j < Num_J_hole_block[i]; j++ ) {
	
        	                delete [] NN[n][i][j];
                
			}

        		delete [] NN[n][i];

		}
	
		delete [] NN[n];

        }

        delete [] NN;

}

//=========================Create Space for T operators==========================
void SubCommon::CreateSpace_T() {

	int counter;

  //----Num_block_T, number of blocks for each hole quantum number, related to J
	Num_block_T = new int [ Block_Num_hole - 1 ];

	for( int i = 0; i < Block_Num_hole - 1; i++ ) 

		Num_block_T[i] = 0;

  //----
	for( int i = 0; i < Block_Num_hole - 1; i++ ) 
	if( Num_hole_block[ i + 1 ] - Num_hole_block[ i ] == 1 ) {

		counter = 0;

		for( int j = 0; j < Num_J_hole_block[i]; j++ )  // J number of the bra vector
		for( int k = 0; k < Num_J_hole_block[ i + 1 ]; k++ )  // J number of the ket vector
		if( abs( Value_J_block[i][j] - Value_J_block[ i + 1 ][k] ) == 1 ) // must be different

			counter++;

		Num_block_T[i] = counter;  // counter denotes the number of the subblocks

	}

  //----Find J_block_T_bra and J_block_T_ket
	J_block_T_bra = new int * [ Block_Num_hole - 1 ];	
	J_block_T_ket = new int * [ Block_Num_hole - 1 ];

	for( int i = 0; i < Block_Num_hole - 1; i++ ) {

		J_block_T_bra[i] = new int [ Num_block_T[i] ];	
		J_block_T_ket[i] = new int [ Num_block_T[i] ];

		for( int j = 0; j < Num_block_T[i]; j++ ) {

			J_block_T_bra[i][j] = 0;  J_block_T_ket[i][j] = 0;

		}

	}

	for( int i = 0; i < Block_Num_hole - 1; i++ )
	if( Num_hole_block[ i + 1 ] - Num_hole_block[ i ] == 1 ) {
 	
		counter = 0;

		for( int j = 0; j < Num_J_hole_block[i]; j++ )
		for( int k = 0; k < Num_J_hole_block[ i + 1 ]; k++ )
		if( abs( Value_J_block[i][j] - Value_J_block[ i + 1 ][k] ) == 1 ) {

			J_block_T_bra[i][counter] = j;  
			J_block_T_ket[i][counter] = k;
			counter++;

		}

	}

}

//==========================Create Space for T operators==========================
void SubCommon::Create_T() {

	int Dim = 0;

	T = new double *** [ operator_number_T ];  // site number

	for( int n = 0; n < operator_number_T; n++ ) {

		T[n] = new double ** [ Block_Num_hole - 1 ];  // hole quantum number 

		for( int i = 0; i < Block_Num_hole - 1; i++ )
 		if( Num_hole_block[ i + 1 ] - Num_hole_block[ i ] == 1 ) {  // hole number difference is 1

			T[n][i] = new double * [ Num_block_T[i] ];  // with different J quantum numbers

			for( int j = 0; j < Num_block_T[i]; j++ ) {

				Dim = Dim_J_block[i][ J_block_T_bra[i][j] ] * Dim_J_block[ i + 1 ][ J_block_T_ket[i][j] ];

				T[n][i][j] = new double [ Dim ];

				for( int k = 0; k < Dim; k++ )

					T[n][i][j][k] = 0.0;

			}

		}

	}

}

//=======================================Free Space for T operators====================================================
inline void SubCommon::FreeSpace_T_operator() {

	for( int n = 0; n < operator_number_T; n++ ) {

		for( int i = 0; i < Block_Num_hole - 1; i++ )
 		if( Num_hole_block[ i + 1 ] - Num_hole_block[ i ] == 1 ) {

			for( int j = 0; j < Num_block_T[i]; j++ ) {

				delete [] T[n][i][j];

			}

			delete [] T[n][i];

		}

		delete [] T[n];

	}
 
	delete [] T;

}

//================================Create Space H====================================
void SubCommon::Create_H() {

	int Dim = 0;

	H = new double ** [ Block_Num_hole ];

        for( int i = 0; i < Block_Num_hole; i++ ) {

		H[i] = new double * [ Num_J_hole_block[i] ];

		for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

	                Dim = Dim_J_block[i][j] * Dim_J_block[i][j];

        	        H[i][j] = new double [ Dim ];

                	for( int k = 0; k < Dim; k++ )

                        	H[i][j][k] = (double) 0;

		}

        }

}

//===============================Free Space H===============================
inline void SubCommon::FreeSpace_H() {

	for( int i = 0; i < Block_Num_hole; i++ ) {

		for( int j = 0; j < Num_J_hole_block[i]; j++ ) {
	
			delete [] H[i][j];
		
		}

		delete [] H[i];
	}

	delete [] H;

}

//============================Create Space 3===================================
//Create Spaces for the eigenvalues and eigenstates of reduced density matrices
//=============================================================================
void SubCommon::CreateSpace3() {

        dm_eig = new double ** [ Block_Num_hole ];
        dm_wave = new double ** [ Block_Num_hole ];

        for( int i = 0; i < Block_Num_hole; i++ ) {

		dm_eig[i] = new double * [ Num_J_hole_block[i] ];
		dm_wave[i] = new double * [ Num_J_hole_block[i] ];

		for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

        	        dm_eig[i][j] = new double [ Dim_J_block[i][j] ];

			for( int k = 0; k < Dim_J_block[i][j]; k++ )
				dm_eig[i][j][k] = (double) 0;
	
	                int Dim = Dim_J_block[i][j] * Dim_J_block[i][j];

	                dm_wave[i][j] = new double [ Dim ];

                	for( int k = 0; k < Dim; k++ ) 
                	        dm_wave[i][j][k] = (double) 0;

		}

        }

}

//==============================Print Space====================================
void SubCommon::Print_space( const int &i, const int &lsys ) {

	FILE *fp = fopen( Combine ( Combine ( "mid/space/", i ), lsys ), "wb" );

        fwrite( &TotSiteNo, sizeof(int), 1, fp );
        fwrite( &Block_Num_hole, sizeof(int), 1, fp );
	fwrite( &operator_number_J, sizeof(int), 1, fp );
        fwrite( Num_hole_block, sizeof(int), Block_Num_hole, fp );
        fwrite( Num_J_hole_block, sizeof(int), Block_Num_hole, fp );

        for( int i = 0; i < Block_Num_hole; i++ )
                fwrite( Value_J_block[i], sizeof(int), Num_J_hole_block[i], fp );

        for( int i = 0; i < Block_Num_hole; i++ )
                fwrite( Dim_J_block[i], sizeof(int), Num_J_hole_block[i], fp );

	fwrite( &operator_number_T, sizeof(int), 1, fp );
        fwrite( Num_block_T, sizeof(int), Block_Num_hole - 1, fp );

        for( int i = 0; i < Block_Num_hole - 1; i++ )
                fwrite( J_block_T_bra[i], sizeof(int), Num_block_T[i], fp );
 
        for( int i = 0; i < Block_Num_hole - 1; i++ )
                fwrite( J_block_T_ket[i], sizeof(int), Num_block_T[i], fp );

        fclose(fp);

}

//==============================Print_operator_H=====================================
void SubCommon::Print_H( const int &i, const int &lsys ) {

	int Dim;

        FILE *fp = fopen( Combine ( Combine ( "mid/operator/H/", i ), lsys ), "wb" );

        for( int i = 0; i < Block_Num_hole; i++ ) 
	for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

                Dim = Dim_J_block[i][j] * Dim_J_block[i][j];

                fwrite( H[i][j], sizeof(double), Dim, fp );

        }

        fclose(fp);

}

//=========================Print_operator_S_Dia===================================
void SubCommon::Print_S_Dia( const int &i, const int &lsys ) {

	int Dim;

        FILE *fp = fopen( Combine ( Combine ( "mid/operator/S_Dia/", i ), lsys ), "wb" );

        for( int n = 0; n < operator_number_J; n++ )
        for( int i = 0; i < Block_Num_hole; i++ ) 
	for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

                Dim = Dim_J_block[i][j] * Dim_J_block[i][j];

                fwrite( S_Dia[n][i][j], sizeof(double), Dim, fp );

        }

        fclose(fp);

}

//===============================Print_operator_S_M_Dia=================================
void SubCommon::Print_S_M_Dia( const int &i, const int &lsys ) {

	int Dim;

        FILE *fp = fopen( Combine ( Combine ( "mid/operator/S_M_Dia/", i ), lsys ), "wb" );

        for( int n = 0; n < operator_number_J; n++ )
        for( int i = 0; i < Block_Num_hole; i++ ) 
	for( int j = 0; j < Num_J_hole_block[i] - 1; j++ ) {

                Dim = Dim_J_block[i][j] * Dim_J_block[i][ j + 1 ];

                fwrite( S_M_Dia[n][i][j], sizeof(double), Dim, fp );

        }

        fclose(fp);

}

//==============================Print_operator_N======================================
void SubCommon::Print_NN( const int &i, const int &lsys ) {

	int Dim;

        FILE *fp = fopen( Combine ( Combine ( "mid/operator/NN/", i ), lsys ), "wb" );

        for( int n = 0; n < operator_number_J; n++ )
        for( int i = 0; i < Block_Num_hole; i++ ) 
	for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

	        Dim = Dim_J_block[i][j] * Dim_J_block[i][j];

                fwrite( NN[n][i][j], sizeof(double), Dim, fp );

	}

        fclose(fp);

}

//===========================================Print T operators========================================
void SubCommon::Print_T( const int &i, const int &lsys ) {

	int Dim;

        FILE *fp = fopen( Combine ( Combine ( "mid/operator/T/", i ), lsys ), "wb" );

        for( int n = 0; n < operator_number_T; n++ )
        for( int i = 0; i < Block_Num_hole - 1; i++ ) 
	for( int j = 0; j < Num_block_T[i]; j++ ) {

		Dim = Dim_J_block[i][ J_block_T_bra[i][j] ] * Dim_J_block[ i + 1 ][ J_block_T_ket[i][j] ];

                fwrite( T[n][i][j], sizeof(double), Dim, fp );

        }

        fclose(fp);

}

//========================Read angular space from hard disk==========================
SubCommon::SubCommon( const int &i, const int &lsys ) {

        IndexNo = 3;      

        FILE *fr = fopen( Combine ( Combine ( "mid/space/", i ), lsys ), "rb" );

        fread( &TotSiteNo, sizeof(int), 1, fr );
        fread( &Block_Num_hole, sizeof(int), 1, fr );
	fread( &operator_number_J, sizeof(int), 1, fr );

        CreateSpace_hole_block();
        fread( Num_hole_block, sizeof(int), Block_Num_hole, fr );
        fread( Num_J_hole_block, sizeof(int), Block_Num_hole, fr );

	CreateSpace_J_block();
        for( int i = 0; i < Block_Num_hole; i++ )
                fread( Value_J_block[i], sizeof(int), Num_J_hole_block[i], fr );

	for( int i = 0; i < Block_Num_hole; i++ )
                fread( Dim_J_block[i], sizeof(int), Num_J_hole_block[i], fr );

	fread( &operator_number_T, sizeof(int), 1, fr );

	CreateSpace_T();
        fread( Num_block_T, sizeof(int), Block_Num_hole - 1, fr );

        for( int i = 0; i < Block_Num_hole - 1; i++ )
                fread( J_block_T_bra[i], sizeof(int), Num_block_T[i], fr );

        for( int i = 0; i < Block_Num_hole - 1; i++ )
                fread( J_block_T_ket[i], sizeof(int), Num_block_T[i], fr );

        fclose(fr);

}

//=============================================Delete the class SubCommon=============================
SubCommon::~SubCommon() {

	if( IndexNo == 0 ) {

		FreeSpace_H();
		FreeSpace_T_operator();
		FreeSpace_T();
		FreeSpace_NN();
		FreeSpace_S_M_Dia();
		FreeSpace_S_Dia();

		FreeSpace_J_block();
		FreeSpace_hole_block();

	}

	else if( IndexNo == 1  ||  IndexNo == 3 ) {

		FreeSpace_T();
		FreeSpace_J_block();
		FreeSpace_hole_block();

	}

	else if( IndexNo == 2 ) {

		if(ope_sign=='H')	FreeSpace_H();

		else if(ope_sign=='T')	FreeSpace_T_operator();

		else if(ope_sign=='N')	FreeSpace_NN();

		else if(ope_sign=='S')	FreeSpace_S_Dia();

		else if(ope_sign=='M')	FreeSpace_S_M_Dia();

		FreeSpace_T();
		FreeSpace_J_block();
		FreeSpace_hole_block();

	}

	else if( IndexNo == 4 ) {

		for(int i=0; i<Block_Num_hole; i++) {

			for(int j=0; j<Num_J_hole_block[i]; j++) {
				delete [] dm_eig[i][j];	delete [] dm_wave[i][j];
			}	
			
			delete [] dm_eig[i];	delete [] dm_wave[i];
		}
		delete [] dm_eig;	delete [] dm_wave;

		FreeSpace_T();
		FreeSpace_J_block();
		FreeSpace_hole_block();

	}

	else if( IndexNo == 5 ) {
	
//		delete [] Old_hole;
//		for(int i=0; i<Block_Num_hole; i++)
//			delete [] Old_J[i];
//		delete [] Old_J;

		FreeSpace_T();
		FreeSpace_J_block();
		FreeSpace_hole_block();

	}

}

//===========================================FreeSpace_hole_block=====================================
inline void SubCommon::FreeSpace_hole_block() {

	delete [] Num_hole_block;	delete [] Num_J_hole_block;

}

//===========================================FreeSpace_J_hole_block===================================
inline void SubCommon::FreeSpace_J_block() {

	if( IndexNo == 1 ) {

		for( int i = 0; i < Block_Num_hole; i++ ) {

			for( int j = 0; j < 3; j++ ) {	

				delete [] Hole_blockOld[i][j];  delete [] J_blockOld[i][j];  delete [] Start[i][j];

			}

		delete [] Hole_blockOld[i];  delete [] J_blockOld[i];  delete [] Start[i];

		}

		delete [] Hole_blockOld;  delete [] J_blockOld;  delete [] Start;

	}

	for( int i = 0; i < Block_Num_hole; i++ ) {

		delete [] Value_J_block[i];	delete [] Dim_J_block[i];

	}

	delete [] Value_J_block;		delete [] Dim_J_block;

}

//=======================================Free Space for T operators==================================
inline void SubCommon::FreeSpace_T() {

  //----Delete J_block_T_bra and J_block_T_ket
	for( int i = 0; i < Block_Num_hole - 1; i++ ) {

		delete [] J_block_T_bra[i];	delete [] J_block_T_ket[i];

	}

	delete [] J_block_T_bra;		delete [] J_block_T_ket;

  //----Delete Num_block_T
	delete [] Num_block_T;

}
//============================================================================================================
