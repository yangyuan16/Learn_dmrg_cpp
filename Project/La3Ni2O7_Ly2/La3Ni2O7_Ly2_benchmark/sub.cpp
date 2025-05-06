#include<iostream>
using namespace std;
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_coupling.h>

#include"common.h"
#include"sub.h"

#define constant sqrt(6.0) * 0.5
#define hoping -sqrt(2.0)

//====================================BLAS ROUTINES===================================
extern "C" {
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
}

//===================================BLAS ROUTINES=====================================
extern "C" {
void dsymm_(char *side, char *uplo, const int *m, const int *n, const double *alpha, double *a, const int *lda, double *b, const int *ldb, const double *beta, double *c, const int *ldc);

void dgemm_(char *transa, char *transb, const int *m, const int *n, const int *k, const double *alpha, double *a, const int *lda, double *b, const int *ldb, const double *beta, double *c, const int *ldc);
}

//=========================Initialize system and environment blocks=====================
//      Initialize the system block: In this code, only a site in the initial step
//======================================================================================
Sub::Sub( Parameter &para ):SubCommon() {
	
	IndexNo = 0;			//control "create space"

	//------
	TotSiteNo = 1;			//Total site of subsystem
	Block_Num_hole = 2;		//Block number of hole quantum number, measured from half-filling

	CreateSpace_hole_block();	//hole number and J number in each block

	for( int i = 0; i < Block_Num_hole; i++ ) {

		Num_hole_block[i] = i;    //hole number: from small to large!!!
		Num_J_hole_block[i] = 1;  //hole=0, J=1/2; hole=1, J=0.

	}
	
	//------
	CreateSpace_J_block();		//J value and dimension in each J block

	Value_J_block[0][0] = 1;	Value_J_block[1][0] = 0;
	Dim_J_block[0][0] = 1;		Dim_J_block[1][0] = 1;		

	//------
	operator_number_J = 1;            //para.Table_sys[0]==para.Table_env[0]=1

	Create_S_Dia();
	S_Dia[0][0][0][0] = constant;

	Create_S_M_Dia();

	//------
	Create_NN();
	NN[0][0][0][0] = - hoping;

	//------
	operator_number_T = 1;
	CreateSpace_T();

	Create_T();
	T[0][0][0][0] = hoping;

	//------
	Create_H();
	
}

//==================================Read in data from hard disks===========================
//-----------Read in space variables---------------------
Sub::Sub( const int &i, const int &lsys ):SubCommon( i, lsys ) {}

//-----------Read in operators---------------------------
Sub::Sub( const int &i, const int &lsys, Sub &space ):SubCommon() {

	IndexNo = 0;

	//----
	TotSiteNo = space.TotSiteNo;
	Block_Num_hole = space.Block_Num_hole;

	CreateSpace_hole_block();       //hole number and J number in each block

        for( int i = 0; i < Block_Num_hole; i++ ) { 

                Num_hole_block[i] = space.Num_hole_block[i];  

		Num_J_hole_block[i] = space.Num_J_hole_block[i];

        }

	//----
        CreateSpace_J_block();          //J value and dimension in each J block

	for( int i = 0; i < Block_Num_hole; i++ )
	for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

		Value_J_block[i][j] = space.Value_J_block[i][j];	

		Dim_J_block[i][j] = space.Dim_J_block[i][j];

	}

	//----
	operator_number_J = space.operator_number_J;

	Create_S_Dia();

        FILE *f = fopen( Combine(Combine("mid/operator/S_Dia/", i), lsys), "rb" );
        for( int n = 0; n < operator_number_J; n++ )
        for( int i = 0; i < Block_Num_hole; i++ )
        for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

                int Dim = Dim_J_block[i][j] * Dim_J_block[i][j];
                fread( S_Dia[n][i][j], sizeof(double), Dim, f );

        }

	fclose(f);


	Create_S_M_Dia();

        FILE *p = fopen( Combine(Combine("mid/operator/S_M_Dia/", i), lsys), "rb" );
        for( int n = 0; n < operator_number_J; n++ )
        for( int i = 0; i < Block_Num_hole; i++ ) 
	for( int j = 0; j < Num_J_hole_block[i]-1; j++ ) {

                int Dim = Dim_J_block[i][j] * Dim_J_block[i][j+1];
                fread( S_M_Dia[n][i][j], sizeof(double), Dim, p );

        }

        fclose(p);

	//----
	Create_NN();

        FILE *v = fopen( Combine(Combine("mid/operator/NN/", i), lsys), "rb" );
        for( int n = 0; n < operator_number_J; n++ )
        for( int i = 0; i < Block_Num_hole; i++ )
        for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

		int Dim = Dim_J_block[i][j] * Dim_J_block[i][j];
                fread( NN[n][i][j], sizeof(double), Dim, v );

	}

	fclose(v);

	//----
	operator_number_T = space.operator_number_T;
        CreateSpace_T();

	Create_T();

        FILE *x = fopen( Combine(Combine("mid/operator/T/", i), lsys), "rb" );
        for( int n = 0; n < operator_number_T; n++ )
        for( int i = 0; i < Block_Num_hole - 1; i++ ) 
	for( int j = 0; j < Num_block_T[i]; j++ ) {

                int Dim = Dim_J_block[i][ J_block_T_bra[i][j] ] * Dim_J_block[ i + 1 ][ J_block_T_ket[i][j] ];
                fread( T[n][i][j], sizeof(double), Dim, x );

        }

        fclose(x);

	//----
	Create_H();

        FILE *h = fopen( Combine(Combine("mid/operator/H/", i), lsys), "rb" );
        for( int i = 0; i < Block_Num_hole; i++ )
        for( int j = 0; j < Num_J_hole_block[i]; j++ ) {

                int Dim = Dim_J_block[i][j] * Dim_J_block[i][j];
                fread( H[i][j], sizeof(double), Dim, h );

        }

	fclose(h);

}

//==========================New space for both sysnew and envnew===========================
Sub::Sub(const int &i, Parameter &para, Sub &old):SubCommon() {

        IndexNo = 1;

	//----
        TotSiteNo = old.TotSiteNo + 1;
	
	operator_number_J = para.Table_J_sys[ old.TotSiteNo ];  //sys.TotSiteNo==TotSiteNo-1
	operator_number_T = para.Table_T_sys[ old.TotSiteNo ];  //sys.TotSiteNo==TotSiteNo-1

	//----
	Find_Block_Num_hole( para, old );

//	cout<<"\n TotSiteNo="<<TotSiteNo<<"\t New Block_Num_hole="<<Block_Num_hole<<endl;

	CreateSpace_hole_block();       	//hole number and J number in each block

	//----
	Find_Num_hole_block( para, old );	//Find Num_hole_block

//	for(int i=0; i<Block_Num_hole; i++)	cout<<"\n Num_hole_block["<<i<<"]="<<Num_hole_block[i]<<endl;

	Find_Num_J_hole_block( para, old );	//Find Num_J_hole_block

	//----
	CreateSpace_T();

  //------Store the angular momentum magnitude for wave function transformation in class superenergy
        FILE *fw = fopen(Combine(Combine("new_block/", i), TotSiteNo), "wb");//"i"=1 (2) stand for sys (env).

        fwrite( &Block_Num_hole, sizeof(int), 1, fw );
	fwrite( Num_hole_block, sizeof(int), Block_Num_hole, fw );
	fwrite( Num_J_hole_block, sizeof(int), Block_Num_hole, fw );

	for(int i=0; i<Block_Num_hole; i++) {
	        fwrite(Value_J_block[i], sizeof(int), Num_J_hole_block[i], fw);
        	fwrite(Dim_J_block[i], sizeof(int), Num_J_hole_block[i], fw);
	}

	for(int i=0; i<Block_Num_hole; i++)
        for(int j=0; j<3; j++) {
                fwrite(Hole_blockOld[i][j], sizeof(int), Num_J_hole_block[i], fw);
                fwrite(J_blockOld[i][j], sizeof(int), Num_J_hole_block[i], fw);
		fwrite(Start[i][j], sizeof(int), Num_J_hole_block[i], fw);
        }

        fclose(fw);

	//----
        FILE *ft = fopen( Combine ( Combine ( "new_block_T/", i ), TotSiteNo ), "wb" );//"i"=1/2 for sys/env

	fwrite( Num_block_T, sizeof(int), Block_Num_hole - 1, ft );

	for( int i = 0; i < Block_Num_hole - 1; i++ ) {

	        fwrite( J_block_T_bra[i], sizeof(int), Num_block_T[i], ft );
        	fwrite( J_block_T_ket[i], sizeof(int), Num_block_T[i], ft );

	}

        fclose(ft);

}

//================================================================
//    New system and environment for finite sweep with only 
//creating new Hilbert space, but not to create and renew operators
//================================================================
Sub::Sub(Parameter &para, char *nooperator, Sub &old) {

        IndexNo = 1;

        TotSiteNo = old.TotSiteNo + 1;
	
	operator_number_J = para.Table_J_sys[ old.TotSiteNo ];  //sys.TotSiteNo==TotSiteNo-1
	operator_number_T = para.Table_T_sys[ old.TotSiteNo ];  //sys.TotSiteNo==TotSiteNo-1

	Find_Block_Num_hole(para, old);

	CreateSpace_hole_block();       	//hole number and J number in each block

	Find_Num_hole_block(para, old);		//Find Num_hole_block

	Find_Num_J_hole_block(para, old);	//Find Num_J_hole_block

	CreateSpace_T();

}

//============================New operators for sysnew========================
Sub::Sub(char &sign, Parameter &para, Sub &sys, Sub &sysnew_space):SubCommon() {

	IndexNo = 2;
	ope_sign = sign;

  //----
	TotSiteNo = sysnew_space.TotSiteNo;
	Block_Num_hole = sysnew_space.Block_Num_hole;

  //----
	CreateSpace_hole_block();       //hole number and J number in each block

        for(int i=0; i<Block_Num_hole; i++) { 

                Num_hole_block[i] = sysnew_space.Num_hole_block[i];  
		Num_J_hole_block[i] = sysnew_space.Num_J_hole_block[i];

        }

  //----
        CreateSpace_J_block();          //J value and dimension in each J block

	for(int i=0; i<Block_Num_hole; i++)
	for(int j=0; j<Num_J_hole_block[i]; j++) {

		Value_J_block[i][j] = sysnew_space.Value_J_block[i][j];  
		Dim_J_block[i][j] = sysnew_space.Dim_J_block[i][j];

	}

  //---- 
	operator_number_J = para.Table_J_sys[sys.TotSiteNo];  //sys.TotSiteNo==TotSiteNo-1
	operator_number_T = para.Table_T_sys[sys.TotSiteNo];  //sys.TotSiteNo==TotSiteNo-1

	CreateSpace_T();

  //----
	if( ope_sign == 'T' ) {

		Create_T();
                New_T(1, para, sysnew_space, sys, para.Table_T_sys_site[sys.TotSiteNo-1], para.Table_T_sys_site[TotSiteNo-1]);

	}

	else if( ope_sign == 'S' ) {

		Create_S_Dia();
                NewS_Dia(1, para, sysnew_space, sys, para.Table_J_sys_site[sys.TotSiteNo-1], para.Table_J_sys_site[TotSiteNo-1]);

	}

	else if( ope_sign == 'M' ) {

		Create_S_M_Dia();
		NewS_M_Dia(1, para, sysnew_space, sys, para.Table_J_sys_site[sys.TotSiteNo-1], para.Table_J_sys_site[TotSiteNo-1]);

	}

	else if( ope_sign == 'N' ) {

		Create_NN();
                New_NN(para, sysnew_space, sys, para.Table_J_sys_site[sys.TotSiteNo-1], para.Table_J_sys_site[TotSiteNo-1]);

	}

	else if( ope_sign == 'H' ) {

		Create_H();
		New_H_Sys(para, sysnew_space, sys);

	}

}

//==========================New operators for envnew==========================
Sub::Sub(Parameter &para, char &sign, Sub &env, Sub &envnew_space):SubCommon() {

	IndexNo = 2;
	ope_sign = sign;

	//----
	TotSiteNo = envnew_space.TotSiteNo;
	Block_Num_hole = envnew_space.Block_Num_hole;

	CreateSpace_hole_block();       //hole number and J number in each block

        for(int i=0; i<Block_Num_hole; i++) { 

                Num_hole_block[i]=envnew_space.Num_hole_block[i];    
		Num_J_hole_block[i]=envnew_space.Num_J_hole_block[i];

        }

	//----
        CreateSpace_J_block();          //J value and dimension in each J block

	for(int i=0; i<Block_Num_hole; i++)
	for(int j=0; j<Num_J_hole_block[i]; j++) {

		Value_J_block[i][j]=envnew_space.Value_J_block[i][j];	
		Dim_J_block[i][j]=envnew_space.Dim_J_block[i][j];

	}

//	StartSite = 2 * envnew_space.TotSiteNo - 1;   //For 1D chain

	int totalsite = 2 * TotSiteNo;  //Total site of sys+ns+ne+env

        for( int i = 1; i < para.N_x; i++ ) {

                if( ( i + 1 ) * para.N_u * para.N_y >= totalsite ) { //(i+1) is the number of column

                        StartSite = ( i + 1 ) * para.N_u * para.N_y - 1; //the first site of envnew!!!

                        break;

                }

        }

	//----
	operator_number_J = para.Table_J_env[env.TotSiteNo]; 
	operator_number_T = para.Table_T_env[env.TotSiteNo];

	CreateSpace_T();

	if(ope_sign=='T') {

		Create_T();
                New_T(2, para, envnew_space, env, para.Table_T_env_site[env.TotSiteNo-1], para.Table_T_env_site[TotSiteNo-1]);

	}

	else if(ope_sign=='S') {

		Create_S_Dia();
                NewS_Dia(2, para, envnew_space, env, para.Table_J_env_site[env.TotSiteNo-1], para.Table_J_env_site[TotSiteNo-1]);

	}

	else if(ope_sign=='M') {

		Create_S_M_Dia();
		NewS_M_Dia(2, para, envnew_space, env, para.Table_J_env_site[env.TotSiteNo-1], para.Table_J_env_site[TotSiteNo-1]);


	}

	else if(ope_sign=='N') {

		Create_NN();
                New_NN(para, envnew_space, env, para.Table_J_env_site[env.TotSiteNo-1], para.Table_J_env_site[TotSiteNo-1]);

	}

	else if(ope_sign=='H') {

		Create_H();
		New_H_Env(para, envnew_space, env);

	}

}

//=========================New operators for envnew===========================
Sub::Sub(Sub &env, Sub &envnew_space, Parameter &para, char &sign):SubCommon() {

	IndexNo=2;
	ope_sign=sign;

	//----
	TotSiteNo=envnew_space.TotSiteNo;
	Block_Num_hole=envnew_space.Block_Num_hole;

	//----
	CreateSpace_hole_block();       //hole number and J number in each block

        for(int i=0; i<Block_Num_hole; i++) { 

                Num_hole_block[i]=envnew_space.Num_hole_block[i];    
		Num_J_hole_block[i]=envnew_space.Num_J_hole_block[i];

        }

        CreateSpace_J_block();          //J value and dimension in each J block

	for(int i=0; i<Block_Num_hole; i++)
	for(int j=0; j<Num_J_hole_block[i]; j++) {
		
		Value_J_block[i][j]=envnew_space.Value_J_block[i][j];	
		Dim_J_block[i][j]=envnew_space.Dim_J_block[i][j];

	}

	//----
	StartSite = para.Total_N - 1; 

	operator_number_J=para.Table_J_env[env.TotSiteNo]; 
	operator_number_T=para.Table_T_env[env.TotSiteNo];

	CreateSpace_T();

	if(ope_sign=='T') {

		Create_T();
                New_T(2, para, envnew_space, env, para.Table_T_env_site[env.TotSiteNo-1], para.Table_T_env_site[TotSiteNo-1]);

	}

	else if(ope_sign=='S') {

		Create_S_Dia();
                NewS_Dia(2, para, envnew_space, env, para.Table_J_env_site[env.TotSiteNo-1], para.Table_J_env_site[TotSiteNo-1]);

	}

	else if(ope_sign=='M') {

		Create_S_M_Dia();
		NewS_M_Dia(2, para, envnew_space, env, para.Table_J_env_site[env.TotSiteNo-1], para.Table_J_env_site[TotSiteNo-1]);

	}

	else if(ope_sign=='N') {

		Create_NN();
                New_NN(para, envnew_space, env, para.Table_J_env_site[env.TotSiteNo-1], para.Table_J_env_site[TotSiteNo-1]);

	}

	else if(ope_sign=='H') {

		Create_H();
		New_H_Env(para, envnew_space, env);

	}

}

//===================================================================================
//                     Reduced density matrix of system and environment
//===================================================================================
Sub::Sub(Sub &space):SubCommon() {

        IndexNo=4;

	TotSiteNo=space.TotSiteNo;
	Block_Num_hole=space.Block_Num_hole;

	CreateSpace_hole_block();       //hole number and J number in each block
        for(int i=0; i<Block_Num_hole; i++) { 
                Num_hole_block[i]=space.Num_hole_block[i];    Num_J_hole_block[i]=space.Num_J_hole_block[i];
        }

        CreateSpace_J_block();          //J value and dimension in each J block
	for(int i=0; i<Block_Num_hole; i++)
	for(int j=0; j<Num_J_hole_block[i]; j++) {
		Value_J_block[i][j]=space.Value_J_block[i][j];	Dim_J_block[i][j]=space.Dim_J_block[i][j];
	}

	operator_number_J=space.operator_number_J;
	operator_number_T=space.operator_number_T;

	CreateSpace_T();

        CreateSpace3();

}

//=====================================Find_Block_Number_hole===================================
inline void Sub::Find_Block_Num_hole( Parameter &para, Sub &old ) {

	int * new_hole_block = new int [ 2 * old.Block_Num_hole ];

	for( int i = 0; i < old.Block_Num_hole; i++ ) {

		new_hole_block[ 2*i ] = old.Num_hole_block[i];	
		new_hole_block[ 2*i + 1 ] = old.Num_hole_block[i] + 1;

	}

	int symbol = new_hole_block[0];

	int counter = 1;

	for( int i = 1; i < 2*old.Block_Num_hole; i++ )
	if( new_hole_block[i] <= para.Total_h ) {	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		if( new_hole_block[i] != symbol ) {

			symbol = new_hole_block[i];	
			counter++;

		}

	}

	Block_Num_hole = counter;
	
	delete [] new_hole_block;

}

//=======================================Find Num_hole_block=============================
inline void Sub::Find_Num_hole_block( Parameter &para, Sub &old ) {

	int * new_hole_block = new int [ 2*old.Block_Num_hole ];

	for( int i = 0; i < old.Block_Num_hole; i++ ) {

		new_hole_block[2*i] = old.Num_hole_block[i];	
		new_hole_block[2*i+1] = old.Num_hole_block[i] + 1;

	}

	Num_hole_block[0] = new_hole_block[0];

	int counter = 0;

	for( int i = 1; i < 2*old.Block_Num_hole; i++ )
	if( new_hole_block[i] <= para.Total_h ) {	//!!!!!!!!!!!!!!!!!!!!!!
		
		if( new_hole_block[i] != Num_hole_block[counter] ) {

			counter++;
			Num_hole_block[counter] = new_hole_block[i];	

		}

	}

	delete [] new_hole_block;

}

//=======================================Find Num_J_hole_block=============================
inline void Sub::Find_Num_J_hole_block( Parameter &para, Sub &old ) {

  //----create space for new_J_hole_block--------------------------------------------------
	int * new_J_hole_block = new int [ 2*old.Block_Num_hole ];

	for( int i = 0; i < 2*old.Block_Num_hole; i++ ) {

		if( i%2 == 0 )		new_J_hole_block[i] = 2*old.Num_J_hole_block[i/2];  //coupled with spin-1/2

		else if( i%2 == 1 )	new_J_hole_block[i] = old.Num_J_hole_block[i/2];  //coupled with spin-0

	}
//*******
//	for(int i=0; i<2*old.Block_Num_hole; i++)	cout<<"\n new_J_hole_block["<<i<<"]="<<new_J_hole_block[i]<<endl;

  //----create new_J_block-------------------------------------------------------------
	int ** new_J_block = new int * [ 2*old.Block_Num_hole ];  //All the J value in each block

	for( int i = 0; i < 2*old.Block_Num_hole; i++ ) {

		new_J_block[i] = new int [new_J_hole_block[i]];

		for( int j = 0; j < new_J_hole_block[i]; j++ )
			new_J_block[i][j] = 0;

	}

	for( int i = 0; i < 2*old.Block_Num_hole; i++ ) {

		if( i%2 == 0 ) {

			for( int j = 0; j < new_J_hole_block[i]/2; j++ )
                                new_J_block[i][j] = old.Value_J_block[i/2][j] + 1;

                        for( int j = new_J_hole_block[i]/2; j < new_J_hole_block[i]; j++ )
                                new_J_block[i][j] = old.Value_J_block[i/2][j-new_J_hole_block[i]/2] - 1;

		}

		else if( i%2 == 1 ) {

			for( int j = 0; j < new_J_hole_block[i]; j++ )
				new_J_block[i][j] = old.Value_J_block[i/2][j];

		}

	}
//******
//	for(int i=0; i<2*old.Block_Num_hole; i++) 
//	for(int j=0; j<new_J_hole_block[i]; j++)
//		cout<<"\n new_J_block["<<i<<"]["<<j<<"]="<<new_J_block[i][j]<<endl;

  //----Create space for total_J_num-----------------------------------------------
	int * total_J_num = new int [Block_Num_hole];

	for( int i = 0; i < Block_Num_hole; i++ )  total_J_num[i]=0;

	for( int i = 0; i < Block_Num_hole; i++ ) {

		for( int j = 0; j < 2*old.Block_Num_hole; j++ )
		if( ( old.Num_hole_block[j/2] + j%2 ) == Num_hole_block[i] )
			total_J_num[i]+=new_J_hole_block[j];

	}
//******
//	for(int i=0; i<Block_Num_hole; i++)	cout<<"\n total_J_num["<<i<<"]="<<total_J_num[i]<<endl;

  //----Create space for total_J_value---------------------------------------------
	int ** total_J_value=new int * [Block_Num_hole];  int ** old_hole=new int * [Block_Num_hole];
	int ** old_J=new int * [Block_Num_hole];  int ** old_index=new int * [Block_Num_hole];
	int ** old_dim=new int * [Block_Num_hole];

	for(int i=0; i<Block_Num_hole; i++) {

		total_J_value[i]=new int [total_J_num[i]];  old_hole[i]=new int [total_J_num[i]];	
		old_J[i]=new int [total_J_num[i]];  old_index[i]=new int [total_J_num[i]];
		old_dim[i]=new int [total_J_num[i]];

		for(int j=0; j<total_J_num[i]; j++) {
			total_J_value[i][j]=0;  old_hole[i][j]=0;  old_J[i][j]=0;  old_index[i][j]=0;  old_dim[i][j]=0;
		}

	}

	for(int i=0; i<Block_Num_hole; i++) {

		int counter=0;

		for(int j=0; j<2*old.Block_Num_hole; j++)
		if((old.Num_hole_block[j/2]+j%2)==Num_hole_block[i]) {

			for(int k=0; k<new_J_hole_block[j]; k++) {

				total_J_value[i][counter]=new_J_block[j][k];

				if(j%2==0) {

                                        if(k/(new_J_hole_block[j]/2)==0) {
						old_hole[i][counter]=j/2;
						old_J[i][counter]=k;
						old_index[i][counter]=1;
                                                old_dim[i][counter]=old.Dim_J_block[j/2][k%(new_J_hole_block[j]/2)];
                                        }

                                        else if(k/(new_J_hole_block[j]/2)==1) {
						old_hole[i][counter]=j/2;
						old_J[i][counter]=k-new_J_hole_block[j]/2;	
						old_index[i][counter]=2;
                                                old_dim[i][counter]=old.Dim_J_block[j/2][k%(new_J_hole_block[j]/2)];
                                        }

				}

				else if(j%2==1) {
					old_hole[i][counter]=j/2;
					old_J[i][counter]=k;
					old_index[i][counter]=0;
					old_dim[i][counter]=old.Dim_J_block[j/2][k];
				}

				counter++;

			}
		
		}

	}
//******
//	for(int i=0; i<Block_Num_hole; i++) 
//	for(int j=0; j<total_J_num[i]; j++)
//		cout<<"\n total_J_value["<<i<<"]["<<j<<"]="<<total_J_value[i][j]<<endl;

  //----Find Num_J_hole_block-------------------------------------------------------
	for(int i=0; i<Block_Num_hole; i++) {

		if(total_J_value[i][0]>=0) {

			int counter=1;
	
			for(int j=1; j<total_J_num[i]; j++) 
			if(total_J_value[i][j]>=0) {

				int index=0;

				for(int k=0; k<j; k++) 
				if(total_J_value[i][k]==total_J_value[i][j]) 
					index++;

				if(index==0)	counter++;

			}

			Num_J_hole_block[i]=counter;	

		}

		else if(total_J_value[i][0]<0) {

			int counter=1;
	
			for(int j=2; j<total_J_num[i]; j++) 
			if(total_J_value[i][j]>=0) {

				int index=0;

				for(int k=1; k<j; k++) 
				if(total_J_value[i][k]==total_J_value[i][j]) 
					index++;

				if(index==0)	counter++;

			}

			Num_J_hole_block[i]=counter;	

		}

	}
//******
//	for(int i=0; i<Block_Num_hole; i++)	cout<<"\n Num_J_hole_block["<<i<<"]="<<Num_J_hole_block[i]<<endl;

  //----create space for J_block--------------------------------------------------
//******
//	cout<<"\n IndexNo="<<IndexNo<<endl;

	CreateSpace_J_block();

  //----Find Value_J_block and Dim_J_block----------------------------------------
	for(int i=0; i<Block_Num_hole; i++) {

		int counter=0, start=0;

		for(int j=0; j<total_J_num[i]; j++)
		if(total_J_value[i][j]>=0) {
			Value_J_block[i][0]=total_J_value[i][j];
			Dim_J_block[i][0]+=old_dim[i][j];
			Hole_blockOld[i][old_index[i][j]][0]=old_hole[i][j];
			J_blockOld[i][old_index[i][j]][0]=old_J[i][j];
			Start[i][old_index[i][j]][0]=0;
			counter++;
			start=j;
			break;
		}

		for(int j=start+1; j<total_J_num[i]; j++) 
		if(total_J_value[i][j]>=0) {

			int index=0;

			for(int k=0; k<counter; k++) 
			if(Value_J_block[i][k]==total_J_value[i][j]) { 
				Start[i][old_index[i][j]][k]=Dim_J_block[i][k];
				Dim_J_block[i][k]+=old_dim[i][j];
				Hole_blockOld[i][old_index[i][j]][k]=old_hole[i][j];
				J_blockOld[i][old_index[i][j]][k]=old_J[i][j];
				index++;
			}

			if(index==0) {
				Value_J_block[i][counter]=total_J_value[i][j];
				Dim_J_block[i][counter]+=old_dim[i][j];
				Hole_blockOld[i][old_index[i][j]][counter]=old_hole[i][j];
	 			J_blockOld[i][old_index[i][j]][counter]=old_J[i][j];
				Start[i][old_index[i][j]][counter]=0;
				counter++;
			}

		}

	}

  //----Rearrange the J from small to large---------------------------------------
	for(int i=0; i<Block_Num_hole; i++) {

		for(int j=0; j<Num_J_hole_block[i]; j++)
		for(int k=j+1; k<Num_J_hole_block[i]; k++) {

			if(Value_J_block[i][j]>Value_J_block[i][k]) {

				int m=Value_J_block[i][j];
				Value_J_block[i][j]=Value_J_block[i][k];
				Value_J_block[i][k]=m;

				m=Dim_J_block[i][j];
				Dim_J_block[i][j]=Dim_J_block[i][k];
				Dim_J_block[i][k]=m;

				for(int n=0; n<3; n++) {
					m=Hole_blockOld[i][n][j]; Hole_blockOld[i][n][j]=Hole_blockOld[i][n][k]; Hole_blockOld[i][n][k]=m;

					m=J_blockOld[i][n][j]; J_blockOld[i][n][j]=J_blockOld[i][n][k]; J_blockOld[i][n][k]=m;
					m=Start[i][n][j];    Start[i][n][j]=Start[i][n][k];       Start[i][n][k]=m;
				}

			}

		}

	}
//******
/*	for(int i=0; i<Block_Num_hole; i++) 
	for(int j=0; j<Num_J_hole_block[i]; j++) 
		cout<<"\n Value_J["<<i<<"]["<<j<<"]="<<Value_J_block[i][j]<<"\t Dim_J["<<i<<"]["<<j<<"]="<<Dim_J_block[i][j]<<endl;

	for(int i=0; i<Block_Num_hole; i++) 
	for(int n=0; n<3; n++)
	for(int j=0; j<Num_J_hole_block[i]; j++) 
		cout<<"\n HoldOld["<<i<<"]["<<n<<"]["<<j<<"]="<<Hole_blockOld[i][n][j]<<"\t J_old["<<i<<"]["<<n<<"]["<<j<<"]="<<J_blockOld[i][n][j]<<"\t Start["<<i<<"]["<<n<<"]["<<j<<"]="<<Start[i][n][j]<<endl;
*/

  //----delete spaces-------------------------------------------------------------
	for(int i=0; i<Block_Num_hole; i++) {
		delete [] total_J_value[i];  delete [] old_hole[i];  delete [] old_J[i];  
		delete [] old_index[i];   delete [] old_dim[i];
	}
	delete [] total_J_value;  delete [] old_hole;  delete [] old_J;  delete [] old_index;  delete [] old_dim;

	delete [] total_J_num;

	for(int i=0; i<2*old.Block_Num_hole; i++)
		delete [] new_J_block[i];
	delete [] new_J_block;

	delete [] new_J_hole_block;

}

//==================================New_T=====================================
inline void Sub::New_T( const int &block, Parameter &para, Sub &space, Sub &old, int *table_old, int *table_new ) {

  //----Variables in the subroutine
	int old_hm, old_hl, old_Jm, old_Jl, position, SubDim, a_new, old_J, sign;
	double factor;

  //----Create 6j_T
	double **** six_j_T;

	six_j_T = new double *** [ Block_Num_hole - 1 ];

	for( int n = 0; n < Block_Num_hole - 1; n ++ ) {

		six_j_T[ n ] = new double ** [ Num_block_T[ n ] ];

		for( int i = 0; i < Num_block_T[ n ]; i ++ ) {

			six_j_T[ n ][ i ] = new double * [ 3 ];

			for( int j = 0; j < 3; j ++ ) {

				six_j_T[ n ][ i ][ j ] = new double [ 3 ];

				for( int k = 0; k < 3; k ++ )

					six_j_T[ n ][ i ][ j ][ k ] = 0.0;
			
			}

		}

	}

	//Initialization
	for( int n = 0; n < Block_Num_hole - 1; n ++ )
	for( int i = 0; i < Num_block_T[ n ]; i ++ )
	for( int m = 0; m < 3; m ++ )
	for( int l = 0; l < 3; l ++ )
	if( ( old_hm = space.Hole_blockOld[ n ][ m ][ J_block_T_bra[ n ][ i ] ] ) != -1  &&  ( old_hl = space.Hole_blockOld[ n + 1 ][ l ][ J_block_T_ket[ n ][ i ] ] ) != -1  &&  ( old_Jm = space.J_blockOld[ n ][ m ][ J_block_T_bra[ n ][ i ] ] ) != -1  &&  ( old_Jl = space.J_blockOld[ n + 1 ][ l ][ J_block_T_ket[ n ][ i ] ] ) != -1  &&  ( old.Num_hole_block[ old_hl ] - old.Num_hole_block[ old_hm ] ) == 1  &&  abs( old.Value_J_block[ old_hm ][ old_Jm ] - old.Value_J_block[ old_hl ][ old_Jl ] ) == 1 ) { 

		sign = 1;

		if( m == 0  &&  l == 0 )  
			sign = 0;

		six_j_T[ n ][ i ][ m ][ l ] = gsl_sf_coupling_6j( old.Value_J_block[ old_hm ][ old_Jm ], Value_J_block[ n ][ J_block_T_bra[ n ][ i ] ], sign, Value_J_block[ n + 1 ][ J_block_T_ket[ n ][ i ] ], old.Value_J_block[ old_hl ][ old_Jl ], 1 );

	}

	//save
	FILE *fp = fopen(Combine(Combine("6j_factor/T_old/", block),  space.TotSiteNo), "wb");

        for( int n = 0; n < Block_Num_hole - 1; n ++ )
	for( int i = 0; i < Num_block_T[ n ]; i ++ )
        for( int m = 0; m < 3; m ++ ) {

                fwrite( six_j_T[ n ][ i ][ m ], sizeof(double), 3, fp );

	}

        fclose(fp);

  //----New T operators
	for( int site = 0; site < operator_number_T; site ++ ) {  // number of new site

		for( int site_old = 0; site_old < old.operator_number_T; site_old ++ )  // number of old site
		if( table_old[ site_old ] == table_new[ site ] ) {

			for( int n = 0; n < Block_Num_hole - 1; n ++ )  // with given hole number
			for( int i = 0; i < Num_block_T[ n ]; i ++ ) {  // with given J numbers

				for( int m = 0; m < 3; m ++ )  // degree of freedom of the new site, bra
				for( int l = 0; l < 3; l ++ )  // degree of freedom of the new site, ket
				if( ( old_hm = space.Hole_blockOld[ n ][ m ][ J_block_T_bra[ n ][ i ] ] ) != -1  &&  ( old_hl = space.Hole_blockOld[ n + 1 ][ l ][ J_block_T_ket[ n ][ i ] ] ) != -1  &&  ( old_Jm = space.J_blockOld[ n ][ m ][ J_block_T_bra[ n ][ i ] ] ) != -1  &&  ( old_Jl = space.J_blockOld[ n + 1 ][ l ][ J_block_T_ket[ n ][ i ] ] ) != -1  &&  ( old.Num_hole_block[ old_hl ] - old.Num_hole_block[ old_hm ] ) == 1  &&  abs( old.Value_J_block[ old_hm ][ old_Jm ] - old.Value_J_block[ old_hl ][ old_Jl ] ) == 1 ) {  // quantum numbers fixed

					//---- J numbers of old T-block
					old_J = 0;

					for( int old_i = 0; old_i < old.Num_block_T[ old_hm ]; old_i ++ )
					if( old.J_block_T_bra[ old_hm ][ old_i ] == old_Jm  &&  old.J_block_T_ket[ old_hm ][ old_i ] == old_Jl )
						old_J = old_i;  // old_J denotes J numbers of old T-block

					//---- sign is the contribution of "1/2 + s" in the "factor"
					sign = 2;

					if( m == 0  &&  l == 0 )  sign = 1;

					//----
					position = Dim_J_block[ n ][ J_block_T_bra[ n ][ i ] ] * space.Start[ n + 1 ][ l ][ J_block_T_ket[ n ][ i ] ] + space.Start[ n ][ m ][ J_block_T_bra[ n ][ i ] ];

					factor = pow( -1.0, ( old.Value_J_block[ old_hm ][ old_Jm ] + Value_J_block[ n + 1 ][ J_block_T_ket[ n ][ i ] ] + sign ) / 2 ) * sqrt( ( Value_J_block[ n ][ J_block_T_bra[ n ][ i ] ] + 1.0 ) * ( Value_J_block[ n + 1 ][ J_block_T_ket[ n ][ i ] ] + 1.0 ) ) * six_j_T[ n ][ i ][ m ][ l ];

					SubDim = old.Dim_J_block[ old_hm ][ old_Jm ] * old.Dim_J_block[ old_hl ][ old_Jl ];

					//----
					for( int a_old = 0; a_old < SubDim; a_old ++ ) {

						a_new = position + Dim_J_block[ n ][ J_block_T_bra[ n ][ i ] ] * ( a_old / old.Dim_J_block[ old_hm ][ old_Jm ] ) + a_old % old.Dim_J_block[ old_hm ][ old_Jm ];

						T[ site ][ n ][ i ][ a_new ] += factor * old.T[ site_old ][ old_hm ][ old_J ][ a_old ];

					}

				}

			}

		}

		if( table_new[ site ] == TotSiteNo - 1 ) {

			for( int n = 0; n < Block_Num_hole - 1; n ++ ) 
			for( int i = 0; i < Num_block_T[ n ]; i ++ ) {

				for( int m = 0; m < 3; m ++ )  // bra label
				for( int l = 0; l < 3; l ++ )  // ket label
				if( m != 0  &&  l == 0  &&  ( old_hm = space.Hole_blockOld[ n ][ m ][ J_block_T_bra[ n ][ i ] ] ) != -1  &&  ( old_hl = space.Hole_blockOld[ n + 1 ][ l ][ J_block_T_ket[ n ][ i ] ] ) != -1  &&  ( old_Jm = space.J_blockOld[ n ][ m ][ J_block_T_bra[ n ][ i ] ] ) != -1  &&  ( old_Jl = space.J_blockOld[ n + 1 ][ l ][ J_block_T_ket[ n ][ i ] ] ) != -1  &&  old_hm == old_hl  &&  old_Jm == old_Jl ) {
	
					//----				
					position = Dim_J_block[ n ][ J_block_T_bra[ n ][ i ] ] * space.Start[ n + 1 ][ l ][ J_block_T_ket[ n ][ i ] ] + space.Start[ n ][ m ][ J_block_T_bra[ n ][ i ] ];

					factor = - sqrt( Value_J_block[ n ][ J_block_T_bra[ n ][ i ] ] + 1.0 );
					//multipled with the hoping term -sqrt(2).

					//--- phase factor
					if( ( old.TotSiteNo - old.Num_hole_block[ old_hl ] ) % 2 == 1 )
					factor = - factor;

					for( int a_old = 0; a_old < old.Dim_J_block[ old_hm ][ old_Jm ]; a_old ++ ) {
						a_new = position + ( Dim_J_block[ n ][ J_block_T_bra[ n ][ i ] ] + 1 ) * a_old;
						T[ site ][ n ][ i ][ a_new ] += factor;
					}

				}

			}

		}

	}

//*******
/*	for(int site=0; site<operator_number_T; site++)
	for(int n=0; n<Block_T; n++)
	for(int i=0; i<Num_block_T[n]; i++)
	for(int j=0; j<Dim_J_block[Index_block_T[n]][J_block_T_bra[n][i]]*Dim_J_block[Index_block_T[n]+1][J_block_T_ket[n][i]]; j++)
		cout<<"\n T["<<site<<"]["<<n<<"]["<<i<<"]["<<j<<"]="<<T[site][n][i][j];*/
//*******

  //----Delete 6j_T
	for(int n=0; n < Block_Num_hole - 1; n++) {

		for(int i=0; i<Num_block_T[n]; i++) {

			for(int j=0; j<3; j++) {

				delete [] six_j_T[n][i][j];

			}

		delete [] six_j_T[n][i];

		}

	delete [] six_j_T[n];

	}

	delete [] six_j_T;

}

//=====================================NewS_Dia========================================
//                      S_Dia matrix is also Hermitian matrix
//=====================================================================================
inline void Sub::NewS_Dia(const int &block, Parameter &para, Sub &space, Sub &old, int *table_old, int *table_new) {

  //----Define variables in the subroutine
	int oldhm, oldhl, oldJm, oldJl, sign, position, SubDim, a_new, a_new_t;
	double factor;

  //----Create_6j_S_Dia
        double ****six_j_S_Dia;

        six_j_S_Dia = new double *** [ Block_Num_hole ];

        for(int n=0; n<Block_Num_hole; n++) {

                six_j_S_Dia[n] = new double ** [Num_J_hole_block[n]];

                for(int i=0; i<Num_J_hole_block[n]; i++) {

                        six_j_S_Dia[n][i]=new double * [3];

                        for(int j=0; j<3; j++) { 

                                six_j_S_Dia[n][i][j]=new double [3];

				for(int k=0; k<3; k++) 

					six_j_S_Dia[n][i][j][k]=0.0;

                        }

                }

        }

	//Initialize
	for(int n=0; n<Block_Num_hole; n++) 
	if(Num_hole_block[n] != TotSiteNo) {

		for(int i=0; i<Num_J_hole_block[n]; i++) 
		if( Value_J_block[n][i] != 0 ) {

	        	for(int j=0; j<3; j++)  
			for(int k=0; k<3; k++) 
			if( ( oldhm = space.Hole_blockOld[n][j][i] ) != -1  &&  space.Hole_blockOld[n][j][i] == space.Hole_blockOld[n][k][i]  &&  ( oldJm = space.J_blockOld[n][j][i] ) != -1  &&  ( oldJl = space.J_blockOld[n][k][i] ) != -1 ) {

				if( j == 0  &&  k == 0 )  //the added site has a hole

					six_j_S_Dia[n][i][j][k] = gsl_sf_coupling_6j( old.Value_J_block[oldhm][oldJm], Value_J_block[n][i], 0, Value_J_block[n][i], old.Value_J_block[oldhm][oldJl], 2 );

				else  //the added site has an electron
	
					six_j_S_Dia[n][i][j][k] = gsl_sf_coupling_6j( old.Value_J_block[oldhm][oldJm], Value_J_block[n][i], 1, Value_J_block[n][i], old.Value_J_block[oldhm][oldJl], 2 );

			}

		}

	}

	//save
        FILE *fp = fopen(Combine(Combine("6j_factor/S_Dia_old/", block),  space.TotSiteNo), "wb");

        for(int n=0; n<Block_Num_hole; n++)
	for(int i=0; i<Num_J_hole_block[n]; i++)
        for(int j=0; j<3; j++) {

                fwrite(six_j_S_Dia[n][i][j], sizeof(double), 3, fp);

	}

        fclose(fp);

//------Create_6j_S_Dia_N
        double ***six_j_S_Dia_N;

        six_j_S_Dia_N=new double ** [Block_Num_hole];

        for(int n=0; n<Block_Num_hole; n++) {

                six_j_S_Dia_N[n]=new double * [Num_J_hole_block[n]];

                for(int i=0; i<Num_J_hole_block[n]; i++) {

			six_j_S_Dia_N[n][i]=new double [3];

			for(int j=0; j<3; j++)

                        	six_j_S_Dia_N[n][i][j]=0.0;

                }

        }

	//Initialize
        for(int n=0; n<Block_Num_hole; n++) 
	if( Num_hole_block[n] != TotSiteNo ) {

		for(int i=0; i<Num_J_hole_block[n]; i++) 
		if( Value_J_block[n][i] != 0 ) {

			for(int j=1; j<3; j++)
			if((oldhm=space.Hole_blockOld[n][j][i])!=-1 && (oldJm=space.J_blockOld[n][j][i])!=-1) {

				six_j_S_Dia_N[n][i][j] = gsl_sf_coupling_6j( Value_J_block[n][i], 2, Value_J_block[n][i], 1, old.Value_J_block[oldhm][oldJm], 1 );

			}

		}

	}

	//save
        FILE *fr=fopen(Combine(Combine("6j_factor/S_Dia_n/", block), space.TotSiteNo), "wb");

	for(int n=0; n<Block_Num_hole; n++)
        for(int i=0; i<Num_J_hole_block[n]; i++) {

                fwrite(six_j_S_Dia_N[n][i], sizeof(double), 3, fr);

	}

        fclose(fr);

  //----New S_Dia operator
	for( int site = 0; site < operator_number_J; site++ ) {

		for( int site_old = 0; site_old < old.operator_number_J; site_old++ )
		if( table_old[ site_old ] == table_new[ site ] ) {

			for( int n = 0; n < Block_Num_hole; n++ )
			if( Num_hole_block[n] != TotSiteNo ) {	//hole number<=site number

				for( int i = 0; i < Num_J_hole_block[n]; i++ )
				if( Value_J_block[n][i] != 0 ) {

					for( int m = 0; m < 3; m++ )
					for( int l = 0; l < 3; l++ )
					if( (oldhm = space.Hole_blockOld[n][m][i]) != -1  &&  (oldhl = space.Hole_blockOld[n][l][i]) != -1  &&  (oldJm = space.J_blockOld[n][m][i]) != -1  &&  (oldJl = space.J_blockOld[n][l][i]) != -1  &&  oldhm == oldhl  &&  old.Num_hole_block[ oldhm ] != old.TotSiteNo ) {

						if( m == l  &&  oldJm == oldJl  &&  old.Value_J_block[oldhm][oldJm] != 0 ) {

							position = (Dim_J_block[n][i] + 1) * space.Start[n][m][i];

							if( m == 0 )  factor = 1.0;			
							else 	      factor = pow( -1.0, ( old.Value_J_block[oldhm][oldJm] + Value_J_block[n][i] + 3 ) / 2 ) * (Value_J_block[n][i] + 1.0 ) * six_j_S_Dia[n][i][m][l];

							SubDim = old.Dim_J_block[oldhm][oldJm] * old.Dim_J_block[oldhl][oldJl];

							for( int a_old = 0; a_old < SubDim; a_old++ ) {

								a_new = position + Dim_J_block[n][i] * ( a_old / old.Dim_J_block[oldhm][oldJm] ) + a_old % old.Dim_J_block[oldhm][oldJm];

								S_Dia[site][n][i][a_new] += factor * old.S_Dia[site_old][oldhm][oldJm][a_old];

							}

						}

						else if( oldJl - oldJm  == 1  &&  old.Value_J_block[oldhl][oldJl] - old.Value_J_block[oldhm][oldJm] == 2 ) {

							position = Dim_J_block[n][i] * space.Start[n][l][i] + space.Start[n][m][i];

							factor = pow(-1.0, ( old.Value_J_block[oldhm][oldJm] + Value_J_block[n][i] + 3 ) / 2 ) * (Value_J_block[n][i] + 1.0) * six_j_S_Dia[n][i][m][l];

							SubDim = old.Dim_J_block[oldhm][oldJm] * old.Dim_J_block[oldhl][oldJl];

							for(int a_old=0; a_old<SubDim; a_old++) {

								a_new = position + Dim_J_block[n][i] * (a_old/old.Dim_J_block[oldhm][oldJm]) + a_old%old.Dim_J_block[oldhm][oldJm];
								a_new_t = Dim_J_block[n][i] * (a_new%Dim_J_block[n][i]) + a_new/Dim_J_block[n][i];

								S_Dia[site][n][i][a_new] += factor * old.S_M_Dia[site_old][oldhm][oldJm][a_old];
								S_Dia[site][n][i][a_new_t] += factor * old.S_M_Dia[site_old][oldhm][oldJm][a_old];

							}

						}

					}

				}

			}

		}

		if( table_new[site] == TotSiteNo - 1 ) {

			for( int n = 0; n < Block_Num_hole; n++ )
			if( Num_hole_block[n] != TotSiteNo ) {	//hole number<=site number

				for( int i = 0; i < Num_J_hole_block[n]; i++ )
				if( Value_J_block[n][i] != 0 ) {

					for( int m = 1; m < 3; m++ ) 
					if( (oldhm=space.Hole_blockOld[n][m][i]) != -1  &&  (oldJm=space.J_blockOld[n][m][i]) != -1 ) {

						position = Dim_J_block[n][i] + 1;

						factor = pow(-1.0, (old.Value_J_block[oldhm][oldJm] + Value_J_block[n][i] + 3) / 2 ) * (Value_J_block[n][i] + 1.0) * constant * six_j_S_Dia_N[n][i][m];

						for(int a_old=0; a_old<old.Dim_J_block[oldhm][oldJm]; a_old++) {

							a_new = position * ( space.Start[n][m][i] + a_old );

							S_Dia[site][n][i][a_new] += factor;

						}

					}

				}

			}

		}

	}

  //----Delete 6j_S_Dia
	for(int n=0; n<Block_Num_hole; n++) {
		for(int i=0; i<Num_J_hole_block[n]; i++) {
			for(int j=0; j<3; j++) {
				delete [] six_j_S_Dia[n][i][j];
			}
		delete [] six_j_S_Dia[n][i];
		}
	delete [] six_j_S_Dia[n];
	}
	delete [] six_j_S_Dia;

  //----Delete 6j_S_Dia_N
	for(int n=0; n<Block_Num_hole; n++) {
		for(int i=0; i<Num_J_hole_block[n]; i++) {
			delete [] six_j_S_Dia_N[n][i];
		}
		delete [] six_j_S_Dia_N[n];
	}
	delete [] six_j_S_Dia_N;

}

//=====================================NewS_M_Dia===================================
inline void Sub::NewS_M_Dia(const int &block, Parameter &para, Sub &space, Sub &old, int *table_old, int *table_new) {

  //----Variables in the subroutine
	int oldhm, oldhl, oldJm, oldJl, position, SubDim, a_new;
	double factor;

  //----Create 6j_S_M_Dia
	double ****six_j_S_M_Dia;
        six_j_S_M_Dia=new double *** [Block_Num_hole];
        for(int n=0; n<Block_Num_hole; n++) {
                six_j_S_M_Dia[n]=new double ** [Num_J_hole_block[n]-1];
                for(int i=0; i<Num_J_hole_block[n]-1; i++) {
                        six_j_S_M_Dia[n][i]=new double * [3];
                        for(int j=0; j<3; j++) { 
                                six_j_S_M_Dia[n][i][j]=new double [3];
				for(int k=0; k<3; k++) 
					six_j_S_M_Dia[n][i][j][k]=0.0;
                        }
                }
        }

	//Initialization
	for(int n=0; n<Block_Num_hole; n++) 
	if(Num_hole_block[n] != TotSiteNo) {

		for(int i=0; i<Num_J_hole_block[n]-1; i++) 
		if((Value_J_block[n][i+1]-Value_J_block[n][i])==2) {

        		for(int j=0; j<3; j++)  
			for(int k=0; k<3; k++) 
			if((oldhm=space.Hole_blockOld[n][j][i])!=-1 && space.Hole_blockOld[n][j][i]==space.Hole_blockOld[n][k][i+1] && (old.Num_hole_block[oldhm] != old.TotSiteNo) && (oldJm=space.J_blockOld[n][j][i])!=-1 && (oldJl=space.J_blockOld[n][k][i+1])!=-1 && (old.Value_J_block[oldhm][oldJl]-old.Value_J_block[oldhm][oldJm])<=2) {

				if(j==0 && k==0)
				six_j_S_M_Dia[n][i][j][k]=gsl_sf_coupling_6j(old.Value_J_block[oldhm][oldJm], Value_J_block[n][i], 0, Value_J_block[n][i+1], old.Value_J_block[oldhm][oldJl], 2);	

				else	six_j_S_M_Dia[n][i][j][k]=gsl_sf_coupling_6j(old.Value_J_block[oldhm][oldJm], Value_J_block[n][i], 1, Value_J_block[n][i+1], old.Value_J_block[oldhm][oldJl], 2);

			}
		}
 	}

	//Save
        FILE *fp=fopen(Combine(Combine("6j_factor/S_M_Dia_old/", block), space.TotSiteNo), "wb");

        for(int n=0; n<Block_Num_hole; n++)
	for(int i=0; i<Num_J_hole_block[n]-1; i++)
        for(int j=0; j<3; j++) {

                fwrite(six_j_S_M_Dia[n][i][j], sizeof(double), 3, fp);

	}

        fclose(fp);

//------Create_6j_S_M_Dia_N
        double ****six_j_S_M_Dia_N;
        six_j_S_M_Dia_N=new double *** [Block_Num_hole];
        for(int n=0; n<Block_Num_hole; n++) {
                six_j_S_M_Dia_N[n]=new double ** [Num_J_hole_block[n]-1];
                for(int i=0; i<Num_J_hole_block[n]-1; i++) {
			six_j_S_M_Dia_N[n][i]=new double * [3];
			for(int j=0; j<3; j++) {
				six_j_S_M_Dia_N[n][i][j]=new double [3];
				for(int k=0; k<3; k++)
 		                       	six_j_S_M_Dia_N[n][i][j][k]=0.0;
			}
                }
        }

	//Initialize
        for(int n=0; n<Block_Num_hole; n++) 
	if(Num_hole_block[n]!=TotSiteNo) {

		for(int i=0; i<Num_J_hole_block[n]-1; i++) 
		if((Value_J_block[n][i+1]-Value_J_block[n][i])==2) {

			for(int j=1; j<3; j++)
			for(int k=1; k<3; k++)
			if((oldhm=space.Hole_blockOld[n][j][i])!=-1 && space.Hole_blockOld[n][j][i]==space.Hole_blockOld[n][k][i+1] && (oldJm=space.J_blockOld[n][j][i])!=-1 && space.J_blockOld[n][j][i]==space.J_blockOld[n][k][i+1]) {

				six_j_S_M_Dia_N[n][i][j][k]=gsl_sf_coupling_6j(1, Value_J_block[n][i], old.Value_J_block[oldhm][oldJm], Value_J_block[n][i+1], 1, 2);

			}

		}
	}

	//save
        FILE *fr=fopen(Combine(Combine("6j_factor/S_M_Dia_n/", block), space.TotSiteNo), "wb");

	for(int n=0; n<Block_Num_hole; n++)
        for(int i=0; i<Num_J_hole_block[n]-1; i++)
	for(int j=0; j<3; j++) {

                fwrite(six_j_S_M_Dia_N[n][i][j], sizeof(double), 3, fr);

	}

        fclose(fr);

  //----New S_M_Dia
	for(int site=0; site<operator_number_J; site++) {

		for(int site_old=0; site_old<old.operator_number_J; site_old++)
		if(table_old[site_old]==table_new[site]) {

			for(int n=0; n<Block_Num_hole; n++)
			if(Num_hole_block[n] != TotSiteNo) {	//hole number<=site number

				for(int i=0; i<Num_J_hole_block[n]-1; i++)
				if((Value_J_block[n][i+1]-Value_J_block[n][i])==2) {

					for(int m=0; m<3; m++)
					for(int l=0; l<3; l++)
					if( (oldhm=space.Hole_blockOld[n][m][i]) != -1  &&  space.Hole_blockOld[n][m][i] == space.Hole_blockOld[n][l][i+1]  &&  ( old.Num_hole_block[oldhm] != old.TotSiteNo )  &&  (oldJm=space.J_blockOld[n][m][i]) != -1  &&  (oldJl=space.J_blockOld[n][l][i+1]) != -1  &&  (old.Value_J_block[oldhm][oldJl]-old.Value_J_block[oldhm][oldJm]) <= 2 ) {

						if((oldJm+1)==oldJl && (old.Value_J_block[oldhm][oldJl]-old.Value_J_block[oldhm][oldJm])==2) {

							position=Dim_J_block[n][i]*space.Start[n][l][i+1]+space.Start[n][m][i];

							if(m==0 && l==0)  factor=1.0;

							else  factor=pow(-1.0, (old.Value_J_block[oldhm][oldJm]+Value_J_block[n][i+1]+3)/2)*sqrt((Value_J_block[n][i]+1.0)*(Value_J_block[n][i+1]+1.0))*six_j_S_M_Dia[n][i][m][l];

							SubDim=old.Dim_J_block[oldhm][oldJm]*old.Dim_J_block[oldhm][oldJl];

							for(int a_old=0; a_old<SubDim; a_old++) {

								a_new=position+Dim_J_block[n][i]*(a_old/old.Dim_J_block[oldhm][oldJm])+a_old%old.Dim_J_block[oldhm][oldJm];

								S_M_Dia[site][n][i][a_new]+=factor*old.S_M_Dia[site_old][oldhm][oldJm][a_old];

							}

						}
	
						else if(oldJm==oldJl && old.Value_J_block[oldhm][oldJm]!=0) {

							position=Dim_J_block[n][i]*space.Start[n][l][i+1]+space.Start[n][m][i];	

							factor=pow(-1.0, (old.Value_J_block[oldhm][oldJm]+Value_J_block[n][i+1]+3)/2)*sqrt((Value_J_block[n][i]+1.0)*(Value_J_block[n][i+1]+1.0))*six_j_S_M_Dia[n][i][m][l];

							SubDim=old.Dim_J_block[oldhm][oldJm]*old.Dim_J_block[oldhm][oldJl];

							for(int a_old=0; a_old<SubDim; a_old++) {

								a_new=position+Dim_J_block[n][i]*(a_old/old.Dim_J_block[oldhm][oldJm])+a_old%old.Dim_J_block[oldhm][oldJm];
								S_M_Dia[site][n][i][a_new]+=factor*old.S_Dia[site_old][oldhm][oldJm][a_old];

							}

						}

					}
				}

			}
	
		}

                if(table_new[site]==TotSiteNo-1) {

			for(int n=0; n<Block_Num_hole; n++)
			if(Num_hole_block[n] != TotSiteNo) {	//hole number<=site number

				for(int i=0; i<Num_J_hole_block[n]-1; i++)
				if((Value_J_block[n][i+1]-Value_J_block[n][i])==2) {

					for(int m=1; m<3; m++)
					for(int l=1; l<3; l++)
					if((oldhm=space.Hole_blockOld[n][m][i])!=-1 && space.Hole_blockOld[n][m][i]==space.Hole_blockOld[n][l][i+1] && (oldJm=space.J_blockOld[n][m][i])!=-1 && space.J_blockOld[n][m][i]==space.J_blockOld[n][l][i+1]) {

						position=Dim_J_block[n][i]*space.Start[n][l][i+1]+space.Start[n][m][i];

						factor=pow(-1.0, (old.Value_J_block[oldhm][oldJm]+Value_J_block[n][i]+3)/2)*sqrt((Value_J_block[n][i]+1.0)*(Value_J_block[n][i+1]+1.0))*six_j_S_M_Dia_N[n][i][m][l]*constant;

						for(int a_old=0; a_old<old.Dim_J_block[oldhm][oldJm]; a_old++) {

							a_new=position+(Dim_J_block[n][i]+1)*a_old;
							S_M_Dia[site][n][i][a_new]+=factor;

						}

					}

				}

			}

		}

	}

  //----Delete 6j_S_M_Dia
	for(int n=0; n<Block_Num_hole; n++) {
		for(int i=0; i<Num_J_hole_block[n]-1; i++) {
			for(int j=0; j<3; j++) {
				delete [] six_j_S_M_Dia[n][i][j];
			}
		delete [] six_j_S_M_Dia[n][i];
		}
	delete [] six_j_S_M_Dia[n];
	}
	delete [] six_j_S_M_Dia;

  //----Delete 6j_S_M_Dia_N
	for(int n=0; n<Block_Num_hole; n++) {
		for(int i=0; i<Num_J_hole_block[n]-1; i++) {
			for(int j=0; j<3; j++) {
				delete [] six_j_S_M_Dia_N[n][i][j];
			}
		delete [] six_j_S_M_Dia_N[n][i];
		}
	delete [] six_j_S_M_Dia_N[n];
	}
	delete [] six_j_S_M_Dia_N;

}

//==============================New_N===================================
inline void Sub::New_NN(Parameter &para, Sub &space, Sub &old, int *table_old, int *table_new) {

	int oldhm, oldJm, position, SubDim, a_new;
	double factor;

  //----New_N_Sys_Old
	for( int site = 0; site < operator_number_J; site++ ) {

		for( int site_old = 0; site_old < old.operator_number_J; site_old++ )
		if( table_old[ site_old ] == table_new[ site ] ) {

			for( int n = 0; n < Block_Num_hole; n++ )
			if( Num_hole_block[ n ] != TotSiteNo ) {

				for( int i = 0; i < Num_J_hole_block[ n ]; i++ )
				for( int m = 0; m < 3; m++ )
				if( ( oldhm = space.Hole_blockOld[ n ][ m ][ i ] ) != -1  &&  ( oldJm = space.J_blockOld[ n ][ m ][ i ] ) != -1  &&  old.Num_hole_block[ oldhm ] != old.TotSiteNo ) {

					//----
					position = ( Dim_J_block[ n ][ i ] + 1 ) * space.Start[ n ][ m ][ i ];

					SubDim = old.Dim_J_block[ oldhm ][ oldJm ] * old.Dim_J_block[ oldhm ][ oldJm ];

					//----
					factor = sqrt( ( Value_J_block[ n ][ i ] + 1.0 ) / ( old.Value_J_block[ oldhm ][ oldJm ] + 1.0 ) );

					//----
					for( int a_old = 0; a_old < SubDim; a_old++ ) {

						a_new = position + Dim_J_block[ n ][ i ] * ( a_old / old.Dim_J_block[ oldhm ][ oldJm ] ) + a_old % old.Dim_J_block[ oldhm ][ oldJm ];

						NN[ site ][ n ][ i ][ a_new ] += factor * old.NN[ site_old ][ oldhm ][ oldJm ][ a_old ];

					}

				}

			}

		}

		if( table_new[ site ] == TotSiteNo - 1 ) {

			for( int n = 0; n < Block_Num_hole; n++ )
			if( Num_hole_block[ n ] != TotSiteNo ) {	//hole number<=site number

				for( int i = 0; i < Num_J_hole_block[ n ]; i++ ) 
				for( int m = 1; m < 3; m++ )  // hole number is zero
				if( ( oldhm = space.Hole_blockOld[ n ][ m ][ i ] ) != -1  &&  ( oldJm = space.J_blockOld[ n ][ m ][ i ] ) != -1 ) {

					//----
					position = Dim_J_block[ n ][ i ] + 1;

					factor = sqrt( Value_J_block[ n ][ i ] + 1.0 );

					//----
					for( int a_old = 0; a_old < old.Dim_J_block[ oldhm ][ oldJm ]; a_old++ ) {

						a_new = position * ( space.Start[ n ][ m ][ i ] + a_old );

						NN[ site ][ n ][ i ][ a_new ] += factor;

					}

				}

			}

		}

	}
//******
//	for(int n=0; n<operator_number_J; n++)
//	for(int i=0; i<Block_Num_hole; i++)
//	for(int j=0; j<Num_J_hole_block[i]; j++)
//	for(int l=0; l<Dim_J_block[i][j]*Dim_J_block[i][j]; l++)
//		cout<<"\n NN["<<n<<"]["<<i<<"]["<<j<<"]["<<l<<"]="<<NN[n][i][j][l];
//******
}

//====================================NewH_Sys====================================
inline void Sub::New_H_Sys(Parameter &para, Sub &space, Sub &sys) {

	int oldhm, oldhl, oldJm, oldJl, position, SubDim, a_new, a_new_t, old_T, old_J;

	int SiteNum = para.Total_N*sys.TotSiteNo;	

	double factor;

  //----Create 6j for H
	double **** six_j_H;
	six_j_H=new double *** [Block_Num_hole];
	for(int n=0; n<Block_Num_hole; n++) {
		six_j_H[n]=new double ** [Num_J_hole_block[n]];
		for(int i=0; i<Num_J_hole_block[n]; i++) {
			six_j_H[n][i]=new double * [3];
			for(int m=0; m<3; m++) {
				six_j_H[n][i][m]=new double [3];
				for(int l=0; l<3; l++) {
					six_j_H[n][i][m][l]=0.0;
				}
			}
		}
	}

  //----Initialization
	for(int n=0; n<Block_Num_hole; n++)
	for(int i=0; i<Num_J_hole_block[n]; i++)
	for(int m=0; m<3; m++)
	for(int l=0; l<3; l++)
	if( (oldhm=space.Hole_blockOld[n][m][i])!=-1 && (oldhl=space.Hole_blockOld[n][l][i])!=-1 && oldhm==oldhl && sys.Num_hole_block[oldhm]==Num_hole_block[n]  &&  (oldJm=space.J_blockOld[n][m][i])!=-1 && (oldJl=space.J_blockOld[n][l][i])!=-1 ) 
		six_j_H[n][i][m][l]=gsl_sf_coupling_6j(Value_J_block[n][i], 1, sys.Value_J_block[oldhm][oldJm], 2, sys.Value_J_block[oldhl][oldJl], 1);

	//save
        FILE *fr=fopen(Combine(Combine("6j_factor/H/", 1), space.TotSiteNo), "wb");

	for(int n=0; n<Block_Num_hole; n++)
        for(int i=0; i<Num_J_hole_block[n]; i++)
	for(int m=0; m<3; m++) {

                fwrite(six_j_H[n][i][m], sizeof(double), 3, fr);

	}

        fclose(fr);

  //----NewH_Sys_Old
	for(int n=0; n<Block_Num_hole; n++)
	if(space.Num_hole_block[n] != space.TotSiteNo) {

		for(int i=0; i<Num_J_hole_block[n]; i++)
		for(int m=0; m<3; m++)
		if((oldhm=space.Hole_blockOld[n][m][i])!=-1 && (oldJm=space.J_blockOld[n][m][i])!=-1) {

			position=(Dim_J_block[n][i]+1)*space.Start[n][m][i];
			SubDim=sys.Dim_J_block[oldhm][oldJm]*sys.Dim_J_block[oldhm][oldJm];

			for(int a_old=0; a_old<SubDim; a_old++) {

				a_new=position+Dim_J_block[n][i]*(a_old/sys.Dim_J_block[oldhm][oldJm])+a_old%sys.Dim_J_block[oldhm][oldJm];
				H[n][i][a_new]+=sys.H[oldhm][oldJm][a_old];

			}

		}

	}

  //----NewH_Sys_Old_N: T couplings
	for( int site = 0; site < sys.operator_number_T; site ++ )
        if( para.Table_T[ SiteNum + para.Table_T_sys_site[ sys.TotSiteNo - 1 ][ site ] ] == 1 ) {

		for( int n = 0; n < Block_Num_hole; n ++ )
		for( int i = 0; i < Num_J_hole_block[n]; i ++ )
		for( int m = 0; m < 3; m ++ )
		for( int l = 0; l < 3; l ++ )
		if( m == 0  &&  l != 0  &&  ( oldhm = space.Hole_blockOld[ n ][ m ][ i ] ) != -1  &&  ( oldhl = space.Hole_blockOld[ n ][ l ][ i ] ) != -1  &&  ( sys.Num_hole_block[ oldhm ] + 1 ) == sys.Num_hole_block[ oldhl ]  &&  ( oldJm = space.J_blockOld[ n ][ m ][ i ] ) != -1  &&  ( oldJl = space.J_blockOld[ n ][ l ][ i ] ) != -1  &&  abs( sys.Value_J_block[ oldhm ][ oldJm ] - sys.Value_J_block[ oldhl ][ oldJl ] ) == 1 ) {

			//----
			position = Dim_J_block[ n ][ i ] * space.Start[ n ][ l ][ i ] + space.Start[ n ][ m ][ i ];

			SubDim = sys.Dim_J_block[ oldhm ][ oldJm ] * sys.Dim_J_block[ oldhl ][ oldJl ];

			factor = para.Interaction_T[ SiteNum + para.Table_T_sys_site[ sys.TotSiteNo - 1 ][ site ] ] / sqrt( Value_J_block[ n ][ i ] + 1.0 );

			//---- fermion phase factor
			if( ( sys.TotSiteNo - sys.Num_hole_block[ oldhl ] ) % 2 == 1 )
				factor = - factor; //sys.TotSite is the particle number at half-filling

			//----
			for( int old_i = 0; old_i < sys.Num_block_T[ oldhm ]; old_i ++ )
			if( sys.J_block_T_bra[ oldhm ][ old_i ] == oldJm  &&  sys.J_block_T_ket[ oldhm ][ old_i ] == oldJl)
				old_J = old_i;

			//----
			for( int a_old = 0; a_old < SubDim; a_old ++ ) {

				a_new = position + Dim_J_block[ n ][ i ] * ( a_old / sys.Dim_J_block[ oldhm ][ oldJm ] ) + a_old % sys.Dim_J_block[oldhm][oldJm];
				a_new_t = Dim_J_block[ n ][ i ] * ( a_new % Dim_J_block[ n ][ i ] ) + a_new / Dim_J_block[ n ][ i ];

				H[ n ][ i ][ a_new ] += factor * sys.T[ site ][ oldhm ][ old_J ][ a_old ];
				H[ n ][ i ][ a_new_t ] += factor * sys.T[ site ][ oldhm ][ old_J ][ a_old ];

			}

		}

	}

  //----NewH_Sys_Old_N: J couplings
	for( int site = 0; site < sys.operator_number_J; site++ )
        if( para.Table_J[ SiteNum + para.Table_J_sys_site[ sys.TotSiteNo - 1 ][ site ] ] == 1 ) {

		for(int n=0; n<Block_Num_hole; n++)
		for(int i=0; i<Num_J_hole_block[n]; i++)
		for(int m=0; m<3; m++)
		for(int l=0; l<3; l++)
		if( (oldhm = space.Hole_blockOld[n][m][i] ) != -1  &&  ( oldhl = space.Hole_blockOld[n][l][i] ) != -1  &&  oldhm == oldhl  &&  sys.Num_hole_block[oldhm] == Num_hole_block[n]  &&  ( oldJm = space.J_blockOld[n][m][i] ) != -1  &&  ( oldJl = space.J_blockOld[n][l][i] ) != -1 ) {

			if( oldJm == oldJl  &&  sys.Value_J_block[oldhm][oldJm] != 0 ) {

                                position = Dim_J_block[n][i] * space.Start[n][l][i] + space.Start[n][m][i];

                                factor = para.Interaction_J[SiteNum+para.Table_J_sys_site[sys.TotSiteNo-1][site]] * pow( -1.0, (Value_J_block[n][i]+sys.Value_J_block[oldhl][oldJl]+1) / 2 ) * constant * six_j_H[n][i][m][l];

                                SubDim = sys.Dim_J_block[oldhm][oldJm] * sys.Dim_J_block[oldhl][oldJl];

                                for(int a_sys=0; a_sys<SubDim; a_sys++) {

                                        a_new = position + Dim_J_block[n][i] * ( a_sys / sys.Dim_J_block[oldhm][oldJm] ) + a_sys % sys.Dim_J_block[oldhm][oldJm];

                                        H[n][i][a_new] += factor * sys.S_Dia[site][oldhm][oldJm][a_sys];

                                }

			}

			else if( (oldJm+1) == oldJl  &&  sys.Value_J_block[oldhl][oldJl] - sys.Value_J_block[oldhm][oldJm] == 2 ) {

				position = Dim_J_block[n][i] * space.Start[n][l][i] + space.Start[n][m][i];

				factor = para.Interaction_J[SiteNum+para.Table_J_sys_site[sys.TotSiteNo-1][site]] * pow( -1.0, ( Value_J_block[n][i]+sys.Value_J_block[oldhl][oldJl]+1 ) / 2 ) * constant * six_j_H[n][i][m][l];

                                SubDim = sys.Dim_J_block[oldhm][oldJm]*sys.Dim_J_block[oldhl][oldJl];

                                for(int a_sys=0; a_sys<SubDim; a_sys++) {

                                        a_new=position+Dim_J_block[n][i]*(a_sys/sys.Dim_J_block[oldhm][oldJm])+a_sys%sys.Dim_J_block[oldhm][oldJm];

					a_new_t=Dim_J_block[n][i]*(a_new%Dim_J_block[n][i])+a_new/Dim_J_block[n][i];

					H[n][i][a_new]+=factor*sys.S_M_Dia[site][oldhm][oldJm][a_sys];

					H[n][i][a_new_t]+=factor*sys.S_M_Dia[site][oldhm][oldJm][a_sys];

				}

			}

		}

	}

  //----NewH_Sys_Old_N: n couplings
	for( int site = 0; site < sys.operator_number_J; site++ )
        if( para.Table_N[ SiteNum + para.Table_N_sys_site[ sys.TotSiteNo - 1 ][ site ] ] == 1 ) {

		for( int n = 0; n < Block_Num_hole; n++ )
		if( Num_hole_block[ n ] != TotSiteNo ) {

			for( int i = 0; i < Num_J_hole_block[ n ]; i++ )
			for( int m = 1; m < 3; m++ )
			if( ( oldhm = space.Hole_blockOld[ n ][ m ][ i ] ) != -1  &&  ( oldJm = space.J_blockOld[ n ][ m ][ i ] ) != -1 ) {

				//----
				position = Dim_J_block[ n ][ i ] * space.Start[ n ][ m ][ i ] + space.Start[ n ][ m ][ i ];

	                        factor = para.Interaction_N[ SiteNum + para.Table_N_sys_site[ sys.TotSiteNo - 1 ][ site ] ] / sqrt( sys.Value_J_block[ oldhm ][ oldJm ] + 1.0 );

        	                SubDim = sys.Dim_J_block[ oldhm ][ oldJm ] * sys.Dim_J_block[ oldhm ][ oldJm ];
	
				//----
        	                for( int a_sys = 0; a_sys < SubDim; a_sys++ ) {

                	        	a_new = position + Dim_J_block[ n ][ i ] * ( a_sys / sys.Dim_J_block[ oldhm ][ oldJm ] ) + a_sys % sys.Dim_J_block[ oldhm ][ oldJm ];

                        	        H[ n ][ i ][ a_new ] += factor * sys.NN[ site ][ oldhm ][ oldJm ][ a_sys ];

                        	}

			}

		}

	}

  //----Delete six_j_H
	for(int n=0; n<Block_Num_hole; n++) {
		for(int i=0; i<Num_J_hole_block[n]; i++) {
			for(int j=0; j<3; j++) {
				delete [] six_j_H[n][i][j];
			}
		delete [] six_j_H[n][i];
		}
	delete [] six_j_H[n];
	}
	delete [] six_j_H;

//*******
//	for(int i=0; i<Block_Num_hole; i++)
//	for(int j=0; j<Num_J_hole_block[i]; j++)
//	for(int k=0; k<Dim_J_block[i][j]*Dim_J_block[i][j]; k++)
//		cout<<"\n H["<<i<<"]["<<j<<"]["<<k<<"]="<<H[i][j][k];
//*******Check the Hamiltonian elements
/*	for(int i=0; i<Block_Num_hole; i++)
	for(int j=0; j<Num_J_hole_block[i]; j++) {

		char jobz='N';       				//only eigenvalues are computed
                char uplo='U';       				//a stores the upper triangular part of A
                int n=Dim_J_block[i][j]; 		        //The order of the matrix A
                int lda=n;  				        //The first dimension of the array a
                int lwork=40*n;     //The dimension of the array work. Constraint:lwork>=max(1, 3n-1), check work(1)!
		int info;

		double * a=new double [lda*n];
                double * work=new double [lwork];
                double * w=new double [n];//contains the eigenvalues of the matrix A in ascending order!!!

		for(int m=0; m<lda*n; m++)
			a[m]=H[i][j][m];

                for(int m=0; m<lwork; m++) 
                        work[m]=0.0;

                for(int m=0; m<n; m++)
                        w[m]=0.0;     

        	dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

		for(int m=0; m<n; m++)
			cout<<"\n Eigenvalues["<<i<<"]["<<j<<"]["<<m<<"]="<<w[m];

		delete [] a;  delete [] work;  delete [] w;

	}
*/
}

//===================================NewH_Env====================================
inline void Sub::New_H_Env(Parameter &para, Sub &space, Sub &env) {

	int oldhm, oldhl, oldJm, oldJl, position, SubDim, a_new, a_new_t, old_T, old_J;
        int SiteNum=para.Total_N*(StartSite-env.TotSiteNo)+StartSite;
	double factor;

  //----Create 6j for H
	double **** six_j_H;
	six_j_H=new double *** [Block_Num_hole];
	for(int n=0; n<Block_Num_hole; n++) {
		six_j_H[n]=new double ** [Num_J_hole_block[n]];
		for(int i=0; i<Num_J_hole_block[n]; i++) {
			six_j_H[n][i]=new double * [3];
			for(int m=0; m<3; m++) {
				six_j_H[n][i][m]=new double [3];
				for(int l=0; l<3; l++) {
					six_j_H[n][i][m][l]=0.0;
				}
			}
		}
	}

  //----Initialization
	for(int n=0; n<Block_Num_hole; n++)
	for(int i=0; i<Num_J_hole_block[n]; i++)
	for(int m=0; m<3; m++)
	for(int l=0; l<3; l++)
	if( (oldhm=space.Hole_blockOld[n][m][i])!=-1 && (oldhl=space.Hole_blockOld[n][l][i])!=-1 && oldhm==oldhl && env.Num_hole_block[oldhm]==Num_hole_block[n]  &&  (oldJm=space.J_blockOld[n][m][i])!=-1 && (oldJl=space.J_blockOld[n][l][i])!=-1) 
		six_j_H[n][i][m][l]=gsl_sf_coupling_6j(Value_J_block[n][i], 1, env.Value_J_block[oldhm][oldJm], 2, env.Value_J_block[oldhl][oldJl], 1);

	//save
        FILE *fr=fopen(Combine(Combine("6j_factor/H/", 2), space.TotSiteNo), "wb");

	for(int n=0; n<Block_Num_hole; n++)
        for(int i=0; i<Num_J_hole_block[n]; i++)
	for(int m=0; m<3; m++) {

                fwrite(six_j_H[n][i][m], sizeof(double), 3, fr);

	}

        fclose(fr);

  //------NewH_Env_Old
	for(int n=0; n<Block_Num_hole; n++)
	if(space.Num_hole_block[n] != space.TotSiteNo) {

		for(int i=0; i<Num_J_hole_block[n]; i++)
		for(int m=0; m<3; m++)
		if((oldhm=space.Hole_blockOld[n][m][i])!=-1 && (oldJm=space.J_blockOld[n][m][i])!=-1) {

			position=(Dim_J_block[n][i]+1)*space.Start[n][m][i];
			SubDim=env.Dim_J_block[oldhm][oldJm]*env.Dim_J_block[oldhm][oldJm];

			for(int a_old=0; a_old<SubDim; a_old++) {

				a_new=position+Dim_J_block[n][i]*(a_old/env.Dim_J_block[oldhm][oldJm])+a_old%env.Dim_J_block[oldhm][oldJm];
				H[n][i][a_new]+=env.H[oldhm][oldJm][a_old];

			}

		}

	}

  //----NewH_Env_Old_N: T couplings
	for(int site=0; site<env.operator_number_T; site++)
        if(para.Table_T[SiteNum-para.Table_T_env_site[env.TotSiteNo-1][site]]==1) {

		for(int n=0; n<Block_Num_hole; n++)
		for(int i=0; i<Num_J_hole_block[n]; i++)
		for(int m=0; m<3; m++)
		for(int l=0; l<3; l++)
		if(m==0 && l!=0 && (oldhm=space.Hole_blockOld[n][m][i])!=-1 && (oldhl=space.Hole_blockOld[n][l][i])!=-1 && (env.Num_hole_block[oldhm]+1)==env.Num_hole_block[oldhl] && (oldJm=space.J_blockOld[n][m][i])!=-1 && (oldJl=space.J_blockOld[n][l][i])!=-1 && abs(env.Value_J_block[oldhm][oldJm]-env.Value_J_block[oldhl][oldJl])==1) {

			//----
			position=Dim_J_block[n][i]*space.Start[n][l][i]+space.Start[n][m][i];

			SubDim=env.Dim_J_block[oldhm][oldJm]*env.Dim_J_block[oldhl][oldJl];

			factor=para.Interaction_T[SiteNum-para.Table_T_env_site[env.TotSiteNo-1][site]]/sqrt(Value_J_block[n][i]+1.0);

			if((env.TotSiteNo-env.Num_hole_block[oldhl])%2==1)
				factor = -factor;	//fermion phase

			for(int old_n=0; old_n < env.Block_Num_hole - 1; old_n++)
			for(int old_i=0; old_i<env.Num_block_T[old_n]; old_i++)
			if(old_n == oldhm  &&  env.J_block_T_bra[old_n][old_i]==oldJm && env.J_block_T_ket[old_n][old_i]==oldJl) {

				old_T=old_n;	old_J=old_i;

			}

			for(int a_old=0; a_old<SubDim; a_old++) {

				a_new=position+Dim_J_block[n][i]*(a_old/env.Dim_J_block[oldhm][oldJm])+a_old%env.Dim_J_block[oldhm][oldJm];

				a_new_t=Dim_J_block[n][i]*(a_new%Dim_J_block[n][i])+a_new/Dim_J_block[n][i];

				H[n][i][a_new]+=factor*env.T[site][old_T][old_J][a_old];

				H[n][i][a_new_t]+=factor*env.T[site][old_T][old_J][a_old];

			}

		}

	}

  //----NewH_Env_Old_N: J couplings
	for(int site=0; site<env.operator_number_J; site++)
        if(para.Table_J[SiteNum-para.Table_J_env_site[env.TotSiteNo-1][site]]==1) {

		for(int n=0; n<Block_Num_hole; n++)
		for(int i=0; i<Num_J_hole_block[n]; i++)
		for(int m=0; m<3; m++)
		for(int l=0; l<3; l++)
		if((oldhm=space.Hole_blockOld[n][m][i])!=-1 && (oldhl=space.Hole_blockOld[n][l][i])!=-1 && oldhm==oldhl && env.Num_hole_block[oldhm]==Num_hole_block[n]  &&  (oldJm=space.J_blockOld[n][m][i])!=-1 && (oldJl=space.J_blockOld[n][l][i])!=-1) {

			if(oldJm==oldJl && env.Value_J_block[oldhm][oldJm]!=0) {

                                position=Dim_J_block[n][i]*space.Start[n][l][i]+space.Start[n][m][i];

                                factor=para.Interaction_J[SiteNum-para.Table_J_env_site[env.TotSiteNo-1][site]]*pow(-1.0,(Value_J_block[n][i]+env.Value_J_block[oldhl][oldJl]+1)/2)*constant*six_j_H[n][i][m][l];

                                SubDim=env.Dim_J_block[oldhm][oldJm]*env.Dim_J_block[oldhl][oldJl];

                                for(int a_env=0; a_env<SubDim; a_env++) {

                                        a_new=position+Dim_J_block[n][i]*(a_env/env.Dim_J_block[oldhm][oldJm])+a_env%env.Dim_J_block[oldhm][oldJm];

                                        H[n][i][a_new]+=factor*env.S_Dia[site][oldhm][oldJm][a_env];

                                }

			}

			else if((oldJm+1)==oldJl && env.Value_J_block[oldhl][oldJl]-env.Value_J_block[oldhm][oldJm]==2) {

				position=Dim_J_block[n][i]*space.Start[n][l][i]+space.Start[n][m][i];

				factor=para.Interaction_J[SiteNum-para.Table_J_env_site[env.TotSiteNo-1][site]]*pow(-1.0,(Value_J_block[n][i]+env.Value_J_block[oldhl][oldJl]+1)/2)*constant*six_j_H[n][i][m][l];

                                SubDim=env.Dim_J_block[oldhm][oldJm]*env.Dim_J_block[oldhl][oldJl];

                                for(int a_env=0; a_env<SubDim; a_env++) {

                                        a_new=position+Dim_J_block[n][i]*(a_env/env.Dim_J_block[oldhm][oldJm])+a_env%env.Dim_J_block[oldhm][oldJm];

					a_new_t=Dim_J_block[n][i]*(a_new%Dim_J_block[n][i])+a_new/Dim_J_block[n][i];

					H[n][i][a_new]+=factor*env.S_M_Dia[site][oldhm][oldJm][a_env];

					H[n][i][a_new_t]+=factor*env.S_M_Dia[site][oldhm][oldJm][a_env];

				}

			}

		}

	}

  //----NewH_Env_Old_N: n couplings
	for(int site=0; site<env.operator_number_J; site++)
        if(para.Table_N[SiteNum-para.Table_N_env_site[env.TotSiteNo-1][site]]==1) {

		for(int n=0; n<Block_Num_hole; n++)
		if( Num_hole_block[ n ] != TotSiteNo ) {

			for(int i=0; i<Num_J_hole_block[n]; i++)
			for(int m=0; m<3; m++)
			if(m!=0 && (oldhm=space.Hole_blockOld[n][m][i])!=-1 && (oldJm=space.J_blockOld[n][m][i])!=-1) {

				position=Dim_J_block[n][i]*space.Start[n][m][i]+space.Start[n][m][i];

	                        factor=para.Interaction_N[SiteNum-para.Table_N_env_site[env.TotSiteNo-1][site]]/sqrt(env.Value_J_block[oldhm][oldJm]+1.0);

        	                SubDim=env.Dim_J_block[oldhm][oldJm]*env.Dim_J_block[oldhm][oldJm];
	
        	                for(int a_env=0; a_env<SubDim; a_env++) {

                	        	a_new=position+Dim_J_block[n][i]*(a_env/env.Dim_J_block[oldhm][oldJm])+a_env%env.Dim_J_block[oldhm][oldJm];

                        	        H[n][i][a_new]+=factor*env.NN[site][oldhm][oldJm][a_env];

                        	}

			}

		}

	}

  //----Delete six_j_H
	for(int n=0; n<Block_Num_hole; n++) {
		for(int i=0; i<Num_J_hole_block[n]; i++) {
			for(int j=0; j<3; j++) {
				delete [] six_j_H[n][i][j];
			}
		delete [] six_j_H[n][i];
		}
	delete [] six_j_H[n];
	}
	delete [] six_j_H;

  //------
/*	for(int i=0; i<Block_Num_hole; i++)
	for(int j=0; j<Num_J_hole_block[i]; j++) {

		char jobz='N';       				//only eigenvalues are computed
                char uplo='U';       				//a stores the upper triangular part of A
                int n=Dim_J_block[i][j]; 		        //The order of the matrix A
                int lda=n;  				        //The first dimension of the array a
                int lwork=40*n;     //The dimension of the array work. Constraint:lwork>=max(1, 3n-1), check work(1)!
		int info;

		double * a=new double [lda*n];
                double * work=new double [lwork];
                double * w=new double [n];//contains the eigenvalues of the matrix A in ascending order!!!

		for(int m=0; m<lda*n; m++)
			a[m]=H[i][j][m];

                for(int m=0; m<lwork; m++) 
                        work[m]=0.0;

                for(int m=0; m<n; m++)
                        w[m]=0.0;     

        	dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

		for(int m=0; m<n; m++)
			cout<<"\n Eigenvalues["<<i<<"]["<<j<<"]["<<m<<"]="<<w[m];

		delete [] a;  delete [] work;  delete [] w;

	}
*/
}

//===============================================================================
//                 Truncate and renew the operators in infinite sweep
//===============================================================================
Sub::Sub(Parameter &para, Sub &density):SubCommon() {

        IndexNo=5;      //only CreateSpace1();

        TotSiteNo=density.TotSiteNo;

        Truncate(para, density, para.StateNoKept);

	operator_number_J=density.operator_number_J;
	operator_number_T=density.operator_number_T;

        if(mmin<para.StateNoKept)  para.untruncated_site++;//count the number of untruncated sites!!!
        cout<<"\n para.untruncated_site="<<para.untruncated_site<<endl;

}

//===============================================================================================================
//===============================================================================================================
//                              Truncate and renew operators in finite sweep
//===============================================================================================================
Sub::Sub(Parameter &para, Sub &density, const int &optimalstate):SubCommon() {

        IndexNo=5;

        TotSiteNo=density.TotSiteNo;

        Truncate(para, density, optimalstate);

	operator_number_J=density.operator_number_J;
	operator_number_T=density.operator_number_T;

}

//===============================================================================================================
Sub::Sub(Parameter &para, char &sign, Sub &density, Sub &trun_space, Sub &old):SubCommon() {

	IndexNo=2;
	ope_sign=sign;

        trans_N='N';    trans_T='T';
        alpha=1.0;      beta=0.0;

  //----Create spaces for variables, read from trun_space
        TotSiteNo=trun_space.TotSiteNo;
        operator_number_J=trun_space.operator_number_J;
        operator_number_T=trun_space.operator_number_T;

	Block_Num_hole=trun_space.Block_Num_hole;

        CreateSpace_hole_block();       //hole number and J number in each block
	Old_hole=new int [Block_Num_hole];
	for(int i=0; i<Block_Num_hole; i++) {
		Num_hole_block[i]=trun_space.Num_hole_block[i];  Num_J_hole_block[i]=trun_space.Num_J_hole_block[i];
		Old_hole[i]=trun_space.Old_hole[i];
	}

	CreateSpace_J_block();		//J value and dimension in each J block
	for(int i=0; i<Block_Num_hole; i++)
	for(int j=0; j<Num_J_hole_block[i]; j++) {
		Value_J_block[i][j]=trun_space.Value_J_block[i][j]; Dim_J_block[i][j]=trun_space.Dim_J_block[i][j];
	}

	Old_J=new int * [Block_Num_hole];
	for(int i=0; i<Block_Num_hole; i++) {
		Old_J[i]=new int [Num_J_hole_block[i]];
		for(int j=0; j<Num_J_hole_block[i]; j++)
			Old_J[i][j]=trun_space.Old_J[i][j];
	}

	CreateSpace_T();

        if(ope_sign=='H') {
                Create_H();
                Truncate_H(density, old);
        }

	else if(ope_sign=='T') {
                Create_T();
		Truncate_T(density, old);
        }

	else if(ope_sign=='S') {
                Create_S_Dia();
		Truncate_S_Dia(density, old);
        }

	else if(ope_sign=='M') {
                Create_S_M_Dia();
		Truncate_S_M_Dia(density, old);
        }

	else if(ope_sign=='N') {
                Create_NN();
		Truncate_NN(density, old);
        }

	delete [] Old_hole;
	for(int i=0; i<Block_Num_hole; i++)
		delete [] Old_J[i];
	delete [] Old_J;

}

//===============================================================================================================
//                              Truncate:Relocate the good quantum numbers
//===============================================================================================================
inline void Sub::Truncate(Parameter &para, Sub &old, const int &optimalstate) {

  //----find the m largest eigenvalues of reduced density matrix
        int m=0, TotStateNo=0;

        for(int i=0; i<old.Block_Num_hole; i++) 
	for(int j=0; j<old.Num_J_hole_block[i]; j++)
                TotStateNo+=old.Dim_J_block[i][j];    //total number of the untruncated basis

	int *N_h=new int [TotStateNo];	//store hole number for each basis
	int *N_J=new int [TotStateNo];  //store J number for each basis
        double *N_e=new double [TotStateNo];	//store reduced density matrix eigenvalue of each basis

        for(int i=0; i<old.Block_Num_hole; i++) 
	for(int j=0; j<old.Num_J_hole_block[i]; j++) {

		for(int k=0; k<old.Dim_J_block[i][j]; k++) {

                       	N_h[m]=i;  N_J[m]=j;  N_e[m++]=old.dm_eig[i][j][k]; 

                }

        }

        for(int i=0; i<TotStateNo; i++)
        for(int j=i+1; j<TotStateNo; j++) {

                if(N_e[i]<N_e[j]) {

                        int n=N_h[i];  N_h[i]=N_h[j];  N_h[j]=n;  n=N_J[i];  N_J[i]=N_J[j];  N_J[j]=n;
                        double e=N_e[i];  N_e[i]=N_e[j];  N_e[j]=e;

                }

        }   // rearrange from big to small

  //----calculate truncation error
        mmin=(m>optimalstate ? optimalstate : m );

	int *Dim1=new int [old.Block_Num_hole];	//store basis number for each hole number
	int **Dim2=new int * [old.Block_Num_hole];//store basis number for each hole and J numbers
	for(int i=0; i<old.Block_Num_hole; i++) {
		Dim1[i]=0;
		Dim2[i]=new int [old.Num_J_hole_block[i]];
		for(int j=0; j<old.Num_J_hole_block[i]; j++)
			Dim2[i][j]=0;
	}

	int kept_state=0;
        double trunc_error=1.0;

        for(int i=0; i<mmin; i++) {

		++(Dim1[N_h[i]]);
                ++(Dim2[N_h[i]][N_J[i]]);       //find the numbers of the new basis with certain h and J number
		kept_state+=old.Value_J_block[N_h[i]][N_J[i]]+1;
                trunc_error-=N_e[i]*(old.Value_J_block[N_h[i]][N_J[i]]+1.0);

        }

        cout<<"\n U(1) equivalent state="<<kept_state<<"\t Truncation Error="<<trunc_error<<endl;

  //----calculate truncated Block_Num_hole
	Block_Num_hole=0;
	for(int i=0; i<old.Block_Num_hole; i++)
	if(Dim1[i]!=0)
		Block_Num_hole++;

  //----calculate truncated Num_hole_block  and Num_J_hole_block
	CreateSpace_hole_block();       //hole number and J number in each block

	Old_hole=new int [Block_Num_hole];

	m=0;
	for(int i=0; i<old.Block_Num_hole; i++) 
	if(Dim1[i]!=0) {

		Num_hole_block[m]=old.Num_hole_block[i];
		Old_hole[m]=i;

		for(int k=0; k<old.Num_J_hole_block[i]; k++)
		if(Dim2[i][k]!=0)
			Num_J_hole_block[m]++;

		m++;

	}

  //----calculate truncated Value_J_block and Dim_J_block
	CreateSpace_J_block();		//J value and dimension in each J block

	Old_J=new int * [Block_Num_hole];
	for(int i=0; i<Block_Num_hole; i++) {
		Old_J[i]=new int [Num_J_hole_block[i]];
		for(int j=0; j<Num_J_hole_block[i]; j++)
			Old_J[i][j]=0;
	}

	for(int n=0; n<Block_Num_hole; n++) {

		m=0;
		for(int i=0; i<old.Num_J_hole_block[Old_hole[n]]; i++)
		if(Dim2[Old_hole[n]][i]!=0) {
			Value_J_block[n][m]=old.Value_J_block[Old_hole[n]][i];
			Dim_J_block[n][m]=Dim2[Old_hole[n]][i];
			Old_J[n][m++]=i;
		}

	}

  //----Create space for T
        CreateSpace_T();

  //----delete memory
	for(int i=0; i<old.Block_Num_hole; i++)
		delete [] Dim2[i];
	delete [] Dim2;

	delete [] Dim1;

	delete [] N_h;  delete [] N_J;  delete [] N_e;

}

//================================================================================================================
//                                      Truncate and renew the T
//================================================================================================================
inline void Sub::Truncate_T(Sub &density, Sub &old) {
//-------------------------------Matrix-Matrix multiplication by BLAS routine dgemm------------------------------
	int old_T, old_N;

        for(int site=0; site<operator_number_T; site++) {

		for(int n=0; n < Block_Num_hole - 1; n++)
                for(int i=0; i<Num_block_T[n]; i++) {

			for(int n_old=0; n_old < old.Block_Num_hole - 1; n_old++)
			for(int i_old=0; i_old<old.Num_block_T[n_old]; i_old++)
			if(Num_hole_block[n]==old.Num_hole_block[n_old] && Value_J_block[n][J_block_T_bra[n][i]]==old.Value_J_block[n_old][old.J_block_T_bra[n_old][i_old]] && Value_J_block[n+1][J_block_T_ket[n][i]]==old.Value_J_block[n_old+1][old.J_block_T_ket[n_old][i_old]]) {
				old_T=n_old;	old_N=i_old;
			}

                        SubDim=old.Dim_J_block[old_T][old.J_block_T_bra[old_T][old_N]]*Dim_J_block[n+1][J_block_T_ket[n][i]];
 
                        double * f_sub=new double [SubDim];
                        for(int j=0; j<SubDim; j++)     f_sub[j]=(double) 0;
 
                        dgemm_(&trans_N, &trans_N, &old.Dim_J_block[old_T][old.J_block_T_bra[old_T][old_N]], &Dim_J_block[n+1][J_block_T_ket[n][i]], &old.Dim_J_block[old_T+1][old.J_block_T_ket[old_T][old_N]], &alpha, old.T[site][old_T][old_N], &old.Dim_J_block[old_T][old.J_block_T_bra[old_T][old_N]], density.dm_wave[old_T+1][old.J_block_T_ket[old_T][old_N]], &old.Dim_J_block[old_T+1][old.J_block_T_ket[old_T][old_N]], &beta, f_sub, &old.Dim_J_block[old_T][old.J_block_T_bra[old_T][old_N]]);
 
                        dgemm_(&trans_T, &trans_N, &Dim_J_block[n][J_block_T_bra[n][i]], &Dim_J_block[n+1][J_block_T_ket[n][i]], &old.Dim_J_block[old_T][old.J_block_T_bra[old_T][old_N]], &alpha, density.dm_wave[old_T][old.J_block_T_bra[old_T][old_N]], &old.Dim_J_block[old_T][old.J_block_T_bra[old_T][old_N]], f_sub, &old.Dim_J_block[old_T][old.J_block_T_bra[old_T][old_N]], &beta, T[site][n][i], &Dim_J_block[n][J_block_T_bra[n][i]]);
 
                        delete [] f_sub;
                }

        }

}

//================================================================================================================
//                                              Truncate and renew H
//================================================================================================================
inline void Sub::Truncate_H(Sub &density, Sub &old) {
//------------------------------------Matrix-Matrix multiplication by BLAS routine dgemm-------------------------
	for(int n=0; n<Block_Num_hole; n++)
        for(int i=0; i<Num_J_hole_block[n]; i++) {

                SubDim=old.Dim_J_block[Old_hole[n]][Old_J[n][i]]*Dim_J_block[n][i];//Dim of middle matrix

                double * f_sub=new double [SubDim];
                for(int j=0; j<SubDim; j++)     f_sub[j]=(double) 0;

                dgemm_(&trans_N, &trans_N, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &Dim_J_block[n][i], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &alpha, old.H[Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], density.dm_wave[Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &beta, f_sub, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]]);

                dgemm_(&trans_T, &trans_N, &Dim_J_block[n][i], &Dim_J_block[n][i], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &alpha, density.dm_wave[Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], f_sub, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &beta, H[n][i], &Dim_J_block[n][i]);

                delete [] f_sub;

        }
}

//================================================================================================================
//                                      Truncate and renew the S_Dia
//================================================================================================================
inline void Sub::Truncate_S_Dia(Sub &density, Sub &old) {
//-------------------------------Matrix-Matrix multiplication by BLAS routine dgemm------------------------------
        for(int site=0; site<operator_number_J; site++) {

		for(int n=0; n<Block_Num_hole; n++)
                for(int i=0; i<Num_J_hole_block[n]; i++) {

                        SubDim=old.Dim_J_block[Old_hole[n]][Old_J[n][i]]*Dim_J_block[n][i];
 
                        double * f_sub=new double [SubDim];
                        for(int j=0; j<SubDim; j++)     f_sub[j]=(double) 0;
 
                        dgemm_(&trans_N, &trans_N, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &Dim_J_block[n][i], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &alpha, old.S_Dia[site][Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], density.dm_wave[Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &beta, f_sub, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]]);
 
                        dgemm_(&trans_T, &trans_N, &Dim_J_block[n][i], &Dim_J_block[n][i], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &alpha, density.dm_wave[Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], f_sub, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &beta, S_Dia[site][n][i], &Dim_J_block[n][i]);
 
                        delete [] f_sub;
                }

        }

}

//================================================================================================================
//                                      Truncate and renew the S_M_Dia
//================================================================================================================
inline void Sub::Truncate_S_M_Dia(Sub &density, Sub &old) {
//-----------------------------Matrix-Matrix multiplication by BLAS routine dgemm--------------------------------
        for(int site=0; site<operator_number_J; site++) {

		for(int n=0; n<Block_Num_hole; n++)
                for(int i=0; i<Num_J_hole_block[n]-1; i++)
                if((Value_J_block[n][i+1]-Value_J_block[n][i])==2 && (old.Value_J_block[Old_hole[n]][Old_J[n][i+1]]-old.Value_J_block[Old_hole[n]][Old_J[n][i]])==2) {

                        SubDim=old.Dim_J_block[Old_hole[n]][Old_J[n][i]]*Dim_J_block[n][i+1];

                        double * f_sub=new double [SubDim];
                        for(int j=0; j<SubDim; j++)     f_sub[j]=(double) 0;

                        dgemm_(&trans_N, &trans_N, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &Dim_J_block[n][i+1], &old.Dim_J_block[Old_hole[n]][Old_J[n][i+1]], &alpha, old.S_M_Dia[site][Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], density.dm_wave[Old_hole[n]][Old_J[n][i+1]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i+1]], &beta, f_sub, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]]);

                        dgemm_(&trans_T, &trans_N, &Dim_J_block[n][i], &Dim_J_block[n][i+1], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &alpha, density.dm_wave[Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], f_sub, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &beta, S_M_Dia[site][n][i], &Dim_J_block[n][i]);

                        delete [] f_sub;

                }

        }

}

//================================================================================================================
//                                      Truncate and renew the NN
//================================================================================================================
inline void Sub::Truncate_NN(Sub &density, Sub &old) {
//-------------------------------Matrix-Matrix multiplication by BLAS routine dgemm------------------------------
        for(int site=0; site<operator_number_J; site++) {

                for(int n=0; n<Block_Num_hole; n++)
                for(int i=0; i<Num_J_hole_block[n]; i++) {

                        SubDim=old.Dim_J_block[Old_hole[n]][Old_J[n][i]]*Dim_J_block[n][i];

                        double * f_sub=new double [SubDim];
                        for(int j=0; j<SubDim; j++)     f_sub[j]=(double) 0;

                        dgemm_(&trans_N, &trans_N, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &Dim_J_block[n][i], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &alpha, old.NN[site][Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], density.dm_wave[Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &beta, f_sub, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]]);

                        dgemm_(&trans_T, &trans_N, &Dim_J_block[n][i], &Dim_J_block[n][i], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &alpha, density.dm_wave[Old_hole[n]][Old_J[n][i]], &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], f_sub, &old.Dim_J_block[Old_hole[n]][Old_J[n][i]], &beta, NN[site][n][i], &Dim_J_block[n][i]);

                        delete [] f_sub;
                }

        }

}

//==========================Delete the Sub class: the delete program is in the class SubCommon===================
Sub::~Sub() {}
//===============================================END=============================================================
