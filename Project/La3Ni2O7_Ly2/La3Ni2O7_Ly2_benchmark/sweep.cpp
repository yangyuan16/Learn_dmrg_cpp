#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
using namespace std;

#include"sub.h"
#include"super.h" 
#include"superenergy.h"
#include"denmat.h"
#include"electrondensity.h"
#include"sqhopping.h"
#include"sqspin.h"
#include"sqdensity.h"
#include"spincorrelation.h"
#include"greenfunction.h"
#include"pairingcorrelation.h"
#include"pairingcorrelation_yy.h"
#include"densitycorrelation.h"
#include"common.h"
#include"sweep.h"

//==============================DMRG process of the 2D t-J lattices==============================
Sweep::Sweep(Parameter &para) {

	if( para.Sign == 'N' ) {

		lsys = 1;         //System contains one site in the beginning
		lns = 1;          //System is added a spin site in each iteraction process
	
		Initial( para );	

		for( ; lsys < para.Total_N/2; )

			Infinite( para, lsys, lns );
//			Infinite_test( para, lsys, lns );


		//------First sweep from Left to Right
	        step = 0;
		num = 0;
		number = 0;
	        initialstate = para.StateNoKept;  //number of kept states for the first sweep from left to right
	        site_start = 1;  //para.untruncated_site/2-1;//left end site of the sweep region
        	site_end = para.Total_N - 2;  //para.untruncated_site/2;//right end site of the sweep region
		sign = 'R';  //from left to right
		precision = 1.e-4;

		Finite_Sweep( para, para.Total_N/2, site_end, sign, initialstate, step, precision );//site_end

                FILE *g = fopen( "control_file/max_num", "r+" );
                fscanf( g, "%d\n", &max_num );
                fclose( g );

	        for( step = 1; step < max_num; step++ ) {

	                if( step % 2 == 1 )      {

        	                sign = 'L';
 
	                	FILE *g = fopen( "control_file/max_num", "r+" );
        	    	    	fscanf( g, "%d\n", &max_num );
                		fclose( g );

                                FILE *f = fopen( "control_file/optimalstate", "r+" );
                                fseek( f, number, SEEK_SET );
                                fscanf( f, "%d\n", &optimalstate_sweep );
                                if( optimalstate_sweep<10 )  number = number + 2;
                                else if( 9 < optimalstate_sweep && optimalstate_sweep < 100 )  number = number + 3;
                                else if( 99 < optimalstate_sweep && optimalstate_sweep < 1000 )  number = number + 4;
                                else if( 999 < optimalstate_sweep && optimalstate_sweep < 10000 )  number = number + 5;
                                fclose(f);

                                FILE *v = fopen( "control_file/tolerance", "r+" );
                                fseek( v, num, SEEK_SET );
                                fscanf( v, "%lf\n", &precision );
                                num = num + 7;
                                fclose(v);

                	        Finite_Sweep( para, site_start, site_end, sign, optimalstate_sweep, step, precision );

                	}

	                else if( step % 2 == 0 ) {

	                        sign = 'R';

 	                	FILE *g = fopen( "control_file/max_num", "r+" );
        	    	    	fscanf( g, "%d\n", &max_num );
                		fclose( g );

                                FILE *f = fopen( "control_file/optimalstate", "r+" );
                                fseek( f, number, SEEK_SET );
                                fscanf( f, "%d\n", &optimalstate_sweep );
                                if( optimalstate_sweep<10 )  number = number + 2;
                                else if( 9 < optimalstate_sweep && optimalstate_sweep < 100 )  number = number + 3;
                                else if( 99 < optimalstate_sweep && optimalstate_sweep < 1000 )  number = number + 4;
                                else if( 999 < optimalstate_sweep && optimalstate_sweep < 10000 )  number = number + 5;
                                fclose(f);

                                FILE *v = fopen( "control_file/tolerance", "r+" );
                                fseek( v, num, SEEK_SET );
                                fscanf( v, "%lf\n", &precision );
                                num = num + 7;
                                fclose(v);

	                        Finite_Sweep( para, site_start, site_end, sign, optimalstate_sweep, step, precision );

        	        }
	
		}
	
	}

	else if( para.Sign == 'R' ) {

		lsys = 1;         //System contains one site in the beginning
		lns = 1;          //System is added a spin site in each iteraction process

		//------First sweep
	        step = 24;	//change this!!!
		num = 0;
		number = 0;
	        initialstate = 10000;  //change this !!! number of kept states for the first sweep from left to right
	        site_start = 1;  //para.untruncated_site/2-1;//left end site of the sweep region
        	site_end = para.Total_N - 2;  //para.untruncated_site/2;//right end site of the sweep region
		sign = 'R';  //change this !!!   from left to right
		precision = 1.e-6;  //change this !!! Like Lanczos precision

		Finite_Sweep( para, site_start, site_end, sign, initialstate, step, precision ); //change 'site_start' to the continue site number !!!

                FILE *g = fopen( "control_file/max_num", "r+" );
                fscanf( g, "%d\n", &max_num );
                fclose( g );

	        for( step = 25; step < max_num; step++ ) {  //change the initial "step" number !!!!!!

	                if( step%2 == 1 )      {

        	                sign = 'L';
 
	                	FILE *g = fopen( "control_file/max_num", "r+" );
        	    	    	fscanf( g, "%d\n", &max_num );
                		fclose( g );

                                FILE *f = fopen( "control_file/optimalstate", "r+" );
                                fseek( f, number, SEEK_SET );
                                fscanf( f, "%d\n", &optimalstate_sweep ); cout<<"\n M="<<optimalstate_sweep;

                                if( optimalstate_sweep<10 )  number = number + 2;
                                else if( 9 < optimalstate_sweep && optimalstate_sweep < 100 )  number = number + 3;
                                else if( 99 < optimalstate_sweep && optimalstate_sweep < 1000 )  number = number + 4;
                                else if( 999 < optimalstate_sweep && optimalstate_sweep < 10000 )  number = number + 5;
                                fclose(f);

                                FILE *v = fopen( "control_file/tolerance", "r+" );
                                fseek( v, num, SEEK_SET );
                                fscanf( v, "%lf\n", &precision ); cout<<"\n Precision="<<precision;
                                num = num + 7;
                                fclose(v);

                	        Finite_Sweep( para, site_start, site_end, sign, optimalstate_sweep, step, precision );

                	}

	                else if( step%2 == 0 ) {

	                        sign = 'R';

 	                	FILE *g = fopen( "control_file/max_num", "r+" );
        	    	    	fscanf( g, "%d\n", &max_num );
                		fclose( g );

                                FILE *f = fopen( "control_file/optimalstate", "r+" );
                                fseek( f, number, SEEK_SET );
                                fscanf( f, "%d\n", &optimalstate_sweep );
                                if( optimalstate_sweep<10 )  number = number + 2;
                                else if( 9 < optimalstate_sweep && optimalstate_sweep < 100 )  number = number + 3;
                                else if( 99 < optimalstate_sweep && optimalstate_sweep < 1000 )  number = number + 4;
                                else if( 999 < optimalstate_sweep && optimalstate_sweep < 10000 )  number = number + 5;
                                fclose(f);

                                FILE *v = fopen( "control_file/tolerance", "r+" );
                                fseek( v, num, SEEK_SET );
                                fscanf( v, "%lf\n", &precision );
                                num = num + 7;
                                fclose(v);

	                        Finite_Sweep( para, site_start, site_end, sign, optimalstate_sweep, step, precision );

        	        }

		}

	}

//	else if( para.Sign == 'M' ) {

		step = 0;
		precision = 1.e-6;

                if( step % 2 == 1 )  sign = 'L';
                else if( step % 2 == 0 )  sign = 'R';

                Measurement( para, step, precision );

//	}

}

//===========================================Initialize sys and env blocks=======================================
inline void Sweep::Initial(Parameter &para) {

        Sub *sys = new Sub(para);
	
  //----Print space----------------------
	sys -> Print_space(1, lsys);
	sys -> Print_space(2, lsys);

  //----Print operators------------------
	sys -> Print_S_Dia(1, lsys);
	sys -> Print_S_M_Dia(1, lsys);
	sys -> Print_NN(1, lsys);
	sys -> Print_T(1, lsys);
	sys -> Print_H(1, lsys);

	sys -> Print_S_Dia(2, lsys);
	sys -> Print_S_M_Dia(2, lsys);
	sys -> Print_NN(2, lsys);
	sys -> Print_T(2, lsys);
	sys -> Print_H(2, lsys);

	delete sys;

}

//==========================Infinite Iteration process========================
inline void Sweep::Infinite( Parameter &para, int &lsys, const int &lns ) {

	cout<<"\n System Length="<<lsys<<endl;

//------Read in the data from the hard disk-----------------------------------
    //--Read in sys_space and env_space 
        Sub *sys_space = new Sub( 1, lsys );
	Sub *env_space = new Sub( 2, lsys );

    //--Read in sys and env operators
	Sub *sys_operator = new Sub( 1, lsys, *sys_space );
	Sub *env_operator = new Sub( 2, lsys, *env_space );

//------Space of blocks of sysnew and envnew----------------------------------
        lsys += lns;

        Sub *sysnew_space = new Sub( 1, para, *sys_space );
	Sub *envnew_space = new Sub( 2, para, *env_space );

//	sysnew_space -> Print_space(1, lsys);
//	envnew_space -> Print_space(2, lsys);

//------Diagonalization of the Hamiltonian-------------------------------------
    //--Allocate space for wave function blocks
	Super *sup_space = new Super(para, sys_space, env_space, sysnew_space, envnew_space);

    //--Allocate space for diagonalization
	sign = 'I';
	Super *sup_diagonalization = new Super(sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space);

    //--Diagonalization  of the Hamiltonian by Matrix-Vector multiplication
	para.total_site++;
	precision = 1.e-4;
	SuperEnergy *supeng = new SuperEnergy(para, *sup_space, *sup_diagonalization, precision);

        ground_state_energy = supeng->eigenvalue;
        excited_state_energy = 0.0;
        WriteEnergy(para, ground_state_energy, excited_state_energy, lsys-1, para.total_site);

	delete supeng;
	delete sup_diagonalization;

//------New operator for sysnew-----------------------------------------------------------------
    //--Create space for eigensystem of reduced density matrix
	Sub *reduced_density_sys = new Sub(*sysnew_space);

    //--Find the eigensystem of reduced density matrix
	DenMat *densys = new DenMat(1, para, *sup_space, *reduced_density_sys);
	delete densys;

    //--Create space for truncation
	Sub *systrun_space = new Sub(para, *reduced_density_sys);
	systrun_space -> Print_space(1, lsys);

    //--New and truncate T operator
	sign = 'T';
	Sub *sysnew_operator_T = new Sub(sign, para, *sys_operator, *sysnew_space);
	Sub *systrun_operator_T = new Sub(para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_T);
//	sysnew_operator_T -> Print_T(1, lsys);
	systrun_operator_T -> Print_T(1, lsys);
	delete systrun_operator_T;
	delete sysnew_operator_T;

	sign = 'S';
	Sub *sysnew_operator_S = new Sub(sign, para, *sys_operator, *sysnew_space);
	Sub *systrun_operator_S = new Sub(para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_S);
//	sysnew_operator_S -> Print_S_Dia(1, lsys);
	systrun_operator_S -> Print_S_Dia(1, lsys);
	delete systrun_operator_S;
	delete sysnew_operator_S;

	sign = 'M';
	Sub *sysnew_operator_M = new Sub(sign, para, *sys_operator, *sysnew_space);
	Sub *systrun_operator_M = new Sub(para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_M);
//	sysnew_operator_M -> Print_S_M_Dia(1, lsys);
	systrun_operator_M -> Print_S_M_Dia(1, lsys);
	delete systrun_operator_M;
	delete sysnew_operator_M;

	sign = 'N';
	Sub *sysnew_operator_NN = new Sub(sign, para, *sys_operator, *sysnew_space);
	Sub *systrun_operator_NN = new Sub(para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_NN);
//	sysnew_operator_NN -> Print_NN(1, lsys);
	systrun_operator_NN -> Print_NN(1, lsys);
	delete systrun_operator_NN;
	delete sysnew_operator_NN;

	sign = 'H';
	Sub *sysnew_operator_H = new Sub(sign, para, *sys_operator, *sysnew_space);
	Sub *systrun_operator_H = new Sub(para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_H);
//	sysnew_operator_H -> Print_H(1, lsys);
	systrun_operator_H -> Print_H(1, lsys);
	delete systrun_operator_H;
	delete sysnew_operator_H;

	delete sys_operator;

        sign = 'r';
        SuperEnergy *trun_sys = new SuperEnergy(para, *sup_space, *reduced_density_sys, *systrun_space, sign);
        delete trun_sys;

	delete systrun_space;
	delete reduced_density_sys;

//------New operator for envnew----------------------------------------------------------------
    //--Create space for eigensystem of reduced density matrix
	Sub *reduced_density_env = new Sub(*envnew_space);

    //--Find the eigensystem of reduced density matrix
	DenMat *denenv = new DenMat(2, para, *sup_space, *reduced_density_env);
	delete denenv;

    //--Create space for truncation
        Sub *envtrun_space = new Sub(para, *reduced_density_env);
	envtrun_space->Print_space(2, lsys);

	sign = 'T';
	Sub *envnew_operator_T = new Sub(para, sign, *env_operator, *envnew_space);
	Sub *envtrun_operator_T = new Sub(para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_T);
//	envnew_operator_T -> Print_T(2, lsys);
	envtrun_operator_T -> Print_T(2, lsys);
	delete envtrun_operator_T;
	delete envnew_operator_T;

	sign = 'S';
	Sub *envnew_operator_S = new Sub(para, sign, *env_operator, *envnew_space);
	Sub *envtrun_operator_S = new Sub(para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_S);
//	envnew_operator_S -> Print_S_Dia(2, lsys);
	envtrun_operator_S -> Print_S_Dia(2, lsys);
	delete envtrun_operator_S;
	delete envnew_operator_S;

	sign = 'M';
	Sub *envnew_operator_M = new Sub(para, sign, *env_operator, *envnew_space);
	Sub *envtrun_operator_M = new Sub(para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_M);
//	envnew_operator_M -> Print_S_M_Dia(2, lsys);
	envtrun_operator_M -> Print_S_M_Dia(2, lsys);
	delete envtrun_operator_M;
	delete envnew_operator_M;

	sign = 'N';
	Sub *envnew_operator_NN = new Sub(para, sign, *env_operator, *envnew_space);
	Sub *envtrun_operator_NN = new Sub(para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_NN);
//	envnew_operator_NN -> Print_NN(2, lsys);
	envtrun_operator_NN -> Print_NN(2, lsys);
	delete envtrun_operator_NN;
	delete envnew_operator_NN;

	sign = 'H';
	Sub *envnew_operator_H = new Sub(para, sign, *env_operator, *envnew_space);
	Sub *envtrun_operator_H = new Sub(para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_H);
//	envnew_operator_H -> Print_H(2, lsys);
	envtrun_operator_H -> Print_H(2, lsys);
	delete envtrun_operator_H;
	delete envnew_operator_H;

	delete env_operator;

        sign = 'l';
        SuperEnergy *trun_env = new SuperEnergy(para, *sup_space, *reduced_density_env, *envtrun_space, sign);
        delete trun_env;

	delete envtrun_space;
	delete reduced_density_env;

//------Delete memories---------------------------------------------------------------
	delete sup_space;
	delete sysnew_space;
	delete envnew_space;
	delete sys_space;
	delete env_space;

}

//==========================Infinite Iteration process========================
inline void Sweep::Infinite_test( Parameter &para, int &lsys, const int &lns ) {

	cout<<"\n System Length="<<lsys<<endl;

//------Read in the data from the hard disk-----------------------------------
    //--Read in sys_space and env_space 
        Sub *sys_space = new Sub( 1, lsys );
	Sub *env_space = new Sub( 2, lsys );

    //--Read in sys and env operators
	Sub *sys_operator = new Sub( 1, lsys, *sys_space );
	Sub *env_operator = new Sub( 2, lsys, *env_space );

//        for(int i = 0; i < sys_operator->operator_number_J; i++ )
//        for(int j = 0; j < sys_operator->Block_Num_hole; j++ )
//        for(int k = 0; k < sys_operator->Num_J_hole_block[j]; k++ )
//        for(int l = 0; l < sys_operator->Dim_J_block[j][k] * sys_operator->Dim_J_block[j][k]; l++ )

//                cout<<"\n"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<l<<"\t"<<sys_operator->S_Dia[i][j][k][l];

//------Space of blocks of sysnew and envnew----------------------------------
        lsys += lns;

        Sub *sysnew_space = new Sub( 1, para, *sys_space );
	Sub *envnew_space = new Sub( 2, para, *env_space );

	sysnew_space -> Print_space(1, lsys);
	envnew_space -> Print_space(2, lsys);

//------Diagonalization of the Hamiltonian-------------------------------------
    //--Allocate space for wave function blocks
	Super *sup_space = new Super(para, sys_space, env_space, sysnew_space, envnew_space);

    //--Allocate space for diagonalization
	sign = 'I';
	Super *sup_diagonalization = new Super(sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space);

    //--Diagonalization  of the Hamiltonian by Matrix-Vector multiplication
	para.total_site++;
	precision = 1.e-12;
	SuperEnergy *supeng = new SuperEnergy(para, *sup_space, *sup_diagonalization, precision);

        ground_state_energy = supeng->eigenvalue;
        excited_state_energy = 0.0;
        WriteEnergy(para, ground_state_energy, excited_state_energy, lsys-1, para.total_site);

	delete supeng;
	delete sup_diagonalization;

//------New operator for sysnew-----------------------------------------------------------------
    //--New and truncate T operator
	sign = 'T';
	Sub *sysnew_operator_T = new Sub(sign, para, *sys_operator, *sysnew_space);
	sysnew_operator_T -> Print_T(1, lsys);
	delete sysnew_operator_T;

	sign = 'S';
	Sub *sysnew_operator_S = new Sub(sign, para, *sys_operator, *sysnew_space);
	sysnew_operator_S -> Print_S_Dia(1, lsys);
	delete sysnew_operator_S;

	sign = 'M';
	Sub *sysnew_operator_M = new Sub(sign, para, *sys_operator, *sysnew_space);
	sysnew_operator_M -> Print_S_M_Dia(1, lsys);
	delete sysnew_operator_M;

	sign = 'N';
	Sub *sysnew_operator_NN = new Sub(sign, para, *sys_operator, *sysnew_space);
	sysnew_operator_NN -> Print_NN(1, lsys);
	delete sysnew_operator_NN;

	sign = 'H';
	Sub *sysnew_operator_H = new Sub(sign, para, *sys_operator, *sysnew_space);
	sysnew_operator_H -> Print_H(1, lsys);
	delete sysnew_operator_H;

	delete sys_operator;

//------New operator for envnew----------------------------------------------------------------
	sign = 'T';
	Sub *envnew_operator_T = new Sub(para, sign, *env_operator, *envnew_space);
	envnew_operator_T -> Print_T(2, lsys);
	delete envnew_operator_T;

	sign = 'S';
	Sub *envnew_operator_S = new Sub(para, sign, *env_operator, *envnew_space);
	envnew_operator_S -> Print_S_Dia(2, lsys);
	delete envnew_operator_S;

	sign = 'M';
	Sub *envnew_operator_M = new Sub(para, sign, *env_operator, *envnew_space);
	envnew_operator_M -> Print_S_M_Dia(2, lsys);
	delete envnew_operator_M;

	sign = 'N';
	Sub *envnew_operator_NN = new Sub(para, sign, *env_operator, *envnew_space);
	envnew_operator_NN -> Print_NN(2, lsys);
	delete envnew_operator_NN;

	sign = 'H';
	Sub *envnew_operator_H = new Sub(para, sign, *env_operator, *envnew_space);
	envnew_operator_H -> Print_H(2, lsys);
	delete envnew_operator_H;

	delete env_operator;

//------Delete memories---------------------------------------------------------------
//	delete sup_space;
	delete sysnew_space;
	delete envnew_space;
	delete sys_space;
	delete env_space;

}

//=====================================================================================
//                                Finite Sweep Iteration 
//=====================================================================================
inline void Sweep::Finite_Sweep( Parameter &para, const int &lsys_start, const int &lsys_end, char &D, const int &optimalstate, const int &step, const double &precision ) {

        if( D == 'R' )      

		Finite_Sweep_LeftToRight( para, lsys_start, lsys_end, optimalstate, step, precision );
//		Finite_Sweep_LeftToRight_test( para, lsys_start, lsys_end, optimalstate, step, precision );

        else if( D == 'L' ) 

		Finite_Sweep_RightToLeft( para, lsys_start, lsys_end, optimalstate, step, precision );
//		Finite_Sweep_RightToLeft_test( para, lsys_start, lsys_end, optimalstate, step, precision );

}

//=====================================================================================
//                       Finite_Sweep: From Left to Right process
//=====================================================================================
inline void Sweep::Finite_Sweep_LeftToRight( Parameter &para, const int &lsys_start, const int &lsys_end, const int &optimalstate, const int &step, const double &precision ) {

	for( lsys = lsys_start; lsys < lsys_end; ) { //lsys denotes the length of system block

                cout<<"\n lsys="<<lsys;

                Sub *sys_space = new Sub( 1, lsys );
                Sub *env_space = new Sub( 2, para.Total_N - lsys - 2 );

                Sub *sys_operator = new Sub( 1, lsys, *sys_space );
                Sub *env_operator = new Sub( 2, para.Total_N - lsys - 2, *env_space );

                lsys += lns;

                Sub *sysnew_space = new Sub( 1, para, *sys_space );
                Sub *envnew_space = new Sub( para, "nonoperator", *env_space );

                Super *sup_space = new Super( para, sys_space, env_space, sysnew_space, envnew_space );

                sign = 'F';
                Super *sup_diagonalization = new Super( sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space );

                if( lsys - 1 > site_start ) {

			sign = 'R';
                        para.total_site++;
                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, sign, precision );
//                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, precision );

                        ground_state_energy = supeng->eigenvalue;
                        excited_state_energy = 0.0;//supeng->eigenvalue_excited;
                        WriteEnergy( para, ground_state_energy, excited_state_energy, lsys-1, para.total_site );
                        WriteEnergy_sweep( para, lsys-1, step, ground_state_energy );
			delete supeng;

		}

                else if( lsys - 1 == site_start ) {

//                      SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization,  precision );
                      SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, 2, precision );
                        delete supeng;

                }

		delete sup_diagonalization;

		delete env_operator;

        //------Adding site and truncate the system block--------------------------------------------------------
                Sub *reduced_density_sys = new Sub( *sysnew_space );

		DenMat *densys = new DenMat( 1, para, *sup_space, *reduced_density_sys );
	        delete densys;

	        Sub *systrun_space=new Sub( para, *reduced_density_sys, optimalstate );
	      	systrun_space->Print_space( 1, lsys );

		sign = 'T';
		Sub *sysnew_operator_T = new Sub( sign, para, *sys_operator, *sysnew_space );
		Sub *systrun_operator_T = new Sub( para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_T );
//		sysnew_operator_T->Print_T( 1, lsys );
		systrun_operator_T->Print_T( 1, lsys );
		delete systrun_operator_T;
		delete sysnew_operator_T;

		sign = 'S';
		Sub *sysnew_operator_S = new Sub( sign, para, *sys_operator, *sysnew_space );
		Sub *systrun_operator_S = new Sub( para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_S );
//		sysnew_operator_S->Print_S_Dia( 1, lsys );
		systrun_operator_S->Print_S_Dia( 1, lsys );
		delete systrun_operator_S;
		delete sysnew_operator_S;

		sign = 'M';
		Sub *sysnew_operator_M = new Sub( sign, para, *sys_operator, *sysnew_space );
		Sub *systrun_operator_M = new Sub( para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_M );
//		sysnew_operator_M->Print_S_M_Dia( 1, lsys );
		systrun_operator_M->Print_S_M_Dia( 1, lsys );
		delete systrun_operator_M;
		delete sysnew_operator_M;

		sign = 'N';
		Sub *sysnew_operator_NN = new Sub( sign, para, *sys_operator, *sysnew_space );
		Sub *systrun_operator_NN = new Sub( para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_NN );
//		sysnew_operator_NN->Print_NN( 1, lsys );
		systrun_operator_NN->Print_NN( 1, lsys );
		delete systrun_operator_NN;
		delete sysnew_operator_NN;

		sign = 'H';
		Sub *sysnew_operator_H = new Sub( sign, para, *sys_operator, *sysnew_space );
		Sub *systrun_operator_H = new Sub( para, sign, *reduced_density_sys, *systrun_space, *sysnew_operator_H );
//		sysnew_operator_H->Print_H( 1, lsys );
		systrun_operator_H->Print_H( 1, lsys );
		delete systrun_operator_H;
		delete sysnew_operator_H;

		delete sys_operator;

                if( lsys < site_end ) {     //site increasing steps

                        sign = 'r';
                        SuperEnergy *trun_func = new SuperEnergy( para, *sup_space, *reduced_density_sys, *systrun_space, sign );
                        delete trun_func;

                }

                else if( lsys == site_end ) {       //reach the end of the sweep in a direction

                        SuperEnergy *trun_func=new SuperEnergy( *sup_space, *systrun_space, 1 );
                        delete trun_func;

                }

		delete systrun_space;
                delete reduced_density_sys;

                delete sup_space;
                delete sysnew_space;	delete envnew_space;
                delete sys_space;       delete env_space;

	}

}

//======Test part======
inline void Sweep::Finite_Sweep_LeftToRight_test( Parameter &para, const int &lsys_start, const int &lsys_end, const int &optimalstate, const int &step, const double &precision ) {

	for( lsys = lsys_start; lsys < lsys_end; ) { //lsys denotes the length of system block

                cout<<"\n lsys="<<lsys;

                Sub *sys_space = new Sub( 1, lsys );
                Sub *env_space = new Sub( 2, para.Total_N - lsys - 2 );

                Sub *sys_operator = new Sub( 1, lsys, *sys_space );
                Sub *env_operator = new Sub( 2, para.Total_N - lsys - 2, *env_space );

                lsys += lns;

                Sub *sysnew_space = new Sub( 1, para, *sys_space );
                Sub *envnew_space = new Sub( para, "nonoperator", *env_space );

		sysnew_space->Print_space( 1, lsys );

                Super *sup_space = new Super( para, sys_space, env_space, sysnew_space, envnew_space );

                sign = 'F';
               Super *sup_diagonalization = new Super( sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space );

                if( lsys - 1 > site_start ) {

			sign = 'R';
                        para.total_site++;
//                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, sign, precision );
                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, precision );

                        ground_state_energy = supeng->eigenvalue;
                        excited_state_energy = 0.0;//supeng->eigenvalue_excited;
                        WriteEnergy( para, ground_state_energy, excited_state_energy, lsys-1, para.total_site );
                        WriteEnergy_sweep( para, lsys-1, step, ground_state_energy );
			delete supeng;

		}

                else if( lsys - 1 == site_start ) {

		//	SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, 2, precision );
                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, precision );
                        delete supeng;

                }

		delete sup_diagonalization;
                delete sup_space;

		delete env_operator;

        //------Adding site and truncate the system block--------------------------------------------------------
		sign = 'T';
		Sub *sysnew_operator_T = new Sub( sign, para, *sys_operator, *sysnew_space );
		sysnew_operator_T->Print_T( 1, lsys );
		delete sysnew_operator_T;

		sign = 'S';
		Sub *sysnew_operator_S = new Sub( sign, para, *sys_operator, *sysnew_space );
		sysnew_operator_S->Print_S_Dia( 1, lsys );
		delete sysnew_operator_S;

		sign = 'M';
		Sub *sysnew_operator_M = new Sub( sign, para, *sys_operator, *sysnew_space );
		sysnew_operator_M->Print_S_M_Dia( 1, lsys );
		delete sysnew_operator_M;

		sign = 'N';
		Sub *sysnew_operator_NN = new Sub( sign, para, *sys_operator, *sysnew_space );
		sysnew_operator_NN->Print_NN( 1, lsys );
		delete sysnew_operator_NN;

		sign = 'H';
		Sub *sysnew_operator_H = new Sub( sign, para, *sys_operator, *sysnew_space );
		sysnew_operator_H->Print_H( 1, lsys );
		delete sysnew_operator_H;

		delete sys_operator;

                delete sysnew_space;	delete envnew_space;
                delete sys_space;       delete env_space;

	}

}

//===============================================================================================================
//                                    Finite_Sweep: From Right to Left process
//===============================================================================================================
inline void Sweep::Finite_Sweep_RightToLeft( Parameter &para, const int &lsys_start, const int &lsys_end, const int &optimalstate, const int &step, const double &precision ) {

        for( lsys = lsys_start; lsys < lsys_end; ) { //Here, lsys denotes the length of environment block

                cout<<"\n lsys="<<lsys;

                Sub *sys_space = new Sub( 1, para.Total_N-lsys-2 );
                Sub *env_space = new Sub( 2, lsys );

                Sub *sys_operator = new Sub( 1, para.Total_N-lsys-2, *sys_space );
                Sub *env_operator = new Sub( 2, lsys, *env_space );

                lsys += lns;

                Sub *sysnew_space = new Sub( para, "nonoperator", *sys_space );
                Sub *envnew_space = new Sub( 2, para, *env_space );

                Super *sup_space = new Super( para, sys_space, env_space, sysnew_space, envnew_space );

                sign = 'F';
                Super *sup_diagonalization = new Super( sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space );

                if( lsys - 1 > site_start ) {

			sign = 'L';
                        para.total_site++;
                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, sign, precision );
//                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, precision );

                        ground_state_energy = supeng->eigenvalue;
                        excited_state_energy = 0.0;//supeng->eigenvalue_excited;
                        WriteEnergy( para, ground_state_energy, excited_state_energy, lsys-1, para.total_site );
                        WriteEnergy_sweep( para, lsys-1, step, ground_state_energy );
			delete supeng;

		}

                else if( lsys - 1 == site_start ) {

//                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization,  precision );
                      SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, 1, precision );
                        delete supeng;

                }

		delete sup_diagonalization;

		delete sys_operator;

        //------Adding site and truncate the system block--------------------------------------------------------
                Sub *reduced_density_env = new Sub( *envnew_space );

		DenMat *denenv = new DenMat( 2, para, *sup_space, *reduced_density_env );
	        delete denenv;

	        Sub *envtrun_space = new Sub( para, *reduced_density_env, optimalstate );
        	envtrun_space->Print_space( 2, lsys );

		sign = 'T';
		Sub *envnew_operator_T = new Sub( *env_operator, *envnew_space, para, sign );
		Sub *envtrun_operator_T = new Sub( para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_T );
//		envnew_operator_T->Print_T( 2, lsys );
		envtrun_operator_T->Print_T( 2, lsys );
		delete envtrun_operator_T;
		delete envnew_operator_T;

		sign = 'S';
		Sub *envnew_operator_S = new Sub( *env_operator, *envnew_space, para, sign );
		Sub *envtrun_operator_S = new Sub( para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_S );
//		envnew_operator_S->Print_S_Dia( 2, lsys );
		envtrun_operator_S->Print_S_Dia( 2, lsys );
		delete envtrun_operator_S;
		delete envnew_operator_S;

		sign = 'M';
		Sub *envnew_operator_M = new Sub( *env_operator, *envnew_space, para, sign );
		Sub *envtrun_operator_M = new Sub( para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_M );
//		envnew_operator_M->Print_S_M_Dia( 2, lsys );
		envtrun_operator_M->Print_S_M_Dia( 2, lsys );
		delete envtrun_operator_M;
		delete envnew_operator_M;

		sign = 'N';
		Sub *envnew_operator_NN = new Sub( *env_operator, *envnew_space, para, sign );
		Sub *envtrun_operator_NN = new Sub( para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_NN );
//		envnew_operator_NN->Print_NN(2, lsys);
		envtrun_operator_NN->Print_NN(2, lsys);
		delete envtrun_operator_NN;
		delete envnew_operator_NN;

		sign = 'H';
		Sub *envnew_operator_H = new Sub( *env_operator, *envnew_space, para, sign );
		Sub *envtrun_operator_H = new Sub( para, sign, *reduced_density_env, *envtrun_space, *envnew_operator_H );
//		envnew_operator_H->Print_H( 2, lsys );
		envtrun_operator_H->Print_H( 2, lsys );
		delete envtrun_operator_H;
		delete envnew_operator_H;

		delete env_operator;

                if( lsys < site_end ) {     //site increasing steps

                        sign = 'l';
                        SuperEnergy *trun_func = new SuperEnergy( para, *sup_space, *reduced_density_env, *envtrun_space, sign );
                        delete trun_func;
                
		}

                else if( lsys == site_end ) {       //reach the end of the sweep in a direction

                        SuperEnergy *trun_func = new SuperEnergy( *sup_space, *envtrun_space, 2 );
                        delete trun_func;

                }

		delete envtrun_space;
                delete reduced_density_env;

                delete sup_space;
                delete sysnew_space;            delete envnew_space;
                delete sys_space;               delete env_space;

	}

}

//======Test part======
inline void Sweep::Finite_Sweep_RightToLeft_test( Parameter &para, const int &lsys_start, const int &lsys_end, const int &optimalstate, const int &step, const double &precision ) {

        for( lsys = lsys_start; lsys < lsys_end; ) { //Here, lsys denotes the length of environment block

                cout<<"\n lsys="<<lsys;

                Sub *sys_space = new Sub( 1, para.Total_N-lsys-2 );
                Sub *env_space = new Sub( 2, lsys );

                Sub *sys_operator = new Sub( 1, para.Total_N-lsys-2, *sys_space );
                Sub *env_operator = new Sub( 2, lsys, *env_space );

                lsys += lns;

                Sub *sysnew_space = new Sub( para, "nonoperator", *sys_space );
                Sub *envnew_space = new Sub( 2, para, *env_space );

		envnew_space->Print_space( 2, lsys );

                Super *sup_space = new Super( para, sys_space, env_space, sysnew_space, envnew_space );

                sign = 'F';
                Super *sup_diagonalization = new Super( sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space );

                if( lsys - 1 > site_start ) {

			sign = 'L';
                        para.total_site++;
//                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, sign, precision );
                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, precision );

                        ground_state_energy = supeng->eigenvalue;
                        excited_state_energy = 0.0;//supeng->eigenvalue_excited;
                        WriteEnergy( para, ground_state_energy, excited_state_energy, lsys-1, para.total_site );
                        WriteEnergy_sweep( para, lsys-1, step, ground_state_energy );
			delete supeng;

		}

                else if( lsys - 1 == site_start ) {

//			SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, 1, precision );
                        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization,  precision );
                        delete supeng;

                }

		delete sup_diagonalization;
		delete sup_space;

		delete sys_operator;

        //------Adding site and truncate the system block--------------------------------------------------------
		sign = 'T';
		Sub *envnew_operator_T = new Sub( *env_operator, *envnew_space, para, sign );
		envnew_operator_T->Print_T( 2, lsys );
		delete envnew_operator_T;

		sign = 'S';
		Sub *envnew_operator_S = new Sub( *env_operator, *envnew_space, para, sign );
		envnew_operator_S->Print_S_Dia( 2, lsys );
		delete envnew_operator_S;

		sign = 'M';
		Sub *envnew_operator_M = new Sub( *env_operator, *envnew_space, para, sign );
		envnew_operator_M->Print_S_M_Dia( 2, lsys );
		delete envnew_operator_M;

		sign = 'N';
		Sub *envnew_operator_NN = new Sub( *env_operator, *envnew_space, para, sign );
		envnew_operator_NN->Print_NN(2, lsys);
		delete envnew_operator_NN;

		sign = 'H';
		Sub *envnew_operator_H = new Sub( *env_operator, *envnew_space, para, sign );
		envnew_operator_H->Print_H( 2, lsys );
		delete envnew_operator_H;

		delete env_operator;

                delete sysnew_space;            delete envnew_space;
                delete sys_space;               delete env_space;

	}

}

//===============================================================================================================
//                                      Measurement
//===============================================================================================================
inline void Sweep::Measurement( Parameter &para, const int &step, const double &precision ) {

        if( step == 0 ) {

		cout << "------[begin measurement]------" <<endl;
		cout << "------[Measurement_LeftToRight]------"<<endl;
		Measurement_LeftToRight( para, precision );
		cout << "------[Measurement_LeftToRight_shifted, measure pairing]------"<<endl;
		Measurement_LeftToRight_shifted( para, precision );
	
	}

        else if( step == 1 ) Measurement_RightToLeft( para, precision );


}

//-------------------------------------------------------------------------------------
inline void Sweep::Measurement_RightToLeft( Parameter &para, const double &precision ) {

        lsys = para.Total_N/2 - 1;
        //lsys = 64;
        cout<<"\n lsys="<<lsys<<endl;

	Sub *sys_space = new Sub( 1, para.Total_N - lsys - 2 );
	Sub *env_space = new Sub( 2, lsys );

        Sub *sys_operator = new Sub( 1, para.Total_N - lsys - 2, *sys_space );
        Sub *env_operator = new Sub( 2, lsys, *env_space );

	Sub *sysnew_space = new Sub( para, "nonoperator", *sys_space );
	Sub *envnew_space = new Sub( 2, para, *env_space );

	Super *sup_space = new Super( para, sys_space, env_space, sysnew_space, envnew_space );

        sign = 'F';
        Super *sup_diagonalization = new Super( sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space );

        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, precision );
        WriteWave_Function( sup_space->BlockNumber_for_TargetBlock, sup_space->Dim_block, sup_space->WaveFunction_block );
	delete supeng;
	delete sup_diagonalization;
	delete sup_space;

	delete sysnew_space;  delete envnew_space;
	delete sys_operator;  delete env_operator;
	delete sys_space;     delete env_space;

	Sub *sys_measure = new Sub( 1, para.Total_N - lsys - 2 );
	Sub *env_measure = new Sub( 2, lsys );

	Sub *sysnew_measure = new Sub( para, "nonoperator", *sys_measure );
	Sub *envnew_measure = new Sub( 2, para, *env_measure );

	Super *sup_measure = new Super( sys_measure, env_measure, sysnew_measure, envnew_measure, para );

	ElectronDensity *electrondensity = new ElectronDensity( para.Total_N - lsys - 2, para, *sup_measure );
	delete electrondensity;

	delete sup_measure;
	delete sysnew_measure;  delete envnew_measure;
	delete sys_measure;     delete env_measure;

}

//-------------------------------------------------------------------------------------
inline void Sweep::Measurement_LeftToRight( Parameter &para, const double &precision ) {

	lsys = para.Total_N/2 - 1;
	cout<<"\n lsys="<<lsys<<endl;

	Sub *sys_space = new Sub( 1, lsys );
	Sub *env_space = new Sub( 2, para.Total_N - lsys - 2 );

        Sub *sys_operator = new Sub( 1, lsys, *sys_space );
        Sub *env_operator = new Sub( 2, para.Total_N - lsys - 2, *env_space );

	Sub *sysnew_space = new Sub( 1, para, *sys_space );
	Sub *envnew_space = new Sub( para, "nonoperator", *env_space );

	Super *sup_space = new Super( para, sys_space, env_space, sysnew_space, envnew_space );

        sign = 'F';
        Super *sup_diagonalization = new Super( sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space );

    	sign = 'R';
        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, sign, precision );
        WriteWave_Function( sup_space->BlockNumber_for_TargetBlock, sup_space->Dim_block, sup_space->WaveFunction_block );
	delete supeng;
	delete sup_diagonalization;
	delete sup_space;

	delete sysnew_space;  delete envnew_space;
	delete sys_operator;  delete env_operator;
	delete sys_space;     delete env_space;

	Sub *sys_measure = new Sub( 1, lsys );
	Sub *env_measure = new Sub( 2, para.Total_N - lsys - 2 );

	Sub *sysnew_measure = new Sub( 1, para, *sys_measure );
	Sub *envnew_measure = new Sub( para, "nonoperator", *env_measure );

	Super *sup_measure = new Super( sys_measure, env_measure, sysnew_measure, envnew_measure, para );

	ElectronDensity *electrondensity = new ElectronDensity( lsys, para, *sup_measure );
	delete electrondensity;

	SpinCorrelation *spin = new SpinCorrelation( lsys, para, *sup_measure );
	delete spin;

	//SqSpin *sqspin = new SqSpin( lsys, para, *sup_measure );
	//delete sqspin;

	GreenFunction *green = new GreenFunction( lsys, para, *sup_measure );
	delete green;

	//SqHopping *sqhopping = new SqHopping( lsys, para, *sup_measure );
	//delete sqhopping;

//    SqDensity *sqdensity = new SqDensity( lsys, para, *sup_measure );
//    delete sqdensity;

        DensityCorrelation *density = new DensityCorrelation( lsys, para, *sup_measure );
        delete density;

	delete sup_measure;
	delete sysnew_measure;  delete envnew_measure;
	delete sys_measure;     delete env_measure;

}

//-------------------------------------------------------------------------------------
inline void Sweep::Measurement_LeftToRight_shifted( Parameter &para, const double &precision ) {

	lsys = 47;
	cout<<"\n lsys="<<lsys<<endl;

	Sub *sys_space = new Sub( 1, lsys );
	Sub *env_space = new Sub( 2, para.Total_N - lsys - 2 );

        Sub *sys_operator = new Sub( 1, lsys, *sys_space );
        Sub *env_operator = new Sub( 2, para.Total_N - lsys - 2, *env_space );

	Sub *sysnew_space = new Sub( 1, para, *sys_space );
	Sub *envnew_space = new Sub( para, "nonoperator", *env_space );

	Super *sup_space = new Super( para, sys_space, env_space, sysnew_space, envnew_space );

        sign = 'F';
        Super *sup_diagonalization = new Super( sign, para, sys_operator, env_operator, sysnew_space, envnew_space, *sup_space );

        sign = 'R';
        SuperEnergy *supeng = new SuperEnergy( para, *sup_space, *sup_diagonalization, sign, precision );
        WriteWave_Function_Shifted( sup_space->BlockNumber_for_TargetBlock, sup_space->Dim_block, sup_space->WaveFunction_block );
	delete supeng;
	delete sup_diagonalization;
	delete sup_space;

	delete sysnew_space;  delete envnew_space;
	delete sys_operator;  delete env_operator;
	delete sys_space;     delete env_space;

	Sub *sys_measure = new Sub( 1, lsys );
	Sub *env_measure = new Sub( 2, para.Total_N - lsys - 2 );

	Sub *sysnew_measure = new Sub( 1, para, *sys_measure );
	Sub *envnew_measure = new Sub( para, "nonoperator", *env_measure );

	Super *sup_measure = new Super( sys_measure, env_measure, sysnew_measure, envnew_measure, para );

	PairingCorrelation *pairing = new PairingCorrelation( lsys, para, *sup_measure );
	delete pairing;

        PairingCorrelation_yy *pairing_yy = new PairingCorrelation_yy( lsys, para, *sup_measure);
        delete pairing_yy;

//	DensityCorrelation *density = new DensityCorrelation( lsys, para, *sup_measure );
//	delete density;

//  	ElectronDensity *electrondensity = new ElectronDensity( lsys, para, *sup_measure );
//  	delete electrondensity;

	delete sup_measure;
	delete sysnew_measure;  delete envnew_measure;
	delete sys_measure;     delete env_measure;

}

//================================================================================================================
//                                      Save the ground state energy to file
//===============================================================================================================
inline void Sweep::WriteEnergy(Parameter &para, const double &ground_state_energy, const double &excited_state_energy, const int &lsys, const int &site) {

        FILE *e=fopen("e/ground_state_energy", "a+");
        fprintf(e, "%d\t", site);
        fprintf(e, "%-20.20f\n", ground_state_energy);
        fclose(e);

        if(lsys==para.Total_N/2-1) {
                FILE *b=fopen("e/ground_state_energy_middle", "a+");
                fprintf(b, "%-20.20f\n", ground_state_energy);
                fclose(b);
        }
}

//===============================================================================================================
//                              Save the ground state energy in each sweep step to file
//===============================================================================================================
inline void Sweep::WriteEnergy_sweep( Parameter &para, const int &site_num, const int &step_num, const double &energy ) {
        FILE *e=fopen(Combine("e/ground_state_energy_sweep-", step_num) , "a+");
        fprintf(e, "%d\t", site_num);
        fprintf(e, "%-20.20f\n", energy);
        fclose(e);
}

//===============================================================================================================
//                                      Save the wave function
//===============================================================================================================
inline void Sweep::WriteWave_Function( const int &BlockNumber, int * Dim_block, double ** WaveFunction ) {

        FILE *g = fopen( "wavefunction/wavefunction_ground", "w+" );

        for( int i = 0; i < BlockNumber; i++ )
	for( int j = 0; j < Dim_block[i]; j++ )

                fprintf(g, "%f\n", WaveFunction[i][j]);

        fclose(g);

}

//===============================================================================================================
inline void Sweep::WriteWave_Function_Shifted( const int &BlockNumber, int * Dim_block, double ** WaveFunction ) {

        FILE *g = fopen( "wavefunction/wavefunction_ground_shifted", "w+" );

        for( int i = 0; i < BlockNumber; i++ )
	for( int j = 0; j < Dim_block[i]; j++ )

                fprintf(g, "%f\n", WaveFunction[i][j]);

        fclose(g);

}



//===============================================================================================================


