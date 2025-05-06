#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include "random_list_generator.h"
using namespace std;

#include"parameter.h"

//==========================Transmit the parameters in this class=========================
Parameter::Parameter( char &sign, const int &spin, const int &total_h, const int &total_j, const int &number_hole, const int &statenokept, const int &n_u, const int &n_x, const int &n_y, const int &total_n, const double &t_n, const double &t_nn, const double &j_n, const double &j_nn, const double &n_n, const double &n_nn, const int &desiredzeros, const int &seedvalue ) {

  //----Initialize parameters
	Sign = sign;  
	Spin = spin;  
	StateNoKept = statenokept;  

	Total_h = total_h;  
	Total_J = total_j; 
	Number_hole = number_hole;	
	
	N_u = n_u;  
	N_x = n_x;  
	N_y = n_y;  
	Total_N = total_n; // total lattice sites                  
	TotalSite = Total_N * Total_N;

	untruncated_site = 2;//denote the number of untruncated sites
        total_site = 0;
        
	T_n = t_n;  
	T_nn = t_nn;  
	J_n = j_n;  
	J_nn = j_nn;  
	N_n = n_n;  
	N_nn = n_nn;

	Desiredzeros = desiredzeros;
	Seedvalue = seedvalue;

  //----Set the tables
	QuantumNumber_hole = new int [ Total_N ];
        QuantumNumber_J = new int [ Total_N ];

        for( int i = 0; i < Total_N; i++ ) {

                QuantumNumber_hole[i] = 0;     QuantumNumber_J[i] = 0;

        }


        Table_T = new int [ TotalSite ];
        Interaction_T = new double [ TotalSite ];

        for( int i = 0; i < TotalSite; i++ ) {

                Table_T[i] = 0;     Interaction_T[i] = 0.0;

        }


        Table_J = new int [ TotalSite ];
        Interaction_J = new double [ TotalSite ];

        for( int i = 0; i < TotalSite; i++ ) {

                Table_J[i] = 0;     Interaction_J[i] = 0.0;

        }


        Table_N = new int [ TotalSite ];
	Interaction_N = new double [ TotalSite ];

        for( int i = 0; i < TotalSite; i++ ) {

                Table_N[i] = 0;     Interaction_N[i] = 0.0;

        }

  //----Set quantum numbers for each system lengths
	QuantumNumber_Initial();

  //----Build the Table of interaction
	Interaction_Table_Initial();

}

//====================================================================
inline void Parameter::QuantumNumber_Initial() {

	for( int i = 0; i < Total_N; i++ )
	if( i > 2 ) {  // "i + 1" is the number of total sites 

		//----
		if( Total_h % 2 == 0 ) {
  
			//----
                	QuantumNumber_hole[ i ] = (int) ( ( i + 1 ) * Total_h / Total_N );

                 	if( QuantumNumber_hole[ i ] % 2 == 1 ) {

                        	QuantumNumber_hole[ i ] = QuantumNumber_hole[ i ] + 1;
			}

			//----
			QuantumNumber_J[ i ] = (int) ( Total_J * ( ( i + 1 ) - QuantumNumber_hole[ i ] ) / ( Total_N - Total_h ) );

 
			if( QuantumNumber_J[ i ] % 2 == 1 ) {

                        	QuantumNumber_J[ i ] = QuantumNumber_J[ i ] + 1;
			}

        	}       

		//----
		else if( Total_h % 2 == 1 ) {

			//----
	        	QuantumNumber_hole[ i ] = (int) ( ( i + 1 ) * Total_h / Total_N );

              		if( QuantumNumber_hole[ i ] % 2 == 0 ) {

                        	QuantumNumber_hole[ i ] = QuantumNumber_hole[ i ] + 1;

			}
			//----
			QuantumNumber_J[ i ] = (int) ( Total_J * ( ( i + 1 ) - QuantumNumber_hole[ i ] ) / ( Total_N - Total_h ) );

	                if( QuantumNumber_J[ i ] % 2 == 0 ) {

                        	QuantumNumber_J[ i ] = QuantumNumber_J[ i ] + 1;

			}

         	}

	}
        //---YY
        /*
        for(int i =0; i<Total_N; i++)
        {
          //cout << "i=: " << i << ", QuantumNumber_hole[ i ]: " << QuantumNumber_hole[ i ] << endl; 
          //cout << "i=: " << i << ", QuantumNumber_J[ i ]: " << QuantumNumber_J[ i ] << endl;
        }
        */
        //---YY

}

//====================================================================
inline void Parameter::Interaction_Table_Initial() {

	int site_1, site_2;

//------ bilayer 2d square lattice

	// random list along zz bonds
	randomList = RandomListGenerator::generateRandomList(N_y * N_x, Desiredzeros, Seedvalue);
	std::cout << "Random list of 0s and 1s with seed " << Seedvalue << ": ";
	for (int i = 0; i < N_y * N_x; ++i){
		std::cout << randomList[i] << " ";
	}
	std::cout << std::endl;


//----couplings between two layers

	//----vertical bonds, interlayer, z direaction 
        cout << " vertical bonds, interlayer, z direaction " << endl;

	int it = 0;
	for( int i = 0; i < N_x; i++ )
	for( int j = 0; j < N_y; j++ ) { 

		site_1 = N_y * N_u * i + j * N_u;  site_2 = site_1 + 1;

                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;

		Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
		Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
		Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

		Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_nn;
		Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_nn * randomList[it];
		Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_nn;

		cout << "vertical bonds interaction, Jnn*randomList["<<it<<"]"<< J_nn * randomList[it] << endl;

		it++;

	}
        // ---------[y 方向 up layer] ------------
         cout << "----[y 方向 up layer] ---" << endl;
        for(int i=0; i<N_x; i++)
        for(int j=0; j<N_y-1; j++){
                site_1 = N_y * N_u * i + j * N_u; site_2 = site_1 + 2;
                 cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;
                
                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
		Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
		Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

		Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_n;
		Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_n;
		Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_n;

        }

        // ---------[y 方向 up layer, periodic boundary]--------
	/*
         cout << "---[y 方向 up layer, periodic boundary]---" << endl;
        for(int i=0; i<N_x; i++)
        {
                site_1 = N_y * N_u * i; site_2 = site_1 + (N_y-1) * N_u;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;
                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
		Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
		Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

		Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_n;
		Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_n;
		Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_n;

        }
	*/

        // --------------[y 方向 down layer]-----------
         cout << "---[y 方向 down layer]---" << endl;
        for(int i=0; i<N_x; i++)
        for(int j=0; j<N_y-1; j++){
                site_1 = N_y * N_u * i + j * N_u + 1; site_2 = site_1 + 2;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;
                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
		Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
		Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

		Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_n;
		Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_n;
		Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_n;


        }
        //-----[y 方向 down layer, periodic boundary]---------
	/*
        cout << "[y 方向 down layer, periodic boundary]" << endl;
        for(int i=0; i<N_x; i++)
        {
                site_1 = N_y * N_u * i + 1; site_2 = site_1 + (N_y-1) * N_u;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;
                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
		Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
		Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

		Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_n;
		Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_n;
		Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_n;
        }
	*/
        // ---------------[x 方向 up layer]--------
         cout << "----[x 方向 up layer]----" << endl;
        for(int i=0; i<N_x-1; i++)
        for(int j=0; j<N_y; j++)
        {
                site_1 = N_y * N_u * i + j * N_u;  site_2 = site_1 + N_y * N_u;
                 cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;
                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
		Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
		Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

		Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_n;
		Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_n;
		Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_n;   
        }
        // -----------[x 方向 down layer]----------
         cout << "---[x 方向 down layer]---" << endl;
        for(int i=0; i<N_x-1; i++)
        for(int j=0; j<N_y; j++)
        {
                site_1 = N_y * N_u * i + j * N_u + 1;  site_2 = site_1 + N_y * N_u;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;
                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
		Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
		Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

		Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_n;
		Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_n;
		Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_n; 
        }
//------Number of the sites whose operators need to be stored: sys block
        int index, number;

        Table_T_sys=new int [Total_N-2];
        Table_T_sys_site=new int * [Total_N-2];

        for(int i=0; i<Total_N-2; i++) {//i+1 is the number of sites in the block, sys block reaches the max length N-2

                index=0;

                for(int j=0; j<=i; j++) {
                        number=Total_N*j;

                        for(int l=i+1; l<Total_N; l++)
                        if(Table_T[number+l] != 0 ) {
                                index++;
                                break;
                        }
                }

                Table_T_sys[i]=index;
                Table_T_sys_site[i]=new int [index];

                index=0;
                for(int j=0; j<=i; j++) {
                        number=Total_N*j;

                        for(int l=i+1; l<Total_N; l++)
                        if( Table_T[number+l] != 0 ) {
                                Table_T_sys_site[i][index++]=j;
                                break;
                        }
                }

        }

//------Number of the sites whose operators need to be stored: env block
        Table_T_env=new int [Total_N-2];
        Table_T_env_site=new int * [Total_N-2];

        for(int i=0; i<Total_N-2; i++) {//i+1 is the length of the environment block
                index=0;

                for(int j=0; j<=i; j++) {
                        number=Total_N*(Total_N-j-1);

                        for(int l=0; l<Total_N-i-1; l++)
                        if( Table_T[number+l] != 0 ) {
                                index++;
                                break;
                        }
                }

                Table_T_env[i]=index;
                Table_T_env_site[i]=new int [index];

                index=0;
                for(int j=0; j<=i; j++) {
                        number=Total_N*(Total_N-j-1);

                        for(int l=0; l<Total_N-i-1; l++)
                        if( Table_T[number+l] != 0 ) {
                                Table_T_env_site[i][index++]=j;//j is the position relative to the StartSite
                                break;
                        }
                }

        }

//------Number of the sites whose operators need to be stored: sys block
        Table_J_sys=new int [Total_N-2];
        Table_J_sys_site=new int * [Total_N-2];

        for(int i=0; i<Total_N-2; i++) {//i+1 is the number of sites in the block, sys block reaches the max length N-2
                index=0;

                for(int j=0; j<=i; j++) {
                        number=Total_N*j;

                        for(int l=i+1; l<Total_N; l++)
                        if(Table_J[number+l] != 0 ) {
                                index++;
                                break;
                        }
                }

                Table_J_sys[i]=index;
                Table_J_sys_site[i]=new int [index];

                index=0;
                for(int j=0; j<=i; j++) {
                        number=Total_N*j;

                        for(int l=i+1; l<Total_N; l++)
                        if(Table_J[number+l] !=0 ) {
                                Table_J_sys_site[i][index++]=j;
                                break;
                        }
                }

        }

//------Number of the sites whose operators need to be stored: env block
        Table_J_env=new int [Total_N-2];
        Table_J_env_site=new int * [Total_N-2];

        for(int i=0; i<Total_N-2; i++) {//i+1 is the length of the environment block
                index=0;

                for(int j=0; j<=i; j++) {
                        number=Total_N*(Total_N-j-1);

                        for(int l=0; l<Total_N-i-1; l++)
                        if(Table_J[number+l]!=0) {
                                index++;
                                break;
                        }
                }

                Table_J_env[i]=index;
                Table_J_env_site[i]=new int [index];

                index=0;
                for(int j=0; j<=i; j++) {
                        number=Total_N*(Total_N-j-1);

                        for(int l=0; l<Total_N-i-1; l++)
                        if(Table_J[number+l]!=0) {
                                Table_J_env_site[i][index++]=j;//j is the position relative to the StartSite
                                break;
                        }
                }

        }

//------Number of the sites whose operators need to be stored: sys block
        Table_N_sys=new int [Total_N-2];
        Table_N_sys_site=new int * [Total_N-2];

        for(int i=0; i<Total_N-2; i++) {//i+1 is the number of sites in the block, sys block reaches the max length N-2
                index=0;

                for(int j=0; j<=i; j++) {
                        number=Total_N*j;

                        for(int l=i+1; l<Total_N; l++)
                        if(Table_N[number+l]!=0) {
                                index++;
                                break;
                        }
                }

                Table_N_sys[i]=index;
                Table_N_sys_site[i]=new int [index];

                index=0;
                for(int j=0; j<=i; j++) {
                        number=Total_N*j;

                        for(int l=i+1; l<Total_N; l++)
                        if(Table_N[number+l]!=0) {
                                Table_N_sys_site[i][index++]=j;
                                break;
                        }
                }

        }

//------Number of the sites whose operators need to be stored: env block
        Table_N_env=new int [Total_N-2];
        Table_N_env_site=new int * [Total_N-2];

        for(int i=0; i<Total_N-2; i++) {//i+1 is the length of the environment block
                index=0;

                for(int j=0; j<=i; j++) {
                        number=Total_N*(Total_N-j-1);

                        for(int l=0; l<Total_N-i-1; l++)
                        if(Table_N[number+l]!=0) {
                                index++;
                                break;
                        }
                }

                Table_N_env[i]=index;
                Table_N_env_site[i]=new int [index];

                index=0;
                for(int j=0; j<=i; j++) {
                        number=Total_N*(Total_N-j-1);

                        for(int l=0; l<Total_N-i-1; l++)
                        if(Table_N[number+l]!=0) {
                                Table_N_env_site[i][index++]=j;//j is the position relative to the StartSite
                                break;
                        }
                }

        }

}

//============================Delete parameter========================
Parameter::~Parameter() {

	delete [] QuantumNumber_hole;	delete [] QuantumNumber_J;
        delete [] randomList;

        delete [] Interaction_T;	delete [] Interaction_J;	delete [] Interaction_N;
        delete [] Table_T;		delete [] Table_J;		delete [] Table_N;
        delete [] Table_T_sys;          delete [] Table_J_sys;		delete [] Table_N_sys;
	delete [] Table_T_env;		delete [] Table_J_env;		delete [] Table_N_env;

        for(int i=0; i<Total_N-2; i++) {
                delete [] Table_T_sys_site[i];	delete [] Table_J_sys_site[i];  delete [] Table_N_sys_site[i];
		delete [] Table_T_env_site[i];	delete [] Table_J_env_site[i];  delete [] Table_N_env_site[i];
        }
        delete [] Table_T_sys_site;		delete [] Table_J_sys_site;	delete [] Table_N_sys_site;
	delete [] Table_T_env_site;		delete [] Table_J_env_site;	delete [] Table_N_env_site;

}
