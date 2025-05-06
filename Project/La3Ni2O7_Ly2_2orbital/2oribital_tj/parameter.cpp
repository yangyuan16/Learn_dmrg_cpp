#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include "random_list_generator.h"
using namespace std;

#include"parameter.h"

//==========================Transmit the parameters in this class=========================
Parameter::Parameter( char &sign, const int &spin, const int &total_h, const int &total_j, const int &number_hole, const int &statenokept, const int &n_u, const int &n_x, const int &n_y, const int &total_n, const double &t_a_xx, const double &j_a_xx, const double &n_a_xx, const double &t_a_zz, const double &j_a_zz, const double &n_a_zz, const double &t_a_xz, const double &j_a_xz, const double &n_a_xz, const double &t_b_xx, const double &j_b_xx, const double &n_b_xx, const double &t_b_zz, const double &j_b_zz, const double &n_b_zz, const double &t_b_xz, const double &j_b_xz, const double &n_b_xz, const double &t_ab_xx, const double &j_ab_xx, const double &n_ab_xx, const double &t_ab_zz, const double &j_ab_zz, const double &n_ab_zz, const int &desiredzeros, const int &seedvalue ) {

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
        
	T_a_xx = t_a_xx;  
	J_a_xx = j_a_xx;  
	N_a_xx = n_a_xx;
	
	T_a_zz = t_a_zz;
	J_a_zz = j_a_zz;
	N_a_zz = n_a_zz;

        T_a_xz = t_a_xz;
	J_a_xz = j_a_xz;
	N_a_xz = n_a_xz;

	T_b_xx = t_b_xx;
	J_b_xx = j_b_xx;
	N_b_xx = n_b_xx;

	T_b_zz = t_b_zz;
	J_b_zz = j_b_zz;
	N_b_zz = n_b_zz;

	T_b_xz = t_b_xz;
	J_b_xz = j_b_xz;
	N_b_xz = n_b_xz;

	T_ab_xx = t_ab_xx;
	J_ab_xx = j_ab_xx;
	N_ab_xx = n_ab_xx;

	T_ab_zz = t_ab_zz;
	J_ab_zz = j_ab_zz;
	N_ab_zz = n_ab_zz;

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

        // -------- vertical bonds, interlayer between a and b, zz orbital
	// 1------5-------9------13
	// |      |       |      |
	// 0------4-------8------12
	
	cout << " vertical bonds, interlayer between a and b, zz orbital " << endl;
	for( int i=0; i < N_x; i++ ){
	         
	       site_1 = N_y * N_u * i; site_2 = site_1 + 1;
               
	       cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl; 

	       Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
               Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
               Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

               Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_ab_zz;
               Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_ab_zz;
               Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_ab_zz;
	}


        // -------- vertical bonds, interlayer between a and b, xx orbital
        // 3------7-------11------15
        // |      |       |      |
        // 2------6-------10------14
        
        cout << " vertical bonds, interlayer between a and b, xx orbital " << endl;
        for( int i=0; i < N_x; i++ ){
        
               site_1 = N_y * N_u * i + N_u; site_2 = site_1 + 1;

               cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;

               Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
               Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
               Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

               Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_ab_xx;
               Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_ab_xx;
               Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_ab_xx;
        }
        //============================================================================================= 
        // ---------[y 方向 up layer a] ------------
	//    3-----7-----11----15
	//   /     /     /     /
	//   1----5------9----13
         cout << "----[up layer a, x orbital to z orbital] ---" << endl;
        for(int i=0; i<N_x; i++){
                site_1 = N_y * N_u * i + 1; site_2 = site_1 + 2;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;
                
                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
		Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
		Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

		Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_a_xz;
		Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_a_xz;
		Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_a_xz;

        }
        
        // ---------[x 方向 up layer a, z orbital] ------------
        //     ----- -----  ----
        //   /     /     /     /
        //   1----5------9----13
         cout << "----[up layer a, z orbital] ---" << endl;
        for(int i=0; i<N_x-1; i++){
                site_1 = N_y * N_u * i + 1; site_2 = site_1 + N_u * N_y;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;

                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
                Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
                Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

                Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_a_zz;
                Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_a_zz;
                Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_a_zz;

        }

        // ---------[x 方向 up layer a, x orbital] ------------
        //    3-----7-----11----15
        //   /     /     /     /
        //   *----*------*----*
         cout << "----[up layer a, x orbital] ---" << endl;
        for(int i=0; i<N_x-1; i++){
                site_1 = N_y * N_u * i + 3; site_2 = site_1 + N_u * N_y;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;

                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
                Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
                Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

                Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_a_xx;
                Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_a_xx;
                Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_a_xx;

        }

        //=======================================================================================================
        // ---------[y 方向 down layer b] ------------
        //    2-----6-----10----14
        //   /     /     /     /
        //   0----4------8----12
         cout << "----[down layer b, x orbital to z orbital] ---" << endl;
        for(int i=0; i<N_x; i++){
                site_1 = N_y * N_u * i; site_2 = site_1 + 2;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;

                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
                Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
                Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

                Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_b_xz;
                Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_b_xz;
                Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_b_xz;

        }

        // ---------[x 方向 down layer b, z orbital] ------------
        //     ----- -----  ----
        //   /     /     /     /
        //   0----4-----8----12
         cout << "----[down layer b, z orbital] ---" << endl;
        for(int i=0; i<N_x-1; i++){
                site_1 = N_y * N_u * i; site_2 = site_1 + N_u * N_y;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;

                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
                Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
                Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

                Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_b_zz;
                Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_b_zz;
                Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_b_zz;

        }
	//---------[x 方向 down layer b, x orbital] ----------------
	//    2-----6-----10----14
        //   /     /     /     /
        //   *----*------*----*
         cout << "----[down layer b, x orbital] ---" << endl;
        for(int i=0; i<N_x-1; i++){
                site_1 = N_y * N_u * i + 2; site_2 = site_1 + N_u * N_y;
                cout << "site_1: " << site_1 << ", site_2: " << site_2 << endl;

                Table_T[ Total_N * site_1 + site_2 ] = Table_T[ Total_N * site_2 + site_1 ] = 1;
                Table_J[ Total_N * site_1 + site_2 ] = Table_J[ Total_N * site_2 + site_1 ] = 1;
                Table_N[ Total_N * site_1 + site_2 ] = Table_N[ Total_N * site_2 + site_1 ] = 1;

                Interaction_T[ Total_N * site_1 + site_2 ] = Interaction_T[ Total_N * site_2 + site_1 ] = T_b_xx;
                Interaction_J[ Total_N * site_1 + site_2 ] = Interaction_J[ Total_N * site_2 + site_1 ] = J_b_xx;
                Interaction_N[ Total_N * site_1 + site_2 ] = Interaction_N[ Total_N * site_2 + site_1 ] = N_b_xx;

        }
        //-------------------------------------------------------------------------------------------

	
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
