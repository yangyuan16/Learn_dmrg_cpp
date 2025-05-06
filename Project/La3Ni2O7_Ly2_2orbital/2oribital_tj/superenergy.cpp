#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<iomanip>
#include<time.h>
#include<math.h>
#include<gsl/gsl_sf_coupling.h>
using namespace std;

#include"sub.h"
#include"super.h"
#include"superenergy.h"
#include"common.h"
//=================================================BLAS ROUTINES=================================================
extern "C" {
void daxpy_(const int *n, double *alpha, double *x, const int *incx, double *y, const int *incy);

void dgemm_(char *transa, char *transb, const int *m, const int *n, const int *k, const double *alpha, double *a, const int *lda, double *b, const int *ldb, const double *beta, double *c, const int *ldc);
}

//==============================Diagonalizing the Hamiltonian in infinite sweep===================================
SuperEnergy::SuperEnergy(Parameter &para, Super &sup_space, Super &sup, const double &lan_precision):Conjugate(sup.Dim) {

        if(sup.Dim<300)  
		SupConjugate(sup_space, sup);

        else {
                Davidson_Parameter(0, sup, lan_precision);
                Davidson(para, sup_space, sup);
        }

}

//================================================================================================================
//      			Diagonalization Method: 1.Conjugate method, 2.Davidson method
//================================================================================================================
//==============================1.Diagonalizing the Hamiltonian by the conjugate method===========================
//================================================================================================================
inline void SuperEnergy::SupConjugate(Super &sup_space, Super &sup) {
        int j;
        for(j=0; j<500; j++) {  //"500" see conjugate::abc_2(iter)
                if(j==0)  sup.H_V(f0, f1, sup.Dim) ; // f1 = H f0 //f1 = 0 ,f0 = sqrt(1./Dim) ;
                if(abc_2(j)) break;
                sup.H_V(f2, f3, sup.Dim); // f3 = H f2
                abc_4();
        }

        eigenvalue=eng;
//      eigenvalue_excited=eng;
        sup.NormalizedCopy(f0, sup_space.WaveFunction);

        for(int i=0; i<sup.Dim; i++)
                sup_space.WaveFunction_block[sup.Table_1to2_Num[i]][sup.Table_1to2_Site[i]]=sup_space.WaveFunction[i];

        cout<<"\n Energy="<<setprecision(10)<<eigenvalue<<endl;

}

//================================================================================================================
//==============================2.Diagonalizing the Hamiltonian by the Davidson method============================
//================================================================================================================
inline void SuperEnergy::Davidson_Parameter(const int &sign, Super &sup, const double &lan_precision) {

        n=sup.Dim;
        eigennum=1;
        input=0;
        maxiteration=800;

        if( sign == 0 ) { minibase=24; totalbase=36; precision=lan_precision; }
        else if( sign == 1 ) { minibase=12;  totalbase=24; precision=lan_precision; }

        EIGS=new double [10];

        X=new double [n];

        if(sign==0) {

                double sum=0.0;

                for(int j=0; j<n; j++) {
                        double s=drand48();
                        sum+=s*s;
                        X[j]=s;
                }

                sum=1.0/sqrt(sum);

                for(int j=0; j<n; j++)
                        X[j]*=sum;

        }

}

//===============================================================================================================
inline void SuperEnergy::Davidson(Parameter &para, Super &sup_space, Super &sup) {

        FILE *e=fopen("e/iteration_step", "a+");
        fprintf(e, "%d\t", para.total_site);
        fclose(e);

        sup.diagonalizeLSMatrix(n, eigennum, precision, totalbase, minibase, maxiteration, input, X, EIGS);

        sup.NormalizedCopy(&X[0], sup_space.WaveFunction);
//      sup.NormalizedCopy(&X[n], sup_space.WaveFunction_excited);

        eigenvalue=EIGS[0];

        cout<<"\n E(0)="<<EIGS[0]<<endl;//"\t E(1)="<<EIGS[1]<<endl;

        for(int i=0; i<sup.Dim; i++) {
                sup_space.WaveFunction_block[sup.Table_1to2_Num[i]][sup.Table_1to2_Site[i]]=sup_space.WaveFunction[i];
//                sup_space.WaveFunction_excited_block[sup.Table_1to2_Num[i]][sup.Table_1to2_Site[i]]=sup_space.WaveFunction_excited[i];                
//		cout<<"\n"<<sup_space.WaveFunction[i]<<endl;
        }

        delete [] EIGS;         delete [] X;

}

//==============================Diagonalizing the Hamiltonian in finite sweep====================================
SuperEnergy::SuperEnergy(Parameter &para, Super &sup_space, Super &sup, char &sign, const double &lan_precision):Conjugate(sup.Dim) {

  //------Wave function transformations and initial guess wave function for diagonalization----------------------
        trans_N = 'N';    trans_T = 'T';
        alpha = 1.0;      beta = 0.0;

        if( sign == 'L' )       Initialtrialfunction_Right_to_Left(para, sup_space);
        else if( sign == 'R' )  Initialtrialfunction_Left_to_Right(para, sup_space);

  //------DPJDREVCOM diagonalization in finite sweep: using initial guess wave function--------------------------
        Davidson_Parameter(1, sup, lan_precision);

        for(int i=0; i<sup.BlockNumber_for_TargetBlock; i++)
        for(int j=0; j<sup.Dim_block[i]; j++) {
                X[sup.Table_2to1[i][j]]=sup_space.WaveFunction_block[i][j];
//              X[sup.Dim+sup.Table_2to1[i][j]]=sup_space.WaveFunction_excited_block[i][j];
	}

//        for(int i=0; i<sup.BlockNumber_for_TargetBlock; i++)
//        for(int j=0; j<sup.Dim_block[i]; j++)
//		cout<<"\n"<<sup.Table_2to1[i][j]<<"\t"<<X[sup.Table_2to1[i][j]];

        Davidson(para, sup_space, sup);

}

//=================Diagonalizing the Hamiltonian in finite sweep for the start point in each sweep===============
SuperEnergy::SuperEnergy(Parameter &para, Super &sup_space, Super &sup, const int &n, const double &lan_precision):Conjugate(sup.Dim) {

        Davidson_Parameter(1, sup, lan_precision);

        if(n==1) {
                FILE *fp=fopen(Combine(Combine("truncated_wave_function/", n), sup.sys->TotSiteNo+1), "rb");
                fread(&X[0], sizeof(double), sup.Dim, fp);
		fclose(fp);
	}

        else if(n==2) {
                FILE *fp=fopen(Combine(Combine("truncated_wave_function/", n), sup.env->TotSiteNo+1), "rb");
                fread(&X[0], sizeof(double), sup.Dim, fp);
                fclose(fp);
        }

        Davidson(para, sup_space, sup);

}

//===========================================
//Wave function transformations and initializing the trial function for diagonalization
//============================================
inline void SuperEnergy::Initialtrialfunction_Left_to_Right(Parameter &para, Super &sup_space) {

//----Second transformation from (sys, ns+envtrun) to (sys+ns, envtrun)
//Transform (Block_sys[n+1], S[n+2]+Block_env[n+3]) to (Block_sys[n+1]+S[n+2], Block_env[n+3])
  //------Transformed wave function(Block_sys[n+1]+S[n+2], Block_env[n+3])
    //------Read envtrun:Block_env[n+3]
        FILE *fp=fopen(Combine(Combine("truncated_density_eigenvector/", 2), sup_space.env_space->TotSiteNo+1), "rb");
  //---------
        fread(&truncated_Block_Num_hole, sizeof(int), 1, fp);
  //---------
        truncated_Num_hole_block=new int [truncated_Block_Num_hole];
        truncated_Num_J_hole_block=new int [truncated_Block_Num_hole];

        fread(truncated_Num_hole_block, sizeof(int), truncated_Block_Num_hole, fp);
        fread(truncated_Num_J_hole_block, sizeof(int), truncated_Block_Num_hole, fp);
  //---------
        truncated_Value_J_block=new int * [truncated_Block_Num_hole];
        truncated_Dim_J_block=new int * [truncated_Block_Num_hole];
        truncated_density_dim=new int * [truncated_Block_Num_hole];

        for(int i=0; i<truncated_Block_Num_hole; i++) {

                truncated_Value_J_block[i]=new int [truncated_Num_J_hole_block[i]];
                truncated_Dim_J_block[i]=new int [truncated_Num_J_hole_block[i]];
                truncated_density_dim[i]=new int [truncated_Num_J_hole_block[i]];

                for(int j=0; j<truncated_Num_J_hole_block[i]; j++) {

                        truncated_Value_J_block[i][j]=0;  
			truncated_Dim_J_block[i][j]=0;  
			truncated_density_dim[i][j]=0;

                }

        }

        for(int i=0; i<truncated_Block_Num_hole; i++) {

                fread(truncated_Value_J_block[i], sizeof(int), truncated_Num_J_hole_block[i], fp);
                fread(truncated_Dim_J_block[i], sizeof(int), truncated_Num_J_hole_block[i], fp);
                fread(truncated_density_dim[i], sizeof(int), truncated_Num_J_hole_block[i], fp);

        }

  //--------
        truncated_density_eigenvector=new double ** [truncated_Block_Num_hole];

        for(int i=0; i<truncated_Block_Num_hole; i++) {

                truncated_density_eigenvector[i]=new double * [truncated_Num_J_hole_block[i]];

                for(int j=0; j<truncated_Num_J_hole_block[i]; j++) {

                        truncated_density_eigenvector[i][j]=new double [truncated_density_dim[i][j]];

                        for(int k=0; k<truncated_density_dim[i][j]; k++)

                                truncated_density_eigenvector[i][j][k]=0.0;

                }

        }

        for(int i=0; i<truncated_Block_Num_hole; i++)
        for(int j=0; j<truncated_Num_J_hole_block[i]; j++)

                fread(truncated_density_eigenvector[i][j], sizeof(double), truncated_density_dim[i][j], fp);

  //--------	//Why we read these??????????????????????????
        truncated_Old_hole=new int [truncated_Block_Num_hole];
        fread(truncated_Old_hole, sizeof(int), truncated_Block_Num_hole, fp);

  //--------
        truncated_Old_J=new int * [truncated_Block_Num_hole];

        for(int i=0; i<truncated_Block_Num_hole; i++) {

                truncated_Old_J[i]=new int [truncated_Num_J_hole_block[i]];

                for(int j=0; j<truncated_Num_J_hole_block[i]; j++)

                        truncated_Old_J[i][j]=0.0;

        }

        for(int i=0; i<truncated_Block_Num_hole; i++)

                fread(truncated_Old_J[i], sizeof(int), truncated_Num_J_hole_block[i], fp);

        fclose(fp);

    //--space for transformed wave function:(Block_sys[n+1]+S[n+2], Block_env[n+3])
    //Block_sys[n+1]+S[n+2]==sup.Block_sysnew
        truncated_block_number=0;

        for(int n_1=0; n_1<sup_space.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<truncated_Block_Num_hole; n_2++)
        if(sup_space.sysnew_space->Num_hole_block[n_1]+truncated_Num_hole_block[n_2]==para.Total_h) {

                for(int j_sysnew=0; j_sysnew<sup_space.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envtrun=0; j_envtrun<truncated_Num_J_hole_block[n_2]; j_envtrun++) {//Block_env[n+3]

                        J_min=abs(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]-truncated_Value_J_block[n_2][j_envtrun]);
                        J_max=(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]+truncated_Value_J_block[n_2][j_envtrun]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if(para.Total_J==(J_min+2*n)) {

                                for(int j_s=0; j_s<3; j_s++)
                                if(sup_space.sysnew_space->Hole_blockOld[n_1][j_s][j_sysnew]!=-1 && sup_space.sysnew_space->J_blockOld[n_1][j_s][j_sysnew]!=-1)

                                        truncated_block_number++;

                        }

                }

        }

    //------Create space for the angular numbers of the transformed blocks
        H_trun=new int [truncated_block_number];
        J_trun=new int [truncated_block_number];//J_trun=Block_env[n+3]:truncated-(Sys_Number_Jn,Sys_Value_Jn,..)

        H_old=new int [truncated_block_number];
        J_old=new int [truncated_block_number];//J_old=Block_sys[n+1]:sup.sys->(Sys_Number_Jn,Sys_Value_Jn,...)

        H_new=new int [truncated_block_number];
        J_new=new int [truncated_block_number];//J_new=Block_sys[n+1]+S[n+2]:sup.sysnew->(Sys_Number_Jn,...)

        truncated_block_dim=new int [truncated_block_number];//truncated_block_dim=Dim.Block_env[n+3]*Dim.Block_sys[n+1]

        for(int i=0; i<truncated_block_number; i++) {

                H_trun[i]=0;  J_trun[i]=0;  
		H_old[i]=0;  J_old[i]=0;  
		H_new[i]=0;  J_new[i]=0;  
		truncated_block_dim[i]=0;

        }

    //------Initialize the values of the above angular numbers
        index=0;

        for(int n_1=0; n_1<sup_space.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<truncated_Block_Num_hole; n_2++)
        if(sup_space.sysnew_space->Num_hole_block[n_1]+truncated_Num_hole_block[n_2]==para.Total_h) {

                for(int j_sysnew=0; j_sysnew<sup_space.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)//Block_sys[n+1]+S[n+2]
                for(int j_envtrun=0; j_envtrun<truncated_Num_J_hole_block[n_2]; j_envtrun++) {//Block_env[n+3]

                        J_min=abs(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]-truncated_Value_J_block[n_2][j_envtrun]);
                        J_max=(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]+truncated_Value_J_block[n_2][j_envtrun]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if(para.Total_J==(J_min+2*n)) {

                                for(int j_s=0; j_s<3; j_s++)
                                if((oldH_sys=sup_space.sysnew_space->Hole_blockOld[n_1][j_s][j_sysnew])!=-1 && (oldJ_sys=sup_space.sysnew_space->J_blockOld[n_1][j_s][j_sysnew])!=-1) {
                                        H_trun[index]=n_2;  J_trun[index]=j_envtrun;

                                        H_old[index]=oldH_sys;  J_old[index]=oldJ_sys;

                                        H_new[index]=n_1;  J_new[index]=j_sysnew;

                                        truncated_block_dim[index]=sup_space.sys_space->Dim_J_block[oldH_sys][oldJ_sys]*truncated_Dim_J_block[n_2][j_envtrun];

					index++;

                                }

			}

		}

	}

    //------Create space for the transformed wave function
        truncated_wave_function=new double * [truncated_block_number];

        for(int i=0; i<truncated_block_number; i++) {

                truncated_wave_function[i]=new double [truncated_block_dim[i]];

                for(int j=0; j<truncated_block_dim[i]; j++)

                        truncated_wave_function[i][j]=(double) 0;

        }

  //----Wave function that to be transformed: read from the last step
   //(Block_sys[n+1], Block_env[n+3]+S[n+2]), obtained from the last step
   //---space for the initial wave function and read it
        FILE *fw=fopen(Combine(Combine("truncated_wave_function/", 1), sup_space.sys_space->TotSiteNo), "rb");

        fread(&untruncated_block_number, sizeof(int), 1, fw);

        H_sys_untrun=new int [untruncated_block_number];
        J_sys_untrun=new int [untruncated_block_number];
        H_env_untrun=new int [untruncated_block_number];
        J_env_untrun=new int [untruncated_block_number];
        H_envnew_untrun=new int [untruncated_block_number];
        J_envnew_untrun=new int [untruncated_block_number];
        untruncated_block_dim=new int [untruncated_block_number];

        fread(H_sys_untrun, sizeof(int), untruncated_block_number, fw);
        fread(J_sys_untrun, sizeof(int), untruncated_block_number, fw);
        fread(H_env_untrun, sizeof(int), untruncated_block_number, fw);
        fread(J_env_untrun, sizeof(int), untruncated_block_number, fw);
        fread(H_envnew_untrun, sizeof(int), untruncated_block_number, fw);
        fread(J_envnew_untrun, sizeof(int), untruncated_block_number, fw);
        fread(untruncated_block_dim, sizeof(int), untruncated_block_number, fw);

        untruncated_wave_function=new double * [untruncated_block_number];

        for(int i=0; i<untruncated_block_number; i++) {

                untruncated_wave_function[i]=new double [untruncated_block_dim[i]];

                for(int j=0; j<untruncated_block_dim[i]; j++)

                        untruncated_wave_function[i][j]=0.0;

        }

        for(int i=0; i<untruncated_block_number; i++)

                fread(untruncated_wave_function[i], sizeof(double), untruncated_block_dim[i], fw);

        fclose(fw);

    //--Read in the block spin value of Block_env[n+3]+S[n+2]
        FILE *fr=fopen(Combine(Combine("new_block/", 2), sup_space.envnew_space->TotSiteNo+1), "rb");

        //---
        fread(&new_Block_Num_hole, sizeof(int), 1, fr);

        //---
        new_Num_hole_block=new int [new_Block_Num_hole];
        new_Num_J_hole_block=new int [new_Block_Num_hole];

        fread(new_Num_hole_block, sizeof(int), new_Block_Num_hole, fr);
        fread(new_Num_J_hole_block, sizeof(int), new_Block_Num_hole, fr);

        //---
        new_Value_J_block=new int * [new_Block_Num_hole];
        new_Dim_J_block=new int * [new_Block_Num_hole];

        for(int i=0; i<new_Block_Num_hole; i++) {

                new_Value_J_block[i]=new int [new_Num_J_hole_block[i]];
                new_Dim_J_block[i]=new int [new_Num_J_hole_block[i]];

                for(int j=0; j<new_Num_J_hole_block[i]; j++) {

                        new_Value_J_block[i][j]=0.0;  new_Dim_J_block[i][j]=0.0;

                }

        }

        for(int i=0; i<new_Block_Num_hole; i++) {

                fread(new_Value_J_block[i], sizeof(int), new_Num_J_hole_block[i], fr);
                fread(new_Dim_J_block[i], sizeof(int), new_Num_J_hole_block[i], fr);

        }

        //---
        new_Hole_blockOld=new int ** [new_Block_Num_hole];
        new_J_blockOld=new int ** [new_Block_Num_hole];
        new_Start=new int ** [new_Block_Num_hole];

        for(int i=0; i<new_Block_Num_hole; i++) {

                new_Hole_blockOld[i]=new int * [3];
                new_J_blockOld[i]=new int * [3];
                new_Start[i]=new int * [3];

                for(int j=0; j<3; j++) {

                        new_Hole_blockOld[i][j]=new int [new_Num_J_hole_block[i]];
                        new_J_blockOld[i][j]=new int [new_Num_J_hole_block[i]];
                        new_Start[i][j]=new int [new_Num_J_hole_block[i]];

                        for(int k=0; k<new_Num_J_hole_block[i]; k++) {

                                new_Hole_blockOld[i][j][k]=0.0;
                                new_J_blockOld[i][j][k]=0.0;
                                new_Start[i][j][k]=0.0;

                        }

                }

        }

        for(int i=0; i<new_Block_Num_hole; i++)
        for(int j=0; j<3; j++) {

                fread(new_Hole_blockOld[i][j], sizeof(int), new_Num_J_hole_block[i], fr);
                fread(new_J_blockOld[i][j], sizeof(int), new_Num_J_hole_block[i], fr);
                fread(new_Start[i][j], sizeof(int), new_Num_J_hole_block[i], fr);

        }

        fclose(fr);
 
    //------The 6-j coefficient for the basis transformation in this step
        six_j_basis_transformation=new double * [truncated_block_number];//truncated

        for(int i=0; i<truncated_block_number; i++) {

                index=0;

                for(int j=0; j<untruncated_block_number; j++)//untruncated
                if(H_sys_untrun[j]==H_old[i] && J_sys_untrun[j]==J_old[i] && H_env_untrun[j]==H_trun[i] && J_env_untrun[j]==J_trun[i])

                        index++;

                six_j_basis_transformation[i]=new double [index];

                for(int n=0; n<index; n++)

                        six_j_basis_transformation[i][n] = 0.0;

                index=0;

                for(int j=0; j<untruncated_block_number; j++)
                if(H_sys_untrun[j]==H_old[i] && J_sys_untrun[j]==J_old[i] && H_env_untrun[j]==H_trun[i] && J_env_untrun[j]==J_trun[i]) {//Block_sys[n+1] and Block_env[n+3]

                        alpha = sqrt( ( sup_space.sysnew_space->Value_J_block[H_new[i]][J_new[i]] + 1.0 ) * ( new_Value_J_block[H_envnew_untrun[j]][J_envnew_untrun[j]] + 1.0 ) );

                        if( (sup_space.sysnew_space->Num_hole_block[H_new[i]] - sup_space.sys_space->Num_hole_block[H_old[i]]) == 1 ) {  // a hole on the site

                                six_j_basis_transformation[i][index] += alpha * gsl_sf_coupling_6j( sup_space.sys_space->Value_J_block[H_old[i]][J_old[i]], 0, sup_space.sysnew_space->Value_J_block[H_new[i]][J_new[i]], truncated_Value_J_block[H_trun[i]][J_trun[i]], para.Total_J, new_Value_J_block[H_envnew_untrun[j]][J_envnew_untrun[j]] );

                                if( ( ( sup_space.sys_space->Value_J_block[H_old[i]][J_old[i]] + para.Total_J + truncated_Value_J_block[H_trun[i]][J_trun[i]] ) / 2 ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

                                if( ( ( truncated_Value_J_block[H_trun[i]][J_trun[i]] - new_Value_J_block[H_envnew_untrun[j]][J_envnew_untrun[j]] ) / 2 ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];


                                index++;

			}

			else if( sup_space.sysnew_space->Num_hole_block[H_new[i]] == sup_space.sys_space->Num_hole_block[H_old[i]] ) { // on n+2 site, we have a spin !!!

                                six_j_basis_transformation[i][index] += alpha * gsl_sf_coupling_6j( sup_space.sys_space->Value_J_block[H_old[i]][J_old[i]], para.Spin, sup_space.sysnew_space->Value_J_block[H_new[i]][J_new[i]], truncated_Value_J_block[H_trun[i]][J_trun[i]], para.Total_J, new_Value_J_block[H_envnew_untrun[j]][J_envnew_untrun[j]] );

                                if( ( ( sup_space.sys_space->Value_J_block[H_old[i]][J_old[i]] + para.Spin + para.Total_J + truncated_Value_J_block[H_trun[i]][J_trun[i]] ) / 2 ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

                                if( ( ( para.Spin + truncated_Value_J_block[H_trun[i]][J_trun[i]] - new_Value_J_block[H_envnew_untrun[j]][J_envnew_untrun[j]] ) / 2 ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

                                if( ( sup_space.env_space->TotSiteNo + 1 - truncated_Num_hole_block[H_trun[i]] ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

                                index++;

			}

		}

	}

  //------Transform wave function from (Block_sys[n+1], Block_env[n+3]+S[n+2]) to (Block_sys[n+1]+S[n+2], Block_env[n+3]): multiply the 6-j coefficient obtained above
 
        int inc=1;              //increament index in BLAS subroutine daxpy_()

        for(int i=0; i<truncated_block_number; i++) {

                index=0;

                for(int j=0; j<untruncated_block_number; j++)
                if(H_sys_untrun[j]==H_old[i] && J_sys_untrun[j]==J_old[i] && H_env_untrun[j]==H_trun[i] && J_env_untrun[j]==J_trun[i]) {

                        daxpy_(&truncated_block_dim[i], &six_j_basis_transformation[i][index], untruncated_wave_function[j], &inc, truncated_wave_function[i], &inc);

                        index++;

                }

        }

//------Delete the untruncated block quantities to reuse them for the next transformation step-------------------
        for(int i=0; i<untruncated_block_number; i++) {
                delete [] untruncated_wave_function[i];
        }
        delete [] untruncated_wave_function;

        delete [] H_envnew_untrun;  delete [] J_envnew_untrun;  delete [] H_sys_untrun;  delete [] J_sys_untrun;
        delete [] H_env_untrun;  delete [] J_env_untrun;  delete [] untruncated_block_dim;

//------Third wave function transformation: from (Block_sys[n+1]+S[n+2], Block_env[n+3]) to (Block_sys[n+1]+S[n+2]
//, S[n+3]+Block_env[n+4]), which is the last step of transformation to provide the initial guess function for 
//diagonalization. The space for the transformed wave function has been created in Class Super!------------------
  //------Create space for the transformed wave function (use untruncated block quantities)----------------------
   //------Find block number------------------------------------------------------------------------------------
        untruncated_block_number=0;

        for(int n_1=0; n_1<sup_space.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup_space.envnew_space->Block_Num_hole; n_2++)
        if(sup_space.sysnew_space->Num_hole_block[n_1]+sup_space.envnew_space->Num_hole_block[n_2]==para.Total_h) {

                for(int j_sysnew=0; j_sysnew<sup_space.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envnew=0; j_envnew<sup_space.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                        J_min=abs(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]-sup_space.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_max=(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]+sup_space.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if(para.Total_J==(J_min+2*n)) {

                                for(int j_s=0; j_s<3; j_s++)
                                if(sup_space.sysnew_space->Hole_blockOld[n_1][j_s][j_sysnew]!=-1 && sup_space.sysnew_space->J_blockOld[n_1][j_s][j_sysnew]!=-1)
                                        untruncated_block_number++;

                        }

                }

        }

    //------Create space for angular momentum numbers------------------------------------------------------------
        H_sysnew_untrun=new int [untruncated_block_number];
        J_sysnew_untrun=new int [untruncated_block_number];
        H_envnew_untrun=new int [untruncated_block_number];
        J_envnew_untrun=new int [untruncated_block_number];
        H_sys_untrun=new int [untruncated_block_number];
        J_sys_untrun=new int [untruncated_block_number];
        untruncated_block_dim=new int [untruncated_block_number];

        for(int i=0; i<untruncated_block_number; i++) {
                H_sysnew_untrun[i]=0;  J_sysnew_untrun[i]=0;   H_envnew_untrun[i]=0;  J_envnew_untrun[i]=0;
                H_sys_untrun[i]=0;  J_sys_untrun[i]=0;  untruncated_block_dim[i]=0;
        }

    //------Initialize the above angular momentum numbers--------------------------------------------------------
        index=0;

        for(int n_1=0; n_1<sup_space.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup_space.envnew_space->Block_Num_hole; n_2++)
        if(sup_space.sysnew_space->Num_hole_block[n_1]+sup_space.envnew_space->Num_hole_block[n_2]==para.Total_h) {

                for(int j_sysnew=0; j_sysnew<sup_space.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envnew=0; j_envnew<sup_space.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                        J_min=abs(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]-sup_space.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_max=(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]+sup_space.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if(para.Total_J==(J_min+2*n)) {

                                for(int j_s=0; j_s<3; j_s++)
                                if((oldH_sys=sup_space.sysnew_space->Hole_blockOld[n_1][j_s][j_sysnew])!=-1 && (oldJ_sys=sup_space.sysnew_space->J_blockOld[n_1][j_s][j_sysnew])!=-1) {
                                        H_sysnew_untrun[index]=n_1;  J_sysnew_untrun[index]=j_sysnew;
                                        H_envnew_untrun[index]=n_2;  J_envnew_untrun[index]=j_envnew;
                                        H_sys_untrun[index]=oldH_sys;  J_sys_untrun[index]=oldJ_sys;
                                        untruncated_block_dim[index++]=sup_space.sys_space->Dim_J_block[oldH_sys][oldJ_sys]*sup_space.envnew_space->Dim_J_block[n_2][j_envnew];
                                }

                        }

                }

	} 

    //------Create space for untruncated_wave_function-----------------------------------------------------------
        untruncated_wave_function=new double * [untruncated_block_number];

        for(int i=0; i<untruncated_block_number; i++) {

                untruncated_wave_function[i]=new double [untruncated_block_dim[i]];

                for(int j=0; j<untruncated_block_dim[i]; j++)
                        untruncated_wave_function[i][j]=(double) 0;

        }

  //------Wave function transformation from (Block_sys[n+1]+S[n+2], Block_env[n+3]) to (Block_sys[n+1]+S[n+2], 
  //S[n+3]+Block_env[n+4]) by matrix-matrix multiplication------------------------------------------------------- 
        for(int i=0; i<untruncated_block_number; i++)
        for(int j=0; j<truncated_block_number; j++)
        if(H_new[j]==H_sysnew_untrun[i] && J_new[j]==J_sysnew_untrun[i] && truncated_Num_hole_block[H_trun[j]]==sup_space.envnew_space->Num_hole_block[H_envnew_untrun[i]] && truncated_Value_J_block[H_trun[j]][J_trun[j]]==sup_space.envnew_space->Value_J_block[H_envnew_untrun[i]][J_envnew_untrun[i]] && H_old[j]==H_sys_untrun[i] && J_old[j]==J_sys_untrun[i]) {
                dgemm_(&trans_N, &trans_T, &sup_space.sys_space->Dim_J_block[H_sys_untrun[i]][J_sys_untrun[i]], &sup_space.envnew_space->Dim_J_block[H_envnew_untrun[i]][J_envnew_untrun[i]], &truncated_Dim_J_block[H_trun[j]][J_trun[j]], &alpha, truncated_wave_function[j], &sup_space.sys_space->Dim_J_block[H_sys_untrun[i]][J_sys_untrun[i]], truncated_density_eigenvector[H_trun[j]][J_trun[j]], &sup_space.envnew_space->Dim_J_block[H_envnew_untrun[i]][J_envnew_untrun[i]], &beta, untruncated_wave_function[i], &sup_space.sys_space->Dim_J_block[H_sys_untrun[i]][J_sys_untrun[i]]);
        }
  
  //------Read from untruncated_wave_function to sup.WaveFunction_block------------------------------------------
        for(int i=0; i<sup_space.BlockNumber_for_TargetBlock; i++)
        for(int j=0; j<untruncated_block_number; j++)
        if(sup_space.H_sysnew[i]==H_sysnew_untrun[j] && sup_space.J_sysnew[i]==J_sysnew_untrun[j] && sup_space.H_envnew[i]==H_envnew_untrun[j] && sup_space.J_envnew[i]==J_envnew_untrun[j] && sup_space.H_sys[i]==H_sys_untrun[j] && sup_space.J_sys[i]==J_sys_untrun[j]) {

                for(int j_e=0; j_e<3; j_e++)
                if((oldH_env=sup_space.envnew_space->Hole_blockOld[H_envnew_untrun[j]][j_e][J_envnew_untrun[j]])!=-1 && (oldJ_env=sup_space.envnew_space->J_blockOld[H_envnew_untrun[j]][j_e][J_envnew_untrun[j]])!=-1 && oldH_env==sup_space.H_env[i] && oldJ_env==sup_space.J_env[i]) {

                        for(int a_env=0; a_env<sup_space.env_space->Dim_J_block[sup_space.H_env[i]][sup_space.J_env[i]]; a_env++)
                        for(int a_sys=0; a_sys<sup_space.sys_space->Dim_J_block[sup_space.H_sys[i]][sup_space.J_sys[i]]; a_sys++) {

                                position_old=(sup_space.envnew_space->Start[H_envnew_untrun[j]][j_e][J_envnew_untrun[j]]+a_env)*sup_space.sys_space->Dim_J_block[sup_space.H_sys[i]][sup_space.J_sys[i]]+a_sys;
                                position_new=a_env*sup_space.sys_space->Dim_J_block[sup_space.H_sys[i]][sup_space.J_sys[i]]+a_sys;
                                sup_space.WaveFunction_block[i][position_new]=untruncated_wave_function[j][position_old];

                        }

                }

        }
 
//------Delete all the created space------------------------------------------------------------------------------
	//-----
        for(int i=0; i<untruncated_block_number; i++) {
                delete [] untruncated_wave_function[i];
        }
        delete [] untruncated_wave_function;

        delete [] H_sysnew_untrun;  delete [] J_sysnew_untrun;  
	delete [] H_envnew_untrun;  delete [] J_envnew_untrun;            
	delete [] H_sys_untrun;     delete [] J_sys_untrun;
        delete [] untruncated_block_dim;

	//------
        for(int i=0; i<truncated_block_number; i++)
                delete [] six_j_basis_transformation[i];
        delete [] six_j_basis_transformation;

	//------
        for(int i=0; i<new_Block_Num_hole; i++) {
                for(int j=0; j<3; j++) {
                        delete [] new_Hole_blockOld[i][j];  delete [] new_J_blockOld[i][j];  delete [] new_Start[i][j];
                }
                delete [] new_Hole_blockOld[i];  delete [] new_J_blockOld[i];  delete [] new_Start[i];
        }
        delete [] new_Hole_blockOld;  delete [] new_J_blockOld;  delete [] new_Start;

        for(int i=0; i<new_Block_Num_hole; i++) {
                delete [] new_Value_J_block[i];  delete [] new_Dim_J_block[i];
        }
        delete [] new_Value_J_block;  delete [] new_Dim_J_block;

        delete [] new_Num_hole_block;  delete [] new_Num_J_hole_block;

	//------
        for(int i=0; i<truncated_block_number; i++)
                delete [] truncated_wave_function[i];
        delete [] truncated_wave_function;

        delete [] H_trun;  delete [] J_trun;  
	delete [] H_old;   delete [] J_old;  
	delete [] H_new;   delete [] J_new;              
	delete [] truncated_block_dim;

	//------
        for(int i=0; i<truncated_Block_Num_hole; i++)
		delete [] truncated_Old_J[i];
	delete [] truncated_Old_J;

	delete [] truncated_Old_hole;	

        for(int i=0; i<truncated_Block_Num_hole; i++) {
                for(int j=0; j<truncated_Num_J_hole_block[i]; j++) {
                        delete [] truncated_density_eigenvector[i][j];
                }
                delete [] truncated_density_eigenvector[i];
	}
        delete [] truncated_density_eigenvector;

        for(int i=0; i<truncated_Block_Num_hole; i++) {
           delete [] truncated_Value_J_block[i];  delete [] truncated_Dim_J_block[i];  delete [] truncated_density_dim[i];
        }
        delete [] truncated_Value_J_block;  delete [] truncated_Dim_J_block;  delete [] truncated_density_dim;

        delete [] truncated_Num_hole_block;  delete [] truncated_Num_J_hole_block;

}

//=============================================================================
inline void SuperEnergy::Initialtrialfunction_Right_to_Left(Parameter &para, Super &sup_space) {

//------Second wave function transformation from (Block_sys[n]+S[n+1], Block_env[n+2]) to (Block_sys[n], S[n+1]+ Block_env[n+2])
  //--Transformed wave function for basis space: (Block_sys[n], S[n+1]+Block_env[n+2])
        //------Read systrun:Block_sys[n]----------------------------
        FILE *fp=fopen(Combine(Combine("truncated_density_eigenvector/", 1), sup_space.sys_space->TotSiteNo+1), "rb");
  //---------
        fread(&truncated_Block_Num_hole, sizeof(int), 1, fp);
  //---------
        truncated_Num_hole_block=new int [truncated_Block_Num_hole];
        truncated_Num_J_hole_block=new int [truncated_Block_Num_hole];

        fread(truncated_Num_hole_block, sizeof(int), truncated_Block_Num_hole, fp);
        fread(truncated_Num_J_hole_block, sizeof(int), truncated_Block_Num_hole, fp);
  //---------
        truncated_Value_J_block=new int * [truncated_Block_Num_hole];
        truncated_Dim_J_block=new int * [truncated_Block_Num_hole];
        truncated_density_dim=new int * [truncated_Block_Num_hole];
        for(int i=0; i<truncated_Block_Num_hole; i++) {
                truncated_Value_J_block[i]=new int [truncated_Num_J_hole_block[i]];
                truncated_Dim_J_block[i]=new int [truncated_Num_J_hole_block[i]];
                truncated_density_dim[i]=new int [truncated_Num_J_hole_block[i]];
                for(int j=0; j<truncated_Num_J_hole_block[i]; j++) {
                        truncated_Value_J_block[i][j]=0;  truncated_Dim_J_block[i][j]=0;  truncated_density_dim[i][j]=0;
                }
        }

        for(int i=0; i<truncated_Block_Num_hole; i++) {
                fread(truncated_Value_J_block[i], sizeof(int), truncated_Num_J_hole_block[i], fp);
                fread(truncated_Dim_J_block[i], sizeof(int), truncated_Num_J_hole_block[i], fp);
                fread(truncated_density_dim[i], sizeof(int), truncated_Num_J_hole_block[i], fp);
        }
  //--------
        truncated_density_eigenvector=new double ** [truncated_Block_Num_hole];
        for(int i=0; i<truncated_Block_Num_hole; i++) {
                truncated_density_eigenvector[i]=new double * [truncated_Num_J_hole_block[i]];
                for(int j=0; j<truncated_Num_J_hole_block[i]; j++) {
                        truncated_density_eigenvector[i][j]=new double [truncated_density_dim[i][j]];
                        for(int k=0; k<truncated_density_dim[i][j]; k++)
                                truncated_density_eigenvector[i][j][k]=0.0;
                }
        }

        for(int i=0; i<truncated_Block_Num_hole; i++)
        for(int j=0; j<truncated_Num_J_hole_block[i]; j++)
                fread(truncated_density_eigenvector[i][j], sizeof(double), truncated_density_dim[i][j], fp);

  //--------
        truncated_Old_hole=new int [truncated_Block_Num_hole];
        fread(truncated_Old_hole, sizeof(int), truncated_Block_Num_hole, fp);
  //--------
        truncated_Old_J=new int * [truncated_Block_Num_hole];
        for(int i=0; i<truncated_Block_Num_hole; i++) {
                truncated_Old_J[i]=new int [truncated_Num_J_hole_block[i]];
                for(int j=0; j<truncated_Num_J_hole_block[i]; j++)
                        truncated_Old_J[i][j]=0.0;
        }

        for(int i=0; i<truncated_Block_Num_hole; i++)
                fread(truncated_Old_J[i], sizeof(int), truncated_Num_J_hole_block[i], fp);

        fclose(fp);

 //Create space for the transformed wave function:(Block_sys[n], S[n+1]+Block_env[n+2])
   //S[n+1]+Block_env[n+2]==sup.Block_envnew
        truncated_block_number=0;

        for(int n_1=0; n_1<truncated_Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup_space.envnew_space->Block_Num_hole; n_2++)
        if(sup_space.envnew_space->Num_hole_block[n_2]+truncated_Num_hole_block[n_1]==para.Total_h) {

                for(int j_systrun=0; j_systrun<truncated_Num_J_hole_block[n_1]; j_systrun++)
                for(int j_envnew=0; j_envnew<sup_space.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {//Block_env[n+3]

                        J_min=abs(sup_space.envnew_space->Value_J_block[n_2][j_envnew]-truncated_Value_J_block[n_1][j_systrun]);
                        J_max=(sup_space.envnew_space->Value_J_block[n_2][j_envnew]+truncated_Value_J_block[n_1][j_systrun]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if(para.Total_J==(J_min+2*n)) {

                                for(int j_e=0; j_e<3; j_e++)
                                if(sup_space.envnew_space->Hole_blockOld[n_2][j_e][j_envnew]!=-1 && sup_space.envnew_space->J_blockOld[n_2][j_e][j_envnew]!=-1)
                                        truncated_block_number++;
                        }

                }

        }

    //------Create space for the angular numbers of the transformed blocks
        H_trun=new int [truncated_block_number];
        J_trun=new int [truncated_block_number];//J_trun=Block_sys[n]:truncated-(Sys_Number_Jn,Sys_Value_Jn,..)
        H_old=new int [truncated_block_number];
        J_old=new int [truncated_block_number];//J_old=Block_env[n+2]:sup.env->(Sys_Number_Jn,Sys_Value_Jn,...)
        H_new=new int [truncated_block_number];
        J_new=new int [truncated_block_number];//J_new=S[n+1]+Block_env[n+2]:sup.envnew->(Sys_Number_Jn,...)
        truncated_block_dim=new int [truncated_block_number];//truncated_block_dim=Dim.Block_env[n+2]*Dim.Block_sys[n]
        for(int i=0; i<truncated_block_number; i++) {
                H_trun[i]=0;  J_trun[i]=0;  H_old[i]=0;  J_old[i]=0;  H_new[i]=0;  J_new[i]=0;  truncated_block_dim[i]=0;
        }

    //------Initialize the values of the above angular numbers
        index=0;

        for(int n_1=0; n_1<truncated_Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup_space.envnew_space->Block_Num_hole; n_2++)
        if(sup_space.envnew_space->Num_hole_block[n_2]+truncated_Num_hole_block[n_1]==para.Total_h) {

                for(int j_systrun=0; j_systrun<truncated_Num_J_hole_block[n_1]; j_systrun++)
                for(int j_envnew=0; j_envnew<sup_space.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {//Block_env[n+3]

                        J_min=abs(sup_space.envnew_space->Value_J_block[n_2][j_envnew]-truncated_Value_J_block[n_1][j_systrun]);
                        J_max=(sup_space.envnew_space->Value_J_block[n_2][j_envnew]+truncated_Value_J_block[n_1][j_systrun]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if(para.Total_J==(J_min+2*n)) {

                                for(int j_e=0; j_e<3; j_e++)
                                if((oldH_env=sup_space.envnew_space->Hole_blockOld[n_2][j_e][j_envnew])!=-1 && (oldJ_env=sup_space.envnew_space->J_blockOld[n_2][j_e][j_envnew])!=-1) {

                                        H_trun[index]=n_1;  J_trun[index]=j_systrun;
                                        H_old[index]=oldH_env;  J_old[index]=oldJ_env;
                                        H_new[index]=n_2;  J_new[index]=j_envnew;
                                        truncated_block_dim[index++]=sup_space.env_space->Dim_J_block[oldH_env][oldJ_env]*truncated_Dim_J_block[n_1][j_systrun];

                                }

                        }

                }
 
	}

    //------Create space for the transformed wave function
        truncated_wave_function=new double * [truncated_block_number];

        for(int i=0; i<truncated_block_number; i++) {

                truncated_wave_function[i]=new double [truncated_block_dim[i]];

                for(int j=0; j<truncated_block_dim[i]; j++)
                        truncated_wave_function[i][j]=(double) 0;

        }

 //Wave function that to be transformed: read from the truncated wave function in the last step
 //(Block_sys[n]+S[n+1], Block_env[n+2]), obtained from the last step
 //Create space for initial wave function and read it from the truncated wave function
        FILE *fw=fopen(Combine(Combine("truncated_wave_function/", 2), sup_space.env_space->TotSiteNo), "rb");

        fread(&untruncated_block_number, sizeof(int), 1, fw);

        H_env_untrun=new int [untruncated_block_number];
        J_env_untrun=new int [untruncated_block_number];
        H_sys_untrun=new int [untruncated_block_number];
        J_sys_untrun=new int [untruncated_block_number];
        H_sysnew_untrun=new int [untruncated_block_number];
        J_sysnew_untrun=new int [untruncated_block_number];
        untruncated_block_dim=new int [untruncated_block_number];

        fread(H_env_untrun, sizeof(int), untruncated_block_number, fw);
        fread(J_env_untrun, sizeof(int), untruncated_block_number, fw);
        fread(H_sys_untrun, sizeof(int), untruncated_block_number, fw);
        fread(J_sys_untrun, sizeof(int), untruncated_block_number, fw);
        fread(H_sysnew_untrun, sizeof(int), untruncated_block_number, fw);
        fread(J_sysnew_untrun, sizeof(int), untruncated_block_number, fw);
        fread(untruncated_block_dim, sizeof(int), untruncated_block_number, fw);

        untruncated_wave_function=new double * [untruncated_block_number];
        for(int i=0; i<untruncated_block_number; i++) {
                untruncated_wave_function[i]=new double [untruncated_block_dim[i]];
                for(int j=0; j<untruncated_block_dim[i]; j++)
                        untruncated_wave_function[i][j]=0.0;
        }

        for(int i=0; i<untruncated_block_number; i++)
                fread(untruncated_wave_function[i], sizeof(double), untruncated_block_dim[i], fw);

        fclose(fw);

    //------Read in the block spin value of Block_sys[n]+S[n+1]
        FILE *fr=fopen(Combine(Combine("new_block/", 1), sup_space.sysnew_space->TotSiteNo+1), "rb");
        //---
        fread(&new_Block_Num_hole, sizeof(int), 1, fr);
        //---
        new_Num_hole_block=new int [new_Block_Num_hole];
        new_Num_J_hole_block=new int [new_Block_Num_hole];

        fread(new_Num_hole_block, sizeof(int), new_Block_Num_hole, fr);
        fread(new_Num_J_hole_block, sizeof(int), new_Block_Num_hole, fr);
        //---
        new_Value_J_block=new int * [new_Block_Num_hole];
        new_Dim_J_block=new int * [new_Block_Num_hole];
        for(int i=0; i<new_Block_Num_hole; i++) {
                new_Value_J_block[i]=new int [new_Num_J_hole_block[i]];
                new_Dim_J_block[i]=new int [new_Num_J_hole_block[i]];
                for(int j=0; j<new_Num_J_hole_block[i]; j++) {
                        new_Value_J_block[i][j]=0.0;  new_Dim_J_block[i][j]=0.0;
                }
        }

        for(int i=0; i<new_Block_Num_hole; i++) {
                fread(new_Value_J_block[i], sizeof(int), new_Num_J_hole_block[i], fr);
                fread(new_Dim_J_block[i], sizeof(int), new_Num_J_hole_block[i], fr);
        }
        //---
        new_Hole_blockOld=new int ** [new_Block_Num_hole];
        new_J_blockOld=new int ** [new_Block_Num_hole];
        new_Start=new int ** [new_Block_Num_hole];
        for(int i=0; i<new_Block_Num_hole; i++) {
                new_Hole_blockOld[i]=new int * [3];
                new_J_blockOld[i]=new int * [3];
                new_Start[i]=new int * [3];
                for(int j=0; j<3; j++) {
                        new_Hole_blockOld[i][j]=new int [new_Num_J_hole_block[i]];
                        new_J_blockOld[i][j]=new int [new_Num_J_hole_block[i]];
                        new_Start[i][j]=new int [new_Num_J_hole_block[i]];
                        for(int k=0; k<new_Num_J_hole_block[i]; k++) {
                                new_Hole_blockOld[i][j][k]=0.0;
                                new_J_blockOld[i][j][k]=0.0;
                                new_Start[i][j][k]=0.0;
                        }
                }
        }

        for(int i=0; i<new_Block_Num_hole; i++)
        for(int j=0; j<3; j++) {
                fread(new_Hole_blockOld[i][j], sizeof(int), new_Num_J_hole_block[i], fr);
                fread(new_J_blockOld[i][j], sizeof(int), new_Num_J_hole_block[i], fr);
                fread(new_Start[i][j], sizeof(int), new_Num_J_hole_block[i], fr);
        }

        fclose(fr);

    //------The 6-j coefficient for the basis transformation in this step
        six_j_basis_transformation=new double * [truncated_block_number];

        for(int i=0; i<truncated_block_number; i++) {

                index=0;

                for(int j=0; j<untruncated_block_number; j++)
                if(H_sys_untrun[j]==H_trun[i] && J_sys_untrun[j]==J_trun[i] && H_env_untrun[j]==H_old[i] && J_env_untrun[j]==J_old[i])

                        index++;

                six_j_basis_transformation[i]=new double [index];
                for(int n=0; n<index; n++)
                        six_j_basis_transformation[i][n]=(double) 0;

                index=0;

                for(int j=0; j<untruncated_block_number; j++)
                if(H_sys_untrun[j]==H_trun[i] && J_sys_untrun[j]==J_trun[i] && H_env_untrun[j]==H_old[i] && J_env_untrun[j]==J_old[i]) {//Block_sys[n+1] and Block_env[n+3]

			alpha = sqrt( ( sup_space.envnew_space->Value_J_block[H_new[i]][J_new[i]] + 1.0 ) * ( new_Value_J_block[H_sysnew_untrun[j]][J_sysnew_untrun[j]] + 1.0 ) );

                        if( (sup_space.envnew_space->Num_hole_block[H_new[i]] - sup_space.env_space->Num_hole_block[H_old[i]]) == 1 ) {// on n+1 site, we have a hole !!!

                                six_j_basis_transformation[i][index] += alpha * gsl_sf_coupling_6j( truncated_Value_J_block[H_trun[i]][J_trun[i]], 0, new_Value_J_block[H_sysnew_untrun[j]][J_sysnew_untrun[j]], sup_space.env_space->Value_J_block[H_old[i]][J_old[i]], para.Total_J, sup_space.envnew_space->Value_J_block[H_new[i]][J_new[i]] );

                                if( ( ( sup_space.env_space->Value_J_block[H_old[i]][J_old[i]] + para.Total_J + truncated_Value_J_block[H_trun[i]][J_trun[i]] ) / 2 ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

                                if( ( ( sup_space.env_space->Value_J_block[H_old[i]][J_old[i]] - sup_space.envnew_space->Value_J_block[H_new[i]][J_new[i]] ) / 2 ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

				index++;

			}

                        else if( sup_space.envnew_space->Num_hole_block[H_new[i]] == sup_space.env_space->Num_hole_block[H_old[i]] ) { // on n+1 site, we have a spin !!!

                                six_j_basis_transformation[i][index] += alpha * gsl_sf_coupling_6j( truncated_Value_J_block[H_trun[i]][J_trun[i]], para.Spin, new_Value_J_block[H_sysnew_untrun[j]][J_sysnew_untrun[j]], sup_space.env_space->Value_J_block[H_old[i]][J_old[i]], para.Total_J, sup_space.envnew_space->Value_J_block[H_new[i]][J_new[i]] );

                                if( ( ( para.Spin + sup_space.env_space->Value_J_block[H_old[i]][J_old[i]] + para.Total_J + truncated_Value_J_block[H_trun[i]][J_trun[i]] ) / 2 ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

                                if( ( ( para.Spin + sup_space.env_space->Value_J_block[H_old[i]][J_old[i]] - sup_space.envnew_space->Value_J_block[H_new[i]][J_new[i]] ) / 2 ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

                                if( ( sup_space.env_space->TotSiteNo - sup_space.env_space->Num_hole_block[H_old[i]] ) % 2 == 1 )

                                        six_j_basis_transformation[i][index] = - six_j_basis_transformation[i][index];

				index++;

			}

		}

	}

  //------Transform the wave function from the basis (Block_sys[n]+S[n+1], Block_env[n+2]) to (Block_sys[n], 
    //S[n+1]+Block_env[n+2]): multiply a factor of the 6-j coefficient obtained above------------------------------
        int inc=1;              //increament index in BLAS subroutine daxpy_()
        for(int i=0; i<truncated_block_number; i++) {

                index=0;
                for(int j=0; j<untruncated_block_number; j++)
                if(H_sys_untrun[j]==H_trun[i] && J_sys_untrun[j]==J_trun[i] && H_env_untrun[j]==H_old[i] && J_env_untrun[j]==J_old[i]) {

                        daxpy_(&truncated_block_dim[i], &six_j_basis_transformation[i][index], untruncated_wave_function[j], &inc, truncated_wave_function[i], &inc);

                        index++;

                }

        }

//------Delete the untruncated block quantities to reuse them for the next transformation step-------------------
        for(int i=0; i<untruncated_block_number; i++) {
                delete [] untruncated_wave_function[i];
        }
        delete [] untruncated_wave_function;

        delete [] H_sysnew_untrun;  delete [] J_sysnew_untrun;  delete [] H_sys_untrun;  delete [] J_sys_untrun;
        delete [] H_env_untrun;  delete [] J_env_untrun;  delete [] untruncated_block_dim;

//------Third wave function transformation: from (Block_sys[n], S[n+1]+Block_env[n+2]) to (Block_sys[n-1]+S[n],
//S[n+1]+Block_env[n+2]), which is the last step of transformation to provide the initial guess function for diag
//-onalization. The space for the transformed wave function has been created in Class Super!----------------------
  //------Create space for the transformed wave function (use untruncated block quantities)----------------------
        //------Find block number------------------------------------------------------------------------------------
        untruncated_block_number=0;

        for(int n_1=0; n_1<sup_space.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup_space.envnew_space->Block_Num_hole; n_2++)
        if(sup_space.sysnew_space->Num_hole_block[n_1]+sup_space.envnew_space->Num_hole_block[n_2]==para.Total_h) {

                for(int j_sysnew=0; j_sysnew<sup_space.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envnew=0; j_envnew<sup_space.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                        J_min=abs(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]-sup_space.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_max=(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]+sup_space.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if(para.Total_J==(J_min+2*n)) {

                                for(int j_e=0; j_e<3; j_e++)
                                if(sup_space.envnew_space->Hole_blockOld[n_2][j_e][j_envnew]!=-1 && sup_space.envnew_space->J_blockOld[n_2][j_e][j_envnew]!=-1)
                                        untruncated_block_number++;

                        }

                }

        }

    //------Create space for angular momentum numbers------------------------------------------------------------
        H_sysnew_untrun=new int [untruncated_block_number];
        J_sysnew_untrun=new int [untruncated_block_number];
        H_envnew_untrun=new int [untruncated_block_number];
        J_envnew_untrun=new int [untruncated_block_number];
        H_env_untrun=new int [untruncated_block_number];
        J_env_untrun=new int [untruncated_block_number];
        untruncated_block_dim=new int [untruncated_block_number];

        for(int i=0; i<untruncated_block_number; i++) {
                H_sysnew_untrun[i]=0;  J_sysnew_untrun[i]=0;   H_envnew_untrun[i]=0;  J_envnew_untrun[i]=0;
                H_env_untrun[i]=0;  J_env_untrun[i]=0;  untruncated_block_dim[i]=0;
        }

    //------Initialize the above angular momentum numbers--------------------------------------------------------
        index=0;

        for(int n_1=0; n_1<sup_space.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup_space.envnew_space->Block_Num_hole; n_2++)
        if(sup_space.sysnew_space->Num_hole_block[n_1]+sup_space.envnew_space->Num_hole_block[n_2]==para.Total_h) {

                for(int j_sysnew=0; j_sysnew<sup_space.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envnew=0; j_envnew<sup_space.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                        J_min=abs(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]-sup_space.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_max=(sup_space.sysnew_space->Value_J_block[n_1][j_sysnew]+sup_space.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if(para.Total_J==(J_min+2*n)) {

                                for(int j_e=0; j_e<3; j_e++)
                                if((oldH_env=sup_space.envnew_space->Hole_blockOld[n_2][j_e][j_envnew])!=-1 && (oldJ_env=sup_space.envnew_space->J_blockOld[n_2][j_e][j_envnew])!=-1) {
                                        H_sysnew_untrun[index]=n_1;  J_sysnew_untrun[index]=j_sysnew;
                                        H_envnew_untrun[index]=n_2;  J_envnew_untrun[index]=j_envnew;
                                        H_env_untrun[index]=oldH_env;  J_env_untrun[index]=oldJ_env;
                                        untruncated_block_dim[index++]=sup_space.env_space->Dim_J_block[oldH_env][oldJ_env]*sup_space.sysnew_space->Dim_J_block[n_1][j_sysnew];
                                }

                        }

                }
    
	}

    //------Create space for untruncated_wave_function-----------------------------------------------------------
        untruncated_wave_function=new double * [untruncated_block_number];

        for(int i=0; i<untruncated_block_number; i++) {

                untruncated_wave_function[i]=new double [untruncated_block_dim[i]];

                for(int j=0; j<untruncated_block_dim[i]; j++)
                        untruncated_wave_function[i][j]=(double) 0;

        }

    //------Wave function transformation from (Block_sys[n], S[n+1]+Block_env[n+2]) to (Block_sys[n-1]+S[n], 
        //S[n+1]+Block_env[n+2]) by matrix-matrix multiplication-----------------------------------------------------
        for(int i=0; i<untruncated_block_number; i++)
        for(int j=0; j<truncated_block_number; j++)
        if(H_new[j]==H_envnew_untrun[i] && J_new[j]==J_envnew_untrun[i] && truncated_Num_hole_block[H_trun[j]]==sup_space.sysnew_space->Num_hole_block[H_sysnew_untrun[i]] && truncated_Value_J_block[H_trun[j]][J_trun[j]]==sup_space.sysnew_space->Value_J_block[H_sysnew_untrun[i]][J_sysnew_untrun[i]] && H_old[j]==H_env_untrun[i] && J_old[j]==J_env_untrun[i]) {
                dgemm_(&trans_N, &trans_N, &sup_space.sysnew_space->Dim_J_block[H_sysnew_untrun[i]][J_sysnew_untrun[i]], &sup_space.env_space->Dim_J_block[H_env_untrun[i]][J_env_untrun[i]], &truncated_Dim_J_block[H_trun[j]][J_trun[j]], &alpha, truncated_density_eigenvector[H_trun[j]][J_trun[j]], &sup_space.sysnew_space->Dim_J_block[H_sysnew_untrun[i]][J_sysnew_untrun[i]], truncated_wave_function[j], &truncated_Dim_J_block[H_trun[j]][J_trun[j]], &beta, untruncated_wave_function[i], &sup_space.sysnew_space->Dim_J_block[H_sysnew_untrun[i]][J_sysnew_untrun[i]]);
        }
 
    //------Read from untruncated_wave_function to sup.WaveFunction_block----------------------------------------
        for(int i=0; i<sup_space.BlockNumber_for_TargetBlock; i++)
        for(int j=0; j<untruncated_block_number; j++)
        if(sup_space.H_sysnew[i]==H_sysnew_untrun[j] && sup_space.J_sysnew[i]==J_sysnew_untrun[j] && sup_space.H_envnew[i]==H_envnew_untrun[j] && sup_space.J_envnew[i]==J_envnew_untrun[j] && sup_space.H_env[i]==H_env_untrun[j] && sup_space.J_env[i]==J_env_untrun[j]) {

                for(int j_s=0; j_s<3; j_s++)
                if( (oldH_sys=sup_space.sysnew_space->Hole_blockOld[H_sysnew_untrun[j]][j_s][J_sysnew_untrun[j]])!=-1 && (oldJ_sys=sup_space.sysnew_space->J_blockOld[H_sysnew_untrun[j]][j_s][J_sysnew_untrun[j]])!=-1 && oldH_sys==sup_space.H_sys[i] && oldJ_sys==sup_space.J_sys[i] ) {

                        for(int a_env=0; a_env<sup_space.env_space->Dim_J_block[sup_space.H_env[i]][sup_space.J_env[i]]; a_env++)
                        for(int a_sys=0; a_sys<sup_space.sys_space->Dim_J_block[sup_space.H_sys[i]][sup_space.J_sys[i]]; a_sys++) {

                                position_old=a_env*sup_space.sysnew_space->Dim_J_block[sup_space.H_sysnew[i]][sup_space.J_sysnew[i]]+sup_space.sysnew_space->Start[H_sysnew_untrun[j]][j_s][J_sysnew_untrun[j]]+a_sys;
                                position_new=a_env*sup_space.sys_space->Dim_J_block[oldH_sys][oldJ_sys]+a_sys;
				sup_space.WaveFunction_block[i][position_new]=untruncated_wave_function[j][position_old];

			}

		}

	}

//------Delete the created space---------------------------------------------------------------------------------
	//------
        for(int i=0; i<untruncated_block_number; i++) {
                delete [] untruncated_wave_function[i];      //delete [] untruncated_wave_function_excited[i];
        }
        delete [] untruncated_wave_function;                 //delete [] untruncated_wave_function_excited;

        delete [] H_sysnew_untrun;  delete [] J_sysnew_untrun;
	delete [] H_envnew_untrun;  delete [] J_envnew_untrun;      
	delete [] H_env_untrun;     delete [] J_env_untrun;
        delete [] untruncated_block_dim;

	//------
        for(int i=0; i<truncated_block_number; i++)
                delete [] six_j_basis_transformation[i];
        delete [] six_j_basis_transformation;

	//------
        for(int i=0; i<new_Block_Num_hole; i++) {
		for(int j=0; j<3; j++) {
			delete [] new_Hole_blockOld[i][j];  delete [] new_J_blockOld[i][j];  delete [] new_Start[i][j];
		}
		delete [] new_Hole_blockOld[i];  delete [] new_J_blockOld[i];  delete [] new_Start[i];
	}
	delete [] new_Hole_blockOld;  delete [] new_J_blockOld;  delete [] new_Start;
	
        for(int i=0; i<new_Block_Num_hole; i++) {
		delete [] new_Value_J_block[i];  delete [] new_Dim_J_block[i];
	}
	delete [] new_Value_J_block;  delete [] new_Dim_J_block;

	delete [] new_Num_hole_block;  delete [] new_Num_J_hole_block;

	//------
        for(int i=0; i<truncated_block_number; i++) {
                delete [] truncated_wave_function[i];        //delete [] truncated_wave_function_excited[i];
        }
        delete [] truncated_wave_function;                   //delete [] truncated_wave_function_excited;

        delete [] H_trun;  delete [] J_trun;
	delete [] H_old;  delete [] J_old;
	delete [] H_new;  delete [] J_new;
        delete [] truncated_block_dim;

	//------
        for(int i=0; i<truncated_Block_Num_hole; i++)
		delete [] truncated_Old_J[i];
	delete [] truncated_Old_J;

	delete [] truncated_Old_hole;

        for(int i=0; i<truncated_Block_Num_hole; i++) {
		for(int j=0; j<truncated_Num_J_hole_block[i]; j++) {
			delete [] truncated_density_eigenvector[i][j];
		}
		delete [] truncated_density_eigenvector[i];
	}
	delete [] truncated_density_eigenvector;	

        for(int i=0; i<truncated_Block_Num_hole; i++) {
		delete [] truncated_Value_J_block[i];  delete [] truncated_Dim_J_block[i];  delete [] truncated_density_dim[i];
	}
	delete [] truncated_Value_J_block;  delete [] truncated_Dim_J_block;  delete [] truncated_density_dim;

	delete [] truncated_Num_hole_block;  delete [] truncated_Num_J_hole_block;

}

//================================================================================================================
//To store the density eigenvectors and to perform the first wave function transformation for the guess initial 
//state function in finite sweep. For sysnew:from (Block_sys[n]+S[n+1], S[n+2]+Block_env[n+3]) to (Block_sys[n+1],
//S[n+2]+Block_env[n+3]); for envnew:from (Block_sys[n]+S[n+1], S[n+2]+Block_env[n+3]) to (Block_sys[n]+S[n+1],
///Block_env[n+2]).
//================================================================================================================
//==============Truncate wave function and store the truncated density eigenvectors for finite sweep==============
SuperEnergy::SuperEnergy( Parameter &para, Super &sup, Sub &old, Sub &trun, char &sign ):Conjugate( sup.Dim ) {

        trans_N = 'N';    trans_T = 'T';    

	alpha = 1.0;      beta = 0.0;

        if( sign == 'r' )

                Truncate_sysnew_density_eigenvector( para, sup, old, trun );

        else if( sign == 'l' )

                Truncate_envnew_density_eigenvector( para, sup, old, trun );

}

//=======Truncate density eigenvector and wave function:sysnew======
inline void SuperEnergy::Truncate_sysnew_density_eigenvector(Parameter &para, Super &sup, Sub &old, Sub &trun) {

//------Create space for truncated_density_eigenvector
        truncated_density_dim = new int * [ trun.Block_Num_hole ];

        for( int i = 0; i < trun.Block_Num_hole; i++ ) {

                truncated_density_dim[i] = new int [ trun.Num_J_hole_block[i] ];

                for( int j = 0; j < trun.Num_J_hole_block[i]; j++ )

                        truncated_density_dim[i][j] = trun.Dim_J_block[i][j] * sup.sysnew_space->Dim_J_block[ trun.Old_hole[i] ][ trun.Old_J[i][j] ];   //total dimension of each truncated block of reduced density matrix

        }

        truncated_density_eigenvector = new double ** [ trun.Block_Num_hole ];

        for( int i = 0; i < trun.Block_Num_hole; i++ ) {

                truncated_density_eigenvector[i] = new double * [ trun.Num_J_hole_block[i] ];

                for( int j = 0; j < trun.Num_J_hole_block[i]; j++ ) {

                        truncated_density_eigenvector[i][j] = new double [ truncated_density_dim[i][j] ];

                        for( int k = 0; k < truncated_density_dim[i][j]; k++ )

                                truncated_density_eigenvector[i][j][k] = old.dm_wave[ trun.Old_hole[i] ][ trun.Old_J[i][j] ][ k ];

                }

        }

//------Truncate the obtained wave function from (sys+ns,env+ne) to (systrun,env+ne)
    //--Allocate space for the truncated wave function
        truncated_block_number = 0;

        for( int n_1 = 0; n_1 < trun.Block_Num_hole; n_1++ )
        for( int n_2 = 0; n_2 < sup.envnew_space->Block_Num_hole; n_2++ )
        if( trun.Num_hole_block[ n_1 ] + sup.envnew_space->Num_hole_block[ n_2 ] == sup.QuantumNumber_hole ) {

                for( int j_systrun = 0; j_systrun < trun.Num_J_hole_block[ n_1 ]; j_systrun++ )
                for( int j_envnew = 0; j_envnew < sup.envnew_space->Num_J_hole_block[ n_2 ]; j_envnew++ ) {

                        J_min = abs( trun.Value_J_block[ n_1 ][ j_systrun ] - sup.envnew_space->Value_J_block[ n_2 ][ j_envnew ] );
                        J_max = ( trun.Value_J_block[ n_1 ][ j_systrun ] + sup.envnew_space->Value_J_block[ n_2 ][ j_envnew ] );
                        J_num = ( J_max - J_min )/2 + 1;

                        for( int n = 0; n < J_num; n++ )
                        if( sup.QuantumNumber_J == ( J_min + 2 * n ) ) {

                                for( int j_e = 0; j_e < 3; j_e++ )
                                if( sup.envnew_space->Hole_blockOld[ n_2 ][ j_e ][ j_envnew ] != -1 && sup.envnew_space->J_blockOld[ n_2 ][ j_e ][ j_envnew ] != -1 )

                                        truncated_block_number++;

                       }

                }

        }

     //------Create space for the angular indices of the blocks
        H_trun = new int [ truncated_block_number ];                //H_trun=H_systrun
        J_trun = new int [ truncated_block_number ];                //J_trun=J_systrun
        H_old = new int [ truncated_block_number ];                 //H_old=H_env
        J_old = new int [ truncated_block_number ];                 //J_old=J_env
        H_new = new int [ truncated_block_number ];                 //H_new=H_envnew
        J_new = new int [ truncated_block_number ];                 //J_new=J_envnew
        truncated_block_dim = new int [ truncated_block_number ];   //truncated_block_dim=Dim_systrun*Dim_env

        for( int i = 0; i < truncated_block_number; i++ ) {
                H_trun[i] = 0;  J_trun[i] = 0;  H_old[i] = 0;  J_old[i] = 0;  H_new[i] = 0;  J_new[i] = 0;  
		truncated_block_dim[i] = 0;
        }

     //------Initialize the values of the above angular indices
        index = 0;

        for( int n_1 = 0; n_1 < trun.Block_Num_hole; n_1++ )
        for( int n_2 = 0; n_2 < sup.envnew_space->Block_Num_hole; n_2++ )
        if( trun.Num_hole_block[n_1] + sup.envnew_space->Num_hole_block[n_2] == sup.QuantumNumber_hole ) {

                for( int j_systrun = 0; j_systrun < trun.Num_J_hole_block[n_1]; j_systrun++ )
                for( int j_envnew = 0; j_envnew < sup.envnew_space->Num_J_hole_block[n_2]; j_envnew++ ) {

                        J_min = abs( trun.Value_J_block[n_1][j_systrun] - sup.envnew_space->Value_J_block[n_2][j_envnew] );
                        J_max = ( trun.Value_J_block[n_1][j_systrun] + sup.envnew_space->Value_J_block[n_2][j_envnew] );
                        J_num = (J_max-J_min)/2 + 1;

                        for( int n = 0; n < J_num; n++ )
                        if( sup.QuantumNumber_J == ( J_min + 2*n ) ) {

                                for( int j_e = 0; j_e < 3; j_e++ )
                                if( ( oldH_env = sup.envnew_space->Hole_blockOld[n_2][j_e][j_envnew] ) != -1 && ( oldJ_env = sup.envnew_space->J_blockOld[n_2][j_e][j_envnew] ) != -1 ) {

                                        H_trun[ index ] = n_1;  J_trun[ index ] = j_systrun;
					H_old[ index ] = oldH_env;  J_old[ index ] = oldJ_env;  
					H_new[ index ] = n_2;  J_new[ index ] = j_envnew;
                                        truncated_block_dim[ index ] = trun.Dim_J_block[ n_1 ][ j_systrun ] * sup.env_space->Dim_J_block[ oldH_env ][ oldJ_env ];

					index++;

                                }

                       }
     
		}

	}

     //------Create space for the truncated wave function
        truncated_wave_function = new double * [ truncated_block_number ];

        for(int i=0; i<truncated_block_number; i++) {

                truncated_wave_function[i] = new double [ truncated_block_dim[i] ];

                for(int j=0; j < truncated_block_dim[i]; j++)

                        truncated_wave_function[i][j] = 0.0;

        }

    //------Allocate space for the untruncated wave function: angular numbers J_sys for a given J_sysnew are combined: in order to perform the matrix-matrix multiplication of calculating the truncated wave function
        untruncated_block_number=0;

        for( int n_1 = 0; n_1 < sup.sysnew_space->Block_Num_hole; n_1++ )
        for( int n_2 = 0; n_2 < sup.envnew_space->Block_Num_hole; n_2++ )
        if( sup.sysnew_space->Num_hole_block[n_1] + sup.envnew_space->Num_hole_block[n_2] == sup.QuantumNumber_hole ) {

                for( int j_sysnew = 0; j_sysnew < sup.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++ )
                for( int j_envnew = 0; j_envnew < sup.envnew_space->Num_J_hole_block[n_2]; j_envnew++ ) {

                        J_min = abs( sup.sysnew_space->Value_J_block[n_1][j_sysnew] - sup.envnew_space->Value_J_block[n_2][j_envnew] );
                        J_max = ( sup.sysnew_space->Value_J_block[n_1][j_sysnew] + sup.envnew_space->Value_J_block[n_2][j_envnew] );
                        J_num = (J_max-J_min)/2 + 1;

                        for( int n = 0; n < J_num; n++ )
                        if( sup.QuantumNumber_J == ( J_min + 2*n) ) {

                                for( int j_e = 0; j_e < 3; j_e++ )
                                if( sup.envnew_space->Hole_blockOld[n_2][j_e][j_envnew] != -1 && sup.envnew_space->J_blockOld[n_2][j_e][j_envnew] != -1 )

                                        untruncated_block_number++;

                        }

                }

        }

   //------Create space for the angular numbers of the untruncated blocks
        H_sysnew_untrun = new int [untruncated_block_number];
        J_sysnew_untrun = new int [untruncated_block_number];
        H_envnew_untrun = new int [untruncated_block_number];
        J_envnew_untrun = new int [untruncated_block_number];
        H_env_untrun = new int [untruncated_block_number];
        J_env_untrun = new int [untruncated_block_number];
        untruncated_block_dim = new int [untruncated_block_number];

        for( int i=0; i<untruncated_block_number; i++ ) {

                H_sysnew_untrun[i]=0;  J_sysnew_untrun[i]=0;  
		
		H_envnew_untrun[i]=0;  J_envnew_untrun[i]=0;  

		H_env_untrun[i]=0;  J_env_untrun[i]=0;  

		untruncated_block_dim[i]=0;

        }

   //------Initialize the above angular numbers
        index = 0;

        for(int n_1=0; n_1<sup.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup.envnew_space->Block_Num_hole; n_2++)
        if(sup.sysnew_space->Num_hole_block[n_1]+sup.envnew_space->Num_hole_block[n_2]== sup.QuantumNumber_hole ) {

                for(int j_sysnew=0; j_sysnew<sup.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envnew=0; j_envnew<sup.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                        J_min=abs(sup.sysnew_space->Value_J_block[n_1][j_sysnew]-sup.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_max=(sup.sysnew_space->Value_J_block[n_1][j_sysnew]+sup.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if( sup.QuantumNumber_J == (J_min+2*n)) {

                                for(int j_e=0; j_e<3; j_e++)
                                if((oldH_env=sup.envnew_space->Hole_blockOld[n_2][j_e][j_envnew])!=-1 && (oldJ_env=sup.envnew_space->J_blockOld[n_2][j_e][j_envnew])!=-1) {

                                        H_sysnew_untrun[index]=n_1;  
					J_sysnew_untrun[index]=j_sysnew;

                                        H_envnew_untrun[index]=n_2;  
					J_envnew_untrun[index]=j_envnew;

                                        H_env_untrun[index]=oldH_env;
                                        J_env_untrun[index]=oldJ_env;

                                        untruncated_block_dim[index]=sup.sysnew_space->Dim_J_block[n_1][j_sysnew]*sup.env_space->Dim_J_block[oldH_env][oldJ_env];

					index++;
	
				}

			}

		}

	}

   //---Create space for the untruncated function and initialize the untruncated wave function
        untruncated_wave_function=new double * [untruncated_block_number];

        for(int i=0; i<untruncated_block_number; i++) {

                untruncated_wave_function[i]=new double [untruncated_block_dim[i]];

                for(int j=0; j<untruncated_block_dim[i]; j++)

                        untruncated_wave_function[i][j] = (double) 0;

        }

        for(int i=0; i<untruncated_block_number; i++) {

                for(int ns=0; ns<3; ns++)
                if((oldH_sys=sup.sysnew_space->Hole_blockOld[H_sysnew_untrun[i]][ns][J_sysnew_untrun[i]])!=-1 && (oldJ_sys=sup.sysnew_space->J_blockOld[H_sysnew_untrun[i]][ns][J_sysnew_untrun[i]])!=-1) {

                        for(int j=0; j<sup.BlockNumber_for_TargetBlock; j++)
                        if(sup.H_sysnew[j]==H_sysnew_untrun[i] && sup.J_sysnew[j]==J_sysnew_untrun[i] && sup.H_envnew[j]==H_envnew_untrun[i] && sup.J_envnew[j]==J_envnew_untrun[i] && sup.H_sys[j]==oldH_sys && sup.J_sys[j]==oldJ_sys && sup.H_env[j]==H_env_untrun[i] && sup.J_env[j]==J_env_untrun[i]) {

                                for(int a_env=0; a_env<sup.env_space->Dim_J_block[H_env_untrun[i]][J_env_untrun[i]]; a_env++)
                                for(int a_sys=0; a_sys<sup.sys_space->Dim_J_block[oldH_sys][oldJ_sys]; a_sys++) {

                                        position_old=a_env*sup.sys_space->Dim_J_block[oldH_sys][oldJ_sys]+a_sys;

                                        position_new=a_env*sup.sysnew_space->Dim_J_block[H_sysnew_untrun[i]][J_sysnew_untrun[i]]+sup.sysnew_space->Start[H_sysnew_untrun[i]][ns][J_sysnew_untrun[i]]+a_sys;

                                        untruncated_wave_function[i][position_new]=sup.WaveFunction_block[j][position_old];

                                }

                        }

                }

        }

    //------Truncate the wave function by matrix-matrix multiplication
       for( int j_trun = 0; j_trun < truncated_block_number; j_trun++ )       
       for( int j_untrun = 0; j_untrun < untruncated_block_number; j_untrun++ )  
       if( H_old[ j_trun ] == H_env_untrun[ j_untrun ] && J_old[ j_trun ] == J_env_untrun[ j_untrun ] && H_new[ j_trun ] == H_envnew_untrun[ j_untrun ] && J_new[ j_trun ] == J_envnew_untrun[ j_untrun ] && trun.Old_hole[ H_trun[j_trun] ] == H_sysnew_untrun[ j_untrun ] && trun.Old_J[ H_trun[j_trun] ][ J_trun[j_trun] ] == J_sysnew_untrun[ j_untrun ] ) {

//		cout<<"\n M="<<trun.Dim_J_block[ H_trun[j_trun] ][ J_trun[j_trun] ]<<"\t N="<<sup.env_space->Dim_J_block[ H_old[j_trun] ][ J_old[j_trun] ]<<"\t K="<<sup.sysnew_space->Dim_J_block[ H_sysnew_untrun[j_untrun] ][ J_sysnew_untrun[j_untrun] ];

                dgemm_( &trans_T, &trans_N, &trun.Dim_J_block[ H_trun[j_trun] ][ J_trun[j_trun] ], &sup.env_space->Dim_J_block[ H_old[j_trun] ][ J_old[j_trun] ], &sup.sysnew_space->Dim_J_block[ H_sysnew_untrun[j_untrun] ][ J_sysnew_untrun[j_untrun] ], &alpha, truncated_density_eigenvector[ H_trun[j_trun] ][ J_trun[j_trun] ], &sup.sysnew_space->Dim_J_block[ H_sysnew_untrun[j_untrun] ][ J_sysnew_untrun[j_untrun] ], untruncated_wave_function[ j_untrun ], &sup.sysnew_space->Dim_J_block[ H_sysnew_untrun[j_untrun] ][ J_sysnew_untrun[j_untrun] ], &beta, truncated_wave_function[j_trun], &trun.Dim_J_block[ H_trun[j_trun] ][ J_trun[j_trun] ] );

	}

    //------Print the truncated wave function
        FILE *fw=fopen(Combine(Combine("truncated_wave_function/", 1), trun.TotSiteNo), "wb");

        fwrite(&truncated_block_number, sizeof(int), 1, fw);
        fwrite(H_trun, sizeof(int), truncated_block_number, fw);
        fwrite(J_trun, sizeof(int), truncated_block_number, fw);
        fwrite(H_old, sizeof(int), truncated_block_number, fw);
        fwrite(J_old, sizeof(int), truncated_block_number, fw);
        fwrite(H_new, sizeof(int), truncated_block_number, fw);
        fwrite(J_new, sizeof(int), truncated_block_number, fw);
        fwrite(truncated_block_dim, sizeof(int), truncated_block_number, fw);

        for(int i=0; i<truncated_block_number; i++)

                fwrite(truncated_wave_function[i], sizeof(double), truncated_block_dim[i], fw);

        fclose(fw);

  //------Print truncated_density_eigenvector: valid for both systrun and envtrun
        FILE *fp=fopen(Combine(Combine("truncated_density_eigenvector/", 1), trun.TotSiteNo), "wb");

        fwrite(&trun.Block_Num_hole, sizeof(int), 1, fp);
        fwrite(trun.Num_hole_block, sizeof(int), trun.Block_Num_hole, fp);
        fwrite(trun.Num_J_hole_block, sizeof(int), trun.Block_Num_hole, fp);

        for(int i=0; i<trun.Block_Num_hole; i++) {

                fwrite(trun.Value_J_block[i], sizeof(int), trun.Num_J_hole_block[i], fp);
                fwrite(trun.Dim_J_block[i], sizeof(int), trun.Num_J_hole_block[i], fp);
                fwrite(truncated_density_dim[i], sizeof(int), trun.Num_J_hole_block[i], fp);

        }

        for(int i=0; i<trun.Block_Num_hole; i++)
        for(int j=0; j<trun.Num_J_hole_block[i]; j++)

                fwrite(truncated_density_eigenvector[i][j], sizeof(double), truncated_density_dim[i][j], fp);

        fwrite(trun.Old_hole, sizeof(int), trun.Block_Num_hole, fp);

        for(int i=0; i<trun.Block_Num_hole; i++)

                fwrite(trun.Old_J[i], sizeof(int), trun.Num_J_hole_block[i], fp);

        fclose(fp);

  //------Delete untruncated and truncated wave function
        for(int i=0; i<untruncated_block_number; i++)
                delete [] untruncated_wave_function[i];         //delete [] untruncated_wave_function_excited[i];
        delete [] untruncated_wave_function;                    //delete [] untruncated_wave_function_excited;

        delete [] H_sysnew_untrun;  delete [] J_sysnew_untrun;  delete [] H_envnew_untrun;  delete [] J_envnew_untrun;
        delete [] H_env_untrun;     delete [] J_env_untrun;     delete [] untruncated_block_dim;

        for(int i=0; i<truncated_block_number; i++)
                delete [] truncated_wave_function[i];           //delete [] truncated_wave_function_excited[i];
        delete [] truncated_wave_function;                      //delete [] truncated_wave_function_excited;

        delete [] H_trun;  delete [] J_trun;  delete [] H_old;  delete [] J_old; delete [] H_new;  delete [] J_new;
        delete [] truncated_block_dim;

  //------Delete truncated_density_eigenvector
        for(int i=0; i<trun.Block_Num_hole; i++) {
                for(int j=0; j<trun.Num_J_hole_block[i]; j++)
                        delete [] truncated_density_eigenvector[i][j];
                delete [] truncated_density_eigenvector[i];
        }
        delete [] truncated_density_eigenvector;

        for(int i=0; i<trun.Block_Num_hole; i++)
                delete [] truncated_density_dim[i];
        delete [] truncated_density_dim;

        delete [] trun.Old_hole;

        for(int i=0; i<trun.Block_Num_hole; i++)
                delete [] trun.Old_J[i];
        delete [] trun.Old_J;

}

//============================Truncated density eigenvector and wave function:envnew==============================
inline void SuperEnergy::Truncate_envnew_density_eigenvector(Parameter &para, Super &sup, Sub &old, Sub &trun) {
//------Create space for truncated_density_eigenvector: valid for both systrun and envtrun------------------------
        truncated_density_dim=new int * [trun.Block_Num_hole];//Dimension of truncated density eigenvector matrix

        for(int i=0; i<trun.Block_Num_hole; i++) {

                truncated_density_dim[i]=new int [trun.Num_J_hole_block[i]];

                for(int j=0; j<trun.Num_J_hole_block[i]; j++)

                        truncated_density_dim[i][j]=trun.Dim_J_block[i][j]*sup.envnew_space->Dim_J_block[trun.Old_hole[i]][trun.Old_J[i][j]];   //total dimension of each truncated block of reduced density matrix

        }

        truncated_density_eigenvector=new double ** [trun.Block_Num_hole];

        for(int i=0; i<trun.Block_Num_hole; i++) {

                truncated_density_eigenvector[i]=new double * [trun.Num_J_hole_block[i]];

                for(int j=0; j<trun.Num_J_hole_block[i]; j++) {

                        truncated_density_eigenvector[i][j]=new double [truncated_density_dim[i][j]];

                        for(int k=0; k<truncated_density_dim[i][j]; k++)

                                truncated_density_eigenvector[i][j][k]=old.dm_wave[trun.Old_hole[i]][trun.Old_J[i][j]][k];

                }

	}

//------Truncate the obtained wave function from (sys+ns,env+ne) to (sys+ns,envtrun)------------------------------
    //--Allocate space for the truncated wave function
        truncated_block_number=0;

        for(int n_1=0; n_1<sup.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<trun.Block_Num_hole; n_2++)
        if(sup.sysnew_space->Num_hole_block[n_1]+trun.Num_hole_block[n_2]== sup.QuantumNumber_hole ) {

                for(int j_sysnew=0; j_sysnew<sup.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envtrun=0; j_envtrun<trun.Num_J_hole_block[n_2]; j_envtrun++) {

                        J_min=abs(trun.Value_J_block[n_2][j_envtrun]-sup.sysnew_space->Value_J_block[n_1][j_sysnew]);
                        J_max=(trun.Value_J_block[n_2][j_envtrun]+sup.sysnew_space->Value_J_block[n_1][j_sysnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if( sup.QuantumNumber_J == (J_min+2*n)) {

                                for(int j_s=0; j_s<3; j_s++)
                                if(sup.sysnew_space->Hole_blockOld[n_1][j_s][j_sysnew]!=-1 && sup.sysnew_space->J_blockOld[n_1][j_s][j_sysnew]!=-1) {
                                        truncated_block_number++;
                                }

                       }

                }

        }

     //------Create space for the angular indices of the blocks--------------------------------------------------
        H_trun=new int [truncated_block_number];                //H_trun=H_envtrun
        J_trun=new int [truncated_block_number];                //J_trun=J_envtrun
        H_old=new int [truncated_block_number];                 //H_old=H_sys
        J_old=new int [truncated_block_number];                 //J_old=J_env
        H_new=new int [truncated_block_number];                 //H_new=H_sysnew
        J_new=new int [truncated_block_number];                 //J_new=J_sysnew
        truncated_block_dim=new int [truncated_block_number];   //truncated_block_dim=Dim_envtrun*Dim_sys

        for(int i=0; i<truncated_block_number; i++) {
                H_trun[i]=0;  J_trun[i]=0;  H_old[i]=0;  J_old[i]=0;  H_new[i]=0;  J_new[i]=0;  truncated_block_dim[i]=0;
        }

     //------Initialize the values of the above angular indices--------------------------------------------------
        index=0;
        for(int n_1=0; n_1<sup.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<trun.Block_Num_hole; n_2++)
        if(sup.sysnew_space->Num_hole_block[n_1]+trun.Num_hole_block[n_2]== sup.QuantumNumber_hole ) {

                for(int j_sysnew=0; j_sysnew<sup.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envtrun=0; j_envtrun<trun.Num_J_hole_block[n_2]; j_envtrun++) {

                        J_min=abs(trun.Value_J_block[n_2][j_envtrun]-sup.sysnew_space->Value_J_block[n_1][j_sysnew]);
                        J_max=(trun.Value_J_block[n_2][j_envtrun]+sup.sysnew_space->Value_J_block[n_1][j_sysnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if( sup.QuantumNumber_J == (J_min+2*n)) {

                                for(int j_s=0; j_s<3; j_s++)
                                if((oldH_sys=sup.sysnew_space->Hole_blockOld[n_1][j_s][j_sysnew])!=-1 && (oldJ_sys=sup.sysnew_space->J_blockOld[n_1][j_s][j_sysnew])!=-1) {

                                        H_trun[index]=n_2;  J_trun[index]=j_envtrun;  H_old[index]=oldH_sys;
                                        J_old[index]=oldJ_sys;  H_new[index]=n_1;  J_new[index]=j_sysnew;
                                        truncated_block_dim[index++]=trun.Dim_J_block[n_2][j_envtrun]*sup.sys_space->Dim_J_block[oldH_sys][oldJ_sys];

                                }

                       }

		}

	}

     //------Create space for the truncated wave function--------------------------------------------------------
        truncated_wave_function=new double * [truncated_block_number];

        for(int i=0; i<truncated_block_number; i++) {

                truncated_wave_function[i]=new double [truncated_block_dim[i]];

                for(int j=0; j<truncated_block_dim[i]; j++)

                        truncated_wave_function[i][j]=(double) 0;
        }

     //------Allocate space for the untruncated wave function: angular numbers J_env for a given J_envnew are combined: in order to perform the matrix-matrix multiplication of calculating the truncated wave function
        untruncated_block_number=0;

        for(int n_1=0; n_1<sup.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup.envnew_space->Block_Num_hole; n_2++)
        if(sup.sysnew_space->Num_hole_block[n_1]+sup.envnew_space->Num_hole_block[n_2]== sup.QuantumNumber_hole ) {

                for(int j_sysnew=0; j_sysnew<sup.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envnew=0; j_envnew<sup.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                        J_min=abs(sup.sysnew_space->Value_J_block[n_1][j_sysnew]-sup.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_max=(sup.sysnew_space->Value_J_block[n_1][j_sysnew]+sup.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if( sup.QuantumNumber_J == (J_min+2*n)) {

                                for(int j_s=0; j_s<3; j_s++)
                                if(sup.sysnew_space->Hole_blockOld[n_1][j_s][j_sysnew]!=-1 && sup.sysnew_space->J_blockOld[n_1][j_s][j_sysnew]!=-1)
                                        untruncated_block_number++;
                        }
                }

        }

   //------Create space for the angular numbers of the untruncated blocks----------------------------------------
        H_sysnew_untrun=new int [untruncated_block_number];
        J_sysnew_untrun=new int [untruncated_block_number];
        H_envnew_untrun=new int [untruncated_block_number];
        J_envnew_untrun=new int [untruncated_block_number];
        H_sys_untrun=new int [untruncated_block_number];
        J_sys_untrun=new int [untruncated_block_number];
        untruncated_block_dim=new int [untruncated_block_number];

        for(int i=0; i<untruncated_block_number; i++) {
                H_sysnew_untrun[i]=0;  J_sysnew_untrun[i]=0;  H_envnew_untrun[i]=0;  J_envnew_untrun[i]=0;  H_sys_untrun[i]=0;  J_sys_untrun[i]=0;  untruncated_block_dim[i]=0;
        }

   //------Initialize the above angular numbers------------------------------------------------------------------
        index=0;
        for(int n_1=0; n_1<sup.sysnew_space->Block_Num_hole; n_1++)
        for(int n_2=0; n_2<sup.envnew_space->Block_Num_hole; n_2++)
        if(sup.sysnew_space->Num_hole_block[n_1]+sup.envnew_space->Num_hole_block[n_2]== sup.QuantumNumber_hole ) {

                for(int j_sysnew=0; j_sysnew<sup.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
                for(int j_envnew=0; j_envnew<sup.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                        J_min=abs(sup.sysnew_space->Value_J_block[n_1][j_sysnew]-sup.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_max=(sup.sysnew_space->Value_J_block[n_1][j_sysnew]+sup.envnew_space->Value_J_block[n_2][j_envnew]);
                        J_num=(J_max-J_min)/2+1;

                        for(int n=0; n<J_num; n++)
                        if( sup.QuantumNumber_J == (J_min+2*n)) {

                                for(int j_s=0; j_s<3; j_s++)
                                if((oldH_sys=sup.sysnew_space->Hole_blockOld[n_1][j_s][j_sysnew])!=-1 && (oldJ_sys=sup.sysnew_space->J_blockOld[n_1][j_s][j_sysnew])!=-1) {

                                        H_sysnew_untrun[index]=n_1;  J_sysnew_untrun[index]=j_sysnew;
                                        H_envnew_untrun[index]=n_2;  J_envnew_untrun[index]=j_envnew;
                                        H_sys_untrun[index]=oldH_sys;
                                        J_sys_untrun[index]=oldJ_sys;
                                        untruncated_block_dim[index++]=sup.envnew_space->Dim_J_block[n_2][j_envnew]*sup.sys_space->Dim_J_block[oldH_sys][oldJ_sys];

				}

			}

		}

	}

   //------Create space for the untruncated function and initialize the untruncated wave function----------------
        untruncated_wave_function=new double * [untruncated_block_number];

        for(int i=0; i<untruncated_block_number; i++) {

                untruncated_wave_function[i]=new double [untruncated_block_dim[i]];

                for(int j=0; j<untruncated_block_dim[i]; j++)

                        untruncated_wave_function[i][j]=(double) 0;

        }

        for(int i=0; i<untruncated_block_number; i++) {

                for(int ne=0; ne<3; ne++)
                if((oldH_env=sup.envnew_space->Hole_blockOld[H_envnew_untrun[i]][ne][J_envnew_untrun[i]])!=-1 && (oldJ_env=sup.envnew_space->J_blockOld[H_envnew_untrun[i]][ne][J_envnew_untrun[i]])!=-1) {

                        for(int j=0; j<sup.BlockNumber_for_TargetBlock; j++)//find corresponding block
                        if(sup.H_sysnew[j]==H_sysnew_untrun[i] && sup.J_sysnew[j]==J_sysnew_untrun[i] && sup.H_envnew[j]==H_envnew_untrun[i] && sup.J_envnew[j]==J_envnew_untrun[i] && sup.H_sys[j]==H_sys_untrun[i] && sup.J_sys[j]==J_sys_untrun[i] && sup.H_env[j]==oldH_env && sup.J_env[j]==oldJ_env) {

                                for(int a_env=0; a_env<sup.env_space->Dim_J_block[oldH_env][oldJ_env]; a_env++)
                                for(int a_sys=0; a_sys<sup.sys_space->Dim_J_block[H_sys_untrun[i]][J_sys_untrun[i]]; a_sys++) {

                                       position_old=a_env*sup.sys_space->Dim_J_block[H_sys_untrun[i]][J_sys_untrun[i]]+a_sys;
                                        position_new=(a_env+sup.envnew_space->Start[H_envnew_untrun[i]][ne][J_envnew_untrun[i]])*sup.sys_space->Dim_J_block[H_sys_untrun[i]][J_sys_untrun[i]]+a_sys;
                                        untruncated_wave_function[i][position_new]=sup.WaveFunction_block[j][position_old];
                                }
                        }

                }

        }

     //------Truncate the wave function by matrix-matrix multiplication------------------------------------------
        for(int j_trun=0; j_trun<truncated_block_number; j_trun++)          //truncated wave function block
        for(int j_untrun=0; j_untrun<untruncated_block_number; j_untrun++)  //untruncated wave function block
        if(H_old[j_trun]==H_sys_untrun[j_untrun] && J_old[j_trun]==J_sys_untrun[j_untrun] && H_new[j_trun]==H_sysnew_untrun[j_untrun] && J_new[j_trun]==J_sysnew_untrun[j_untrun] && trun.Old_hole[H_trun[j_trun]]==H_envnew_untrun[j_untrun] && trun.Old_J[H_trun[j_trun]][J_trun[j_trun]]==J_envnew_untrun[j_untrun])

                dgemm_(&trans_N, &trans_N, &sup.sys_space->Dim_J_block[H_old[j_trun]][J_old[j_trun]], &trun.Dim_J_block[H_trun[j_trun]][J_trun[j_trun]], &sup.envnew_space->Dim_J_block[H_envnew_untrun[j_untrun]][J_envnew_untrun[j_untrun]], &alpha, untruncated_wave_function[j_untrun], &sup.sys_space->Dim_J_block[H_old[j_trun]][J_old[j_trun]], truncated_density_eigenvector[H_trun[j_trun]][J_trun[j_trun]], &sup.envnew_space->Dim_J_block[H_envnew_untrun[j_untrun]][J_envnew_untrun[j_untrun]], &beta, truncated_wave_function[j_trun], &sup.sys_space->Dim_J_block[H_old[j_trun]][J_old[j_trun]]);

    //------Print the truncated wave function--------------------------------------------------------------------
        FILE *fw=fopen(Combine(Combine("truncated_wave_function/", 2), trun.TotSiteNo), "wb");

        fwrite(&truncated_block_number, sizeof(int), 1, fw);
        fwrite(H_trun, sizeof(int), truncated_block_number, fw);
        fwrite(J_trun, sizeof(int), truncated_block_number, fw);
        fwrite(H_old, sizeof(int), truncated_block_number, fw);
        fwrite(J_old, sizeof(int), truncated_block_number, fw);
        fwrite(H_new, sizeof(int), truncated_block_number, fw);
        fwrite(J_new, sizeof(int), truncated_block_number, fw);
        fwrite(truncated_block_dim, sizeof(int), truncated_block_number, fw);

        for(int i=0; i<truncated_block_number; i++)
                fwrite(truncated_wave_function[i], sizeof(double), truncated_block_dim[i], fw);

        fclose(fw);

    //------Print truncated_density_eigenvector: valid for both systrun and envtrun------------------------------
        FILE *fp=fopen(Combine(Combine("truncated_density_eigenvector/", 2), trun.TotSiteNo), "wb");

        fwrite(&trun.Block_Num_hole, sizeof(int), 1, fp);
        fwrite(trun.Num_hole_block, sizeof(int), trun.Block_Num_hole, fp);
        fwrite(trun.Num_J_hole_block, sizeof(int), trun.Block_Num_hole, fp);

        for(int i=0; i<trun.Block_Num_hole; i++) {
                fwrite(trun.Value_J_block[i], sizeof(int), trun.Num_J_hole_block[i], fp);
                fwrite(trun.Dim_J_block[i], sizeof(int), trun.Num_J_hole_block[i], fp);
                fwrite(truncated_density_dim[i], sizeof(int), trun.Num_J_hole_block[i], fp);
        }

        for(int i=0; i<trun.Block_Num_hole; i++)
        for(int j=0; j<trun.Num_J_hole_block[i]; j++)
                fwrite(truncated_density_eigenvector[i][j], sizeof(double), truncated_density_dim[i][j], fp);

        fwrite(trun.Old_hole, sizeof(int), trun.Block_Num_hole, fp);

        for(int i=0; i<trun.Block_Num_hole; i++)
                fwrite(trun.Old_J[i], sizeof(int), trun.Num_J_hole_block[i], fp);

        fclose(fp);

//------Delete untruncated and truncated wave function-----------------------------------------------------------
        for(int i=0; i<untruncated_block_number; i++)
                delete [] untruncated_wave_function[i];         //delete [] untruncated_wave_function_excited[i];
        delete [] untruncated_wave_function;                    //delete [] untruncated_wave_function_excited;

        delete [] H_sysnew_untrun;  delete [] J_sysnew_untrun;  delete [] H_envnew_untrun;  delete [] J_envnew_untrun;
        delete [] H_sys_untrun;     delete [] J_sys_untrun;     delete [] untruncated_block_dim;

        for(int i=0; i<truncated_block_number; i++)
                delete [] truncated_wave_function[i];           //delete [] truncated_wave_function_excited[i];
        delete [] truncated_wave_function;                      //delete [] truncated_wave_function_excited;

        delete [] H_trun;  delete [] J_trun;  delete [] H_old;  delete [] J_old; delete [] H_new;  delete [] J_new;
        delete [] truncated_block_dim;

  //------Delete truncated_density_eigenvector-------------------------------------------------------------------
        for(int i=0; i<trun.Block_Num_hole; i++) {
                for(int j=0; j<trun.Num_J_hole_block[i]; j++)
                        delete [] truncated_density_eigenvector[i][j];
                delete [] truncated_density_eigenvector[i];
        }
        delete [] truncated_density_eigenvector;

        for(int i=0; i<trun.Block_Num_hole; i++)
                delete [] truncated_density_dim[i];
        delete [] truncated_density_dim;

        delete [] trun.Old_hole;

        for(int i=0; i<trun.Block_Num_hole; i++)
                delete [] trun.Old_J[i];
        delete [] trun.Old_J;

}

//===============================================================================================================
//In the end of each finite sweep, the wave function does not need to be truncated, and the truncated density 
//eigenvectors may not be stored, which would not be used
//===============================================================================================================
SuperEnergy::SuperEnergy(Super &sup, Sub &trun, const int &n):Conjugate(sup.Dim) {

        FILE *fp=fopen(Combine(Combine("truncated_wave_function/", n), trun.TotSiteNo), "wb");
        fwrite(sup.WaveFunction, sizeof(double), sup.Dim, fp);
	fclose(fp);

        delete [] trun.Old_hole;

        for(int i=0; i<trun.Block_Num_hole; i++)
                delete [] trun.Old_J[i];
        delete [] trun.Old_J;

}

//=================================================Delete superenergy=============================================
SuperEnergy::~SuperEnergy() {}
//=====================================================END========================================================
