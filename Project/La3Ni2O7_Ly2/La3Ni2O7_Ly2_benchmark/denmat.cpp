#include<iostream>
#include<iomanip>
#include<fstream>
using namespace std;
#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

#include"common.h"
#include"sub.h"
#include"super.h"
#include"denmat.h"

//=================================================BLAS ROUTINES=================================================
extern "C" {
void dgemm_(char *transa, char *transb, const int *m, const int *n, const int *k, const double *alpha, double *a, const int *lda, double *b, const int *ldb, const double *beta, double *c, const int *ldc);
//===============================================================================================================
//==============================================LAPACK SUBROUTINE================================================
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
}
//===============================================================================================================

//================================Reduced density matrix of blocks===============================================
DenMat::DenMat(const int &f, Parameter &para, Super &sup, Sub &sysnew) {

        CreateSpace(para, sup);

        if(f==1) {
                DenMat_Jacobi_sysnew(sup, sysnew);
        }
        else if(f==2) {
                DenMat_Jacobi_envnew(sup, sysnew);
        }

        FreeSpace();

}

//========================================Create Space for wave functions========================================
inline void DenMat::CreateSpace(Parameter &para, Super &sup) {
  //------Create block space-------------------------------------------------------------------------------------
        int index, J_min, J_max, J_num, oldH_sys, oldH_env, oldJ_sys, oldJ_env, position_old, position_new;

        trans_N='N';    trans_T='T';
        index=0;        BlockNumber=0;
        beta=1.0;

	QuantumNumber_hole = sup.QuantumNumber_hole;
	QuantumNumber_J = sup.QuantumNumber_J;

	for(int n_1=0; n_1<sup.sysnew_space->Block_Num_hole; n_1++)
	for(int n_2=0; n_2<sup.envnew_space->Block_Num_hole; n_2++)
        if(sup.sysnew_space->Num_hole_block[n_1]+sup.envnew_space->Num_hole_block[n_2]==QuantumNumber_hole) {

	        for(int j_sysnew=0; j_sysnew<sup.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
        	for(int j_envnew=0; j_envnew<sup.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                	J_min=abs(sup.sysnew_space->Value_J_block[n_1][j_sysnew]-sup.envnew_space->Value_J_block[n_2][j_envnew]);
	                J_max=(sup.sysnew_space->Value_J_block[n_1][j_sysnew]+sup.envnew_space->Value_J_block[n_2][j_envnew]);
        	        J_num=(J_max-J_min)/2+1;

	                for(int n=0; n<J_num; n++)
        	        if(QuantumNumber_J==(J_min+2*n)) 
                	        BlockNumber++;

        	}

	}

	H_sysnew=new int [BlockNumber];	    H_envnew=new int [BlockNumber];
        J_sysnew=new int [BlockNumber];     J_envnew=new int [BlockNumber];     Dim_block=new int [BlockNumber];

        for(int i=0; i<BlockNumber; i++) {
                H_sysnew[i]=0;  H_envnew[i]=0;  J_sysnew[i]=0;  J_envnew[i]=0;  Dim_block[i]=0;
        }

  //------Initialize the values of above angular numbers---------------------------------------------------------
	for(int n_1=0; n_1<sup.sysnew_space->Block_Num_hole; n_1++)
	for(int n_2=0; n_2<sup.envnew_space->Block_Num_hole; n_2++)
        if(sup.sysnew_space->Num_hole_block[n_1]+sup.envnew_space->Num_hole_block[n_2]==QuantumNumber_hole) {

		for(int j_sysnew=0; j_sysnew<sup.sysnew_space->Num_J_hole_block[n_1]; j_sysnew++)
        	for(int j_envnew=0; j_envnew<sup.envnew_space->Num_J_hole_block[n_2]; j_envnew++) {

                	J_min=abs(sup.sysnew_space->Value_J_block[n_1][j_sysnew]-sup.envnew_space->Value_J_block[n_2][j_envnew]);
	                J_max=(sup.sysnew_space->Value_J_block[n_1][j_sysnew]+sup.envnew_space->Value_J_block[n_2][j_envnew]);
        	        J_num=(J_max-J_min)/2+1;

	                for(int n=0; n<J_num; n++)
        	        if(QuantumNumber_J==(J_min+2*n)) {

				H_sysnew[index]=n_1;		H_envnew[index]=n_2;
                	        J_sysnew[index]=j_sysnew;       J_envnew[index]=j_envnew;
                        	Dim_block[index++]=sup.sysnew_space->Dim_J_block[n_1][j_sysnew]*sup.envnew_space->Dim_J_block[n_2][j_envnew];

                	}
        	}

	}

        wavefunc_new=new double * [BlockNumber];

        for(int i=0; i<BlockNumber; i++) {

		wavefunc_new[i]=new double [Dim_block[i]];

                for(int j=0; j<Dim_block[i]; j++)
	                wavefunc_new[i][j]=(double) 0;

	}

        for(int i=0; i<BlockNumber; i++) {

		for(int j_n=0; j_n<3; j_n++)
                for(int j_e=0; j_e<3; j_e++)
		if((oldH_sys=sup.sysnew_space->Hole_blockOld[H_sysnew[i]][j_n][J_sysnew[i]])!=-1 && (oldH_env=sup.envnew_space->Hole_blockOld[H_envnew[i]][j_e][J_envnew[i]])!=-1 && (oldJ_sys=sup.sysnew_space->J_blockOld[H_sysnew[i]][j_n][J_sysnew[i]])!=-1 && (oldJ_env=sup.envnew_space->J_blockOld[H_envnew[i]][j_e][J_envnew[i]])!=-1) {

			for(int j=0; j<sup.BlockNumber_for_TargetBlock; j++)
                        if(sup.H_sysnew[j]==H_sysnew[i] && sup.H_envnew[j]==H_envnew[i] && sup.J_sysnew[j]==J_sysnew[i] && sup.J_envnew[j]==J_envnew[i] && sup.H_sys[j]==oldH_sys && sup.H_env[j]==oldH_env && sup.J_sys[j]==oldJ_sys && sup.J_env[j]==oldJ_env) {

				for(int a_env=0; a_env<sup.env_space->Dim_J_block[oldH_env][oldJ_env]; a_env++)
                                for(int a_sys=0; a_sys<sup.sys_space->Dim_J_block[oldH_sys][oldJ_sys]; a_sys++) {

	                                position_old=a_env*sup.sys_space->Dim_J_block[oldH_sys][oldJ_sys]+a_sys;
                                        position_new=(sup.envnew_space->Start[H_envnew[i]][j_e][J_envnew[i]]+a_env)*sup.sysnew_space->Dim_J_block[H_sysnew[i]][J_sysnew[i]]+sup.sysnew_space->Start[H_sysnew[i]][j_n][J_sysnew[i]]+a_sys;
                                        wavefunc_new[i][position_new]=sup.WaveFunction_block[j][position_old];

                                }

                        }

                }

	}

}

//================================================================================================================
//================================================================================================================
//============================================DenMat_Jacobi_sysnew===============================================
inline void DenMat::DenMat_Jacobi_sysnew(Super &sup, Sub &sysnew) {
  //------Diagonalization by Lapack subroutine dsyev-------------------------------------------------------------

	double entropy=0.0;

        FILE *g=fopen(Combine("entanglement/spectrum/sys-ground_state-", sysnew.TotSiteNo), "w+");

	for(int jn=0; jn<sysnew.Block_Num_hole; jn++)
        for(int js=0; js<sysnew.Num_J_hole_block[jn]; js++)
        if(sysnew.Dim_J_block[jn][js]!=0) {

		double constant=(sysnew.Value_J_block[jn][js]+1.0);

                Initial_Lapack(sysnew.Dim_J_block[jn][js]);
                Find_aa_sysnew(sup, jn, js);
                dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

                for(int i=(sysnew.Dim_J_block[jn][js]-1), it=0; i>=0; i--, it++) {

			sysnew.dm_eig[jn][js][it]=fabs(w[i]);//w:from small to large, dm_eig: from large to small

			if(fabs(w[i])>1.e-15)
                       		entropy+=-fabs(w[i])*log(fabs(w[i]))*constant;

			if(fabs(w[i])>1.e-15) {
				fprintf(g, "%d\t", sysnew.Num_hole_block[jn]);
				fprintf(g, "%d\t", sysnew.Value_J_block[jn][js]);
        	                fprintf(g, "%f\n", -log(fabs(w[i])));
			}
	
			for(int j=0; j<sysnew.Dim_J_block[jn][js]; j++)
                        	sysnew.dm_wave[jn][js][it*sysnew.Dim_J_block[jn][js]+j]=a[i*sysnew.Dim_J_block[jn][js]+j];

		}

	                delete [] a;         delete [] work;      delete [] w;

	}

        fclose(g);

        FILE *f=fopen("entanglement/entropy/sys-ground_state", "a+");
        fprintf(f, "%d\t", sysnew.TotSiteNo);
        fprintf(f, "%f\n", entropy);
        fclose(f);

}

//============================================DenMat_Jacobi_envnew===============================================
inline void DenMat::DenMat_Jacobi_envnew(Super &sup, Sub &envnew) {

	double entropy=0.0;

        FILE *g=fopen(Combine("entanglement/spectrum/env-ground_state-", envnew.TotSiteNo), "w+");

	for(int jn=0; jn<envnew.Block_Num_hole; jn++)
        for(int je=0; je<envnew.Num_J_hole_block[jn]; je++)
        if(envnew.Dim_J_block[jn][je]!=0) {

		double constant=(envnew.Value_J_block[jn][je]+1.0); 

                Initial_Lapack(envnew.Dim_J_block[jn][je]);
                Find_aa_envnew(sup, jn, je);
		dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

                for(int i=(envnew.Dim_J_block[jn][je]-1), it=0; i>=0; i--, it++) {

	                envnew.dm_eig[jn][je][it]=fabs(w[i]);//w: from small to large, dm_eig: from large to small

                        if(fabs(w[i])>1.e-15)
          	              entropy+=-fabs(w[i])*log(fabs(w[i]))*constant;

                        if(fabs(w[i])>1.e-15) {
		                fprintf(g, "%d\t", envnew.Num_hole_block[jn]);
        	                fprintf(g, "%d\t", envnew.Value_J_block[jn][je]);
                	        fprintf(g, "%f\n", -log(fabs(w[i])));
			}

                        for(int j=0; j<envnew.Dim_J_block[jn][je]; j++)
		  	      envnew.dm_wave[jn][je][it*envnew.Dim_J_block[jn][je]+j]=a[i*envnew.Dim_J_block[jn][je]+j];
                        
		}

                delete [] a;    delete [] work;      delete [] w;

	}

        fclose(g);

        FILE *f=fopen("entanglement/entropy/env-ground_state", "a+");
        fprintf(f, "%d\t", envnew.TotSiteNo);
        fprintf(f, "%f\n", entropy);
        fclose(f);

}

//================================Initialize parameters for Lapack subroutine dsyev==============================
inline void DenMat::Initial_Lapack(int &Dim) {

	jobz='V';       //eigenvalues and eigenvectors are computed
        uplo='U';       //a stores the upper triangular part of A
        n=Dim;          //The order of the matrix A
        lda=Dim;        //The first dimension of the array a
        lwork=40*n;     //The dimension of the array work. Constraint:lwork>=max(1, 3n-1), check work(1)!
        a=new double [lda*n];

        work=new double [lwork];

        w=new double [n];//contains the eigenvalues of the matrix A in ascending order!!!

        for(int i=0; i<lda*n; i++)
        	a[i]=0.0;

        for(int i=0; i<lwork; i++)
                work[i]=0.0;

        for(int i=0; i<n; i++)
                w[i]=0.0;

}

//================Find the matrix elements for each subblocks of the sysnew reduced density matrices=============
inline void DenMat::Find_aa_sysnew(Super &sup, int &jn, int &js) {

	factor=1.0/(sup.sysnew_space->Value_J_block[jn][js]+1.0);
        for(int i=0; i<BlockNumber; i++)
        if(H_sysnew[i]==jn && J_sysnew[i]==js)
        	dgemm_(&trans_N, &trans_T, &sup.sysnew_space->Dim_J_block[jn][js], &sup.sysnew_space->Dim_J_block[jn][js], &sup.envnew_space->Dim_J_block[H_envnew[i]][J_envnew[i]], &factor, wavefunc_new[i], &sup.sysnew_space->Dim_J_block[jn][js], wavefunc_new[i], &sup.sysnew_space->Dim_J_block[jn][js], &beta, a, &sup.sysnew_space->Dim_J_block[jn][js]);

}

//================Find the matrix elements for each subblocks of the envnew reduced density matrices=============
inline void DenMat::Find_aa_envnew(Super &sup, int &jn, int &je) {

	factor=1.0/(sup.envnew_space->Value_J_block[jn][je]+1.0);
        for(int i=0; i<BlockNumber; i++)
        if(H_envnew[i]==jn && J_envnew[i]==je)
		dgemm_(&trans_T, &trans_N, &sup.envnew_space->Dim_J_block[jn][je], &sup.envnew_space->Dim_J_block[jn][je], &sup.sysnew_space->Dim_J_block[H_sysnew[i]][J_sysnew[i]], &factor, wavefunc_new[i], &sup.sysnew_space->Dim_J_block[H_sysnew[i]][J_sysnew[i]], wavefunc_new[i], &sup.sysnew_space->Dim_J_block[H_sysnew[i]][J_sysnew[i]], &beta, a, &sup.envnew_space->Dim_J_block[jn][je]);

}

//===================================================Free Space==================================================
inline void DenMat::FreeSpace() {

	for(int i=0; i<BlockNumber; i++)
		delete [] wavefunc_new[i];
        delete [] wavefunc_new;

        delete [] H_sysnew;  delete [] H_envnew;  delete [] J_sysnew;  delete [] J_envnew;  delete [] Dim_block;

}

//===============================================================================================================
DenMat::~DenMat() {}
//===================================================END==========================================================

