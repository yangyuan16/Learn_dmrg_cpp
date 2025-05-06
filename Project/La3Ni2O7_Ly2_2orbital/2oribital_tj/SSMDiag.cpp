#include <iostream>
#include <math.h>
#include <assert.h>

#ifdef DEBUG_ED_TIME
#include <unistd.h>
#include <time.h>
#include <sys/times.h>
#endif

#include "SSMDiag.h"

using namespace std;

extern "C"
{
   void zheev_(const char*jobz,const char*uplo,const int*n,Complex*a,const int*lda,double*w,Complex*work,int*lwork,
                     double*rwok,int*info);

   void dsyev_(const char*jobz,const char*uplo,const int*n,double*a,const int*lda,double*w,double*work,int*lwork,int*info);
}

int SSMDiag::SSMED(double* Matrix, int Dim, double* Eigenvalue, char Order)
{
#ifdef DEBUG_ED_TIME
   int click=sysconf(_SC_CLK_TCK);
   static struct tms tms_begin_ED;
   int time_begin_ED=times(&tms_begin_ED);
#endif

   char jobz='V';
   char uplo='U';
   const int n=Dim;
   const int lda=n;
   int info;

   int lwork=3*Dim;

   double*work=new double[lwork];
   assert(work);

   if(Order=='D')for(int i=0;i<Dim*Dim;i++)Matrix[i]=-Matrix[i];

   dsyev_(&jobz,&uplo,&n,Matrix,&lda,Eigenvalue,work,&lwork,&info);

   if(info==0)
   {
      if(Order=='D')for(int i=0;i<Dim;i++)Eigenvalue[i]=-Eigenvalue[i];

      delete []work;

#ifdef DEBUG_ED_TIME
      static struct tms tms_end_ED;
      int time_end_ED=times(&tms_end_ED);
      cout.precision(10);
      cout<<"Time for ExactDiagonalization Matrix with Dim "<<n<<" is :"<<(time_end_ED-time_begin_ED)/(1.0*click)<<endl;
#endif

      return 1;
   }
   else if(info>0)
   {
      delete []work;

      cout<<"The algorithm failed to converge!"<<endl;
      return 0;
   }
   else
   {
      delete []work;

      cout<<"The "<<info<<" parameter is illegal!"<<endl;
      return 0;
   }
}
//Matrix must arranged by col;
int SSMDiag::SSMED(Complex* Matrix,int Dim,double* Eigenvalue, char Order)
{
#ifdef DEBUG_ED_TIME
   int click=sysconf(_SC_CLK_TCK);
   static struct tms tms_begin_ED;
   int time_begin_ED=times(&tms_begin_ED);
#endif

   assert(Dim>0);

   char jobz='V';
   char uplo='U';
   const int n=Dim;
   const int lda=n;
   int info;
   int lwork=2*Dim;

   Complex*work=new Complex[lwork];
   assert(work);

   double*rwork=new double[3*Dim];
   assert(rwork);

/*Exchange row and column of Matrix,
*for fortran and c in the order of matrix are different
*/

   for(int i=0;i<Dim;i++)
   {
      for(int j=0;j<i;j++)exchange(Matrix[i*Dim+j],Matrix[j*Dim+i]);
   }

   if(Order=='D')for(int i=0;i<Dim*Dim;i++)Matrix[i]=-Matrix[i];

   zheev_(&jobz,&uplo,&n,Matrix,&lda,Eigenvalue,work,&lwork,rwork,&info);

   if(info==0)
   {
      if(Order=='D')for(int i=0;i<Dim;i++)Eigenvalue[i]=-Eigenvalue[i];

      delete []rwork;
      delete []work;

#ifdef DEBUG_ED_TIME
   static struct tms tms_end_ED;
   int time_end_ED=times(&tms_end_ED);
   cout.precision(10);
   cout<<"Time for ExactDiagonalization Matrix with Dim "<<n<<" is :"<<(time_end_ED-time_begin_ED)/(1.0*click)<<endl;
#endif

      return 1;
   }
   else if(info>0)
   {
      delete []rwork;
      delete []work;
      cout<<"The algorithm failed to converge!"<<endl;
      return 0;
   }
   else
   {
      delete []rwork;
      delete []work;
      cout<<"The "<<info<<" parameter is illegal!"<<endl;
      return 0;
   }
}
