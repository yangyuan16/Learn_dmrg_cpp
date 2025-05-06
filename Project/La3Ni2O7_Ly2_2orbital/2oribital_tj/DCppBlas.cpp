/*C++ interface of double precision fortran blas library
 *This library is distributed WITHOUT ANY WARRANTY.
 *Author:  Zhao Jize <dmrglist@gmail.com>
 */

#include <assert.h>
#include "DCppBlas.h"

double C_dasum(const int&n, const double*dx, const int&incx)
{
   return dasum_(&n,dx,&incx);
}

void C_daxpy(const int&n, const double&da, const double*dx, const int&incx, double*dy, const int&incy)
{
   daxpy_(&n,&da,dx,&incx,dy,&incy);
}

void C_dcopy(const int&n, const double*dx, const int&incx, double*dy, const int&incy)
{
   dcopy_(&n,dx,&incx,dy,&incy); 
}

double C_ddot(const int&n, const double*dx, const int&incx, const double*dy,const int&incy)
{
   return ddot_(&n,dx,&incx,dy,&incy);
}

void C_dgemm(char TRANSAC, char TRANSBC, const int MC, const int NC, const int KC,double ALPHAC, const double*AC,
             const double*BC, const double BETAC, double* CC)
{
/* CC(i,j)=opt(AC) * opt(BC) 
 * opt(AC) is MC X KC Matrix
 * opt(BC) is KC X NC Matrix
 * CC is MC X NC Matrix
 */

   assert(TRANSAC=='N'||TRANSAC=='T');
   assert(TRANSBC=='N'||TRANSBC=='T');

   const int M=NC;
   const int N=MC;
   const int K=KC;

   const int LDC=M;

   const double* A=BC;
   const double* B=AC;   //A and B exchange positions here

   const char TRANSA=TRANSBC;
   const char TRANSB=TRANSAC;  //A and B exchange transpose state

   if(TRANSAC=='N')
   {
      if(TRANSBC=='N')
      {
         const int LDA=M;
         const int LDB=K;
         dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHAC, A, &LDA, B, &LDB, &BETAC, CC, &LDC);
      }
      else 
      {
         const int LDA=K;
         const int LDB=K;
         dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHAC, A, &LDA, B, &LDB, &BETAC, CC, &LDC);
      } 
   }
   else
   {
      if(TRANSBC=='N')
      {
         const int LDA=M;
         const int LDB=N;
         dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHAC, A, &LDA, B, &LDB, &BETAC, CC, &LDC);         
      }
      else
      {
         const int LDA=K;
         const int LDB=N;
         dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHAC, A, &LDA, B, &LDB, &BETAC, CC, &LDC); 
      }
   }
}

void C_dgemm2(const char&TRANSAC, const char&TRANSBC, const int&MC, const int&NC, const int&KC,const double&ALPHAC,
              const double*AC, const int&LDAC, const double*BC, const int&LDBC, const double&BETAC, double* CC, 
              const int&LDCC)
{
/* CC(i,j)=opt(AC) * opt(BC) 
 * opt(AC) is MC X KC Matrix
 * opt(BC) is KC X NC Matrix
 * CC is MC X NC Matrix
 *
 * TRANSAC=='N', 'T'
 * TRANSBC=='N', 'T'
 */
/*
   const int M=NC;
   const int N=MC;
   const int K=KC;

   const int LDC=LDCC;

   const double* A=BC;
   const double* B=AC;   //A and B exchange positions here

   const char TRANSA=TRANSBC;
   const char TRANSB=TRANSAC;  //A and B exchange transpose state

   const int LDA=LDBC;
   const int LDB=LDAC;
*/
   dgemm_(&TRANSBC, &TRANSAC, &NC, &MC, &KC, &ALPHAC, BC, &LDBC, AC, &LDAC, &BETAC, CC, &LDCC);
}

void C_dgemv(const char&TRANSC, const int&MC, const int&NC, const double&ALPHAC, const double*AC, const int&LDAC,
             const double*XC, const int&INCXC, const double&BETAC, double*YC, const int&INCYC)
{
   assert(TRANSC=='N'||TRANSC=='T');

   const int M=NC;
   const int N=MC;
   assert(LDAC==M);

   if(TRANSC=='N')
   {
      const char TRANS='T';
      dgemv_(&TRANS,&M,&N,&ALPHAC,AC,&LDAC,XC,&INCXC,&BETAC,YC,&INCYC);
   }
   else
   {
      const char TRANS='N';
      dgemv_(&TRANS,&M,&N,&ALPHAC,AC,&LDAC,XC,&INCXC,&BETAC,YC,&INCYC);
   }
}

double C_dnrm2(const int&NC, const double*XC, const int&INCXC)
{
   return dnrm2_(&NC,XC,&INCXC); 
}

void C_dscal(const int&nc, const double&dac, double*xc, const int&incxc)
{
   dscal_(&nc,&dac,xc,&incxc);
}

void C_dspmv(const char&UPLOC, const int&NC, const double&ALPHAC, const double*APC, const double*XC, const int&INCXC,
             const double&BETAC, const double*YC, const int&INCYC)
{
   assert(UPLOC=='U'||UPLOC=='L');
   if(UPLOC=='U')
   {
      const char UPLO='L';
      dspmv_(&UPLO,&NC,&ALPHAC,APC,XC,&INCXC,&BETAC,YC,&INCYC);
   }
   else
   {
      const char UPLO='U';
      dspmv_(&UPLO,&NC,&ALPHAC,APC,XC,&INCXC,&BETAC,YC,&INCYC);
   }
}

