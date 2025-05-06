/*C++ interface of complex fortran blas library
 * Author:  Zhao Jize <dmrglist@gmail.com>
 */

#include <assert.h>
#include "ZCppBlas.h"

double C_dznrm2(const int&NC, const Complex*XC, const int&INCXC)
{
   return dznrm2_(&NC,XC,&INCXC);
}

void C_zaxpy(const int&n, const Complex&za, const Complex*zx, const int&incx, Complex*zy, const int&incy)
{
   zaxpy_(&n,&za,zx,&incx,zy,&incy);
}

Complex C_zdotu(const int&n, const Complex*zx, const int&incx, const Complex*zy, const int&incy)
{
   return zdotu_(&n,zx,&incx,zy,&incy);
}

void C_zgemm(char TRANSAC, char TRANSBC, const int MC, const int NC, const int KC,Complex ALPHAC, const Complex*AC,
             const Complex*BC, const Complex BETAC, Complex* CC)
{
/* CC(i,j)=opt(AC) * opt(BC) 
 * opt(AC) is MC X KC Matrix
 * opt(BC) is KC X NC Matrix
 * CC is MC X NC Matrix
 */

   assert(TRANSAC=='N'||TRANSAC=='T'||TRANSAC=='C');
   assert(TRANSBC=='N'||TRANSBC=='T'||TRANSBC=='C');
   assert(MC>=1&&NC>=1&&KC>=1);

   const int m=NC;
   const int n=MC;
   const int k=KC;

   Complex alpha=ALPHAC;
   Complex beta=BETAC;

   const int ldc=m;

   const Complex* A=BC;
   const Complex* B=AC;   //A and B exchange positions here
   Complex* C=CC;

   char ta=TRANSBC;
   char tb=TRANSAC;  //A and B exchange transpose state

   if(TRANSAC=='N')
   {
      if(TRANSBC=='N')
      {
         const int lda=m;
         const int ldb=k;
         zgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
      }
      else
      {
         const int lda=k;
         const int ldb=k;
         zgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
      }
   }
   else
   {
      if(TRANSBC=='N')
      {
         const int lda=m;
         const int ldb=n;
         zgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
      }
      else
      {
         const int lda=k;
         const int ldb=n;
         zgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
      }
   }
}

void C_zscal(const int&nc, const Complex&zac, Complex*xc, const int&incxc)
{
   zscal_(&nc,&zac,xc,&incxc);
}

