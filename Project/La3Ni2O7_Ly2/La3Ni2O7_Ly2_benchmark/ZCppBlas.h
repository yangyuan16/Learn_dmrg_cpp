/*C++ interface of complex fortran blas library
 * Author:  Zhao Jize <dmrglist@gmail.com>
 */

#ifndef ZCPPBLAS_H
#define ZCPPBLAS_H

#include "Complex.h"

extern "C"
{
   double dznrm2_(const int*N, const Complex*X, const int*INCX);

   void zaxpy_(const int*n, const Complex*da, const Complex*dx, const int*incx, Complex*dy, const int*incy);

   Complex zdotu_(const int*n,const Complex*zx, const int*incx, const Complex*zy, const int*incy);

   void zgemm_( const char*TRANSA, const char*TRANSB, const int*M, const int*N, const int*K, Complex*ALPHA, const Complex*A, 
                const int*LDA, const Complex*B, const int*LDB, const Complex*BETA, Complex*C, const int*LDC );

   void zscal_(const int*n, const Complex*za, Complex*x, const int*incx);
}

/* Here I write C++ interfaces to the function of the above fortran function one to one correspondence
 * 1. I assume all the matrix are stored continuously in memory, so I don't need "LDA", "LDB" ...
 * 2. Now the column and row storage in the parameters list is transformed from fortran to C convention	
 */

double C_dznrm2(const int&N, const Complex*X, const int&INCX);

void C_zaxpy(const int&n, const Complex&za, const Complex*zx, const int&incx, Complex*zy, const int&incy);

Complex C_zdotu(const int&n, const Complex*zx, const int&incx, const Complex*zy, const int&incy);

void C_zgemm(char TRANSAC, char TRANSBC, const int MC, const int NC, const int KC,Complex ALPHAC, const Complex*AC,
             const Complex*BC, const Complex BETAC, Complex*CC);

void C_zscal(const int&n, const Complex&za, Complex*x, const int&incx);
#endif
