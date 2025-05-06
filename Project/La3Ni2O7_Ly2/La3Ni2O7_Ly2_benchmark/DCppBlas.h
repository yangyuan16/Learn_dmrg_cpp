/*C++ interface of double precision fortran blas library
 *This library is distributed WITHOUT ANY WARRANTY.
 *Author:  Zhao Jize <dmrglist@gmail.com>
 */

#ifndef DCPPBLAS_H
#define DCPPBLAS_H

extern "C"
{
   double dasum_(const int*n, const double*dx, const int*incx);

   void daxpy_(const int*n, const double*da, const double*dx, const int*incx, double*dy, const int*incy);

   void dcopy_(const int*n, const double*dx, const int*incx, double*dy, const int*incy);

   double ddot_(const int*n, const double*dx, const int*incx, const double*dy,const int*incy);

   void dgemm_(const char*TRANSA, const char*TRANSB, const int*M, const int*N, const int*K,const double*ALPHA, 
               const double*A, const int*LDA, const double*B, const int*LDB, const double*BETA, double*C,  const int*LDC);

   void dgemv_(const char*TRANS, const int*M, const int*N, const double*ALPHA, const double*A, const int*LDA,
               const double*X, const int*INCX, const double*BETA, double*Y, const int*INCY);

   double dnrm2_(const int*N, const double*X, const int*INCX);

   void dscal_(const int*n, const double*da, double*x, const int*incx);

   void dspmv_(const char*UPLO, const int*N, const double*ALPHA, const double*AP, const double*X, const int*INCX, 
               const double*BETA, const double*Y, const int*INCY);
}

/* Here I write C++ interfaces to the function of the above fortran function one to one correspondence
 * 1. I assume all the matrix are stored continuously in memory, so I don't need "LDA", "LDB" ...
 * 2. Now the column and row storage in the parameters list is transformed from fortran to C convention	
 */

double C_dasum(const int&n, const double*dx, const int&incx);

void C_daxpy(const int&n, const double&da, const double*dx, const int&incx, double*dy, const int&incy);

void C_dcopy(const int&n, const double*dx, const int&incx, double*dy, const int&incy);

double C_ddot(const int&n, const double*dx, const int&incx, const double*dy,const int&incy);

void C_dgemm(char TRANSAC, char TRANSBC, const int MC, const int NC, const int KC,double ALPHAC, const double*AC,
             const double*BC, const double BETAC, double* CC);

void C_dgemm2(const char&TRANSAC, const char&TRANSBC, const int&MC, const int&NC, const int&KC,const double&ALPHAC, 
              const double*AC, const int&LDAC, const double*BC, const int&LDBC, const double&BETAC, double* CC, 
              const int&LDCC);

void C_dgemv(const char&TRANSC, const int&MC, const int&NC, const double&ALPHAC, const double*AC, const int&LDAC,
             const double*XC, const int&INCXC, const double&BETAC, double*YC, const int&INCYC);

double C_dnrm2(const int&N, const double*X, const int&INCX);

void C_dscal(const int&n, const double&da, double*x, const int&incx);

void C_dspmv(const char&UPLO, const int&N, const double&ALPHA, const double*AP, const double*X, const int&INCX,
             const double&BETA, const double*Y, const int&INCY);
#endif
