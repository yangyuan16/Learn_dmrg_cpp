#include <math.h>

#include "DCppBlas.h"
#include "ZCppBlas.h"
#include "AuxiliaryLSMatrixEigenv.h"

DATATYPE AuxiliaryLSMatrixEigenv::innerProduct(const DATATYPE*LV,const DATATYPE*RV, const int&dim)
{
#ifdef DOUBLEDATATYPE
   return C_ddot(dim, LV, 1, RV, 1);
#else
   return C_zdotu(dim,LV, 1, RV, 1);
#endif
}

/*Function normalize(DATATYPE*,int) normalizes the vector */

double AuxiliaryLSMatrixEigenv::normalize(DATATYPE*vector, const int&dim)
{
#ifdef DOUBLEDATATYPE
   double sum=C_dnrm2(dim,vector,1);
   double invsum=1.0/sum;
   C_dscal(dim,invsum,vector,1);
#else
   double sum=C_dznrm2(dim,vector,1);
   double invsum=1.0/sum;
   C_zscal(dim,invdim,vector,1);
#endif

   return sum;
}

/*orthogonalization orthogonalizes the vector "base" to V. and V is a normal vectors*/

void AuxiliaryLSMatrixEigenv::orthogonalization(const DATATYPE*V, DATATYPE*base, const int&dim)
{
   DATATYPE sum=innerProduct(V, base, dim);
#ifdef DOUBLEDATATYPE
   C_daxpy(dim, -sum, V, 1, base, 1);
#else
   C_zaxpy(dim, -sum, V, 1, base, 1);
#endif
}

void AuxiliaryLSMatrixEigenv::constructEigenvector(const DATATYPE*V, const int&dim, const DATATYPE*sv, const int&svdim, 
                                                   DATATYPE*eigenvector)
{
#ifdef DOUBLEDATATYPE
   C_dgemv('T', svdim, dim, 1.0, V, dim, sv, 1, 0.0, eigenvector, 1);
#else
   for(int i=0;i<dim;i++)
   {
      eigenvector[i]=0.0;
      for(int j=0;j<svdim;j++)
      {
          eigenvector[i]+=sv[j]*V[j*dim+i];
      }
   }
#endif
}


