#if !defined LSMATRIXDIAGONALIZATION_H
#define LSMATRIXDIAGONALIZATION_H

#include "LargeMatrixEigenvector.h"

class LSMatrixDiagonalization:LargeMatrixEigenvector
{
   public:
      int diagonalizeLSMatrix(const int&dim, const int&eigennum, const double&precision, const int&totalbasenum,
                              const int&minibasenum, const int&maxiteration, int&inputinfo, DATATYPE*eigenfunction,
                              double*eigenvalue);
	    
      virtual void H_V(const DATATYPE*, DATATYPE*,const int&dim)=0;

      virtual void getMatrixDiagElement(DATATYPE*,const int&dim)=0;
};
#endif
