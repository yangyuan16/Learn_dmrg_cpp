/*This program define several simple functions that will be 
 *needed in Davidson or Lanczos. 
*/

#if !defined AUXILIARYLSMATRIXEIGENV_H
#define AUXILIARYLSMATRIXEIGENV_H

#include "globalconfig.h"

class AuxiliaryLSMatrixEigenv
{
   public:
      DATATYPE innerProduct(const DATATYPE*LV, const DATATYPE*RV, const int&dim);

      double normalize(DATATYPE*V, const int&dim);
    
      void orthogonalization(const DATATYPE*, DATATYPE*, const int&dim);
    
      void constructEigenvector(const DATATYPE*,const int&, const DATATYPE*, const int&, DATATYPE*);
};
#endif
