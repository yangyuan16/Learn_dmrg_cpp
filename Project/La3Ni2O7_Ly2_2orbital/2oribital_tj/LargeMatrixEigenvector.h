/*Block Jacobi-Davidson   Copyright 2001,All right reserved
*This class is designed to get some lowest eigenvalues 
*and corresponding eigenvectors of a large scale matrix
*/
 
#ifndef LARGEMATRIXEIGENVECTOR_H
#define LARGEMATRIXEIGENVECTOR_H

#include "globalconfig.h"
#include "SSMDiag.h"
#include "AuxiliaryLSMatrixEigenv.h"

class LargeMatrixEigenvector : private AuxiliaryLSMatrixEigenv, virtual public SSMDiag 
{
   public:
#ifdef DEBUG_DAVIDSON_DIAG
      int totalProcessNum;
      int nonDavidsonNum; 
#endif  
   public:
      virtual void H_V(const DATATYPE*,DATATYPE*, const int&)=0;
    
      virtual void getMatrixDiagElement(DATATYPE*, const int&)=0;
    
      int Davidson(DATATYPE*ritzvector,const int&dim, const int&N_Min, const int&N_Max,const int&N_Wanted, 
                            double*ritzvalue, const double&tol, const int&maxiteration);

   private:
      int initializeBaseVector(DATATYPE*, DATATYPE*, const int&, const int&, DATATYPE*,const int&, const int&);

      void constructProjectSubSpace(DATATYPE*V, const int&dim, DATATYPE*sv, const int&svdim, const int&shift, 
                                    const int&svnum, const int&NREQ, DATATYPE*HV, DATATYPE*ritzvector, const int&num);
};
#endif
