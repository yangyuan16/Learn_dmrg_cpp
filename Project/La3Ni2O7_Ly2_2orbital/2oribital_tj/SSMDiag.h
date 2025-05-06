#ifndef SSMDIAG_H
#define SSMDIAG_H

#include "Complex.h"

class SSMDiag
{
   public:
      int SSMED(double*  matrix, int dim, double* eigenvalue, char order='D');
      int SSMED(Complex* matrix, int dim, double* eigenvalue, char order='D');
   private:
      inline void exchange(int& m1,int& m2){int s=m1;m1=m2;m2=s;}
      inline void exchange(double &m1,double &m2){double s=m1;m1=m2;m2=s;}
      inline void exchange(Complex&m1,Complex&m2){Complex s=m1;m1=m2;m2=s;}
};

#endif

