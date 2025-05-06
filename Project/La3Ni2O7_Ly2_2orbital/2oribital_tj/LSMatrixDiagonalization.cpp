#include <assert.h>
#include <iostream>

#include "LSMatrixDiagonalization.h"

using namespace std;

int LSMatrixDiagonalization::diagonalizeLSMatrix(const int&dim, const int&eigennum, const double&precision, 
                             const int&totalbasenum, const int&minibasenum, const int&maxiteration, int&inputinfo,
                             DATATYPE*eigenfunction, double*eigenvalue)
{
   int num_min=minibasenum>totalbasenum/2?totalbasenum/2:minibasenum;   //  ??
   
//   if(totalbasenum>=dim)totalbasenum=dim-1;
 
   return Davidson(eigenfunction,dim,num_min,totalbasenum,eigennum,eigenvalue,precision,maxiteration); 
 
//   cout<<"Before arnoldi"<<endl;
 
//   return Arnoldi(eigenfunction,dim,eigennum,totalbasenum,maxiteration,precision,eigenvalue);
}


    
    
        
