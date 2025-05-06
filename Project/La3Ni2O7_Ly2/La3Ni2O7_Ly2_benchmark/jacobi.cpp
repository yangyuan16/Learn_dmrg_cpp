#include<iostream> 
#include<fstream> 
#include<iomanip>
using namespace std;
#include<math.h> 
#include<stdlib.h> 
#include<assert.h> 
 
#include "jacobi.h" 
//================
Jacobi::Jacobi () {}
/*c   
     computes all eigenvalues and eigenvactors of a real 
     symmetric matrix aa, witch is of size N by N, stored in a 
     phycal NP by NP array.  On output, elements of aa above the 
     diagonal are destroyed.  Eig returns the eigenvalues of aa in 
     its first N elements.  V is a matrix with the same logical 
     and physical dimensions as aa whose columns contain, on 
     output, the normalized eigenvectors of aa.  NROT retures the 
     number of Jacobi rotations which were required. 
*/  
//====================
void Jacobi::JacobiDia( ) {
  JacobiDia2() ;
  OrderEig() ;
  PrintEig() ;
} 
//=========================
inline int Jacobi::JacobiDia2( ) {
  int ITR_MAX = 50 ; // maximum number of iterations 
  int i, j, ip, iq, NROT = 0 ; 
  double tresh, t, g, h, s, ttau, theta, sm, c ; 
               // initialize V to the identity matrix 
  for (i = 0; i < DIM; i++) 
  for (j = 0; j < DIM; j++)
     V[i][j] = (j == i ? 1.0 : 0.0) ;
               // initialize B and Eig to the diagonal of aa 
  for (i = 0; i < DIM; i++) { 
     B[i] = aa[i][i] ;
     Eig[i] = B[i] ;
     Z[i] = 0. ;
  } 
 for ( i = 0; i < ITR_MAX; i++ ) {   
               // Sum off-diagonal elements 
    sm = 0. ;
    for ( ip = 0; ip < DIM-1; ip++) 
    for ( iq = ip + 1; iq < DIM; iq++)
        sm += fabs( aa[iq][ip] ) ; 
    if( sm == 0. ) return 0; 
    if ( i < 3 )  
      tresh = 0.2 * sm / (DIM * DIM) ; 
    else 
      tresh = 0. ; 
    for ( ip = 0; ip < DIM-1; ip++) 
    for ( iq = ip + 1; iq < DIM; iq++) {
        g = 100. * fabs( aa[iq][ip] ) ; 
        if ( (i > 3) && ( fabs(Eig[ip])+g == fabs(Eig[ip]) ) 
           && ( fabs(Eig[iq]) + g == fabs(Eig[iq]) ) )  
          aa[iq][ip] = 0. ; 
        else if( fabs(aa[iq][ip]) > tresh ) { 
          h = Eig[iq] - Eig[ip] ; 
          if( fabs(h) + g == fabs(h) ) t = aa[iq][ip] / h ; 
          else { 
            theta = 0.5 * h / aa[iq][ip] ; 
            t = 1. / (fabs(theta) + sqrt(1. + theta*theta)) ; 
            if( theta < 0. ) t = - t ; 
          } 
          c = 1. / sqrt(1. + t * t) ; 
          s = t * c ; 
          ttau = s / (1. + c) ; 
          h = t * aa[iq][ip] ; 
          Z[ip] -= h ; 
          Z[iq] += h ; 
          Eig[ip] -= h ; 
          Eig[iq] += h ; 
          aa[iq][ip] = 0. ; 
          for (j = 0; j < ip; j++)  
            ForJacobi( s, ttau, ip, j, iq, j) ; 
          for ( j=ip+1; j < iq; j++ )   
            ForJacobi( s, ttau, j, ip, iq, j) ; 
          for ( j = iq+1; j < DIM; j++ )  
            ForJacobi( s, ttau, j, ip, j, iq) ; 
          for ( j = 0; j < DIM; j++)  
            ForJacobi2( s, ttau, ip, j, iq, j) ; 
          NROT++ ; 
      }   
    } 
    for (ip = 0; ip < DIM; ip++) { 
      B[ip] += Z[ip] ; 
      Eig[ip] = B[ip] ; 
      Z[ip] = 0. ; 
    } 
  }  
  return NROT ; 
}  
//===================================================
inline void Jacobi::ForJacobi(double &s, double &ttau, int &i1,
                         int &i2, int &j1, int &j2) {
  double g = aa[i1][i2] ; 
  double h = aa[j1][j2] ; 
  aa[i1][i2] -= s * (h + g * ttau) ; 
  aa[j1][j2] += s * (g - h * ttau) ; 
} 
//===================================================
inline void Jacobi::ForJacobi2(double &s, double &ttau, int &i1,  
                           int &i2, int &j1, int &j2) {
  double g = V[i1][i2] ; 
  double h = V[j1][j2] ; 
  V[i1][i2] -= s * (h + g * ttau) ; 
  V[j1][j2] += s * (g - h * ttau) ; 
} 
//=========================
inline void Jacobi::OrderEig( ) {
  for ( long i = 0 ; i < DIM-1 ; i++ )
  for ( long j = i + 1 ; j < DIM ; j++ )
      if ( Eig[i] > Eig[j] ) {
         double e = Eig[i] ;
         Eig[i] = Eig[j] ;
         Eig[j] = e ;
         for ( long k = 0 ; k < DIM ; k++ ) {
              e = V[i][k] ;
              V[i][k] = V[j][k] ;
              V[j][k] = e ;
         }
      }
} 
//========================
inline void Jacobi::PrintEig( ) {
   for ( long i = 0 ; i < DIM ; i++ ) {
       trace += Eig[i] ;
   }
} 
//=====================
void Jacobi::Create_aa( ) {
  aa = new double*[DIM] ;
  V = new double*[DIM] ;
  for ( int i = 0 ; i < DIM ; i++ ) {
    aa[i] = new double[DIM] ;
    V[i] = new double[DIM] ;
    for ( int j = 0 ; j < DIM ; j++ )
        aa[i][j] = (double) 0 ;
  }
  Eig = new double[DIM] ; assert( Eig != NULL ) ;
  B = new double[DIM] ; assert( B != NULL ) ;
  Z = new double[DIM] ; assert( Z != NULL ) ;
} 
//==================
void Jacobi::Free_aa( ) {
  for ( long i = 0 ; i < DIM ; i++ ) {
    free( aa[i] ) ;
    free( V[i] ) ;
  }
  free( aa ) ;
  free( V ) ;
  delete [] Eig ;
  delete [] B ;
  delete [] Z ;
}
//======================end
