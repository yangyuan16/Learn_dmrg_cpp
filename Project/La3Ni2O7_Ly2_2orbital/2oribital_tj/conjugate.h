/*  
   Conjugate gradient method for diagonalizing  
   a sparse symmetric matrix 
*/ 

class Conjugate { 
public:  
  Conjugate(const int &Dim) ;
  ~Conjugate( ) ;

  void NormTo1 ( double *f ) ; // normalize <f|f> to 1
  double f1timesf2(const double *f, const double *g) ;
  int abc_2(const int &iter) ;
  void abc_4() ;
  double *f0 ;
  double *f1 ;
  double *f2 ;
  double *f3 ;
  int Dim ;
  double eng ;   // EigenValue of Hamiltonian
private: 
  double y00 ;
  double y01 ;
  double y02 ;
  double y22 ;
  double x33_old ;
  inline void CreateSpace() ;
  inline void FreeSpace() ;
  inline void Ini_f0() ;
} ; 

