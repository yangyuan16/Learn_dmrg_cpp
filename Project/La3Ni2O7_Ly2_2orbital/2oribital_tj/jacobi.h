// Jacobi Method 
class Jacobi { 
public : 
  Jacobi () ; 
  void Create_aa() ;
  void Free_aa() ;
  void JacobiDia() ;

  double **aa ;
  double **V ;
  double *Eig ;
  int DIM ;  // Dimension of matrix A[][]
  double trace ;
private :   
  double **copy_aa ;	// useless
  int JacobiDia2() ;
  inline void ForJacobi( double &s, double &ttau, int &i1,
                      int &i2, int &j1, int &j2) ;
  inline void ForJacobi2( double &s, double &ttau, int &i1,
                       int &i2, int &j1, int &j2) ;

  inline void OrderEig() ;
  inline void PrintEig() ;
  double *B ;	// local variable
  double *Z ; // local variable
}; 

