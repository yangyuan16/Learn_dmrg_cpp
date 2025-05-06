/*We define two data type here, be sure to choose 
* one of them. 
*/

/*
#ifndef DATATYPE
 #define COMPLEXDATATYPE
 #define DATATYPE  complex
 #include "complex.h"
#endif 
*/

#ifndef DATATYPE
 #define DOUBLEDATATYPE
 #define DATATYPE double
#endif

#if !defined FERMION
#define FERMION   978013 
#endif

#if !defined BOSON 
#define BOSON 101762
#endif

#if !defined STATISTICTYPE      
#define STATISTICTYPE  FERMION
//#define STATISTICTYPE  BOSON
#endif

