#include<string.h>
#include<iostream>
using namespace std;

char *Combine(char *s1, long ks) { 

  long km = ks ;
  int len = strlen( s1 ) + 1  ;
  while( ks /= 10 ) len++ ;
 
  char *str = new char[len + 1] ;
  str[len] = '\0' ;
  for ( unsigned int i = 0 ; i < strlen(s1) ; i++ )
     str[i] = s1[i] ;
  do { str[--len] = (km % 10) + '0' ;
  } while (km /= 10) ;
  return str ;

  delete [] str;

} 

char *Combine(char *s1, char *s2) { 
  int len = strlen( s1 ) + strlen( s2 )  ;
  char *str = new char[len + 1] ;
  unsigned int i = 0, j = 0 ;
  while ( i < strlen(s1) ) str[i] = s1[i++] ;
  while ( j < strlen(s2) ) str[i++] = s2[j++] ;
       str[len] = '\0' ;
  return str ;

  delete [] str;

} 
