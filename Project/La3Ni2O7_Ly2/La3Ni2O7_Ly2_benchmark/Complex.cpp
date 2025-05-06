#include "Complex.h"

Complex operator+(const double& __x, const Complex& __y)
{
   Complex __r = __y;
   __r += __x;
   return __r;
}

Complex operator-(const double& __x, const Complex& __y)
{
   Complex __r(__x-__y.real(), -__y.imag());
   return __r; 
}

Complex operator*(const double& __x, const Complex& __y)
{
   Complex __r = __y;
   __r *= __x;
   return __r;
}

Complex operator/(const double& __x, const Complex& __y)
{
   Complex __r = __x;
   __r /= __y;
   return __r;
}

bool operator==(const double& __x, const Complex& __y)
{  
   return __x == __y.real() && double() == __y.imag(); 
}

bool operator!=(const double& __x, const Complex& __y)
{
   return __x != __y.real() || double() != __y.imag();
}

Complex operator+(const Complex& __x)
{
   return __x;
}

Complex operator-(const Complex& __x)
{
   return Complex(-__x.real(), -__x.imag());
}

double norm(const Complex&__z)
{
   const double __x = __z.real();
   const double __y = __z.imag();
   return __x * __x + __y * __y;
}

double fabs(const Complex& __z)
{
   double __x = __z.real();
   double __y = __z.imag();
   const double __s = std::max(abs(__x), abs(__y));
   if (__s == double())return __s;
   __x /= __s;
   __y /= __s;
   return __s * sqrt(__x * __x + __y * __y);
}

Complex conj(const Complex& __z)
{
   return Complex(__z.real(), -__z.imag());
}

ostream & operator << ( ostream & __os, const Complex& __x )
{
   __os << '(' << __x.real() << ',' << __x.imag() << ')';
   return __os;
}

