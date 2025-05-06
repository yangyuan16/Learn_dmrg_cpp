// The template and inlines for the -*- C++ -*- complex number classes.
//
// // Copyright (C) 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005,
// // 2006, 2007, 2008, 2009
// // Free Software Foundation, Inc.
// //
// // This file is part of the GNU ISO C++ Library.  This library is free
// // software; you can redistribute it and/or modify it under the
// // terms of the GNU General Public License as published by the
// // Free Software Foundation; either version 3, or (at your option)
// // any later version.
//
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// // GNU General Public License for more details.
//
// // Under Section 7 of GPL version 3, you are granted additional
// // permissions described in the GCC Runtime Library Exception, version
// // 3.1, as published by the Free Software Foundation.
//
// // You should have received a copy of the GNU General Public License and
// // a copy of the GCC Runtime Library Exception along with this program;
// // see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// // <http://www.gnu.org/licenses/>.

#ifndef COMPLEX_H
#define COMPLEX_H

using namespace std;

#include <iostream>
#include <cmath>

class Complex
{
   private:
      double _M_real, _M_imag;

   public:
      inline Complex(const double& __r= 0.0, const double& __i = 0.0):_M_real(__r), _M_imag(__i) { }

      inline Complex(const Complex&__z) : _M_real(__z.real()), _M_imag(__z.imag()) { } 

      inline double real() const { return _M_real; }

      inline double imag() const { return _M_imag; }

      inline Complex&operator+=(const double& __t){ _M_real += __t; return *this; }

      inline Complex& operator-=(const double& __t){ _M_real -= __t; return *this; }

      inline Complex& operator=(const double& __t) { _M_real = __t; _M_imag = double(); return *this; }

      inline Complex&operator*=(const double& __t){ _M_real *= __t; _M_imag *= __t; return *this; }

      inline Complex&operator/=(const double& __t){ _M_real /= __t; _M_imag /= __t; return *this; }

      inline Complex&operator=(const Complex& __z){ _M_real = __z.real(); _M_imag = __z.imag(); return *this; }
 
      inline Complex&operator+=(const Complex& __z){ _M_real += __z.real(); _M_imag += __z.imag(); return *this; }

      inline Complex&operator-=(const Complex& __z){ _M_real -= __z.real(); _M_imag -= __z.imag(); return *this; }

      inline Complex&operator*=(const Complex& __z)
      { 
         const double __r = _M_real * __z.real() - _M_imag * __z.imag();
         _M_imag = _M_real * __z.imag() + _M_imag * __z.real();
         _M_real = __r;
         return *this;
      }
          
      inline Complex&operator/=(const Complex& __z)
      {
         const double __r = _M_real * __z.real() + _M_imag * __z.imag();
         const double __n = norm(__z);
         _M_imag = (_M_imag * __z.real() - _M_real * __z.imag()) / __n;
         _M_real = __r / __n;
         return *this;
      }

      inline Complex operator+(const Complex& __y){ Complex __r = *this; __r += __y; return  __r; }

      inline Complex operator+(const double& __y){ Complex __r = *this; __r += __y; return __r; } 
 
      inline Complex operator-(const Complex& __y){ Complex __r = *this; __r -= __y; return __r; }

      inline Complex operator-( const double& __y){ Complex __r = *this; __r -= __y; return __r; }

      inline Complex operator*(const Complex& __y)const{ Complex __r = *this; __r *= __y; return __r; }

      inline Complex operator*(const double& __y)const{ Complex __r = *this; __r *= __y; return __r;}

      inline Complex operator/(const Complex& __y){ Complex __r = *this; __r/= __y; return __r;}

      inline Complex operator/(const double& __y){ Complex __r = *this; __r/= __y; return __r; }

      inline bool operator==(const Complex& __y){ return _M_real == __y.real() && _M_imag == __y.imag(); }

      inline bool operator==(const double& __y){ return _M_real == __y && _M_imag == double(); }

      inline bool operator!=(const Complex& __y){ return _M_real != __y.real() || _M_imag != __y.imag(); }

      inline bool operator!=(const double& __y) { return _M_real != __y || _M_imag != double();}

 
      friend Complex operator+(const double& __x, const Complex& __y);

      friend Complex operator-(const double& __x, const Complex& __y);

      friend Complex operator*(const double& __x, const Complex& __y);

      friend Complex operator/(const double& __x, const Complex& __y);

      friend bool operator==(const double& __x, const Complex& __y);

      friend bool operator!=(const double& __x, const Complex& __y);

      friend Complex operator+(const Complex& __x);

      friend Complex operator-(const Complex& __x);

      friend double norm(const Complex& __z);

      friend double fabs(const Complex& __z);

      friend Complex conj(const Complex& __z);

      friend ostream & operator << ( ostream & output, const Complex& __z );
};
#endif
