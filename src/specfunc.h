#ifndef MPULSE_SPECFUNC_H
#define MPULSE_SPECFUNC_H

#include <complex>


#ifndef M_2PI
#define M_2PI       6.2831853071795864769252867665590   // 2Pi
#endif


/*-----------------------------------------------------------------------------*\
| Matpack special functions - Faddeeva(z)                             cwofz.cpp |
\*-----------------------------------------------------------------------------*/

// complex<double> Faddeeva (const complex<double>& z)
//
// Given a complex number z = (x,y), this subroutine computes 
// the value of the Faddeeva function w(z) = exp(-z^2)*erfc(-i*z), 
// where erfc is the complex complementary error function and i means sqrt(-1). 

std::complex<double> Faddeeva_2 (const std::complex<double>& z);

#endif
