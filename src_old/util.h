#include <cmath>
#ifndef UTIL_H
#define UTIL_H

/** @file util.h
  * @brief helper functions
  *
  * Defines some helper function templates for mathematical operations
  */

#if ! defined (__GNUG__) && ! defined (__HP_aCC)
/** returns maximum of two T's */
template<class T>
inline T max(T a, T b)
{
  return a > b ? a:b;
}
/** returns maximum of two T's */
template<class T>
inline T min(T a, T b)
{
  return a < b ? a:b;
}
/** returns absolute value */
template<class T>
inline T abs(T a)
{
  return (a > 0) ? a : -a;
}
#endif

/**returns square root */
template<class T>
inline T sqr(T a)
{
  return a*a;
}

/** returns the norm of a*/
inline double norm(double a)
{
  return a*a;
}

/** returns norm of a */
inline float norm(float a)
{
  return a*a;
}

/** returns sign of a*/
template<class T>
inline T sign(T a, T b)
{
  return (b >= 0) ? abs(a) : -abs(a);
}

/** returns logarithm to base 2 */
template<class T>
inline T log2(T a)
{ 
  const double inv_log2 = 1/log(2.);
  return T(inv_log2*log(a));
}

template<class T>
    inline T ipow(T x, unsigned int e)
{
  T r(1);
  for (unsigned int i=0; i<e; ++i)
    r *= x;
  return r;
}

#endif
