/*
 * types.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_TYPES_H
#define MPULSE_TYPES_H

#include "../huerto/types.hpp"

inline bool doDiag(int i, int j, int k) {
  return false && (i==395) && (j==25) && (k==25); //(j>=13) && (j<=16);
}

#endif // MPULSE_TYPES_H
