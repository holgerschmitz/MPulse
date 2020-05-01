/*
 * types.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_TYPES_H
#define MPULSE_TYPES_H

#include "../huerto/types.hpp"

static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true );
static const Stagger bxStaggerYee(false, true,  true );
static const Stagger byStaggerYee(true,  false, true );
static const Stagger bzStaggerYee(true,  true,  false);

inline bool doDiag(int i, int j, int k) {
  return false && (i==395) && (j==25) && (k==25); //(j>=13) && (j<=16);
}

#endif // MPULSE_TYPES_H
