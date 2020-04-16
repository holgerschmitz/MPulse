/*
 * types.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef VELLAMO_TYPES_H
#define VELLAMO_TYPES_H

#include <schnek/grid.hpp>
#include <cstddef>

#ifdef NDEBUG
#define MPulseGridChecker schnek::GridNoArgCheck
#else
#define MPulseGridChecker schnek::GridAssertCheck
#endif

static const size_t DIMENSION = 3;

typedef schnek::Array<int, DIMENSION> Index;
typedef schnek::Array<double, DIMENSION> Vector;
typedef schnek::Grid<double, DIMENSION, MPulseGridChecker> Grid;
typedef boost::shared_ptr<Grid> pGrid;
typedef schnek::Field<double, DIMENSION, MPulseGridChecker> Field;
typedef boost::shared_ptr<Field> pField;
typedef schnek::Range<int, DIMENSION> Range;
typedef schnek::Array<bool, DIMENSION> Stagger;


typedef schnek::Grid<double, 1, MPulseGridChecker> DataLine;
typedef boost::shared_ptr<DataLine> pDataLine;

static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true );
static const Stagger bxStaggerYee(false, true,  true );
static const Stagger byStaggerYee(true,  false, true );
static const Stagger bzStaggerYee(true,  true,  false);


enum Direction {north, south, west, east, up, down};
inline bool doDiag(int i, int j, int k) {
  return false && (i==395) && (j==25) && (k==25); //(j>=13) && (j<=16);
}

#endif // MPULSE_TYPES_H
