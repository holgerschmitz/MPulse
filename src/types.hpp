/*
 * types.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_TYPES_H
#define MPULSE_TYPES_H

#include <schnek/grid.hpp>
#include <cstddef>

#ifdef NDEBUG
#define MPulseGridChecker schnek::GridNoArgCheck
#else
#define MPulseGridChecker schnek::GridAssertCheck
#endif

static const std::size_t DIMENSION = 3;

typedef schnek::Array<int, DIMENSION> Index;
typedef schnek::Array<double, DIMENSION> Vector;

typedef schnek::Grid<double, DIMENSION> Grid;
typedef boost::shared_ptr<Grid> pGrid;

typedef schnek::Field<double, DIMENSION> Field;
typedef boost::shared_ptr<Field> pField;

typedef schnek::Grid<double, 1> DataLine;
typedef boost::shared_ptr<DataLine> pDataLine;

typedef schnek::Range<int, DIMENSION> Range;
typedef schnek::Array<bool, DIMENSION> Stagger;

static const double clight = 299792458;
static const double clight2 = clight*clight;

static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true );

static const Stagger bxStaggerYee(false, true,  true );
static const Stagger byStaggerYee(true,  false, true );
static const Stagger bzStaggerYee(true,  true,  false);

enum Direction {north, south, west, east, up, down};

#endif
