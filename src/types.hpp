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

static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true );

static const Stagger bxStaggerYee(false, true,  true );
static const Stagger byStaggerYee(true,  false, true );
static const Stagger bzStaggerYee(true,  true,  false);

enum Direction {north, south, west, east, up, down};

static const double PI     = 3.141592653589793238462643383279502884L;
static const double TWO_PI = 6.283185307179586476925286766559005768L;
static const double clight = 299792458;
static const double clight2 = clight*clight;
static const double mass_e = 9.10938291e-31;
static const double mass_p = 1.672621777e-27;
static const double unit_charge = 1.602176565e-19;
static const double mu_0 = 4e-7*PI;
static const double eps_0 = 1/(mu_0*clight2);
static const double eps_0_inv = (mu_0*clight2);

#endif
