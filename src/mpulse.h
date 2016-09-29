#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include <schnek/grid.hpp>
#include <schnek/array.hpp>

#ifdef NDEBUG
#define MPulseGridChecker schnek::GridNoArgCheck
#else
#define MPulseGridChecker schnek::GridAssertCheck
#endif

typedef schnek::Grid<double, 1, MPulseGridChecker> DataGrid1d;
typedef schnek::Grid<double, 2, MPulseGridChecker> DataGrid2d;
typedef schnek::Grid<double, 3, MPulseGridChecker> DataGrid3d;

typedef schnek::Array<double, 3> Vector;

typedef DataGrid1d::IndexType GridIndex1d;
typedef DataGrid2d::IndexType GridIndex2d;
typedef DataGrid3d::IndexType GridIndex3d;

static const size_t Dimension = 3;
typedef DataGrid1d DataLine;
typedef DataGrid3d DataGrid;
typedef GridIndex3d GridIndex;

enum Direction {north, south, west, east, up, down};

#endif // MPULSE_MPULSE_H
