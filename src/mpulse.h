#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include <schnek/matrix.h>
#include <schnek/fixedarray.h>

#ifdef NDEBUG
#define MPulseGridChecker schnek::MatrixNoArgCheck
#else
#define MPulseGridChecker schnek::MatrixAssertCheck
#endif

//typedef double REAL;
typedef float REAL;

typedef schnek::Matrix<REAL, 1, MPulseGridChecker> DataGrid1d;
typedef schnek::Matrix<REAL, 2, MPulseGridChecker> DataGrid2d;
typedef schnek::Matrix<REAL, 3, MPulseGridChecker> DataGrid3d;

typedef schnek::FixedArray<REAL, 3> Vector;

typedef DataGrid1d::IndexType GridIndex1d;
typedef DataGrid2d::IndexType GridIndex2d;
typedef DataGrid3d::IndexType GridIndex3d;

static const size_t Dimension = 3;
typedef DataGrid1d DataLine;
typedef DataGrid3d DataGrid;
typedef GridIndex3d GridIndex;

enum Direction {north, south, west, east, up, down};

#endif // MPULSE_MPULSE_H
