#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include <schnek/matrix.h>
#include <schnek/fixedarray.h>

typedef schnek::Matrix<double, 1> DataGrid1d;
typedef schnek::Matrix<double, 2> DataGrid2d;
typedef schnek::Matrix<double, 3> DataGrid3d;

//typedef schnek::Matrix<double, 1, schnek::MatrixAssertCheck> DataGrid1d;
//typedef schnek::Matrix<double, 2, schnek::MatrixAssertCheck> DataGrid2d;
//typedef schnek::Matrix<double, 3, schnek::MatrixAssertCheck> DataGrid3d;

typedef schnek::FixedArray<double, 3> Vector;

typedef DataGrid1d::IndexType GridIndex1d;
typedef DataGrid2d::IndexType GridIndex2d;
typedef DataGrid3d::IndexType GridIndex3d;

static const size_t Dimension = 3;
typedef DataGrid1d DataLine;
typedef DataGrid3d DataGrid;
typedef GridIndex3d GridIndex;

enum Direction {north, south, west, east, up, down};

#endif // MPULSE_MPULSE_H
