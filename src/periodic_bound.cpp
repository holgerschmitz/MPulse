#include "periodic_bound.h"

#include "mpulse.h"
#include "globals.h"

/* **************************************************************
 *                  SinglePeriodicBoundary                      *
 ****************************************************************/

void SinglePeriodicBoundary::exchangeX(DataGrid3d &field)
{
  const GridIndex3d &UBound = field.getHigh();
  const GridIndex3d &LBound = field.getLow();
  
  int yi;
  int zi;

  int mx0=UBound[0], mx1=mx0-1;
  int lx0=LBound[0], lx1=lx0+1;
  
  for (yi = LBound[1]; yi <= UBound[1]; ++yi)
  {
    for (zi = LBound[2]; zi <= UBound[2]; ++zi)
    {
      field(lx0, yi, zi) = field(mx1, yi, zi);
      field(mx0, yi, zi) = field(lx1, yi, zi);
    }
  }
}

void SinglePeriodicBoundary::exchangeY(DataGrid3d &field) 
{
  const GridIndex3d &UBound = field.getHigh();
  const GridIndex3d &LBound = field.getLow();
  
  int xi;
  int zi;

  int my0=UBound[1], my1=my0-1;
  int ly0=LBound[1], ly1=ly0+1;
  
  for (xi = LBound[0]; xi <= UBound[0]; ++xi)
  {
    for (zi = LBound[2]; zi <= UBound[2]; ++zi)
    {
      field(xi, ly0, zi) = field(xi, my1, zi);
      field(xi, my0, zi) = field(xi, ly1, zi);
    }
  }
}

void SinglePeriodicBoundary::exchangeZ(DataGrid3d &field) 
{
  const GridIndex3d &UBound = field.getHigh();
  const GridIndex3d &LBound = field.getLow();
  
  int xi;
  int yi;

  int mz0=UBound[2], mz1=mz0-1;
  int lz0=LBound[2], lz1=lz0+1;
  
  for (xi = LBound[0]; xi <= UBound[0]; ++xi)
  {
    for (yi = LBound[1]; yi <= UBound[1]; ++yi)
    {
      field(xi, yi, lz0) = field(xi, yi, mz1);
      field(xi, yi, mz0) = field(xi, yi, lz1);
    }
  }
}


const GridIndex &SinglePeriodicBoundary::RegionLow() const {
  return Globals::instance().gridLow();
}

const GridIndex &SinglePeriodicBoundary::RegionHigh() const {
  return Globals::instance().gridHigh();
}


/* **************************************************************
 *                  SingleXYPeriodicBoundary                      *
 ****************************************************************/

void SingleXYPeriodicBoundary::exchangeX(DataGrid3d &field)
{
  const GridIndex3d &UBound = field.getHigh();
  const GridIndex3d &LBound = field.getLow();
  
  int yi;
  int zi;

  int mx0=UBound[0], mx1=mx0-1;
  int lx0=LBound[0], lx1=lx0+1;
  
  for (yi = LBound[1]; yi <= UBound[1]; ++yi)
  {
    for (zi = LBound[2]; zi <= UBound[2]; ++zi)
    {
      field(lx0, yi, zi) = field(mx1, yi, zi);
      field(mx0, yi, zi) = field(lx1, yi, zi);
    }
  }
}

void SingleXYPeriodicBoundary::exchangeY(DataGrid3d &field) 
{
  const GridIndex3d &UBound = field.getHigh();
  const GridIndex3d &LBound = field.getLow();
  
  int xi;
  int zi;

  int my0=UBound[1], my1=my0-1;
  int ly0=LBound[1], ly1=ly0+1;
  
  for (xi = LBound[0]; xi <= UBound[0]; ++xi)
  {
    for (zi = LBound[2]; zi <= UBound[2]; ++zi)
    {
      field(xi, ly0, zi) = field(xi, my1, zi);
      field(xi, my0, zi) = field(xi, ly1, zi);
    }
  }
}

void SingleXYPeriodicBoundary::exchangeZ(DataGrid3d &field) 
{}


const GridIndex &SingleXYPeriodicBoundary::RegionLow() const {
  return Globals::instance().gridLow();
}

const GridIndex &SingleXYPeriodicBoundary::RegionHigh() const {
  return Globals::instance().gridHigh();
}
