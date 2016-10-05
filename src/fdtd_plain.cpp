/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "fdtd_plain.hpp"

void FDTD_Plain::initStorage(Storage *storage_)
{  
  storage = storage_;
  
  pEx = storage->addGrid("Ex");
  pEy = storage->addGrid("Ey");
  pEz = storage->addGrid("Ez");

  pBx = storage->addGrid("Bx");
  pBy = storage->addGrid("By");
  pBz = storage->addGrid("Bz");
  
  storage->addToGroup("E", "Ex");
  storage->addToGroup("E", "Ey");
  storage->addToGroup("E", "Ez");

  storage->addToGroup("B", "Bx");
  storage->addToGroup("B", "By");
  storage->addToGroup("B", "Bz");
}

void FDTD_Plain::stepSchemeInit(double dt)
{
  stepB(0.5*dt);
}

void FDTD_Plain::stepScheme(double dt)
{
  stepD(dt);
  stepB(dt);
}

void FDTD_Plain::stepD(double dt)
{
  DataGrid &Ex = *pEx;
  DataGrid &Ey = *pEy;
  DataGrid &Ez = *pEz;
  
  DataGrid &Bx = *pBx;
  DataGrid &By = *pBy;
  DataGrid &Bz = *pBz;
  
  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();

  double dx = storage->getDx();
  double dy = storage->getDy();
  double dz = storage->getDz();

  for (int i=low[0]+1; i<=high[0]; ++i)
    for (int j=low[1]+1; j<=high[1]; ++j)
      for (int k=low[2]+1; k<=high[2]; ++k)
  {
    Ex(i,j,k) = Ex(i,j,k) 
      + dt*(
          (Bz(i,j,k) - Bz(i,j-1,k))/dy
        - (By(i,j,k) - By(i,j,k-1))/dz
      );
      
    Ey(i,j,k) = Ey(i,j,k) 
      + dt*(
          (Bx(i,j,k) - Bx(i,j,k-1))/dz
        - (Bz(i,j,k) - Bz(i-1,j,k))/dx
      );
 
    Ez(i,j,k) = Ez(i,j,k) 
      + dt*(
          (By(i,j,k) - By(i-1,j,k))/dx
        - (Bx(i,j,k) - Bx(i,j-1,k))/dy
      );
  }
  storage->applyBoundary("E");
}


void FDTD_Plain::stepB(double dt)
{
  DataGrid &Ex = *pEx;
  DataGrid &Ey = *pEy;
  DataGrid &Ez = *pEz;
  
  DataGrid &Bx = *pBx;
  DataGrid &By = *pBy;
  DataGrid &Bz = *pBz;

  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();

  double dx = storage->getDx();
  double dy = storage->getDy();
  double dz = storage->getDz();

  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      for (int k=low[2]; k<high[2]; ++k)
  {
    Bx(i,j,k) = Bx(i,j,k) 
      + dt*(
          (Ey(i,j,k+1) - Ey(i,j,k))/dz
        - (Ez(i,j+1,k) - Ez(i,j,k))/dy
      );

    By(i,j,k) = By(i,j,k) 
      + dt*(
          (Ez(i+1,j,k) - Ez(i,j,k))/dx
        - (Ex(i,j,k+1) - Ex(i,j,k))/dz
      );

    Bz(i,j,k) = Bz(i,j,k) 
      + dt*( 
          (Ex(i,j+1,k) - Ex(i,j,k))/dx
        - (Ey(i+1,j,k) - Ey(i,j,k))/dx
      );
  }
  storage->applyBoundary("B");
}
