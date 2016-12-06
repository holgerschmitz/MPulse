/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "fdtd_plain.hpp"

#include <schnek/grid.hpp>

void FDTD_Plain::init()
{
  retrieveData("Ex", pEx);
  retrieveData("Ey", pEy);
  retrieveData("Ez", pEz);

  retrieveData("Bx", pBx);
  retrieveData("By", pBy);
  retrieveData("Bz", pBz);
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
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;
  
  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;
  
  IndexType low = Ex.getInnerLo();
  IndexType high = Ex.getInnerHi();

  Vector dx = MPulse::getDx();

  for (int i=low[0]+1; i<=high[0]; ++i)
    for (int j=low[1]+1; j<=high[1]; ++j)
      for (int k=low[2]+1; k<=high[2]; ++k)
  {
    Ex(i,j,k) = Ex(i,j,k) + dt*clight2*( (Bz(i,j,k) - Bz(i,j-1,k))/dx[1] - (By(i,j,k) - By(i,j,k-1))/dx[2] );
    Ey(i,j,k) = Ey(i,j,k) + dt*clight2*( (Bx(i,j,k) - Bx(i,j,k-1))/dx[2] - (Bz(i,j,k) - Bz(i-1,j,k))/dx[0] );
    Ez(i,j,k) = Ez(i,j,k) + dt*clight2*( (By(i,j,k) - By(i-1,j,k))/dx[0] - (Bx(i,j,k) - Bx(i,j-1,k))/dx[1] );
  }

  schnek::DomainSubdivision<Field> &sub = MPulse::getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_Plain::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;
  
  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  IndexType low = storage->getLow();
  IndexType high = storage->getHigh();

  Vector dx = storage->getDx();

  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      for (int k=low[2]; k<high[2]; ++k)
  {
    Bx(i,j,k) = Bx(i,j,k) 
      + dt*(
          (Ey(i,j,k+1) - Ey(i,j,k))/dx[2]
        - (Ez(i,j+1,k) - Ez(i,j,k))/dx[1]
      );

    By(i,j,k) = By(i,j,k) 
      + dt*(
          (Ez(i+1,j,k) - Ez(i,j,k))/dx[0]
        - (Ex(i,j,k+1) - Ex(i,j,k))/dx[2]
      );

    Bz(i,j,k) = Bz(i,j,k) 
      + dt*( 
          (Ex(i,j+1,k) - Ex(i,j,k))/dx[1]
        - (Ey(i+1,j,k) - Ey(i,j,k))/dx[0]
      );
  }

  schnek::DomainSubdivision<Field> &sub = MPulse::getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
