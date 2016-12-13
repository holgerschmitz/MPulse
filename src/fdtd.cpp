/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "fdtd.hpp"

#include <schnek/grid.hpp>

void FDTDSolver::registerData()
{
//  addData("Ex", Ex);
//  addData("Ey", Ey);
//  addData("Ez", Ez);
//
//  addData("Bx", Bx);
//  addData("By", By);
//  addData("Bz", Bz);
}


void FDTDSolver::init()
{
  schnek::ChildBlock<FDTDSolver>::init();

  retrieveData("Ex", pEx);
  retrieveData("Ey", pEy);
  retrieveData("Ez", pEz);

  retrieveData("Bx", pBx);
  retrieveData("By", pBy);
  retrieveData("Bz", pBz);
}
void FDTDSolver::stepScheme(double dt)
{
  schnek::DomainSubdivision<Field> &subdivision = Simulation::getSubdivision();
  stepE(dt);
  subdivision.exchange(*pEx);
  subdivision.exchange(*pEy);
  subdivision.exchange(*pEz);

  stepB(dt);
  subdivision.exchange(*pBx);
  subdivision.exchange(*pBy);
  subdivision.exchange(*pBz);
}

void FDTDSolver::stepE(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  IndexType low = Ex.getInnerLo();
  IndexType high = Ex.getInnerHi();

  Vector dx = Simulation::getDx();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    if ((i==25) && (j==25) && (k==25)) std::cout << "E: " << Ey(i,j,k) << " (" << Bz(i,j,k) << " " << Bz(i,j-1,k);
    Ex(i,j,k) = Ex(i,j,k) + dt*clight2*( (Bz(i,j,k) - Bz(i,j-1,k))/dx[1] - (By(i,j,k) - By(i,j,k-1))/dx[2] );
    Ey(i,j,k) = Ey(i,j,k) + dt*clight2*( (Bx(i,j,k) - Bx(i,j,k-1))/dx[2] - (Bz(i,j,k) - Bz(i-1,j,k))/dx[0] );
    Ez(i,j,k) = Ez(i,j,k) + dt*clight2*( (By(i,j,k) - By(i-1,j,k))/dx[0] - (Bx(i,j,k) - Bx(i,j-1,k))/dx[1] );
    if ((i==25) && (j==25) && (k==25)) std::cout << ")  -->  " << Ey(i,j,k) << std::endl;
  }
}

void FDTDSolver::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  IndexType low = Bx.getInnerLo();
  IndexType high = Bx.getInnerHi();

  Vector dx = Simulation::getDx();

  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      for (int k=low[2]; k<high[2]; ++k)
  {
    if ((i==25) && (j==25) && (k==25)) std::cout << "B: " << Bz(i,j,k);
    Bx(i,j,k) = Bx(i,j,k) + dt*((Ey(i,j,k+1) - Ey(i,j,k))/dx[2] - (Ez(i,j+1,k) - Ez(i,j,k))/dx[1]);
    By(i,j,k) = By(i,j,k) + dt*((Ez(i+1,j,k) - Ez(i,j,k))/dx[0] - (Ex(i,j,k+1) - Ex(i,j,k))/dx[2]);
    Bz(i,j,k) = Bz(i,j,k) + dt*((Ex(i,j+1,k) - Ex(i,j,k))/dx[1] - (Ey(i+1,j,k) - Ey(i,j,k))/dx[0]);
    if ((i==25) && (j==25) && (k==25)) std::cout << "  -->  " << Bz(i,j,k) << std::endl;
  }
}
