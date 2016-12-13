/*
 * fdtd.cpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#include "fdtd.hpp"

void FieldSolver::registerData() {
  addData("Ex", Ex);
  addData("Ey", Ey);
  addData("Ez", Ez);

  addData("Bx", Bx);
  addData("By", By);
  addData("Bz", Bz);
}

void FieldSolver::stepSchemeInit(double dt) {
  schnek::DomainSubdivision<Field> &subdivision = Simulation::getSubdivision();
  stepB(0.5*dt);
  subdivision.exchange(Bx);
  subdivision.exchange(By);
  subdivision.exchange(Bz);
}

void FieldSolver::stepScheme(double dt) {
  schnek::DomainSubdivision<Field> &subdivision = Simulation::getSubdivision();
  stepD(dt);
  subdivision.exchange(Ex);
  subdivision.exchange(Ey);
  subdivision.exchange(Ez);
  stepB(dt);
  subdivision.exchange(Bx);
  subdivision.exchange(By);
  subdivision.exchange(Bz);
}

void FieldSolver::stepD(double dt) {
  Vector dx = Simulation::getDx();
  IndexType low = Ex.getInnerLo();
  IndexType high = Ex.getInnerHi();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k) {
        Ex(i,j,k) += dt*clight2*( (Bz(i,j,k) - Bz(i,j-1,k))/dx[1] - (By(i,j,k) - By(i,j,k-1))/dx[2] );
        Ey(i,j,k) += dt*clight2*( (Bx(i,j,k) - Bx(i,j,k-1))/dx[2] - (Bz(i,j,k) - Bz(i-1,j,k))/dx[0] );
        Ez(i,j,k) += dt*clight2*( (By(i,j,k) - By(i-1,j,k))/dx[0] - (Bx(i,j,k) - Bx(i,j-1,k))/dx[1] );
      }

}


void FieldSolver::stepB(double dt) {
  Vector dx = Simulation::getDx();
  IndexType low = Bx.getInnerLo();
  IndexType high = Bx.getInnerHi();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k) {
        Bx(i,j,k) += dt*( (Ey(i,j,k+1) - Ey(i,j,k))/dx[2] - (Ez(i,j+1,k) - Ez(i,j,k))/dx[1] );
        By(i,j,k) += dt*( (Ez(i+1,j,k) - Ez(i,j,k))/dx[0] - (Ex(i,j,k+1) - Ex(i,j,k))/dx[2] );
        Bz(i,j,k) += dt*( (Ex(i,j+1,k) - Ex(i,j,k))/dx[1] - (Ey(i+1,j,k) - Ey(i,j,k))/dx[0] );
      }

}
