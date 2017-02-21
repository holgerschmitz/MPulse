/*
 * fdtd.cpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#include "fdtd.hpp"
#include <schnek/tools/fieldtools.hpp>
#include <schnek/tools/literature.hpp>

void FieldSolver::initParameters(BlockParameters& parameters) {
  parameters.addParameter("eps_rel", &eps_rel);

  x_par = parameters.addArrayParameter("", x, BlockParameters::readonly);
  E_par = parameters.addArrayParameter("E", initE, 0.0);
  B_par = parameters.addArrayParameter("B", initB, 0.0);
}

void FieldSolver::registerData() {
  addData("Ex", Ex);
  addData("Ey", Ey);
  addData("Ez", Ez);

  addData("Bx", Bx);
  addData("By", By);
  addData("Bz", Bz);
}

void FieldSolver::init() {
  pBlockVariables blockVars = getVariables();
  pDependencyMap depMap(new DependencyMap(blockVars));
  DependencyUpdater updater(depMap);

  pParametersGroup spaceVars = pParametersGroup(new ParametersGroup());
  spaceVars->addArray(x_par);

  updater.addIndependentArray(x_par);

  fill_field(Ex, x, initE[0], updater, E_par[0]);
  fill_field(Ey, x, initE[1], updater, E_par[1]);
  fill_field(Ez, x, initE[2], updater, E_par[2]);

  fill_field(Bx, x, initB[0], updater, B_par[0]);
  fill_field(By, x, initB[1], updater, B_par[1]);
  fill_field(Bz, x, initB[2], updater, B_par[2]);

  LiteratureArticle Yee1966("Yee1966", "Yee, K",
      "Numerical solution of initial boundary value problems"
      "involving Maxwell's equations in isotropic media.",
      "IEEE Transactions on Antennas and Propagation", "1966", "AP-14", "302--307");

  LiteratureManager::instance().addReference(
      "Integration of electrodynamic fields uses the Finite Difference Time Domain method.",
      Yee1966);
}

void FieldSolver::stepSchemeInit(double dt) {
  stepB(0.5*dt);
  subdivision->exchange(Bx);
  subdivision->exchange(By);
  subdivision->exchange(Bz);
}

void FieldSolver::stepScheme(double dt) {
  stepD(dt);
  subdivision->exchange(Ex);
  subdivision->exchange(Ey);
  subdivision->exchange(Ez);
  stepB(dt);
  subdivision->exchange(Bx);
  subdivision->exchange(By);
  subdivision->exchange(Bz);
}

void FieldSolver::stepD(double dt) {
  IndexType low = Ex.getInnerLo();
  IndexType high = Ex.getInnerHi();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j) {
      Ex(i,j) += dt*clight2/eps_rel*(   (Bz(i,j) - Bz(i,j-1))/dx[1] );
      Ey(i,j) += dt*clight2/eps_rel*( - (Bz(i,j) - Bz(i-1,j))/dx[0] );
      Ez(i,j) += dt*clight2/eps_rel*(   (By(i,j) - By(i-1,j))/dx[0] - (Bx(i,j) - Bx(i,j-1))/dx[1] );
    }

}


void FieldSolver::stepB(double dt) {
  IndexType low = Bx.getInnerLo();
  IndexType high = Bx.getInnerHi();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j) {
      Bx(i,j) += dt*( - (Ez(i,j+1) - Ez(i,j))/dx[1] );
      By(i,j) += dt*(   (Ez(i+1,j) - Ez(i,j))/dx[0] );
      Bz(i,j) += dt*(   (Ex(i,j+1) - Ex(i,j))/dx[1] - (Ey(i+1,j) - Ey(i,j))/dx[0] );
    }

}
