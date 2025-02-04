/*
 * fieldionization.cpp
 *
 * Created on: 06/01/2025
 * Author: hrnair
 */

#include "fieldionization.hpp"

#include "../huerto/maths/functions/core.hpp"

void FieldIonization::initParameters(schnek::BlockParameters &blockPars)
{
  blockPars.addParameter("Z", &Z, 1.0);
  blockPars.addParameter("Uion", &Uion, 1.0);
  blockPars.addParameter("Wion", &Wion, 1.0);
  blockPars.addParameter("neff", &neff, 1.0);
  blockPars.addParameter("rat", &rat, 1.0);
}

void FieldIonization::init()
{
  SimulationEntity::init(this);
  retrieveData("Ex", Ex);
  retrieveData("Ey", Ey);
  retrieveData("Ez", Ez);

  retrieveData("Rho", Electrons);
  retrieveData("NeutralDensity", Neutrals);
}

double FieldIonization::computeIonizationRate(double E) 
{
  const double Eh = 5.13e11; // Atomic unit of electric field (V/m)
  const double e_pi = 6.6e16 / M_PI; // Pre-factor for ADK rate

  double neff = Z / sqrt(Uion / 13.6); // Effective quantum number
  double En = E / Eh; // Normalized electric field

  double factor1 = e_pi * (Z * Z) / pow(neff, 4.5);
  double factor2 = pow(10.87 * Eh / En * pow(Z, 3) / pow(neff, 4), 2 * neff - 1.5);
  double factor3 = exp(-2.0 / 3.0 * Eh / En * pow(Z, 3) / pow(neff, 3));

  return factor1 * factor2 * factor3; // Ionization rate (s^-1)
}

void FieldIonization::execute()
{
  double dt = getContext().getDt();
  Index low = Electrons.getLo();
  Index high = Electrons.getHi();

#ifdef HUERTO_ONE_DIM
  for (int i=low[0]; i<=high[0]; ++i)
  {
    double ex = Ex(i);
    double ey = Ey(i);
    double ez = Ez(i);
    double &rho = Electrons(i);
    double &neut = Neutrals(i);
#endif

#ifdef HUERTO_TWO_DIM
  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
  {
    double ex = Ex(i,j);
    double ey = Ey(i,j);
    double ez = Ez(i,j);
    double &rho = Electrons(i,j);
    double &neut = Neutrals(i,j);
#endif

#ifdef HUERTO_THREE_DIM
  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double ex = Ex(i,j,k);
    double ey = Ey(i,j,k);
    double ez = Ez(i,j,k);
    double &rho = Electrons(i,j,k);
    double &neut = Neutrals(i,j,k);
#endif
    double E = sqrt(ex * ex + ey * ey + ez * ez); // Electric field magnitude

    double Ri = computeIonizationRate(E); // Ionization rate
    double R = Ri * 0.5 * dt;

    rho = (rho * (1 - R) + Ri * rat * dt) / (1 + R);
  }
}
