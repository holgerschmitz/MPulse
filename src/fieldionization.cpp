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
  std::cout << "FieldIonization::init called" << std::endl;

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
  std::cout << "FieldIonization::execute called" << std::endl;

  double dt = getContext().getDt();
  Index low = Electrons.getLo();
  Index high = Electrons.getHi();
  
  #ifdef HUERTO_ONE_DIM
  Range1d::LimitType low1D(low[0]);
  Range1d::LimitType high1D(high[0]);
  Range1d range(low1D, high1D);
  
  std::cout<<"FieldIonization: 1D Range"<<std::endl;

  for (auto &pos : range)
  {
    double ex = Ex[pos];
    double ey = Ey[pos];
    double ez = Ez[pos];
    double &rho = Electrons[pos];
    double &neut = Neutrals[pos];
#endif

#ifdef HUERTO_TWO_DIM
  Range2d::LimitType low2D(low[0], low[1]);
  Range2d::LimitType high2D(high[0], high[1]);
  Range2d range(low2D, high2D);  
  
  std::cout<<"2D Range"<<std::endl;

  for (auto &pos : range)
  {
    double ex = Ex[pos];
    double ey = Ey[pos];
    double ez = Ez[pos];
    double &rho = Electrons[pos];
    double &neut = Neutrals[pos];
#endif

#ifdef HUERTO_THREE_DIM
  Range3d::LimitType low3D(low[0], low[1], low[2]);
  Range3d::LimitType high3D(high[0], high[1], high[2]);
  Range3d range(low3D, high3D);  
  
  std::cout<<"3D Range"<<std::endl;

  for (auto &pos : range)
  {
    double ex = Ex[pos];
    double ey = Ey[pos];
    double ez = Ez[pos];
    double &rho = Electrons[pos];
    double &neut = Neutrals[pos];
#endif
    double E = sqrt(ex * ex + ey * ey + ez * ez); // Electric field magnitude

    double Ri = computeIonizationRate(E); // Ionization rate
    double R = Ri * 0.5 * dt;

    rho = (rho * (1 - R) + Ri * rat * dt) / (1 + R);
  }
}
