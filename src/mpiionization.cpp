/*
 * mpiionization.cpp
 *
 * Created on: 07/01/2025
 * Author: hrnair
 */

#include "mpiionization.hpp"

#include "../huerto/maths/functions/core.hpp"

void MPIIonization::initParameters(schnek::BlockParameters &blockPars)
{
  blockPars.addParameter("rat", &rat, 1.0);
  blockPars.addParameter("mpa", &mpa, 1.0);
  blockPars.addParameter("Wion", &Wion, 1.0);
  blockPars.addParameter("K", &K, 5);
}

void MPIIonization::init()
{
  SimulationEntity::init(this);
  retrieveData("Ex", Ex);
  retrieveData("Ey", Ey);
  retrieveData("Ez", Ez);

  retrieveData("Rho", Electrons);
  retrieveData("NeutralDensity", Neutrals);
}

void MPIIonization::execute()
{
  double dt = getContext().getDt();
  Index low = Electrons.getLo();
  Index high = Electrons.getHi();
  
  const double m = dt*mpa*0.5;
  const double w = mpa*Wion;
  
  const unsigned int Km =  std::max(K - 1, 0);
  
  double sigmaMax = 0;
  
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
    // double &sigma = Sigma(i,j,k); // global constant
    
    double I   = ex*ex + ey*ey + ez*ez;
    double IKm = pow(I, Km);
    double IK  = I*IKm;
    
    rho = ( rho*(1-m*IK) + mpa*IK*rat) / (1 + m*IK);
    if (w*IKm > sigmaMax) sigmaMax = w*IKm;
  }
      
  if (sigmaMax*dt>0.2)
  {
    std::cerr << "Plasma absorption warning! sigma = " << sigmaMax << std::endl; 
  }
}
