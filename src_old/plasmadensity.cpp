#include "plasmadensity.h"
#include "storage.h"
#include "util.h"

void PlasmaDensity::initStorage(Storage *storage_)
{
  storage = storage_;
  
  pEx = storage->addGrid("Ex");
  pEy = storage->addGrid("Ey");
  pEz = storage->addGrid("Ez");
  
  pSigma = storage->addGrid("Sigma");
    
  pRho = storage->addGrid("PlasmaDensity");
}
    
void PlasmaDensity::stepScheme(double dt)
{
  DataGrid &Ex = *pEx;
  DataGrid &Ey = *pEy;
  DataGrid &Ez = *pEz;
  
  DataGrid &Rho = *pRho;
  DataGrid &Sigma = *pSigma;

  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();
  
  const double n = 0.5*dt*nu;
  const double m = dt*mpa;
  const double w = mpa*Wion;
  
  const int Km = K - 1;
  
  double sigmaMax = 0;
  
  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double ex = Ex(i,j,k);
    double ey = Ey(i,j,k);
    double ez = Ez(i,j,k);
    double &rho = Rho(i,j,k);
    double &sigma = Sigma(i,j,k);
    
    double I   = ex*ex + ey*ey + ez*ez;
    double IKm = ipow(I,Km);
    double IK  = I*IKm;
    
    rho = ( rho*(1+n*I) + m*IK ) / (1 - n*I);
    sigma += w*IKm;
    if (sigma > sigmaMax) sigmaMax = sigma;
  }
  
  if (sigmaMax*dt>0.2)
  {
    std::cerr << "Plasma absorption warning! sigma = " << sigmaMax << std::endl; 
  }
}


ParameterMap* PlasmaDensity::MakeParamMap (ParameterMap* pm)
{
  pm = OptField::MakeParamMap(pm);

  (*pm)["nu"] = WParameter(new ParameterValue<double>(&nu,1.0));
  (*pm)["mpa"] = WParameter(new ParameterValue<double>(&mpa,1.0));
  (*pm)["Wion"] = WParameter(new ParameterValue<double>(&Wion,1.0));
  (*pm)["K"] = WParameter(new ParameterValue<int>(&K,5));
  
  return pm;
}
