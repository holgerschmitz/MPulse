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

void ConstantPlasmaDensity::initStorage(Storage *storage_)
{
  storage = storage_;
  pRho = storage->addGrid("PlasmaDensity");
  initialized = false;
}

void ConstantPlasmaDensity::stepScheme(double dt)
{
  if (initialized) return;
  
  DataGrid &Rho = *pRho;
  
  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();

  double dx = Globals::instance().gridDX();
  double dy = Globals::instance().gridDY();
  double dz = Globals::instance().gridDZ();
  double x,y,z;

  int nxh = Globals::instance().gridX()/2;
  int nyh = Globals::instance().gridY()/2;
  
  std::cout << "Intializing: (" << low[0] << ","<< low[1] <<","<< low[2] << "), (" << high[0] << ","<< high[1] <<","<< high[2] << ")\n";
  
  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    //Rho(i,j,k) = minDensity*exp(k*dz/length);
    x = (i-nxh)*dx;
    y = (j-nyh)*dy;
    double pert = pertAmp*sin(2*M_PI*(pertKx*x + pertPhaseX))*sin(2*M_PI*(pertKy*y + pertPhaseY));
    z = k*dz - pert;
    if (z<pos)
      Rho(i,j,k) = minDensity;
    else if (z>=(pos+length))
      Rho(i,j,k) = maxDensity;
    else
      Rho(i,j,k) = minDensity + (maxDensity-minDensity)*(z-pos)/length;
  }
  initialized = true;
}


ParameterMap* ConstantPlasmaDensity::MakeParamMap (ParameterMap* pm)
{
  pm = OptField::MakeParamMap(pm);

  (*pm)["minDensity"] = WParameter(new ParameterValue<double>(&minDensity,1.0));
  (*pm)["maxDensity"] = WParameter(new ParameterValue<double>(&maxDensity,1.0));
  (*pm)["pos"] = WParameter(new ParameterValue<double>(&pos,1.0));
  (*pm)["length"] = WParameter(new ParameterValue<double>(&length,1.0));

  (*pm)["pertKx"] = WParameter(new ParameterValue<double>(&pertKx,1.0));
  (*pm)["pertKy"] = WParameter(new ParameterValue<double>(&pertKy,1.0));
  (*pm)["pertAmp"] = WParameter(new ParameterValue<double>(&pertAmp,0.0));
  (*pm)["pertPhaseX"] = WParameter(new ParameterValue<double>(&pertPhaseX,0.0));
  (*pm)["pertPhaseY"] = WParameter(new ParameterValue<double>(&pertPhaseY,0.0));
  
  return pm;
}
