#include "plasmacurrent.h"
#include "storage.h"
#include "fieldsolver.h"


ParameterMap* PlasmaCurrentFactory::MakeParamMap (ParameterMap* pm)
{
  pm = CurrentFactory::MakeParamMap(pm);

  (*pm)["em"] = WParameter(new ParameterValue<double>(&em,1.0));
  (*pm)["gamma"] = WParameter(new ParameterValue<double>(&gamma,0.1));
  (*pm)["rho_bg"] = WParameter(new ParameterValue<double>(&rho_bg,0.0));
  
  return pm;
}


void PlasmaCurrentFactory::initCurrents(Storage *storage_, FieldSolver *solver)
{
  solver->addCurrent(new PlasmaCurrent(em, gamma, rho_bg));
}

PlasmaCurrent::PlasmaCurrent(double em_, double gamma_, double rho_bg_) 
  : em(em_), gamma(gamma_), rho_bg(rho_bg_)
{}

void PlasmaCurrent::initStorage(Storage *storage_)
{
  storage = storage_;
  
  pEx = storage->addGrid("Ex");
  pEy = storage->addGrid("Ey");
  pEz = storage->addGrid("Ez");
  
  pJx = storage->addGrid("PlasmaJx");
  pJy = storage->addGrid("PlasmaJy");
  pJz = storage->addGrid("PlasmaJz");
  
  pRho = storage->addGrid("PlasmaDensity");
}
    
void PlasmaCurrent::stepScheme(double dt)
{
  DataGrid &Ex = *pEx;
  DataGrid &Ey = *pEy;
  DataGrid &Ez = *pEz;
  
  DataGrid &Jx = *pJx;
  DataGrid &Jy = *pJy;
  DataGrid &Jz = *pJz;
  DataGrid &Rho = *pRho;

  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();
  
  const double gdtn = 1-0.5*gamma*dt;
  const double gdtd = 1+0.5*gamma*dt;
  const double emdt = em*dt;
  
  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      for (int k=low[2]; k<high[2]; ++k)
  {
    double &jx = Jx(i,j,k);
    double &jy = Jy(i,j,k);
    double &jz = Jz(i,j,k);
    double rho = Rho(i,j,k) + rho_bg;
    
    jx = (jx*gdtn - emdt*Ex(i,j,k)*rho)/gdtd;
    jy = (jy*gdtn - emdt*Ey(i,j,k)*rho)/gdtd;
    jz = (jz*gdtn - emdt*Ez(i,j,k)*rho)/gdtd;
  }

}

