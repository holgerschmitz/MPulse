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





ParameterMap* MetalCurrentFactory::MakeParamMap (ParameterMap* pm)
{
  pm = CurrentFactory::MakeParamMap(pm);

  (*pm)["pos"] = WParameter(new ParameterValue<int>(&pos,100));
  (*pm)["amp"] = WParameter(new ParameterValue<double>(&amp,1.0));
  
  return pm;
}


void MetalCurrentFactory::initCurrents(Storage *storage_, FieldSolver *solver)
{
  solver->addCurrent(new MetalCurrent(pos, amp));
}

MetalCurrent::MetalCurrent(int pos_, double amp_) 
  : pos(pos_), amp(amp_)
{}

void MetalCurrent::initStorage(Storage *storage_)
{
  storage = storage_;
  
  pEx = storage->addGrid("Ex");
  pEy = storage->addGrid("Ey");
  pEz = storage->addGrid("Ez");
  
  pBx = storage->addGrid("Bx");
  pBy = storage->addGrid("By");
  pBz = storage->addGrid("Bz");
  
  pJx = storage->addGrid("MetalJx");
  pJy = storage->addGrid("MetalJy");
  pJz = storage->addGrid("MetalJz");
  
  pRho = storage->addGrid("PlasmaDensity");
  
  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();
  
  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    (*pRho)(i,j,k) = (k<pos)?-1.0:1.0;
  }
  
}
    
void MetalCurrent::stepScheme(double dt)
{
  DataGrid &Bx = *pBx;
  DataGrid &By = *pBy;
  DataGrid &Bz = *pBz;
  
  DataGrid &Ex = *pEx;
  DataGrid &Ey = *pEy;
  DataGrid &Ez = *pEz;
  
  DataGrid &Jx = *pJx;
  DataGrid &Jy = *pJy;
  DataGrid &Jz = *pJz;
  DataGrid &Rho = *pRho;

  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();
  
  if ((pos<low[2]) || (pos>high[2])) return;
  int lbound = std::max(pos, low[2]);
  
//   for (int i=low[0]; i<high[0]; ++i)
//     for (int j=low[1]; j<high[1]; ++j)
//       for (int k=lbound; k<high[2]; ++k)
//   {
//     Jx(i,j,k) = -amp*By(i,j,k);
//     Jy(i,j,k) = amp*Bx(i,j,k);
//   }

  for (int i=low[0]; i<high[0]; ++i)
     for (int j=low[1]; j<high[1]; ++j)
     {
       Ex(i,j,pos) = 0;
       Ey(i,j,pos) = 0;
       Ex(i,j,pos+1) = 0;
       Ey(i,j,pos+1) = 0;
       //Bx(i,j,pos-1) = Bx(i,j,pos);
     }
}

