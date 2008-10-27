
#include "fdtd_nl.h"
#include "storage.h"
#include <cmath>

void FDTD_Nonlinear::initStorage(Storage *storage_)
{  
  storage = storage_;
  
  pEx = storage->addGrid("Ex");
  pEy = storage->addGrid("Ey");
  pEz = storage->addGrid("Ez");

  pBx = storage->addGrid("Bx");
  pBy = storage->addGrid("By");
  pBz = storage->addGrid("Bz");
  
  storage->addToGroup("E", "Ex");
  storage->addToGroup("E", "Ey");
  storage->addToGroup("E", "Ez");

  storage->addToGroup("B", "Bx");
  storage->addToGroup("B", "By");
  storage->addToGroup("B", "Bz");
}

void FDTD_Nonlinear::stepSchemeInit(double dt)
{
  stepB(0.5*dt);
}

void FDTD_Nonlinear::stepScheme(double dt)
{
  stepD(dt);
  stepB(dt);
}

void FDTD_Nonlinear::stepD(double dt)
{
  DataGrid &Ex = *pEx;
  DataGrid &Ey = *pEy;
  DataGrid &Ez = *pEz;
  
  DataGrid &Bx = *pBx;
  DataGrid &By = *pBy;
  DataGrid &Bz = *pBz;

  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();

  double dx = storage->getDx();
  double dy = storage->getDy();
  double dz = storage->getDz();

  for (int i=low[0]+1; i<=high[0]; ++i)
    for (int j=low[1]+1; j<=high[1]; ++j)
      for (int k=low[2]+1; k<=high[2]; ++k)
  {
    double ex = Ex(i,j,k);
    double ey = Ey(i,j,k);
    double ez = Ez(i,j,k);

    double E2 = ex*ex + ey*ey + ez*ez;
    
    double Dx = (eps + chi*E2)*ex
      + dt*(
          (Bz(i,j,k) - Bz(i,j-1,k))/dy
        - (By(i,j,k) - By(i,j,k-1))/dz
      );
      
    double Dy = (eps + chi*E2)*ey 
      + dt*(
          (Bx(i,j,k) - Bx(i,j,k-1))/dz
        - (Bz(i,j,k) - Bz(i-1,j,k))/dx
      );
 
    double Dz = (eps + chi*E2)*ez 
      + dt*(
          (By(i,j,k) - By(i-1,j,k))/dx
        - (Bx(i,j,k) - Bx(i,j-1,k))/dy
      );
      
    double E = sqrt(E2);
    double D = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
    
    // Newton iteration
    double Eold;
    do {
      Eold = E;
      E2 = E*E;
      E = (2*chi*E2*E + D) / (3*chi*E2 + eps);
    } while (fabs(E-Eold) > 1e-9);
    
    Ex(i,j,k) = E*Dx/D;
    Ey(i,j,k) = E*Dy/D;
    Ez(i,j,k) = E*Dz/D;
    
  }
  storage->applyBoundary("E");
}


void FDTD_Nonlinear::stepB(double dt)
{
  DataGrid &Ex = *pEx;
  DataGrid &Ey = *pEy;
  DataGrid &Ez = *pEz;
  
  DataGrid &Bx = *pBx;
  DataGrid &By = *pBy;
  DataGrid &Bz = *pBz;

  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();

  double dx = storage->getDx();
  double dy = storage->getDy();
  double dz = storage->getDz();

  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      for (int k=low[2]; k<high[2]; ++k)
  {
    Bx(i,j,k) = Bx(i,j,k) 
      + dt*(
          (Ey(i,j,k+1) - Ey(i,j,k))/dz
        - (Ez(i,j+1,k) - Ez(i,j,k))/dy
      );

    By(i,j,k) = By(i,j,k) 
      + dt*(
          (Ez(i+1,j,k) - Ez(i,j,k))/dx
        - (Ex(i,j,k+1) - Ex(i,j,k))/dz
      );

    Bz(i,j,k) = Bz(i,j,k) 
      + dt*( 
          (Ex(i,j+1,k) - Ex(i,j,k))/dx
        - (Ey(i+1,j,k) - Ey(i,j,k))/dx
      );
  }
  storage->applyBoundary("B");
}


ParameterMap* FDTD_Nonlinear::MakeParamMap (ParameterMap* pm)
{
  pm = FieldSolver::MakeParamMap(pm);

  (*pm)["eps"] = WParameter(new ParameterValue<double>(&eps,1.0));
  (*pm)["chi"] = WParameter(new ParameterValue<double>(&chi,0.1));
  return pm;
}
