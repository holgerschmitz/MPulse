
#include "fdtd_disp.h"
#include "storage.h"
#include "plasmacurrent.h"
#include <cmath>

void FDTD_Dispersion::initStorage(Storage *storage_)
{  
  storage = storage_;
  
  pEx = storage->addGrid("Ex");
  pEy = storage->addGrid("Ey");
  pEz = storage->addGrid("Ez");

  pBx = storage->addGrid("Bx");
  pBy = storage->addGrid("By");
  pBz = storage->addGrid("Bz");
  
  pPx[0] = storage->addGrid("P1x");
  pPy[0] = storage->addGrid("P1y");
  pPz[0] = storage->addGrid("P1z");
  
  pPx[1] = storage->addGrid("P2x");
  pPy[1] = storage->addGrid("P2y");
  pPz[1] = storage->addGrid("P2z");
  
  pPx[2] = storage->addGrid("P3x");
  pPy[2] = storage->addGrid("P3y");
  pPz[2] = storage->addGrid("P3z");
  
  pPxp[0] = storage->addGrid("P1xp");
  pPyp[0] = storage->addGrid("P1yp");
  pPzp[0] = storage->addGrid("P1zp");
  
  pPxp[1] = storage->addGrid("P2xp");
  pPyp[1] = storage->addGrid("P2yp");
  pPzp[1] = storage->addGrid("P2zp");
  
  pPxp[2] = storage->addGrid("P3xp");
  pPyp[2] = storage->addGrid("P3yp");
  pPzp[2] = storage->addGrid("P3zp"); 
  
  storage->addToGroup("E", "Ex");
  storage->addToGroup("E", "Ey");
  storage->addToGroup("E", "Ez");

  storage->addToGroup("B", "Bx");
  storage->addToGroup("B", "By");
  storage->addToGroup("B", "Bz");
  
  for 
  (
    CurrentList::iterator it = currents.begin(); 
    it != currents.end(); 
    ++it
  )
  {
    (*it)->initStorage(storage);
  }

}

void FDTD_Dispersion::stepSchemeInit(double dt)
{
  stepB(0.5*dt);
  
  for 
  (
    CurrentList::iterator it = currents.begin(); 
    it != currents.end(); 
    ++it
  )
  {
    (*it)->stepSchemeInit(dt);
  }
}

void FDTD_Dispersion::stepScheme(double dt)
{
  stepD(dt);
  stepB(dt);
  
  for 
  (
    CurrentList::iterator it = currents.begin(); 
    it != currents.end(); 
    ++it
  )
  {
    (*it)->stepScheme(dt);
  }
}

void FDTD_Dispersion::stepD(double dt)
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

  double dt2 = dt*dt;

  for (int i=low[0]+1; i<=high[0]; ++i)
    for (int j=low[1]+1; j<=high[1]; ++j)
      for (int k=low[2]+1; k<=high[2]; ++k)
  {
    double ex = Ex(i,j,k);
    double ey = Ey(i,j,k);
    double ez = Ez(i,j,k);
    
    double Px = 0;
    double Py = 0;
    double Pz = 0;
    
    double Pxp = 0;
    double Pyp = 0;
    double Pzp = 0;
    
    for (int n=0;n<3;++n)
    {
      REAL &px = pPx[n]->operator()(i,j,k);
      REAL &py = pPy[n]->operator()(i,j,k);
      REAL &pz = pPz[n]->operator()(i,j,k);
      
      REAL &pxp = pPxp[n]->operator()(i,j,k);
      REAL &pyp = pPyp[n]->operator()(i,j,k);
      REAL &pzp = pPzp[n]->operator()(i,j,k);
      
      double a = dt2*LOm2[n];
      double b = a*LEps2[n];
      double pxn = (2-a)*px - pxp + b*ex;
      double pyn = (2-a)*py - pyp + b*ey;
      double pzn = (2-a)*pz - pzp + b*ez;
      
      pxp = px;
      pyp = py;
      pzp = pz;
      
      px = pxn;
      py = pyn;
      pz = pzn;
      
      Px += px;
      Py += py;
      Pz += pz;

      Pxp += pxp;
      Pyp += pyp;
      Pzp += pzp;
    }

//    double E2 = ex*ex + ey*ey + ez*ez;
    
    // after this Dx, Dy and Dz actually contain D - P_L
    
    double Dx = /* (eps + chi*E2)* */ ex + Pxp - Px
      + dt*(
          (Bz(i,j,k) - Bz(i,j-1,k))/dy
        - (By(i,j,k) - By(i,j,k-1))/dz
      );
      
    double Dy = /* (eps + chi*E2)* */ ey + Pyp - Py
      + dt*(
          (Bx(i,j,k) - Bx(i,j,k-1))/dz
        - (Bz(i,j,k) - Bz(i-1,j,k))/dx
      );
 
    double Dz = /* (eps + chi*E2)* */ ez + Pzp - Pz
      + dt*(
          (By(i,j,k) - By(i-1,j,k))/dx
        - (Bx(i,j,k) - Bx(i,j-1,k))/dy
      );
    
    
//    double E = sqrt(E2);
//    double D = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
//    
//    // Scalar Newton iteration because (D-P) parallel to E
//    double Eold;
//    do {
//      Eold = E;
//      E2 = E*E;
//      E = (2*chi*E2*E + D) / (3*chi*E2 + eps);
//    } while (fabs(E-Eold) > 1e-9);
//    
//    Ex(i,j,k) = E*Dx/D;
//    Ey(i,j,k) = E*Dy/D;
//    Ez(i,j,k) = E*Dz/D;

    Ex(i,j,k) = Dx;
    Ey(i,j,k) = Dy;
    Ez(i,j,k) = Dz;
   
  }
  storage->applyBoundary("E");
}


void FDTD_Dispersion::stepB(double dt)
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


ParameterMap* FDTD_Dispersion::MakeParamMap (ParameterMap* pm)
{
  pm = FieldSolver::MakeParamMap(pm);

  (*pm)["eps"] = WParameter(new ParameterValue<double>(&eps,1.0));
  (*pm)["chi"] = WParameter(new ParameterValue<double>(&chi,0.1));

  (*pm)["L1eps2"] = WParameter(new ParameterValue<double>(&LEps2[0],0.0));
  (*pm)["L2eps2"] = WParameter(new ParameterValue<double>(&LEps2[1],0.0));
  (*pm)["L3eps2"] = WParameter(new ParameterValue<double>(&LEps2[2],0.0));
  
  (*pm)["L1Om2"] = WParameter(new ParameterValue<double>(&LOm2[0],0.0));
  (*pm)["L2Om2"] = WParameter(new ParameterValue<double>(&LOm2[1],0.0));
  (*pm)["L3Om2"] = WParameter(new ParameterValue<double>(&LOm2[2],0.0));
  
//  (*pm)["plasma"] = WParameter(
//      new ParameterRebuild<PlasmaCurrent, Current>(&currents)
//  );
  
  return pm;
}
