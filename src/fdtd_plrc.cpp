
#include "fdtd_plrc.h"
#include "storage.h"
#include <cmath>
#include <complex>
#include <algorithm>
#include <limits>

//===============================================================
//==========  FDTD_PLRCCore
//===============================================================

void FDTD_PLRCCore::coreInitStorage(Storage *storage_)
{  
  storage = storage_;
  
  pEx = storage->addGrid("Ex");
  pEy = storage->addGrid("Ey");
  pEz = storage->addGrid("Ez");

  pBx = storage->addGrid("Bx");
  pBy = storage->addGrid("By");
  pBz = storage->addGrid("Bz");
  
  pSigma = storage->addGrid("Sigma");
  
  pKappaEdx = storage->addLine("KappaEdx", 0, 1.0);
  pKappaEdy = storage->addLine("KappaEdy", 1, 1.0);
  pKappaEdz = storage->addLine("KappaEdz", 2, 1.0);

  pKappaHdx = storage->addLine("KappaHdx", 0, 1.0);
  pKappaHdy = storage->addLine("KappaHdy", 1, 1.0);
  pKappaHdz = storage->addLine("KappaHdz", 2, 1.0);

  pPsiRx[0] = storage->addGrid("PsiR1x");
  pPsiRy[0] = storage->addGrid("PsiR1y");
  pPsiRz[0] = storage->addGrid("PsiR1z");
  
  pPsiRx[1] = storage->addGrid("PsiR2x");
  pPsiRy[1] = storage->addGrid("PsiR2y");
  pPsiRz[1] = storage->addGrid("PsiR2z");
  
  pPsiRx[2] = storage->addGrid("PsiR3x");
  pPsiRy[2] = storage->addGrid("PsiR3y");
  pPsiRz[2] = storage->addGrid("PsiR3z");
  
  pPsiIx[0] = storage->addGrid("PsiI1x");
  pPsiIy[0] = storage->addGrid("PsiI1y");
  pPsiIz[0] = storage->addGrid("PsiI1z");
  
  pPsiIx[1] = storage->addGrid("PsiI2x");
  pPsiIy[1] = storage->addGrid("PsiI2y");
  pPsiIz[1] = storage->addGrid("PsiI2z");
  
  pPsiIx[2] = storage->addGrid("PsiI3x");
  pPsiIy[2] = storage->addGrid("PsiI3y");
  pPsiIz[2] = storage->addGrid("PsiI3z"); 
  
  storage->addToGroup("E", "Ex");
  storage->addToGroup("E", "Ey");
  storage->addToGroup("E", "Ez");

  storage->addToGroup("B", "Bx");
  storage->addToGroup("B", "By");
  storage->addToGroup("B", "Bz");
  
  std::for_each(currents.begin(), currents.end(),  InitStorageFunctor<Current>(storage));
  currents.erase(
    std::remove_if(currents.begin(), currents.end(), Current::CurrentInvalidPredicate()),
    currents.end()
  );

  for (CurrentList::iterator it = currents.begin(); it !=currents.end(); ++it)
  {
    if ((*it)->getJx() == 0) {
	std::cerr << "Invalid Current\n";
    }
  }
  
  std::for_each(magCurrents.begin(), magCurrents.end(),  InitStorageFunctor<Current>(storage));
  magCurrents.erase(
    std::remove_if(magCurrents.begin(), magCurrents.end(), Current::CurrentInvalidPredicate()),
     magCurrents.end()
  );
  
}

//===============================================================
//==========  FDTD_PLRCLinCore
//===============================================================

#include <iostream>

void FDTD_PLRCLinCore::plrcStepD(double dt, 
                                 int i, int j, int k, 
                                 double dx, double dy, double dz,
                                 double Jx, double Jy, double Jz)
{
  
//  for (int l=pKappaEdx->getLow()[0]; l<=pKappaEdx->getHigh()[0]; ++l)
//  {
//    std::cout << l << " " << (*pKappaEdx)(l) << std::endl;
//  }
  
  double &ex = (*pEx)(i,j,k);
  double &ey = (*pEy)(i,j,k);
  double &ez = (*pEz)(i,j,k);
  
#ifndef NDEBUG
  {
      double test = ex*ey*ez;
      if ( !((test>0) || (test<1)) )
      {
        std::cerr << "NaN in Field (1) " << ex << " " << ey << " " << ez << "\n";
        exit(-1);
      }
  }
#endif

  double kappaEdx = (*pKappaEdx)(i)*dx;
  double kappaEdy = (*pKappaEdy)(j)*dy;
  double kappaEdz = (*pKappaEdz)(k)*dz;
  
//  std::cerr << i << " " << j << " " << k << " " << kappaEdx << " " << kappaEdy << " " << kappaEdz << " " << std::endl;
         
#ifndef NDEBUG
  {
      double test = kappaEdx*kappaEdy*kappaEdz;
      if ( !((test>0) || (test<1)) || (test==0) )
      {
        std::cerr << "Error in Kappa " << kappaEdx << " " << kappaEdy << " " << kappaEdz << "\n";
        exit(-1);
      }
  }
#endif


  double Psix = 0;
  double Psiy = 0;
  double Psiz = 0;

  for (int n=0;n<3;++n)
  {
    double &pxr = pPsiRx[n]->operator()(i,j,k);
    double &pyr = pPsiRy[n]->operator()(i,j,k);
    double &pzr = pPsiRz[n]->operator()(i,j,k);

    double &pxi = pPsiIx[n]->operator()(i,j,k);
    double &pyi = pPsiIy[n]->operator()(i,j,k);
    double &pzi = pPsiIz[n]->operator()(i,j,k);

    std::complex<double> D = plrcData.dchi0[n]-plrcData.dxi0[n];

    std::complex<double> px = D*ex + std::complex<double>(pxr,pxi);
    std::complex<double> py = D*ey + std::complex<double>(pyr,pyi);
    std::complex<double> pz = D*ez + std::complex<double>(pzr,pzi);

    Psix += std::real(px);
    Psiy += std::real(py);
    Psiz += std::real(pz);

    px = plrcData.dxi0[n]*ex + plrcData.Crec[n]*px;
    py = plrcData.dxi0[n]*ey + plrcData.Crec[n]*py;
    pz = plrcData.dxi0[n]*ez + plrcData.Crec[n]*pz;

    pxr = std::real(px);
    pyr = std::real(py);
    pzr = std::real(pz);

    pxi = std::imag(px);
    pyi = std::imag(py);
    pzi = std::imag(pz);

  }
  
#ifndef NDEBUG
  {
      double test = Psix*Psiy*Psiz;
      if ( !((test>0) || (test<1)) )
      {
	  std::cerr << "NaN in Psi " << Psix << " " << Psiy << " " << Psiz << "\n";
	  exit(-1);
      }
  }
#endif


  double sigma = 0.5*(*pSigma)(i,j,k)*dt;

  double denomX = eps + (plrcData.sumChi0 - plrcData.sumXi0) + sigma;
  double denomY = eps + (plrcData.sumChi0 - plrcData.sumXi0) + sigma;
  double denomZ = eps + (plrcData.sumChi0 - plrcData.sumXi0) + sigma;
  
  double numerX = eps - plrcData.sumXi0 - sigma;
  double numerY = eps - plrcData.sumXi0 - sigma;
  double numerZ = eps - plrcData.sumXi0 - sigma;
  
//    double E2 = ex*ex + ey*ey + ez*ez;

  // after this Dx, Dy and Dz actually contain D - P_L
  
//  if ((Jx!=0.0) || (Jy!=0.0) || (Jz!=0.0))
//    std::cerr << i << ' ' << j << ' ' << k << ' ' << Jx << ' ' << Jy << ' ' << Jz << '\n';
  
//  if (fabs(Jx+Jy+Jz) > 1e-40) std::cerr << "Current " << Jx << " " << Jy << " " << Jz << std::endl;

  double exn =
    ( 
      numerX*ex
      + (
          dt*(
            ((*pBz)(i,j,k) - (*pBz)(i,j-1,k))/kappaEdy
          - ((*pBy)(i,j,k) - (*pBy)(i,j,k-1))/kappaEdz
          + Jx
        )
        + Psix
      )
    ) / denomX;

  double eyn =
    ( 
      numerY*ey
      + (
          dt*(
            ((*pBx)(i,j,k) - (*pBx)(i,j,k-1))/kappaEdz
          - ((*pBz)(i,j,k) - (*pBz)(i-1,j,k))/kappaEdx
          + Jy
        )
        + Psiy  
      )
    ) / denomY;

  double ezn =
    ( 
      numerZ*ez
      + (
          dt*(
            ((*pBy)(i,j,k) - (*pBy)(i-1,j,k))/kappaEdx
          - ((*pBx)(i,j,k) - (*pBx)(i,j-1,k))/kappaEdy
          + Jz
        )
        + Psiz  
      )
    ) / denomZ;

  ex = exn;
  ey = eyn;
  ez = ezn;
  
#ifndef NDEBUG
 {
     double test = ex*ey*ez;
     if ( !((test>0) || (test<1)) )
     {
	 std::cerr << "NaN in Field (2) " << ex << " " << ey << " " << ez << "\n";
	 exit(-1);
     }
 }
#endif

}

void FDTD_PLRCLinCore::plrcStepB(double dt, 
                                 int i, int j, int k, 
                                 double dx, double dy, double dz,
                                 double Jx, double Jy, double Jz)
{
  
  double kappaHdx = (*pKappaHdx)(i)*dx;
  double kappaHdy = (*pKappaHdy)(j)*dy;
  double kappaHdz = (*pKappaHdz)(k)*dz;
  
  (*pBx)(i,j,k) = (*pBx)(i,j,k) 
    + dt*(
        ((*pEy)(i,j,k+1) - (*pEy)(i,j,k))/kappaHdz
      - ((*pEz)(i,j+1,k) - (*pEz)(i,j,k))/kappaHdy
     + Jx
    );

  (*pBy)(i,j,k) = (*pBy)(i,j,k) 
    + dt*(
        ((*pEz)(i+1,j,k) - (*pEz)(i,j,k))/kappaHdx
      - ((*pEx)(i,j,k+1) - (*pEx)(i,j,k))/kappaHdz
     + Jy    
    );

  (*pBz)(i,j,k) = (*pBz)(i,j,k) 
    + dt*( 
        ((*pEx)(i,j+1,k) - (*pEx)(i,j,k))/kappaHdy
      - ((*pEy)(i+1,j,k) - (*pEy)(i,j,k))/kappaHdx
     + Jz     
    );
}

ParameterMap* FDTD_PLRCLinCore::CustomParamMap (ParameterMap* pm)
{
  return pm;
}

//===============================================================
//==========  FDTD_PLRCNonlinCore
//===============================================================

void FDTD_PLRCNonlinCore::plrcStepD(double dt, 
                                    int i, int j, int k, 
                                    double dx, double dy, double dz,
                                    double Jx, double Jy, double Jz)
{
  double &ex = (*pEx)(i,j,k);
  double &ey = (*pEy)(i,j,k);
  double &ez = (*pEz)(i,j,k);

  double kappaEdx = (*pKappaEdx)(i)*dx;
  double kappaEdy = (*pKappaEdy)(j)*dy;
  double kappaEdz = (*pKappaEdz)(k)*dz;

  double Psix = 0;
  double Psiy = 0;
  double Psiz = 0;

  for (int n=0;n<3;++n)
  {
    double &pxr = pPsiRx[n]->operator()(i,j,k);
    double &pyr = pPsiRy[n]->operator()(i,j,k);
    double &pzr = pPsiRz[n]->operator()(i,j,k);

    double &pxi = pPsiIx[n]->operator()(i,j,k);
    double &pyi = pPsiIy[n]->operator()(i,j,k);
    double &pzi = pPsiIz[n]->operator()(i,j,k);

    std::complex<double> D = plrcData.dchi0[n]-plrcData.dxi0[n];

    std::complex<double> px = D*ex + std::complex<double>(pxr,pxi);
    std::complex<double> py = D*ey + std::complex<double>(pyr,pyi);
    std::complex<double> pz = D*ez + std::complex<double>(pzr,pzi);

    Psix += std::real(px);
    Psiy += std::real(py);
    Psiz += std::real(pz);

    px = plrcData.dxi0[n]*ex + plrcData.Crec[n]*px;
    py = plrcData.dxi0[n]*ey + plrcData.Crec[n]*py;
    pz = plrcData.dxi0[n]*ez + plrcData.Crec[n]*pz;

    pxr = std::real(px);
    pyr = std::real(py);
    pzr = std::real(pz);

    pxi = std::imag(px);
    pyi = std::imag(py);
    pzi = std::imag(pz);

  }


  double sigma = 0.5*(*pSigma)(i,j,k)*dt;

  double denom = eps + plrcData.sumChi0 - plrcData.sumXi0 + sigma;
  double numer = eps - plrcData.sumXi0 - sigma;

  double E2 = ex*ex + ey*ey + ez*ez;

//    double E2 = ex*ex + ey*ey + ez*ez;

  // after this Dx, Dy and Dz actually contain D - P_L


  double cx =
    ( 
      (numer + chi*E2)*ex

      + dt*(
          ((*pBz)(i,j,k) - (*pBz)(i,j-1,k))/kappaEdy
        - ((*pBy)(i,j,k) - (*pBy)(i,j,k-1))/kappaEdz
        + Jx
      )
      + Psix
    ) / denom;

  double cy =
    ( 
      (numer + chi*E2)*ey
      + dt*(
          ((*pBx)(i,j,k) - (*pBx)(i,j,k-1))/kappaEdz
        - ((*pBz)(i,j,k) - (*pBz)(i-1,j,k))/kappaEdx
        + Jy
      )
      + Psiy
    ) / denom;

  double cz =
    ( 
      (numer + chi*E2)*ez
      + dt*(
          ((*pBy)(i,j,k) - (*pBy)(i-1,j,k))/kappaEdx
        - ((*pBx)(i,j,k) - (*pBx)(i,j-1,k))/kappaEdy
         + Jz
      )
      + Psiz
    ) / denom;

  double E = sqrt(E2);
  double A = chi/denom;
  double C = sqrt(cx*cx + cy*cy + cz*cz);
      
  if (C>0)
  {
    // Newton iteration
    double Eold;
    do {
      Eold = E;
      E2 = E*E;
      E = (2*A*E2*E + C) / (3*A*E2 + 1);
    } while (fabs(E-Eold) > 1e-9);

    ex = E*cx/C;
    ey = E*cy/C;
    ez = E*cz/C;
  }
  else
  {
    ex = 0.0;
    ey = 0.0;
    ez = 0.0;
  }
   
}

void FDTD_PLRCNonlinCore::plrcStepB(double dt, 
                                    int i, int j, int k, 
                                    double dx, double dy, double dz,
                                    double Jx, double Jy, double Jz)
{
  
  double kappaHdx = (*pKappaHdx)(i)*dx;
  double kappaHdy = (*pKappaHdy)(j)*dy;
  double kappaHdz = (*pKappaHdz)(k)*dz;

  (*pBx)(i,j,k) = (*pBx)(i,j,k) 
    + dt*(
        ((*pEy)(i,j,k+1) - (*pEy)(i,j,k))/kappaHdz
      - ((*pEz)(i,j+1,k) - (*pEz)(i,j,k))/kappaHdy
      + Jx
    );

  (*pBy)(i,j,k) = (*pBy)(i,j,k) 
    + dt*(
        ((*pEz)(i+1,j,k) - (*pEz)(i,j,k))/kappaHdx
      - ((*pEx)(i,j,k+1) - (*pEx)(i,j,k))/kappaHdz
       + Jy
    );

  (*pBz)(i,j,k) = (*pBz)(i,j,k) 
    + dt*( 
        ((*pEx)(i,j+1,k) - (*pEx)(i,j,k))/kappaHdy
      - ((*pEy)(i+1,j,k) - (*pEy)(i,j,k))/kappaHdx
       + Jz
    );
}

ParameterMap* FDTD_PLRCNonlinCore::CustomParamMap (ParameterMap* pm)
{
  (*pm)["chi"] = WParameter(new ParameterValue<double>(&chi,0.1));
  return pm;
}

