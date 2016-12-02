
#include "fdtd_plrc.hpp"

#include <schnek/tools/literature.hpp>

#include <boost/make_shared.hpp>

#include <cmath>
#include <complex>
#include <algorithm>
#include <limits>

//===============================================================
//==========  FDTD_PLRCCore
//===============================================================


void FDTD_PLRCCore::registerData()
{
  schnek::DomainSubdivision<Field> &subdivision = MPulse::getSubdivision();

  schnek::Range<double, DIMENSION> domainSize(schnek::Array<double, DIMENSION>(0,0,0), MPulse::getSize());
  schnek::Array<bool, DIMENSION> stagger;
  stagger = false;


  schnek::Range<double, 1> domainSizeX(schnek::Array<double, 1>(0), schnek::Array<double, 1>(MPulse::getSize()[0]));
  schnek::Range<double, 1> domainSizeY(schnek::Array<double, 1>(0), schnek::Array<double, 1>(MPulse::getSize()[1]));
  schnek::Range<double, 1> domainSizeZ(schnek::Array<double, 1>(0), schnek::Array<double, 1>(MPulse::getSize()[2]));

  schnek::Array<bool, 1> stagger1d;
  stagger1d = false;

  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  pSigma = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  
  pKappaEdx = boost::make_shared<DataLine>(schnek::Array<int, 1>(lowIn[0]), schnek::Array<int, 1>(highIn[0]), domainSizeX, stagger1d, 2);
  pKappaEdy = boost::make_shared<DataLine>(schnek::Array<int, 1>(lowIn[1]), schnek::Array<int, 1>(highIn[1]), domainSizeY, stagger1d, 2);
  pKappaEdz = boost::make_shared<DataLine>(schnek::Array<int, 1>(lowIn[2]), schnek::Array<int, 1>(highIn[2]), domainSizeZ, stagger1d, 2);

  pKappaHdx = boost::make_shared<DataLine>(schnek::Array<int, 1>(lowIn[0]), schnek::Array<int, 1>(highIn[0]), domainSizeX, stagger1d, 2);
  pKappaHdy = boost::make_shared<DataLine>(schnek::Array<int, 1>(lowIn[1]), schnek::Array<int, 1>(highIn[1]), domainSizeY, stagger1d, 2);
  pKappaHdz = boost::make_shared<DataLine>(schnek::Array<int, 1>(lowIn[2]), schnek::Array<int, 1>(highIn[2]), domainSizeZ, stagger1d, 2);

  pPsiRx[0] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiRy[0] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiRz[0] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  
  pPsiRx[1] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiRy[1] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiRz[1] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  
  pPsiRx[2] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiRy[2] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiRz[2] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);

  pPsiIx[0] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiIy[0] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiIz[0] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  
  pPsiIx[1] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiIy[1] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiIz[1] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);

  pPsiIx[2] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiIy[2] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  pPsiIz[2] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);

  addData("KappaEdx", pKappaEdx);
  addData("KappaEdy", pKappaEdy);
  addData("KappaEdz", pKappaEdz);

  addData("KappaHdx", pKappaHdx);
  addData("KappaHdy", pKappaHdy);
  addData("KappaHdz", pKappaHdz);
}


void FDTD_PLRCCore::init()
{

  retrieveData("Ex", pEx);
  retrieveData("Ey", pEy);
  retrieveData("Ez", pEz);

  retrieveData("Bx", pBx);
  retrieveData("By", pBy);
  retrieveData("Bz", pBz);

  
  BOOST_FOREACH(pCurrentBlock current, schnek::BlockContainer<CurrentBlock>::childBlocks())
  {
    current->initCurrents(*this);
  }

  schnek::LiteratureArticle Kelley1996("Kelley1996", "D. F. Kelley and R. J. Luebbers",
      "Piecewise linear recursive convolution for dispersive media using FDTD",
      "IEEE Transactions on Antennas and Propagation", "1996", "44", "792--797");

  schnek::LiteratureManager::instance().addReference(
      "Piecewise Linear Recursive Convolution for dispersive media",
      Kelley1996);

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

void FDTD_PLRCNonlinCore::initParameters(schnek::BlockParameters &blockPars)
{
  FDTD_PLRCCore::initParameters(blockPars);
  blockPars.addParameter("chi", &chi,0.1);
}

