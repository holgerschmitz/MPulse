
#include "fdtd_plrc.hpp"

#include "../huerto/constants.hpp"

#include <schnek/tools/literature.hpp>

#include <cmath>
#include <complex>
#include <algorithm>
#include <limits>
#include <memory>

//===============================================================
//==========  FDTD_PLRCCore
//===============================================================


void FDTD_PLRCCore::registerData()
{
  schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();

  Domain domainSize = subdivision.getInnerExtent(getContext().getSize());
  schnek::Array<bool, DIMENSION> stagger;
  stagger = false;

  Index low  = subdivision.getLo();
  Index high = subdivision.getHi();

  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  pSigma = std::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);

  pKappaEdx = std::make_shared<Grid1d>(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  (*pKappaEdx) = 1.0;
#ifndef HUERTO_ONE_DIM
  pKappaEdy = std::make_shared<Grid1d>(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  (*pKappaEdy) = 1.0;
#endif
#ifdef HUERTO_THREE_DIM
  pKappaEdz = std::make_shared<Grid1d>(schnek::Array<int, 1>(low[2]), schnek::Array<int, 1>(high[2]));
  (*pKappaEdz) = 1.0;
#endif

  pKappaHdx = std::make_shared<Grid1d>(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  (*pKappaHdx) = 1.0;
#ifndef HUERTO_ONE_DIM
  pKappaHdy = std::make_shared<Grid1d>(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  (*pKappaHdy) = 1.0;
#endif
#ifdef HUERTO_THREE_DIM
  pKappaHdz = std::make_shared<Grid1d>(schnek::Array<int, 1>(low[2]), schnek::Array<int, 1>(high[2]));
  (*pKappaHdz) = 1.0;
#endif

  for (int d=0; d<3; d++)
  {
    pPsiRx[d] = std::make_unique<Field>(lowIn, highIn, domainSize, stagger, 2);
    pPsiRy[d] = std::make_unique<Field>(lowIn, highIn, domainSize, stagger, 2);
    pPsiRz[d] = std::make_unique<Field>(lowIn, highIn, domainSize, stagger, 2);

    pPsiIx[d] = std::make_unique<Field>(lowIn, highIn, domainSize, stagger, 2);
    pPsiIy[d] = std::make_unique<Field>(lowIn, highIn, domainSize, stagger, 2);
    pPsiIz[d] = std::make_unique<Field>(lowIn, highIn, domainSize, stagger, 2);
  }

  addData("KappaEdx", pKappaEdx);
#ifndef HUERTO_ONE_DIM
  addData("KappaEdy", pKappaEdy);
#endif
#ifdef HUERTO_THREE_DIM
  addData("KappaEdz", pKappaEdz);
#endif

  addData("KappaHdx", pKappaHdx);
#ifndef HUERTO_ONE_DIM
  addData("KappaHdy", pKappaHdy);
#endif
#ifdef HUERTO_THREE_DIM
  addData("KappaHdz", pKappaHdz);
#endif
}


void FDTD_PLRCCore::init()
{
  retrieveData("Ex", pEx);
  retrieveData("Ey", pEy);
  retrieveData("Ez", pEz);

  retrieveData("Bx", pBx);
  retrieveData("By", pBy);
  retrieveData("Bz", pBz);


  for(pCurrentBlock current: schnek::BlockContainer<CurrentBlock>::childBlocks())
  {
    current->initCurrents(*this);
  }

  CurrentContainer::init(getContext());

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
                                 Index pos,
                                 Vector dx,
                                 double Jx, double Jy, double Jz) {
  double &ex = (*pEx)[pos];
  double &ey = (*pEy)[pos];
  double &ez = (*pEz)[pos];

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

  double kappaEdx = (*pKappaEdx)(pos[0])*dx[0];
#ifndef HUERTO_ONE_DIM
  double kappaEdy = (*pKappaEdy)(pos[1])*dx[1];
#endif
#ifdef HUERTO_THREE_DIM
  double kappaEdz = (*pKappaEdz)(pos[2])*dx[2];
#endif

  double Psix = 0;
  double Psiy = 0;
  double Psiz = 0;

  for (int n=0;n<3;++n)
  {
    double &pxr = (*pPsiRx[n])[pos];
    double &pyr = (*pPsiRy[n])[pos];
    double &pzr = (*pPsiRz[n])[pos];

    double &pxi = (*pPsiIx[n])[pos];
    double &pyi = (*pPsiIy[n])[pos];
    double &pzi = (*pPsiIz[n])[pos];

    // The following lines contain Kelley & Luebbers Eq (28)
    // split into two parts
    std::complex<double> D = plrcData.dchi0[n]-plrcData.dxi0[n];

    std::complex<double> px = D*ex + std::complex<double>(pxr, pxi);
    std::complex<double> py = D*ey + std::complex<double>(pyr, pyi);
    std::complex<double> pz = D*ez + std::complex<double>(pzr, pzi);

    // This intermediate result is needed to update the fields
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

  double sigma = 0.5*(*pSigma)[pos]*dt;

  double denomX = eps + (plrcData.sumChi0 - plrcData.sumXi0) + sigma;
  double denomY = eps + (plrcData.sumChi0 - plrcData.sumXi0) + sigma;
  double denomZ = eps + (plrcData.sumChi0 - plrcData.sumXi0) + sigma;

  double numerX = eps - plrcData.sumXi0 - sigma;
  double numerY = eps - plrcData.sumXi0 - sigma;
  double numerZ = eps - plrcData.sumXi0 - sigma;

#ifdef HUERTO_ONE_DIM
  double exn =
    (numerX*ex + dt*Jx + Psix) / denomX;

  double eyn =
    (
      numerY*ey
      +
        dt*(
          clight2*(
            - ((*pBz)[pos] - (*pBz)(pos[0]-1))/kappaEdx
          )
          - Jy/eps_0
        )
        + Psiy

    ) / denomY;

  double ezn =
    (
      numerZ*ez
      + (
          dt*(
            clight2*(
              ((*pBy)[pos] - (*pBy)(pos[0]-1))/kappaEdx
            )
          - Jz/eps_0
        )
        + Psiz
      )
    ) / denomZ;
#endif

#ifdef HUERTO_TWO_DIM
  double exn =
    (
      numerX*ex
      + (
          dt*(
            clight2*(
              ((*pBz)[pos] - (*pBz)(pos[0], pos[1]-1))/kappaEdy
            )
          - Jx/eps_0
        )
        + Psix
      )
    ) / denomX;

  double eyn =
    (
      numerY*ey
      + (
          dt*(
            clight2*(
              - ((*pBz)[pos] - (*pBz)(pos[0]-1, pos[1]))/kappaEdx
            )
          - Jy/eps_0
        )
        + Psiy
      )
    ) / denomY;

  double ezn =
    (
      numerZ*ez
      + (
          dt*(
              clight2*(
                  ((*pBy)[pos] - (*pBy)(pos[0]-1, pos[1]  ))/kappaEdx
                - ((*pBx)[pos] - (*pBx)(pos[0]  , pos[1]-1))/kappaEdy
              )
          - Jz/eps_0
        )
        + Psiz
      )
    ) / denomZ;
#endif

#ifdef HUERTO_THREE_DIM
  double exn =
    (
      numerX*ex
      + (
          dt*(
            clight2*(
                ((*pBz)[pos] - (*pBz)(pos[0], pos[1]-1, pos[2]  ))/kappaEdy
              - ((*pBy)[pos] - (*pBy)(pos[0], pos[1]  , pos[2]-1))/kappaEdz
            )
          - Jx/eps_0
        )
        + Psix
      )
    ) / denomX;

  double eyn =
    (
      numerY*ey
      + (
          dt*(
            clight2*(
                ((*pBx)[pos] - (*pBx)(pos[0],   pos[1], pos[2]-1))/kappaEdz
              - ((*pBz)[pos] - (*pBz)(pos[0]-1, pos[1], pos[2]  ))/kappaEdx
            )
          - Jy/eps_0
        )
        + Psiy
      )
    ) / denomY;

  double ezn =
    (
      numerZ*ez
      + (
          dt*(
            clight2*(
                ((*pBy)[pos] - (*pBy)(pos[0]-1, pos[1],   pos[2]))/kappaEdx
              - ((*pBx)[pos] - (*pBx)(pos[0]  , pos[1]-1, pos[2]))/kappaEdy
            )
          - Jz/eps_0
        )
        + Psiz
      )
    ) / denomZ;
#endif
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
                                 Index pos,
                                 Vector dx,
                                 double Jx, double Jy, double Jz)
{
  double kappaHdx = (*pKappaHdx)(pos[0])*dx[0];
#ifndef HUERTO_ONE_DIM
  double kappaHdy = (*pKappaHdy)(pos[1])*dx[1];
#endif
#ifdef HUERTO_THREE_DIM
  double kappaHdz = (*pKappaHdz)(pos[2])*dx[2];
#endif

#ifdef HUERTO_ONE_DIM
  (*pBx)[pos] = (*pBx)[pos]
    + dt*(
     + Jx
    );

  (*pBy)[pos] = (*pBy)[pos]
    + dt*(
        ((*pEz)(pos[0]+1) - (*pEz)[pos])/kappaHdx
     + Jy
    );

  (*pBz)[pos] = (*pBz)[pos]
    + dt*(
      - ((*pEy)(pos[0]+1) - (*pEy)[pos])/kappaHdx
     + Jz
    );
#endif

#ifdef HUERTO_TWO_DIM
  (*pBx)[pos] = (*pBx)[pos]
    + dt*(
      - ((*pEz)(pos[0], pos[1]+1) - (*pEz)[pos])/kappaHdy
     + Jx
    );

  (*pBy)[pos] = (*pBy)[pos]
    + dt*(
        ((*pEz)(pos[0]+1,pos[1]) - (*pEz)[pos])/kappaHdx
     + Jy
    );

  (*pBz)[pos] = (*pBz)[pos]
    + dt*(
        ((*pEx)(pos[0],pos[1]+1) - (*pEx)[pos])/kappaHdy
      - ((*pEy)(pos[0]+1,pos[1]) - (*pEy)[pos])/kappaHdx
     + Jz
    );
#endif

#ifdef HUERTO_THREE_DIM
  (*pBx)[pos] = (*pBx)[pos]
    + dt*(
        ((*pEy)(pos[0], pos[1],  pos[2]+1) - (*pEy)[pos])/kappaHdz
      - ((*pEz)(pos[0], pos[1]+1,pos[2]) - (*pEz)[pos])/kappaHdy
     + Jx
    );

  (*pBy)[pos] = (*pBy)[pos]
    + dt*(
        ((*pEz)(pos[0]+1,pos[1],pos[2]) - (*pEz)[pos])/kappaHdx
      - ((*pEx)(pos[0],pos[1],pos[2]+1) - (*pEx)[pos])/kappaHdz
     + Jy
    );

  (*pBz)[pos] = (*pBz)[pos]
    + dt*(
        ((*pEx)(pos[0],pos[1]+1,pos[2]) - (*pEx)[pos])/kappaHdy
      - ((*pEy)(pos[0]+1,pos[1],pos[2]) - (*pEy)[pos])/kappaHdx
     + Jz
    );
#endif
}


//===============================================================
//==========  FDTD_PLRCNonlinCore
//===============================================================

void FDTD_PLRCNonlinCore::plrcStepD(double dt,
                                    Index pos,
                                    Vector dx,
                                    double Jx, double Jy, double Jz) {
  double &ex = (*pEx)[pos];
  double &ey = (*pEy)[pos];
  double &ez = (*pEz)[pos];

  double kappaEdx = (*pKappaEdx)(pos[0])*dx[0];
#ifndef HUERTO_ONE_DIM
  double kappaEdy = (*pKappaEdy)(pos[1])*dx[1];
#endif
#ifdef HUERTO_THREE_DIM
  double kappaEdz = (*pKappaEdz)(pos[2])*dx[2];
#endif

  double Psix = 0;
  double Psiy = 0;
  double Psiz = 0;

  for (int n=0;n<3;++n)
  {
    double &pxr = (*pPsiRx[n])[pos];
    double &pyr = (*pPsiRy[n])[pos];
    double &pzr = (*pPsiRz[n])[pos];

    double &pxi = (*pPsiIx[n])[pos];
    double &pyi = (*pPsiIy[n])[pos];
    double &pzi = (*pPsiIz[n])[pos];

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


  double sigma = 0.5*(*pSigma)[pos]*dt;

  double denom = eps + plrcData.sumChi0 - plrcData.sumXi0 + sigma;
  double numer = eps - plrcData.sumXi0 - sigma;

  double E2 = ex*ex + ey*ey + ez*ez;

  // after this Dx, Dy and Dz actually contain D - P_L


#ifdef HUERTO_ONE_DIM
  double cx =
    (
      (numer + chi*E2)*ex

      + dt*(
        - Jx/eps_0
      )
      + Psix
    ) / denom;

  double cy =
    (
      (numer + chi*E2)*ey
      + dt*(
        clight2*(
          - ((*pBz)[pos] - (*pBz)(pos[0]-1))/kappaEdx
        )
        - Jy/eps_0
      )
      + Psiy
    ) / denom;

  double cz =
    (
      (numer + chi*E2)*ez
      + dt*(
        clight2*(
          ((*pBy)[pos] - (*pBy)(pos[0]-1))/kappaEdx
        )
        - Jz/eps_0
      )
      + Psiz
    ) / denom;
#endif

#ifdef HUERTO_TWO_DIM
  double cx =
    (
      (numer + chi*E2)*ex

      + dt*(
        clight2*(
          ((*pBz)[pos] - (*pBz)(pos[0],pos[1]-1))/kappaEdy
        )
        - Jx/eps_0
      )
      + Psix
    ) / denom;

  double cy =
    (
      (numer + chi*E2)*ey
      + dt*(
        clight2*(
          - ((*pBz)[pos] - (*pBz)(pos[0]-1,pos[1]))/kappaEdx
        )
        - Jy/eps_0
      )
      + Psiy
    ) / denom;

  double cz =
    (
      (numer + chi*E2)*ez
      + dt*(
        clight2*(
            ((*pBy)[pos] - (*pBy)(pos[0]-1,pos[1]))/kappaEdx
          - ((*pBx)[pos] - (*pBx)(pos[0],pos[1]-1))/kappaEdy
        )
        - Jz/eps_0
      )
      + Psiz
    ) / denom;
#endif

#ifdef HUERTO_THREE_DIM
  double cx =
    (
      (numer + chi*E2)*ex

      + dt*(
        clight2*(
            ((*pBz)[pos] - (*pBz)(pos[0],pos[1]-1,pos[2]))/kappaEdy
          - ((*pBy)[pos] - (*pBy)(pos[0],pos[1],pos[2]-1))/kappaEdz
        )
        - Jx/eps_0
      )
      + Psix
    ) / denom;

  double cy =
    (
      (numer + chi*E2)*ey
      + dt*(
        clight2*(
            ((*pBx)[pos] - (*pBx)(pos[0],pos[1],pos[2]-1))/kappaEdz
          - ((*pBz)[pos] - (*pBz)(pos[0]-1,pos[1],pos[2]))/kappaEdx
        )
        - Jy/eps_0
      )
      + Psiy
    ) / denom;

  double cz =
    (
      (numer + chi*E2)*ez
      + dt*(
        clight2*(
            ((*pBy)[pos] - (*pBy)(pos[0]-1,pos[1],pos[2]))/kappaEdx
          - ((*pBx)[pos] - (*pBx)(pos[0],pos[1]-1,pos[2]))/kappaEdy
        )
        - Jz /eps_0
      )
      + Psiz
    ) / denom;
#endif

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
                                    Index pos,
                                    Vector dx,
                                    double Jx, double Jy, double Jz)
{
  double kappaHdx = (*pKappaHdx)(pos[0])*dx[0];
#ifndef HUERTO_ONE_DIM
  double kappaHdy = (*pKappaHdy)(pos[1])*dx[1];
#endif
#ifdef HUERTO_THREE_DIM
  double kappaHdz = (*pKappaHdz)(pos[2])*dx[2];
#endif

#ifdef HUERTO_ONE_DIM
  (*pBx)[pos] = (*pBx)[pos]
    + dt*(
      + Jx
    );

  (*pBy)[pos] = (*pBy)[pos]
    + dt*(
        ((*pEz)(pos[0]+1) - (*pEz)[pos])/kappaHdx
       + Jy
    );

  (*pBz)[pos] = (*pBz)[pos]
    + dt*(
      - ((*pEy)(pos[0]+1) - (*pEy)[pos])/kappaHdx
       + Jz
    );
#endif

#ifdef HUERTO_TWO_DIM
  (*pBx)[pos] = (*pBx)[pos]
    + dt*(
      - ((*pEz)(pos[0],pos[1]+1) - (*pEz)[pos])/kappaHdy
      + Jx
    );

  (*pBy)[pos] = (*pBy)[pos]
    + dt*(
        ((*pEz)(pos[0]+1,pos[1]) - (*pEz)[pos])/kappaHdx
       + Jy
    );

  (*pBz)[pos] = (*pBz)[pos]
    + dt*(
        ((*pEx)(pos[0],pos[1]+1) - (*pEx)[pos])/kappaHdy
      - ((*pEy)(pos[0]+1,pos[1]) - (*pEy)[pos])/kappaHdx
       + Jz
    );
#endif

#ifdef HUERTO_THREE_DIM
  (*pBx)[pos] = (*pBx)[pos]
    + dt*(
        ((*pEy)(pos[0],pos[1],pos[2]+1) - (*pEy)[pos])/kappaHdz
      - ((*pEz)(pos[0],pos[1]+1,pos[2]) - (*pEz)[pos])/kappaHdy
      + Jx
    );

  (*pBy)[pos] = (*pBy)[pos]
    + dt*(
        ((*pEz)(pos[0]+1,pos[1],pos[2]) - (*pEz)[pos])/kappaHdx
      - ((*pEx)(pos[0],pos[1],pos[2]+1) - (*pEx)[pos])/kappaHdz
       + Jy
    );

  (*pBz)[pos] = (*pBz)[pos]
    + dt*(
        ((*pEx)(pos[0],pos[1]+1,pos[2]) - (*pEx)[pos])/kappaHdy
      - ((*pEy)(pos[0]+1,pos[1],pos[2]) - (*pEy)[pos])/kappaHdx
       + Jz
    );
#endif
}

void FDTD_PLRCNonlinCore::initParameters(schnek::BlockParameters &blockPars)
{
  FDTD_PLRCCore::initParameters(blockPars);
  blockPars.addParameter("chi", &chi,0.1);
}

void FDTD_PLRCNonlinCore::init()
{
  FDTD_PLRCCore::init();

  schnek::LiteratureArticle Schmitz2012("Schmitz2012", "Holger Schmitz and Vladimir Mezentsev",
      "Full-vectorial modeling of femtosecond pulses for laser inscription of photonic structures",
      "Journal of the Optical Society of America B", "2012", "29", "1208-1217");

  schnek::LiteratureManager::instance().addReference(
      "Nonlinear Piecewise Linear Recursive Convolution for dispersive media",
      Schmitz2012);
}

