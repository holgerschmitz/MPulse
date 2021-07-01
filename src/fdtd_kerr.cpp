/*
 * fdtd_kerr.cpp
 *
 *  Created on: 21 Dec 2020
 *      Author: Holger Schmitz
 */

#include "../huerto/constants.hpp"

#include <schnek/grid.hpp>
#include <schnek/tools/literature.hpp>

#include <boost/make_shared.hpp>

#include "fdtd_kerr.hpp"


//===============================================================
//==========  FDTD_Kerr (instantaneous nonlinearity)
//===============================================================

void FDTD_Kerr::initParameters(schnek::BlockParameters &blockPars)
{
  blockPars.addParameter("chi", &chi, 0.0);
  blockPars.addParameter("eps", &eps, 1.0);
}

void FDTD_Kerr::registerData() {
#ifdef HUERTO_ONE_DIM
  pKappaEdx = std::make_shared<Grid1d>();
  pKappaHdx = std::make_shared<Grid1d>();

  addData("KappaEdx", pKappaEdx);
  addData("KappaHdx", pKappaHdx);
#endif

#ifdef HUERTO_TWO_DIM
  pKappaEdx = std::make_shared<Grid1d>();
  pKappaEdy = std::make_shared<Grid1d>();

  pKappaHdx = std::make_shared<Grid1d>();
  pKappaHdy = std::make_shared<Grid1d>();

  addData("KappaEdx", pKappaEdx);
  addData("KappaEdy", pKappaEdy);

  addData("KappaHdx", pKappaHdx);
  addData("KappaHdy", pKappaHdy);
#endif

#ifdef HUERTO_THREE_DIM
  pKappaEdx = std::make_shared<Grid1d>();
  pKappaEdy = std::make_shared<Grid1d>();
  pKappaEdz = std::make_shared<Grid1d>();

  pKappaHdx = std::make_shared<Grid1d>();
  pKappaHdy = std::make_shared<Grid1d>();
  pKappaHdz = std::make_shared<Grid1d>();

  addData("KappaEdx", pKappaEdx);
  addData("KappaEdy", pKappaEdy);
  addData("KappaEdz", pKappaEdz);

  addData("KappaHdx", pKappaHdx);
  addData("KappaHdy", pKappaHdy);
  addData("KappaHdz", pKappaHdz);
#endif
}

void FDTD_Kerr::init() {
  SimulationEntity::init(this);

  schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();
  Index low  = subdivision.getLo();
  Index high = subdivision.getHi();

#ifdef HUERTO_ONE_DIM
  pKappaEdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaHdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  (*pKappaEdx) = 1.0;
  (*pKappaHdx) = 1.0;
#endif

#ifdef HUERTO_TWO_DIM
  pKappaEdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaEdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  (*pKappaEdx) = 1.0;
  (*pKappaEdy) = 1.0;

  pKappaHdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaHdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  (*pKappaHdx) = 1.0;
  (*pKappaHdy) = 1.0;
#endif

#ifdef HUERTO_THREE_DIM
  pKappaEdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaEdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  pKappaEdz->resize(schnek::Array<int, 1>(low[2]), schnek::Array<int, 1>(high[2]));
  (*pKappaEdx) = 1.0;
  (*pKappaEdy) = 1.0;
  (*pKappaEdz) = 1.0;

  pKappaHdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaHdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  pKappaHdz->resize(schnek::Array<int, 1>(low[2]), schnek::Array<int, 1>(high[2]));
  (*pKappaHdx) = 1.0;
  (*pKappaHdy) = 1.0;
  (*pKappaHdz) = 1.0;
#endif


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

  CurrentContainer::init(getContext());

  schnek::LiteratureArticle Yee1966("Yee1966", "Yee, K",
      "Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media.",
      "IEEE Transactions on Antennas and Propagation", "1966", "AP-14", "302--307");

  schnek::LiteratureManager::instance().addReference(
      "Integration of electrodynamic fields uses the Finite Difference Time Domain method.",
      Yee1966);
}

void FDTD_Kerr::stepSchemeInit(double dt) {
  stepB(0.5*dt);

  BOOST_FOREACH(pCurrent current, this->currents) {
    current->stepSchemeInit(dt);
  }

  BOOST_FOREACH(pCurrent current, this->magCurrents) {
    current->stepSchemeInit(dt);
  }
}

void FDTD_Kerr::stepScheme(double dt)
{
  BOOST_FOREACH(pCurrent current, this->currents) {
    current->stepScheme(dt);
  }

  stepD(dt);

  BOOST_FOREACH(pCurrent current, this->magCurrents) {
    current->stepScheme(dt);
  }

  stepB(dt);
}

#ifdef HUERTO_ONE_DIM
void FDTD_Kerr::stepD(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaEdx= *pKappaEdx;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
  {
    double jx = (*this->pJx)(i);
    double jy = (*this->pJy)(i);
    double jz = (*this->pJz)(i);

    double kappaEdx = rKappaEdx(i)*dx[0];

    double ex = Ex(i);
    double ey = Ey(i);
    double ez = Ez(i);

    double E2 = ex*ex + ey*ey + ez*ez;

    double Dx = (eps + chi*E2)*ex + dt*jx/eps_0;

    double Dy = (eps + chi*E2)*ey
      + dt*(
          clight2*(
          - (Bz(i) - Bz(i-1))/kappaEdx
          )
        - jy/eps_0
      );

    double Dz = (eps + chi*E2)*ez
      + dt*(
          clight2*(
            (By(i) - By(i-1))/kappaEdx
          )
        - jz/eps_0
      );
    double E = sqrt(E2);
    double D = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

    if (D == 0.0) {
      Ex(i) = 0.0;
      Ey(i) = 0.0;
      Ez(i) = 0.0;
      continue; 
    }

    // Newton iteration
    double Eold;
    do {
      Eold = E;
      E2 = E*E;
      E = (2*chi*E2*E + D) / (3*chi*E2 + eps);
    } while (fabs(E-Eold) > 1e-9);

    Ex(i) = E*Dx/D;
    Ey(i) = E*Dy/D;
    Ez(i) = E*Dz/D;
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_Kerr::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaHdx= *pKappaHdx;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i) {
    double jy = (*this->pMy)(i);
    double jz = (*this->pMz)(i);

    double kappaHdx = rKappaHdx(i)*dx[0];

    By(i) = By(i)
      + dt*(
          (Ez(i+1) - Ez(i))/kappaHdx
        + jy
      );

    Bz(i) = Bz(i)
      + dt*(
        - (Ey(i+1) - Ey(i))/kappaHdx
        + jz
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif

#ifdef HUERTO_TWO_DIM
void FDTD_Kerr::stepD(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaEdx= *pKappaEdx;
  Grid1d &rKappaEdy= *pKappaEdy;


  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
  {
    double jx = (*this->pJx)(i,j);
    double jy = (*this->pJy)(i,j);
    double jz = (*this->pJz)(i,j);

    double kappaEdx = rKappaEdx(i)*dx[0];
    double kappaEdy = rKappaEdy(j)*dx[1];

  //  if (jx != 0.0 || jy != 0.0 || jz != 0.0) {
  //    std::cout << "FDTD stepD " << i << " " << j << " "
  //        << "(" << Ex(i, j) << ", " << Ey(i, j) << ", " << Ez(i, j) << ") "
  //        << "(" << Bx(i, j) << ", " << By(i, j) << ", " << Bz(i, j) << ") "
  //        << "(" << jx << ", " << jy << ", " << jz << ")  "
  //        << kappaEdx << std::endl;
  //  }

    double ex = Ex(i,j);
    double ey = Ey(i,j);
    double ez = Ez(i,j);

    double E2 = ex*ex + ey*ey + ez*ez;

    double Dx = (eps + chi*E2)*ex
      + dt*(
          clight2*(
            (Bz(i,j) - Bz(i,j-1))/kappaEdy
          )
        - jx/eps_0
      );

    double Dy = (eps + chi*E2)*ey
      + dt*(
          clight2*(
          - (Bz(i,j) - Bz(i-1,j))/kappaEdx
          )
        - jy/eps_0
      );

    double Dz = (eps + chi*E2)*ez
      + dt*(
          clight2*(
            (By(i,j) - By(i-1,j))/kappaEdx
          - (Bx(i,j) - Bx(i,j-1))/kappaEdy
          )
        - jz/eps_0
      );

    double E = sqrt(E2);
    double D = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

    if (D == 0.0) {
      Ex(i,j) = 0.0;
      Ey(i,j) = 0.0;
      Ez(i,j) = 0.0;
      continue; 
    }
    // Newton iteration
    double Eold;
    do {
      Eold = E;
      E2 = E*E;
      E = (2*chi*E2*E + D) / (3*chi*E2 + eps);
    } while (fabs(E-Eold) > 1e-9);

    Ex(i,j) = E*Dx/D;
    Ey(i,j) = E*Dy/D;
    Ez(i,j) = E*Dz/D;
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_Kerr::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaHdx= *pKappaHdx;
  Grid1d &rKappaHdy= *pKappaHdy;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
  {
    double jx = (*this->pMx)(i,j);
    double jy = (*this->pMy)(i,j);
    double jz = (*this->pMz)(i,j);

    double kappaHdx = rKappaHdx(i)*dx[0];
    double kappaHdy = rKappaHdy(j)*dx[1];

//    if (jx != 0.0 || jy != 0.0 || jz != 0.0) {
//      std::cout << "FDTD stepB " << i << " " << j << " "
//          << "(" << Ex(i, j) << ", " << Ey(i, j) << ", " << Ez(i, j) << ") "
//          << "(" << Bx(i, j) << ", " << By(i, j) << ", " << Bz(i, j) << ") "
//          << "(" << jx << ", " << jy << ", " << jz << ")  "
//          << kappaHdx << std::endl;
//    }


    Bx(i,j) = Bx(i,j)
      + dt*(
        - (Ez(i,j+1) - Ez(i,j))/kappaHdy
        + jx
      );

    By(i,j) = By(i,j)
      + dt*(
          (Ez(i+1,j) - Ez(i,j))/kappaHdx
        + jy
      );

    Bz(i,j) = Bz(i,j)
      + dt*(
          (Ex(i,j+1) - Ex(i,j))/kappaHdy
        - (Ey(i+1,j) - Ey(i,j))/kappaHdx
        + jz
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif

#ifdef HUERTO_THREE_DIM
void FDTD_Kerr::stepD(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaEdx= *pKappaEdx;
  Grid1d &rKappaEdy= *pKappaEdy;
  Grid1d &rKappaEdz= *pKappaEdz;


  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double jx = (*this->pJx)(i,j,k);
    double jy = (*this->pJy)(i,j,k);
    double jz = (*this->pJz)(i,j,k);

    double kappaEdx = rKappaEdx(i)*dx[0];
    double kappaEdy = rKappaEdy(j)*dx[1];
    double kappaEdz = rKappaEdz(k)*dx[2];

    double ex = Ex(i,j,k);
    double ey = Ey(i,j,k);
    double ez = Ez(i,j,k);

    double E2 = ex*ex + ey*ey + ez*ez;

    double Dx = (eps + chi*E2)*ex
      + dt*(
          clight2*(
            (Bz(i,j,k) - Bz(i,j-1,k))/kappaEdy
          - (By(i,j,k) - By(i,j,k-1))/kappaEdz
          )
        - jx/eps_0
      );

    double Dy = (eps + chi*E2)*ey
      + dt*(
          clight2*(
            (Bx(i,j,k) - Bx(i,j,k-1))/kappaEdz
          - (Bz(i,j,k) - Bz(i-1,j,k))/kappaEdx
          )
        - jy/eps_0
      );

    double Dz = (eps + chi*E2)*ez
      + dt*(
          clight2*(
            (By(i,j,k) - By(i-1,j,k))/kappaEdx
          - (Bx(i,j,k) - Bx(i,j-1,k))/kappaEdy
          )
        - jz/eps_0
      );

    double E = sqrt(E2);
    double D = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

    if (D == 0.0) {
      Ex(i,j,k) = 0.0;
      Ey(i,j,k) = 0.0;
      Ez(i,j,k) = 0.0;
      continue; 
    }

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

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_Kerr::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaHdx= *pKappaHdx;
  Grid1d &rKappaHdy= *pKappaHdy;
  Grid1d &rKappaHdz= *pKappaHdz;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double jx = (*this->pMx)(i,j,k);
    double jy = (*this->pMy)(i,j,k);
    double jz = (*this->pMz)(i,j,k);

    double kappaHdx = rKappaHdx(i)*dx[0];
    double kappaHdy = rKappaHdy(j)*dx[1];
    double kappaHdz = rKappaHdz(k)*dx[2];

    Bx(i,j,k) = Bx(i,j,k)
      + dt*(
          (Ey(i,j,k+1) - Ey(i,j,k))/kappaHdz
        - (Ez(i,j+1,k) - Ez(i,j,k))/kappaHdy
        + jx
      );

    By(i,j,k) = By(i,j,k)
      + dt*(
          (Ez(i+1,j,k) - Ez(i,j,k))/kappaHdx
        - (Ex(i,j,k+1) - Ex(i,j,k))/kappaHdz
        + jy
      );

    Bz(i,j,k) = Bz(i,j,k)
      + dt*(
          (Ex(i,j+1,k) - Ex(i,j,k))/kappaHdy
        - (Ey(i+1,j,k) - Ey(i,j,k))/kappaHdx
        + jz
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif

//===============================================================
//==========  FDTD_KerrAverage (averaged nonlinearity)
//===============================================================

void FDTD_KerrAverage::initParameters(schnek::BlockParameters &blockPars)
{
  blockPars.addParameter("T", &tAverage);
  blockPars.addParameter("chi", &chi, 0.0);
  blockPars.addParameter("eps", &eps, 1.0);
}

void FDTD_KerrAverage::registerData() {
#ifdef HUERTO_ONE_DIM
  pKappaEdx = std::make_shared<Grid1d>();
  pKappaHdx = std::make_shared<Grid1d>();

  addData("KappaEdx", pKappaEdx);
  addData("KappaHdx", pKappaHdx);
#endif

#ifdef HUERTO_TWO_DIM
  pKappaEdx = std::make_shared<Grid1d>();
  pKappaEdy = std::make_shared<Grid1d>();

  pKappaHdx = std::make_shared<Grid1d>();
  pKappaHdy = std::make_shared<Grid1d>();

  addData("KappaEdx", pKappaEdx);
  addData("KappaEdy", pKappaEdy);

  addData("KappaHdx", pKappaHdx);
  addData("KappaHdy", pKappaHdy);
#endif

#ifdef HUERTO_THREE_DIM
  pKappaEdx = std::make_shared<Grid1d>();
  pKappaEdy = std::make_shared<Grid1d>();
  pKappaEdz = std::make_shared<Grid1d>();

  pKappaHdx = std::make_shared<Grid1d>();
  pKappaHdy = std::make_shared<Grid1d>();
  pKappaHdz = std::make_shared<Grid1d>();

  addData("KappaEdx", pKappaEdx);
  addData("KappaEdy", pKappaEdy);
  addData("KappaEdz", pKappaEdz);

  addData("KappaHdx", pKappaHdx);
  addData("KappaHdy", pKappaHdy);
  addData("KappaHdz", pKappaHdz);
#endif

  pE2xAverage = std::make_shared<Field>();
  pE2yAverage = std::make_shared<Field>();
  pE2zAverage = std::make_shared<Field>();
}

void FDTD_KerrAverage::init() {
  SimulationEntity::init(this);

  schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();
  Index low  = subdivision.getLo();
  Index high = subdivision.getHi();
  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  schnek::Range<double, DIMENSION> domainSize = subdivision.getInnerExtent(getContext().getSize());

  pE2xAverage->resize(lowIn, highIn, domainSize, exStaggerYee, 2);
  pE2yAverage->resize(lowIn, highIn, domainSize, eyStaggerYee, 2);
  pE2zAverage->resize(lowIn, highIn, domainSize, ezStaggerYee, 2);

#ifdef HUERTO_ONE_DIM
  pKappaEdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaHdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  (*pKappaEdx) = 1.0;
  (*pKappaHdx) = 1.0;
#endif

#ifdef HUERTO_TWO_DIM
  pKappaEdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaEdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  (*pKappaEdx) = 1.0;
  (*pKappaEdy) = 1.0;

  pKappaHdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaHdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  (*pKappaHdx) = 1.0;
  (*pKappaHdy) = 1.0;
#endif

#ifdef HUERTO_THREE_DIM
  pKappaEdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaEdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  pKappaEdz->resize(schnek::Array<int, 1>(low[2]), schnek::Array<int, 1>(high[2]));
  (*pKappaEdx) = 1.0;
  (*pKappaEdy) = 1.0;
  (*pKappaEdz) = 1.0;

  pKappaHdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaHdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  pKappaHdz->resize(schnek::Array<int, 1>(low[2]), schnek::Array<int, 1>(high[2]));
  (*pKappaHdx) = 1.0;
  (*pKappaHdy) = 1.0;
  (*pKappaHdz) = 1.0;
#endif


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

  CurrentContainer::init(getContext());

  schnek::LiteratureArticle Yee1966("Yee1966", "Yee, K",
      "Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media.",
      "IEEE Transactions on Antennas and Propagation", "1966", "AP-14", "302--307");

  schnek::LiteratureManager::instance().addReference(
      "Integration of electrodynamic fields uses the Finite Difference Time Domain method.",
      Yee1966);
}

void FDTD_KerrAverage::stepSchemeInit(double dt) {
  stepB(0.5*dt);

  BOOST_FOREACH(pCurrent current, this->currents) {
    current->stepSchemeInit(dt);
  }

  BOOST_FOREACH(pCurrent current, this->magCurrents) {
    current->stepSchemeInit(dt);
  }
}

void FDTD_KerrAverage::stepScheme(double dt)
{
  BOOST_FOREACH(pCurrent current, this->currents) {
    current->stepScheme(dt);
  }

  stepD(dt);

  BOOST_FOREACH(pCurrent current, this->magCurrents) {
    current->stepScheme(dt);
  }

  stepB(dt);
}

#ifdef HUERTO_ONE_DIM
void FDTD_KerrAverage::stepD(double dt) {
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &E2xAvg = *pE2xAverage;
  Field &E2yAvg = *pE2yAverage;
  Field &E2zAvg = *pE2zAverage;

  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaEdx= *pKappaEdx;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();
  
  const double eta = dt/tAverage;

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
  {
    double jx = (*this->pJx)(i);
    double jy = (*this->pJy)(i);
    double jz = (*this->pJz)(i);

    double kappaEdx = rKappaEdx(i)*dx[0];

    double &ex = Ex(i);
    double &ey = Ey(i);
    double &ez = Ez(i);

    double &e2xAvg = E2xAvg(i);
    double &e2yAvg = E2yAvg(i);
    double &e2zAvg = E2zAvg(i);

    ex = ex + dt*jx/(eps_0*(eps + chi*e2xAvg));

    ey = ey + dt*(
        clight2*(
        - (Bz(i) - Bz(i-1))/kappaEdx
        )
      - jy/eps_0
    )/(eps + chi*e2yAvg);

    ez = ez + dt*(
        clight2*(
          (By(i) - By(i-1))/kappaEdx
        )
      - jz/eps_0
    )/(eps + chi*e2zAvg);


    double eyp = 0.5*(ey + Ey(i + 1));
    double ezp = 0.5*(ez + Ez(i + 1));

    double exm = 0.5*(ex + Ex(i - 1));

    const double E2x = ex*ex   + eyp*eyp + ezp*ezp;
    const double E2y = exm*exm + ey*ey   + ez*ez;
    const double E2z = exm*exm + ey*ey   + ez*ez;

    e2xAvg = (1-eta)*e2xAvg + eta*E2x;
    e2yAvg = (1-eta)*e2yAvg + eta*E2y;
    e2zAvg = (1-eta)*e2zAvg + eta*E2z;
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_KerrAverage::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaHdx= *pKappaHdx;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i) {
    double jy = (*this->pMy)(i);
    double jz = (*this->pMz)(i);

    double kappaHdx = rKappaHdx(i)*dx[0];

    By(i) = By(i)
      + dt*(
          (Ez(i+1) - Ez(i))/kappaHdx
        + jy
      );

    Bz(i) = Bz(i)
      + dt*(
        - (Ey(i+1) - Ey(i))/kappaHdx
        + jz
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif

#ifdef HUERTO_TWO_DIM
void FDTD_KerrAverage::stepD(double dt) {
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &E2xAvg = *pE2xAverage;
  Field &E2yAvg = *pE2yAverage;
  Field &E2zAvg = *pE2zAverage;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaEdx= *pKappaEdx;
  Grid1d &rKappaEdy= *pKappaEdy;


  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();
  
  const double eta = dt/tAverage;

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
  {
    double jx = (*this->pJx)(i, j);
    double jy = (*this->pJy)(i, j);
    double jz = (*this->pJz)(i, j);

    double kappaEdx = rKappaEdx(i)*dx[0];
    double kappaEdy = rKappaEdy(j)*dx[1];

    double &ex = Ex(i, j);
    double &ey = Ey(i, j);
    double &ez = Ez(i, j);

    double &e2xAvg = E2xAvg(i, j);
    double &e2yAvg = E2yAvg(i, j);
    double &e2zAvg = E2zAvg(i, j);

    ex = ex + dt*(
        clight2*(
          (Bz(i,j) - Bz(i,j-1))/kappaEdy
        )
      - jx/eps_0
    )/(eps + chi*e2xAvg);

    ey = ey + dt*(
        clight2*(
        - (Bz(i,j) - Bz(i-1,j))/kappaEdx
        )
      - jy/eps_0
    )/(eps + chi*e2yAvg);

    ez = ez + dt*(
        clight2*(
          (By(i,j) - By(i-1,j))/kappaEdx
        - (Bx(i,j) - Bx(i,j-1))/kappaEdy
        )
      - jz/eps_0
    )/(eps + chi*e2zAvg);

    double eyax = 0.25*(ey + Ey(i+1, j) + Ey(i, j-1) + Ey(i, j-1));
    double ezax = 0.5*(ez + Ez(i+1, j));

    double exay = 0.25*(ex + Ex(i-1, j) + Ex(i, j+1) + Ex(i-1, j+1));
    double ezay = 0.5*(ez + Ez(i, j+1));

    double exaz = 0.5*(ex + Ex(i-1, j));
    double eyaz = 0.5*(ey + Ey(i, j-1));

    const double E2x = ex*ex     + eyax*eyax + ezax*ezax;
    const double E2y = exay*exay + ey*ey     + ezay*ezay;
    const double E2z = exaz*exaz + eyaz*eyaz + ez*ez;

    e2xAvg = (1-eta)*e2xAvg + eta*E2x;
    e2yAvg = (1-eta)*e2yAvg + eta*E2y;
    e2zAvg = (1-eta)*e2zAvg + eta*E2z;
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_KerrAverage::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaHdx= *pKappaHdx;
  Grid1d &rKappaHdy= *pKappaHdy;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
  {
    double jx = (*this->pMx)(i,j);
    double jy = (*this->pMy)(i,j);
    double jz = (*this->pMz)(i,j);

    double kappaHdx = rKappaHdx(i)*dx[0];
    double kappaHdy = rKappaHdy(j)*dx[1];

//    if (jx != 0.0 || jy != 0.0 || jz != 0.0) {
//      std::cout << "FDTD stepB " << i << " " << j << " "
//          << "(" << Ex(i, j) << ", " << Ey(i, j) << ", " << Ez(i, j) << ") "
//          << "(" << Bx(i, j) << ", " << By(i, j) << ", " << Bz(i, j) << ") "
//          << "(" << jx << ", " << jy << ", " << jz << ")  "
//          << kappaHdx << std::endl;
//    }


    Bx(i,j) = Bx(i,j)
      + dt*(
        - (Ez(i,j+1) - Ez(i,j))/kappaHdy
        + jx
      );

    By(i,j) = By(i,j)
      + dt*(
          (Ez(i+1,j) - Ez(i,j))/kappaHdx
        + jy
      );

    Bz(i,j) = Bz(i,j)
      + dt*(
          (Ex(i,j+1) - Ex(i,j))/kappaHdy
        - (Ey(i+1,j) - Ey(i,j))/kappaHdx
        + jz
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif

#ifdef HUERTO_THREE_DIM
void FDTD_KerrAverage::stepD(double dt) {
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &E2xAvg = *pE2xAverage;
  Field &E2yAvg = *pE2yAverage;
  Field &E2zAvg = *pE2zAverage;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaEdx= *pKappaEdx;
  Grid1d &rKappaEdy= *pKappaEdy;
  Grid1d &rKappaEdz= *pKappaEdz;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  const double eta = dt/tAverage;

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double jx = (*this->pJx)(i,j,k);
    double jy = (*this->pJy)(i,j,k);
    double jz = (*this->pJz)(i,j,k);

    double kappaEdx = rKappaEdx(i)*dx[0];
    double kappaEdy = rKappaEdy(j)*dx[1];
    double kappaEdz = rKappaEdz(k)*dx[2];

    double &ex = Ex(i,j,k);
    double &ey = Ey(i,j,k);
    double &ez = Ez(i,j,k);

    double &e2xAvg = E2xAvg(i, j, k);
    double &e2yAvg = E2yAvg(i, j, k);
    double &e2zAvg = E2zAvg(i, j, k);

    ex = ex
      + dt*(
          clight2*(
            (Bz(i,j,k) - Bz(i,j-1,k))/kappaEdy
          - (By(i,j,k) - By(i,j,k-1))/kappaEdz
          )
        - jx/eps_0
      )/(eps + chi*e2xAvg);

    ey = ey
      + dt*(
          clight2*(
            (Bx(i,j,k) - Bx(i,j,k-1))/kappaEdz
          - (Bz(i,j,k) - Bz(i-1,j,k))/kappaEdx
          )
        - jy/eps_0
      )/(eps + chi*e2yAvg);

    ez = ez
      + dt*(
          clight2*(
            (By(i,j,k) - By(i-1,j,k))/kappaEdx
          - (Bx(i,j,k) - Bx(i,j-1,k))/kappaEdy
          )
        - jz/eps_0
      )/(eps + chi*e2zAvg);

    double eyax = 0.25*(ey + Ey(i+1, j, k) + Ey(i, j-1, k) + Ey(i+1, j-1, k));
    double ezax = 0.25*(ez + Ez(i+1, j, k) + Ez(i, j, k-1) + Ez(i+1, j, k-1));

    double exay = 0.25*(ex + Ex(i-1, j, k) + Ex(i, j+1, k) + Ex(i-1, j+1, k));
    double ezay = 0.25*(ez + Ez(i, j+1, k) + Ez(i, j, k-1) + Ez(i, j+1, k-1));

    double exaz = 0.25*(ex + Ex(i-1, j, k) + Ex(i, j, k+1) + Ex(i-1, j, k+1));
    double eyaz = 0.25*(ey + Ey(i, j-1, k) + Ey(i, j, k+1) + Ey(i, j-1, k+1));

    const double E2x = ex*ex     + eyax*eyax + ezax*ezax;
    const double E2y = exay*exay + ey*ey     + ezay*ezay;
    const double E2z = exaz*exaz + eyaz*eyaz + ez*ez;

    e2xAvg = (1-eta)*e2xAvg + eta*E2x;
    e2yAvg = (1-eta)*e2yAvg + eta*E2y;
    e2zAvg = (1-eta)*e2zAvg + eta*E2z;
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_KerrAverage::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  Grid1d &rKappaHdx= *pKappaHdx;
  Grid1d &rKappaHdy= *pKappaHdy;
  Grid1d &rKappaHdz= *pKappaHdz;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double jx = (*this->pMx)(i,j,k);
    double jy = (*this->pMy)(i,j,k);
    double jz = (*this->pMz)(i,j,k);

    double kappaHdx = rKappaHdx(i)*dx[0];
    double kappaHdy = rKappaHdy(j)*dx[1];
    double kappaHdz = rKappaHdz(k)*dx[2];

    Bx(i,j,k) = Bx(i,j,k)
      + dt*(
          (Ey(i,j,k+1) - Ey(i,j,k))/kappaHdz
        - (Ez(i,j+1,k) - Ez(i,j,k))/kappaHdy
        + jx
      );

    By(i,j,k) = By(i,j,k)
      + dt*(
          (Ez(i+1,j,k) - Ez(i,j,k))/kappaHdx
        - (Ex(i,j,k+1) - Ex(i,j,k))/kappaHdz
        + jy
      );

    Bz(i,j,k) = Bz(i,j,k)
      + dt*(
          (Ex(i,j+1,k) - Ex(i,j,k))/kappaHdy
        - (Ey(i+1,j,k) - Ey(i,j,k))/kappaHdx
        + jz
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif
