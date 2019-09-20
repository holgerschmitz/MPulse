/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "fdtd_plain.hpp"

#include <schnek/grid.hpp>
#include <schnek/tools/literature.hpp>

#include <boost/make_shared.hpp>

void FDTD_Plain::registerData()
{
  pKappaEdx = boost::make_shared<DataLine>();
  pKappaEdy = boost::make_shared<DataLine>();
  pKappaEdz = boost::make_shared<DataLine>();

  pKappaHdx = boost::make_shared<DataLine>();
  pKappaHdy = boost::make_shared<DataLine>();
  pKappaHdz = boost::make_shared<DataLine>();

  addData("KappaEdx", pKappaEdx);
  addData("KappaEdy", pKappaEdy);
  addData("KappaEdz", pKappaEdz);

  addData("KappaHdx", pKappaHdx);
  addData("KappaHdy", pKappaHdy);
  addData("KappaHdz", pKappaHdz);

}

void FDTD_Plain::init()
{
  schnek::DomainSubdivision<Field> &subdivision = MPulse::getSubdivision();
  Index low  = subdivision.getLo();
  Index high = subdivision.getHi();

  pKappaEdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaEdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  pKappaEdz->resize(schnek::Array<int, 1>(low[2]), schnek::Array<int, 1>(high[2]));

  pKappaHdx->resize(schnek::Array<int, 1>(low[0]), schnek::Array<int, 1>(high[0]));
  pKappaHdy->resize(schnek::Array<int, 1>(low[1]), schnek::Array<int, 1>(high[1]));
  pKappaHdz->resize(schnek::Array<int, 1>(low[2]), schnek::Array<int, 1>(high[2]));


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

  CurrentContainer::init();

  schnek::LiteratureArticle Yee1966("Yee1966", "Yee, K",
      "Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media.",
      "IEEE Transactions on Antennas and Propagation", "1966", "AP-14", "302--307");

  schnek::LiteratureManager::instance().addReference(
      "Integration of electrodynamic fields uses the Finite Difference Time Domain method.",
      Yee1966);
}

void FDTD_Plain::stepSchemeInit(double dt)
{
  stepB(0.5*dt);

  BOOST_FOREACH(pCurrent current, this->currents)
  {
    current->stepSchemeInit(dt);
  }

  BOOST_FOREACH(pCurrent current, this->magCurrents)
  {
    current->stepSchemeInit(dt);
  }
}

void FDTD_Plain::stepScheme(double dt)
{
  BOOST_FOREACH(pCurrent current, this->currents)
  {
    current->stepScheme(dt);
  }

  stepD(dt);


  BOOST_FOREACH(pCurrent current, this->magCurrents)
  {
    current->stepScheme(dt);
  }

  stepB(dt);
}

void FDTD_Plain::stepD(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;
  
  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  DataLine &rKappaEdx= *pKappaEdx;
  DataLine &rKappaEdy= *pKappaEdy;
  DataLine &rKappaEdz= *pKappaEdz;

  
  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = MPulse::getDx();

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

    if (doDiag(i,j,k)) {
      std::cout << "Before stepD" << std::endl;
      std::cout << "Pos " << i << " " << j << " " << k << std::endl;
      std::cout << "J " << jx << " " << jy << " " << jz << std::endl;
      std::cout << "B " << Bx(i,j,k) << " " << By(i,j,k) << " " << Bz(i,j,k) << std::endl;
      std::cout << "E " << Ex(i,j,k) << " " << Ey(i,j,k) << " " << Ez(i,j,k) << std::endl;
      std::cout << "dt*c^2 " << dt*clight2 << std::endl;
      std::cout << "kappaEdx " << kappaEdx << std::endl;
    }

    Ex(i,j,k) = Ex(i,j,k) 
      + dt*(
          clight2*(
            (Bz(i,j,k) - Bz(i,j-1,k))/kappaEdy
          - (By(i,j,k) - By(i,j,k-1))/kappaEdz
          )
        + jx/eps_0
      );
      
    Ey(i,j,k) = Ey(i,j,k) 
      + dt*(
          clight2*(
            (Bx(i,j,k) - Bx(i,j,k-1))/kappaEdz
          - (Bz(i,j,k) - Bz(i-1,j,k))/kappaEdx
          )
        + jy/eps_0
      );
 
    Ez(i,j,k) = Ez(i,j,k) 
      + dt*(
          clight2*(
            (By(i,j,k) - By(i-1,j,k))/kappaEdx
          - (Bx(i,j,k) - Bx(i,j-1,k))/kappaEdy
          )
        + jz/eps_0
      );

    if (doDiag(i,j,k)) {
      std::cout << "After stepD" << std::endl;
      std::cout << "E " << Ex(i,j,k) << " " << Ey(i,j,k) << " " << Ez(i,j,k) << std::endl;
    }

  }

  schnek::DomainSubdivision<Field> &sub = MPulse::getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_Plain::stepB(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;
  
  Field &Bx = *pBx;
  Field &By = *pBy;
  Field &Bz = *pBz;

  DataLine &rKappaHdx= *pKappaHdx;
  DataLine &rKappaHdy= *pKappaHdy;
  DataLine &rKappaHdz= *pKappaHdz;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = MPulse::getDx();

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

    if (doDiag(i,j,k)) {
      std::cout << "Before stepB" << std::endl;
      std::cout << "Pos " << i << " " << j << " " << k << std::endl;
      std::cout << "J " << jx << " " << jy << " " << jz << std::endl;
      std::cout << "B " << Bx(i,j,k) << " " << By(i,j,k) << " " << Bz(i,j,k) << std::endl;
      std::cout << "E " << Ex(i,j,k) << " " << Ey(i,j,k) << " " << Ez(i,j,k) << std::endl;
      std::cout << "dt " << dt << std::endl;
      std::cout << "kappaHdx " << kappaHdx << std::endl;
    }


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

    if (doDiag(i,j,k)) {
      std::cout << "After stepB" << std::endl;
      std::cout << "B " << Bx(i,j,k) << " " << By(i,j,k) << " " << Bz(i,j,k) << std::endl;
    }
  }

  schnek::DomainSubdivision<Field> &sub = MPulse::getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
