#include "cpml_border.hpp"

#include "mpulse.hpp"

#include "../huerto/constants.hpp"
#include "../huerto/electromagnetics/fieldsolver.hpp"
#include "../huerto/electromagnetics/source/border.hpp"

#include <schnek/tools/literature.hpp>

#include <memory>

#include <vector>

//===============================================================
//==========  CPMLBorder
//===============================================================

void CPMLBorder::init()
{
  schnek::ChildBlock<CurrentBlock>::init();

  schnek::LiteratureArticle Roden2000("Roden2000", "Roden, J. A. and Gedney, S. D.",
      "An efficient fdtd implementation of the cfs-pml for arbitrary media",
      "Microwave and Optical Technology Letters", "2000", "27", "334--339");

  schnek::LiteratureManager::instance().addReference(
      "Implementation of the Convolution Perfectly Matched Layer (CPML) boundary condition", Roden2000);
}

void CPMLBorder::initCurrents(CurrentContainer &container)
{
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, east, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, west, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
#ifndef HUERTO_ONE_DIM
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, north, this->kappaMax, this->aMax, this->sigmaMax, 1.0, boost::ref(*this))
  );
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, south, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
#endif
#ifdef HUERTO_THREE_DIM
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
#endif

  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, east, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, west, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
#ifndef HUERTO_ONE_DIM
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, north, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, south, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
#endif
#ifdef HUERTO_THREE_DIM
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
#endif

  initCoefficients();

}

void CPMLBorder::initCoefficients()
{
  schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();
  // initialize Kappas here

  Range gdomain = subdivision.getGlobalDomain();
  Index glow  = gdomain.getLo();
  Index ghigh = gdomain.getHi();

  Index low  = subdivision.getInnerLo();
  Index high = subdivision.getInnerHi();

  std::vector<pGrid1d> pKappaEdk(DIMENSION);
  std::vector<pGrid1d> pKappaHdk(DIMENSION);

  std::vector<pGrid1d> pCpmlSigmaE(DIMENSION);
  std::vector<pGrid1d> pCpmlSigmaH(DIMENSION);

  retrieveData("KappaEdx", pKappaEdk[0]);
#ifndef HUERTO_ONE_DIM
  retrieveData("KappaEdy", pKappaEdk[1]);
#endif
#ifdef HUERTO_THREE_DIM
  retrieveData("KappaEdz", pKappaEdk[2]);
#endif

  retrieveData("KappaHdx", pKappaHdk[0]);
#ifndef HUERTO_ONE_DIM
  retrieveData("KappaHdy", pKappaHdk[1]);
#endif
#ifdef HUERTO_THREE_DIM
  retrieveData("KappaHdz", pKappaHdk[2]);
#endif

//  std::cerr << "Field Size Edx " << pKappaEdk[0]->getLo()[0] << ", "
//                                 << pKappaEdk[0]->getHi()[0] <<std::endl;
//
//  std::cerr << "Field Size Edy " << pKappaEdk[1]->getLo()[0] << ", "
//                                 << pKappaEdk[1]->getHi()[0] <<std::endl;
//
//  std::cerr << "Field Size Edz " << pKappaEdk[2]->getLo()[0] << ", "
//                                 << pKappaEdk[2]->getHi()[0] <<std::endl;
//
//  std::cerr << "Field Size Hdx " << pKappaHdk[0]->getLo()[0] << ", "
//                                 << pKappaHdk[0]->getHi()[0] <<std::endl;
//
//  std::cerr << "Field Size Hdy " << pKappaHdk[1]->getLo()[0] << ", "
//                                 << pKappaHdk[1]->getHi()[0] <<std::endl;
//
//  std::cerr << "Field Size Hdz " << pKappaHdk[2]->getLo()[0] << ", "
//                                 << pKappaHdk[2]->getHi()[0] <<std::endl;

  for (size_t dim = 0; dim<DIMENSION; ++dim)
  {
    std::cerr << "Dim " << dim << std::endl;

    (*pKappaEdk[dim]) = 1.0;
    (*pKappaHdk[dim]) = 1.0;


    Index blow, bhigh;
    Direction dir;

    switch (dim)
    {
      case 0: dir = west; break;
#ifndef HUERTO_ONE_DIM
      case 1: dir = south; break;
#endif
#ifdef HUERTO_THREE_DIM
      case 2: dir = down; break;
#endif
    }

    if (getBorderExtent(dir, thickness, 1, blow, bhigh, false, getContext()))
    {
      int lowk  = blow[dim];
      int highk = bhigh[dim];
      int kLimit = highk-lowk + 1;

      for (int k=0; k<kLimit; ++k)
      {
        double x = 1 - double(k)/double(thickness);
        double x3 = x*x*x;

        (*pKappaEdk[dim])(lowk+k) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    if (getBorderExtent(dir, thickness, 1, blow, bhigh, true, getContext()))
    {
      int lowk  = blow[dim];
      int highk = bhigh[dim];
      int kLimit = highk-lowk + 1;

      for (int k=0; k<kLimit; ++k)
      {
        double x = 1 - (double(k) + 0.5)/double(thickness);
        double x3 = x*x*x;

        (*pKappaHdk[dim])(lowk+k) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    switch (dim)
    {
      case 0: dir = east; break;
#ifndef HUERTO_ONE_DIM
      case 1: dir = north; break;
#endif
#ifdef HUERTO_THREE_DIM
      case 2: dir = up; break;
#endif
    }

    if (getBorderExtent(dir, thickness, 1, blow, bhigh, false, getContext()))
    {
      int lowk  = blow[dim];
      int highk = bhigh[dim];
      int kLimit = highk-lowk + 1;

      for (int k=0; k<kLimit; ++k)
      {
        double x = 1 - double(k)/double(thickness);
        double x3 = x*x*x;

        (*pKappaEdk[dim])(highk-k) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    if (getBorderExtent(dir, thickness, 1, blow, bhigh, true, getContext()))
    {
      int lowk  = blow[dim];
      int highk = bhigh[dim];
      int kLimit = highk-lowk + 1;

      for (int k=0; k<kLimit; ++k)
      {
        double x = 1 - (double(k) - 0.5)/double(thickness);
        double x3 = x*x*x;

        (*pKappaHdk[dim])(highk-k) = 1 + (this->kappaMax - 1)*x3;
      }
    }

//
//
//    if (low[dim]<glow[dim]+thickness)
//    {
//
//      for (int i=0; i<=thickness; ++i)
//      {
//        double x  = 1 - double(i)/double(thickness);
//        double x3 = x*x*x;
//        (*pKappaEdk[dim])(low[dim]+i) = 1 + (this->kappaMax - 1)*x3;
//      }
//      for (int i=0; i<thickness; ++i)
//      {
//        double x  = 1 - (double(i)+0.5)/double(thickness);
//        double x3 = x*x*x;
//        (*pKappaHdk[dim])(low[dim]+i) = 1 + (this->kappaMax - 1)*x3;
//      }
//    }
//
//    if (high[dim]>ghigh[dim]-thickness)
//    {
//      for (int i=0; i<=thickness; ++i)
//      {
//        double x  = 1 - double(i)/double(thickness);
//        double x3 = x*x*x;
//        (*pKappaHdk[dim])(high[dim]-i) = 1 + (this->kappaMax - 1)*x3;
//      }
//      for (int i=0; i<thickness; ++i)
//
//      {
//        double x  = 1 - (double(i)+0.5)/double(thickness);
//        double x3 = x*x*x;
//        (*pKappaEdk[dim])(high[dim]-i) = 1 + (this->kappaMax - 1)*x3;
//      }
//    }
  }
}

void CPMLBorder::initParameters(schnek::BlockParameters &blockPars)
{
  CurrentBlock::initParameters(blockPars);

  blockPars.addParameter("d", &this->thickness, 8);
  blockPars.addParameter("kappaMax", &this->kappaMax, 15.0);
  blockPars.addParameter("aMax", &this->aMax, 0.25);
  blockPars.addParameter("sigmaMax", &this->sigmaMax, 3.0);
}

//===============================================================
//==========  CPMLBorderCurrent
//===============================================================


CPMLBorderCurrent::CPMLBorderCurrent(int thickness, Direction dir, bool isH,
                                     double kappaMax, double aMax, double sigmaMax, double eps,
                                     CurrentBlock &borderBlock)
  : reverse(false), thickness(thickness), dir(dir), isH(isH), lowOffset(0), highOffset(0), zerolayer(0),
    kappaMax(kappaMax), aMax(aMax), sigmaMax(sigmaMax), eps(eps), borderBlock(borderBlock)
{
  switch (dir)
  {
    case east:
    case west:  dim = 0;
                transverse1 = 1;
                transverse2 = 2;
                break;
#ifndef HUERTO_ONE_DIM
    case north:
    case south: dim = 1;
                transverse1 = 0;
                transverse2 = 2;
                break;
#endif
#ifdef HUERTO_THREE_DIM
    case up:
    case down:  dim = 2;
                transverse1 = 0;
                transverse2 = 1;
                break;
#endif
  }

//  double eta = sqrt(mu_0/eps_0);
  sigmaMax = sigmaMax *clight* 0.8*4 / borderBlock.getContext().getDx()[dim];
}

void CPMLBorderCurrent::makeCoeff()
{
  double dt = borderBlock.getContext().getDt();

  Index low  = pJx->getLo();
  Index high = pJx->getHi();

  switch (dir)
  {
    case east:
#ifndef HUERTO_ONE_DIM
    case north:
#endif
#ifdef HUERTO_THREE_DIM
    case up:
#endif
      reverse = false;
      break;
    case west:
#ifndef HUERTO_ONE_DIM
    case south:
#endif
#ifdef HUERTO_THREE_DIM
    case down:
#endif
      reverse = true;
      break;
  }


  int lowk  = low[dim];
  int highk = high[dim];

  bCoeff.resize(lowk, highk);
  cCoeff.resize(lowk, highk);

  double offset = 0.0;
  lowOffset = 1;

  if (isH)
  {
    offset = 0.5;
    lowOffset = 0;
  }


  int kLimit = highk-lowk + 1;

  for (int k=0; k<kLimit; ++k)
  {
    double x = 1 - (double(k)-offset)/double(thickness);
    double x3 = x*x*x;

    int pos = reverse ? (lowk+k) : (highk-k);

    double sigma = x3*sigmaMax;
    double kappa = 1 + (kappaMax - 1)*x3;
    double a = aMax*(1-x);

    double b = exp(-(sigma/kappa + a)*dt);
    double c = sigma*(b-1)/(kappa*(sigma+kappa*a));

    bCoeff(pos) = b;
    cCoeff(pos) = c;
  }

}

//===============================================================
//==========  CPMLBorderECurrent
//===============================================================


CPMLBorderECurrent::CPMLBorderECurrent(int thickness,
                                       Direction dir,
                                       double kappaMax,
                                       double aMax,
                                       double sigmaMax,
                                       double eps,
                                       CurrentBlock &borderBlock)
  : CPMLBorderCurrent(thickness,dir,false,kappaMax,aMax,sigmaMax,eps,borderBlock),
    dx(0)
{}

void CPMLBorderECurrent::init()
{
  Index blow, bhigh;

  if (!getBorderExtent(dir, thickness, 1, blow, bhigh, false, borderBlock.getContext())) return;

  pJx = std::make_shared<Grid>(blow, bhigh);
  pJy = std::make_shared<Grid>(blow, bhigh);
  pJz = std::make_shared<Grid>(blow, bhigh);

  switch (dir)
  {
    case east:
    case west:
      pPsi[0] = pJy;
      pPsi[1] = pJz;
      borderBlock.retrieveData("By", pB[0]);
      borderBlock.retrieveData("Bz", pB[1]);
      borderBlock.retrieveData("Bx", pB[2]);

      dx = borderBlock.getContext().getDx()[0];
      break;
#ifndef HUERTO_ONE_DIM
    case north:
    case south:
      pPsi[0] = pJz;
      pPsi[1] = pJx;
      borderBlock.retrieveData("Bz", pB[0]);
      borderBlock.retrieveData("Bx", pB[1]);
      borderBlock.retrieveData("By", pB[2]);

      dx = borderBlock.getContext().getDx()[1];
      break;
#endif
#ifdef HUERTO_THREE_DIM
    case up:
    case down:
      pPsi[0] = pJx;
      pPsi[1] = pJy;
      borderBlock.retrieveData("Bx", pB[0]);
      borderBlock.retrieveData("By", pB[1]);
      borderBlock.retrieveData("Bz", pB[2]);

      dx = borderBlock.getContext().getDx()[2];
      break;
#endif
  }

  makeCoeff();
}

void CPMLBorderECurrent::stepSchemeInit(double dt)
{}

void CPMLBorderECurrent::stepScheme(double dt)
{
  Index low  = pPsi[0]->getLo();
  Index high = pPsi[0]->getHi();

  Grid &Psi0 = *pPsi[0];
  Grid &Psi1 = *pPsi[1];
  Field &B0 = *pB[0];
  Field &B1 = *pB[1];

  Index ind, indn;

  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0]) {
#ifndef HUERTO_ONE_DIM
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1]) {
#endif
#ifdef HUERTO_THREE_DIM
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2]) {
#endif

        int j = ind[dim];
        Index indm(ind);
        --indm[dim];

        Psi0[ind]
          = bCoeff(j)*Psi0[ind]
            - cCoeff(j)*(B1[ind]-B1[indm])/(mu_0*dx);
        Psi1[ind]
          = bCoeff(j)*Psi1[ind]
            + cCoeff(j)*(B0[ind]-B0[indm])/(mu_0*dx);
#ifdef HUERTO_THREE_DIM
      }
#endif
#ifndef HUERTO_ONE_DIM
    }
#endif
  }
}


//===============================================================
//==========  CPMLBorderHCurrent
//===============================================================

CPMLBorderHCurrent::CPMLBorderHCurrent(int thickness,
                                       Direction dir,
                                       double kappaMax,
                                       double aMax,
                                       double sigmaMax,
                                       double eps,
                                       CurrentBlock &borderBlock)
  : CPMLBorderCurrent(thickness,dir,true,kappaMax,aMax,sigmaMax,eps,borderBlock),
    dx(0)
{}

void CPMLBorderHCurrent::init()
{
  Index blow, bhigh;

  if (!getBorderExtent(dir, thickness, 1, blow, bhigh, true, borderBlock.getContext())) return;

  pJx = std::make_shared<Grid>(blow, bhigh);
  pJy = std::make_shared<Grid>(blow, bhigh);
  pJz = std::make_shared<Grid>(blow, bhigh);

  switch (dir)
  {
    case east:
    case west:
      pPsi[0] = pJy;
      pPsi[1] = pJz;
      borderBlock.retrieveData("Ey", pE[0]);
      borderBlock.retrieveData("Ez", pE[1]);
      borderBlock.retrieveData("Ex", pE[2]);

      dx = borderBlock.getContext().getDx()[0];
      break;
#ifndef HUERTO_ONE_DIM
    case north:
    case south:
      pPsi[0] = pJz;
      pPsi[1] = pJx;
      borderBlock.retrieveData("Ez", pE[0]);
      borderBlock.retrieveData("Ex", pE[1]);
      borderBlock.retrieveData("Ey", pE[2]);

      dx = borderBlock.getContext().getDx()[1];
      break;
#endif
#ifdef HUERTO_THREE_DIM
    case up:
    case down:
      pPsi[0] = pJx;
      pPsi[1] = pJy;
      borderBlock.retrieveData("Ex", pE[0]);
      borderBlock.retrieveData("Ey", pE[1]);
      borderBlock.retrieveData("Ez", pE[2]);

      dx = borderBlock.getContext().getDx()[2];
      break;
#endif
  }

  makeCoeff();

}

void CPMLBorderHCurrent::stepSchemeInit(double dt)
{
  stepScheme(0.5*dt);
}

void CPMLBorderHCurrent::stepScheme(double dt)
{
  Index low  = pPsi[0]->getLo();
  Index high = pPsi[0]->getHi();

  Grid &Psi0 = *pPsi[0];
  Grid &Psi1 = *pPsi[1];
  Field &E0 = *pE[0];
  Field &E1 = *pE[1];

  Index ind, indn;

  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0]) {
#ifndef HUERTO_ONE_DIM
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1]) {
#endif
#ifdef HUERTO_THREE_DIM
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2]) {
#endif
        int j = ind[dim];
        Index indp(ind);
        ++indp[dim];

        Psi0[ind]
          = bCoeff(j)*Psi0[ind]
            + cCoeff(j)*(E1[indp]-E1[ind])/dx;

        Psi1[ind]
          = bCoeff(j)*Psi1[ind]
            - cCoeff(j)*(E0[indp]-E0[ind])/dx;
#ifdef HUERTO_THREE_DIM
      }
#endif
#ifndef HUERTO_ONE_DIM
    }
#endif
  }
}


