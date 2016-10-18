#include "cpml_border.hpp"

#include "mpulse.hpp"
#include "border.hpp"
#include "fieldsolver.hpp"

#include <boost/make_shared.hpp>

#include <vector>

//===============================================================
//==========  CPMLBorder
//===============================================================


void CPMLBorder::initCurrents(CurrentContainer &container)
{ 
  container.addCurrent(
    boost::make_shared<CPMLBorderECurrent>(thickness, north, this->kappaMax, this->aMax, this->sigmaMax, 1.0, boost::ref(*this))
  );
  container.addCurrent(
      boost::make_shared<CPMLBorderECurrent>(thickness, south, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addCurrent(
      boost::make_shared<CPMLBorderECurrent>(thickness, east, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addCurrent(
      boost::make_shared<CPMLBorderECurrent>(thickness, west, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addCurrent(
      boost::make_shared<CPMLBorderECurrent>(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addCurrent(
      boost::make_shared<CPMLBorderECurrent>(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  
  container.addMagCurrent(
      boost::make_shared<CPMLBorderHCurrent>(thickness, north, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
      boost::make_shared<CPMLBorderHCurrent>(thickness, south, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
      boost::make_shared<CPMLBorderHCurrent>(thickness, east, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
      boost::make_shared<CPMLBorderHCurrent>(thickness, west, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
      boost::make_shared<CPMLBorderHCurrent>(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
      boost::make_shared<CPMLBorderHCurrent>(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  
  initCoefficients();

}

void CPMLBorder::initCoefficients()
{
  schnek::DomainSubdivision<Field> &subdivision = MPulse::getSubdivision();
  // initialize Kappas here
  
  Index glow  = Index(0);
  Index ghigh = MPulse::getGlobalMax();

  Index low  = subdivision.getInnerLo();
  Index high = subdivision.getInnerHi();
  
  std::vector<pDataLine> pKappaEdk(3);
  std::vector<pDataLine> pKappaHdk(3);

  std::vector<pDataLine> pCpmlSigmaE(3);
  std::vector<pDataLine> pCpmlSigmaH(3);

  retrieveData("KappaEdx", pKappaEdk[0]);
  retrieveData("KappaEdy", pKappaEdk[1]);
  retrieveData("KappaEdz", pKappaEdk[2]);

  retrieveData("KappaHdx", pKappaHdk[0]);
  retrieveData("KappaHdy", pKappaHdk[1]);
  retrieveData("KappaHdz", pKappaHdk[2]);
  
  std::cerr << "Field Size Edx " << pKappaEdk[0]->getLo()[0] << ", "
                                 << pKappaEdk[0]->getHi()[0] <<std::endl;
  
  std::cerr << "Field Size Edy " << pKappaEdk[1]->getLo()[0] << ", "
                                 << pKappaEdk[1]->getHi()[0] <<std::endl;
  
  std::cerr << "Field Size Edz " << pKappaEdk[2]->getLo()[0] << ", "
                                 << pKappaEdk[2]->getHi()[0] <<std::endl;
  
  std::cerr << "Field Size Hdx " << pKappaHdk[0]->getLo()[0] << ", "
                                 << pKappaHdk[0]->getHi()[0] <<std::endl;
  
  std::cerr << "Field Size Hdy " << pKappaHdk[1]->getLo()[0] << ", "
                                 << pKappaHdk[1]->getHi()[0] <<std::endl;
  
  std::cerr << "Field Size Hdz " << pKappaHdk[2]->getLo()[0] << ", "
                                 << pKappaHdk[2]->getHi()[0] <<std::endl;
  
  for (int dim = 0; dim<3; ++dim)
  {
    std::cerr << "Dim " << dim << std::endl;
    if (low[dim]<glow[dim]+thickness)
    {
      (*pKappaEdk[dim]) = 1.0;
      
      (*pKappaHdk[dim]) = 1.0;
      
      for (int i=0; i<=thickness; ++i)
      {
        double x  = 1 - double(i)/double(thickness);
        double x3 = x*x*x;
        (*pKappaEdk[dim])(low[dim]+i) = 1 + (this->kappaMax - 1)*x3;
      }
      for (int i=0; i<thickness; ++i)
      {
        double x  = 1 - (double(i)+0.5)/double(thickness);
        double x3 = x*x*x;
        (*pKappaHdk[dim])(low[dim]+i) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    if (high[dim]>ghigh[dim]-thickness)
    {
      for (int i=0; i<=thickness; ++i)
      {
        double x  = 1 - double(i)/double(thickness);
        double x3 = x*x*x;
        (*pKappaHdk[dim])(high[dim]-i) = 1 + (this->kappaMax - 1)*x3;
      }
      for (int i=0; i<thickness; ++i)
      
      {
        double x  = 1 - (double(i)+0.5)/double(thickness);
        double x3 = x*x*x;
        (*pKappaEdk[dim])(high[dim]-i) = 1 + (this->kappaMax - 1)*x3;
      }
    }
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
//==========  CPMLBorderOneD
//===============================================================


void CPMLBorderOneD::initCurrents(CurrentContainer &container)
{ 

  container.addCurrent(
      boost::make_shared<CPMLBorderECurrent>(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addCurrent(
      boost::make_shared<CPMLBorderECurrent>(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  
  container.addMagCurrent(
      boost::make_shared<CPMLBorderHCurrent>(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  container.addMagCurrent(
      boost::make_shared<CPMLBorderHCurrent>(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1, boost::ref(*this))
  );
  
  initCoefficients();

}

void CPMLBorderOneD::initCoefficients()
{
  // initialize Kappas here

  schnek::DomainSubdivision<Field> &subdivision = MPulse::getSubdivision();
  Index glow  = Index(0);
  Index ghigh = MPulse::getGlobalMax();

  Index low  = subdivision.getInnerLo();
  Index high = subdivision.getInnerHi();

  std::vector<pDataLine> pKappaEdk(3);
  std::vector<pDataLine> pKappaHdk(3);

  std::vector<pDataLine> pCpmlSigmaE(3);
  std::vector<pDataLine> pCpmlSigmaH(3);

  retrieveData("KappaEdx", pKappaEdk[0]);
  retrieveData("KappaEdy", pKappaEdk[1]);
  retrieveData("KappaEdz", pKappaEdk[2]);

  retrieveData("KappaHdx", pKappaHdk[0]);
  retrieveData("KappaHdy", pKappaHdk[1]);
  retrieveData("KappaHdz", pKappaHdk[2]);
    
  *(pKappaEdk[0]) = 1;
  *(pKappaEdk[1]) = 1;

  *(pKappaHdk[0]) = 1;
  *(pKappaHdk[1]) = 1;
    
  std::cerr << "Field Size Edx " << pKappaEdk[2]->getLo()[0] << ", "
                                 << pKappaEdk[2]->getHi()[0] <<std::endl;
  
  std::cerr << "Field Size Hdx " << pKappaHdk[2]->getLo()[0] << ", "
                                 << pKappaHdk[2]->getHi()[0] <<std::endl;
  
  
    if (low[2]<glow[2]+thickness)
    {
      (*pKappaEdk[2]) = 1.0;
      
      (*pKappaHdk[2]) = 1.0;
      
      for (int i=0; i<thickness; ++i)
      {
        double x  = 1 - double(i)/double(thickness);
        double x3 = x*x*x;
        (*pKappaEdk[2])(low[2]+i+1) = 1 + (this->kappaMax - 1)*x3;
      }
      
      for (int i=0; i<thickness; ++i)
      {
        double x  = 1 - (double(i)+0.5)/double(thickness);
        double x3 = x*x*x;
        (*pKappaHdk[2])(low[2]+i) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    if (high[2]>ghigh[2]-thickness)
    {
      for (int i=0; i<thickness; ++i)      
      {
        double x  = 1 - (double(i))/double(thickness);
        double x3 = x*x*x;
        (*pKappaEdk[2])(high[2]-i-1) = 1 + (this->kappaMax - 1)*x3;
      }
      for (int i=0; i<thickness; ++i)
      {
        double x  = 1 - double(i+0.5)/double(thickness);
        double x3 = x*x*x;
        (*pKappaHdk[2])(high[2]-i-1) = 1 + (this->kappaMax - 1)*x3;
      }
    }
}

void CPMLBorderOneD::initParameters(schnek::BlockParameters &blockPars)
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


CPMLBorderCurrent::CPMLBorderCurrent( int thickness_, Direction dir_, bool isH_,
                                      double kappaMax_, double aMax_, double sigmaMax_, double eps_,
                                      CurrentBlock &borderBlock_)
  : thickness(thickness_), dir(dir_), isH(isH_),
    kappaMax(kappaMax_), aMax(aMax_), sigmaMax(sigmaMax_), eps(eps_), borderBlock(borderBlock_)
{
  switch (dir)
  {
    case east:  
    case west:  dim = 0;
                transverse1 = 1;
                transverse2 = 2;
                break;
    case north: 
    case south: dim = 1;
                transverse1 = 0;
                transverse2 = 2;
                break;
    case up:    
    case down:  dim = 2;
                transverse1 = 0;
                transverse2 = 1;
                break;
  }
}

void CPMLBorderCurrent::makeCoeff()
{
  double dt = MPulse::getDt();
  
  Index low  = pJx->getLo();
  Index high = pJx->getHi();
  
  switch (dir)
  {
    case east:  
    case north: 
    case up:    reverse = false; break;
    case west:  
    case south: 
    case down:  reverse = true;  break;
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
  
    
  for (int k=0; k<thickness; ++k)
  {
    double x = 1 - (double(k)-offset)/double(thickness);
    double x3 = x*x*x;
    
    int pos = reverse ? (lowk+k) : (highk-k);
    
    double sigma = x3*sigmaMax;
    double kappa = 1 + (kappaMax - 1)*x3;
    double a = aMax*(1-x);
    
    double b = exp(-(sigma/kappa + a)*dt/eps);
    double c = sigma*(b-1)/(kappa*(sigma+kappa*a));
        
    bCoeff(pos) = b;
    cCoeff(pos) = c;
  }
  
}

//===============================================================
//==========  CPMLBorderECurrent
//===============================================================


CPMLBorderECurrent::CPMLBorderECurrent( int thickness_, Direction dir_, 
                                        double kappaMax_, double aMax_, double sigmaMax_, double eps_,
                                        CurrentBlock &borderBlock_)
  : CPMLBorderCurrent(thickness_,dir_,false,kappaMax_,aMax_,sigmaMax_,eps_,borderBlock_)
{}

void CPMLBorderECurrent::init()
{
  Index blow, bhigh;

  if (!getBorderExtent(dir, thickness, 1, blow, bhigh)) return;

  pJx = boost::make_shared<Grid>(blow, bhigh);
  pJy = boost::make_shared<Grid>(blow, bhigh);
  pJz = boost::make_shared<Grid>(blow, bhigh);

  switch (dir)
  {
    case east:  
    case west:
      pPsi[0] = pJy;
      pPsi[1] = pJz;
      borderBlock.retrieveData("By", pB[0]);
      borderBlock.retrieveData("Bz", pB[1]);
      borderBlock.retrieveData("Bx", pB[2]);

      dx = MPulse::getDx()[0];
      break;
    case north: 
    case south:
      pPsi[0] = pJz;
      pPsi[1] = pJx;
      borderBlock.retrieveData("Bz", pB[0]);
      borderBlock.retrieveData("Bx", pB[1]);
      borderBlock.retrieveData("By", pB[2]);

      dx = MPulse::getDx()[1];
      break;
    case up:    
    case down:
      pPsi[0] = pJx;
      pPsi[1] = pJy;
      borderBlock.retrieveData("Bx", pB[0]);
      borderBlock.retrieveData("By", pB[1]);
      borderBlock.retrieveData("Bz", pB[2]);

      dx = MPulse::getDx()[2];
      break;
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
      
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2])
      {
        int j = ind[dim];
        Index indm(ind);
        --indm[dim];
        
        Psi0(ind[0], ind[1], ind[2]) 
          = bCoeff(j)*Psi0(ind[0], ind[1], ind[2]) 
            - cCoeff(j)*(B1(ind[0], ind[1], ind[2])-B1(indm[0], indm[1], indm[2]))/dx;
        Psi1(ind[0], ind[1], ind[2]) 
          = bCoeff(j)*Psi1(ind[0], ind[1], ind[2]) 
            + cCoeff(j)*(B0(ind[0], ind[1], ind[2])-B0(indm[0], indm[1], indm[2]))/dx;
      }
}


//===============================================================
//==========  CPMLBorderHCurrent
//===============================================================

CPMLBorderHCurrent::CPMLBorderHCurrent( int thickness_, Direction dir_, 
                                        double kappaMax_, double aMax_, double sigmaMax_, double eps_,
                                        CurrentBlock &borderBlock_)
  : CPMLBorderCurrent(thickness_,dir_,true,kappaMax_,aMax_,sigmaMax_,eps_,borderBlock_)
{}

void CPMLBorderHCurrent::init()
{
  int distance;
  
  switch (dir)
  {
    case east:  
    case north: 
    case up:    distance = 1; break;
    case west:  
    case south: 
    case down:  distance = 0; break;
    default:    distance = 0; break;
  }

  Index blow, bhigh;

  if (!getBorderExtent(dir, thickness, distance, blow, bhigh)) return;

  pJx = boost::make_shared<Grid>(blow, bhigh);
  pJy = boost::make_shared<Grid>(blow, bhigh);
  pJz = boost::make_shared<Grid>(blow, bhigh);

  switch (dir)
  {
    case east:  
    case west:
      pPsi[0] = pJy;
      pPsi[1] = pJz;
      borderBlock.retrieveData("Ey", pE[0]);
      borderBlock.retrieveData("Ez", pE[1]);
      borderBlock.retrieveData("Ex", pE[2]);

      dx = MPulse::getDx()[0];
      break;
    case north: 
    case south:
      pPsi[0] = pJz;
      pPsi[1] = pJx;
      borderBlock.retrieveData("Ez", pE[0]);
      borderBlock.retrieveData("Ex", pE[1]);
      borderBlock.retrieveData("Ey", pE[2]);

      dx = MPulse::getDx()[1];
      break;
    case up:    
    case down:
      pPsi[0] = pJx;
      pPsi[1] = pJy;
      borderBlock.retrieveData("Ex", pE[0]);
      borderBlock.retrieveData("Ey", pE[1]);
      borderBlock.retrieveData("Ez", pE[2]);

      dx = MPulse::getDx()[2];
      break;
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
      
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2])
      {
        int j = ind[dim];
        Index indp(ind);
        ++indp[dim];
        
        Psi0(ind[0], ind[1], ind[2]) 
          = bCoeff(j)*Psi0(ind[0], ind[1], ind[2]) 
            + cCoeff(j)*(E1(indp[0], indp[1], indp[2])-E1(ind[0], ind[1], ind[2]))/dx;
            
        Psi1(ind[0], ind[1], ind[2]) 
          = bCoeff(j)*Psi1(ind[0], ind[1], ind[2]) 
            - cCoeff(j)*(E0(indp[0], indp[1], indp[2])-E0(ind[0], ind[1], ind[2]))/dx;
      }
}


