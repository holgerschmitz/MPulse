#include "mpulse.h"
#include "cpml_border.h"
#include "fieldsolver.h"
#include "globals.h"
#include "storage.h"

#include <vector>

//===============================================================
//==========  CPMLBorder
//===============================================================


void CPMLBorder::initCurrents(Storage *storage, FieldSolver *solver)
{ 
  solver->addCurrent(
    new CPMLBorderECurrent(thickness, north, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addCurrent(
    new CPMLBorderECurrent(thickness, south, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addCurrent(
    new CPMLBorderECurrent(thickness, east, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addCurrent(
    new CPMLBorderECurrent(thickness, west, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addCurrent(
    new CPMLBorderECurrent(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addCurrent(
    new CPMLBorderECurrent(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  
  solver->addMagCurrent(
    new CPMLBorderHCurrent(thickness, north, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addMagCurrent(
    new CPMLBorderHCurrent(thickness, south, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addMagCurrent(
    new CPMLBorderHCurrent(thickness, east, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addMagCurrent(
    new CPMLBorderHCurrent(thickness, west, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addMagCurrent(
    new CPMLBorderHCurrent(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addMagCurrent(
    new CPMLBorderHCurrent(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  
  initCoefficients(storage);

}

void CPMLBorder::initCoefficients(Storage *storage)
{
  // initialize Kappas here
  
  GridIndex glow  = Globals::instance().gridLow();
  GridIndex ghigh = Globals::instance().gridHigh();

  GridIndex low  = storage->getLow();
  GridIndex high = storage->getHigh();
  
  std::vector<DataLine *> pKappaEdk(3);
  std::vector<DataLine *> pKappaHdk(3);

  std::vector<DataLine *> pCpmlSigmaE(3);
  std::vector<DataLine *> pCpmlSigmaH(3);

  pKappaEdk[0] = storage->addLine("KappaEdx", 0, 1.0);
  pKappaEdk[1] = storage->addLine("KappaEdy", 1, 1.0);
  pKappaEdk[2] = storage->addLine("KappaEdz", 2, 1.0);

  pKappaHdk[0] = storage->addLine("KappaHdx", 0, 1.0);
  pKappaHdk[1] = storage->addLine("KappaHdy", 1, 1.0);
  pKappaHdk[2] = storage->addLine("KappaHdz", 2, 1.0);
  
  std::cerr << "Field Size Edx " << pKappaEdk[0]->getLow()[0] << ", "  
                                 << pKappaEdk[0]->getHigh()[0] <<std::endl;
  
  std::cerr << "Field Size Edy " << pKappaEdk[1]->getLow()[0] << ", "  
                                 << pKappaEdk[1]->getHigh()[0] <<std::endl;
  
  std::cerr << "Field Size Edz " << pKappaEdk[2]->getLow()[0] << ", "  
                                 << pKappaEdk[2]->getHigh()[0] <<std::endl;
  
  std::cerr << "Field Size Hdx " << pKappaHdk[0]->getLow()[0] << ", "  
                                 << pKappaHdk[0]->getHigh()[0] <<std::endl;
  
  std::cerr << "Field Size Hdy " << pKappaHdk[1]->getLow()[0] << ", "  
                                 << pKappaHdk[1]->getHigh()[0] <<std::endl;
  
  std::cerr << "Field Size Hdz " << pKappaHdk[2]->getLow()[0] << ", "  
                                 << pKappaHdk[2]->getHigh()[0] <<std::endl;
  
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

ParameterMap* CPMLBorder::MakeParamMap (ParameterMap* pm)
{
  pm = CurrentFactory::MakeParamMap(pm);
  
  (*pm)["d"] = WParameter(new ParameterValue<int>(&this->thickness,8));
  (*pm)["kappaMax"] = WParameter(new ParameterValue<double>(&this->kappaMax,15));
  (*pm)["aMax"] = WParameter(new ParameterValue<double>(&this->aMax,0.25));
  (*pm)["sigmaMax"] = WParameter(new ParameterValue<double>(&this->sigmaMax,3));
  return pm;
}

//===============================================================
//==========  CPMLBorderOneD
//===============================================================


void CPMLBorderOneD::initCurrents(Storage *storage, FieldSolver *solver)
{ 

  solver->addCurrent(
    new CPMLBorderECurrent(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addCurrent(
    new CPMLBorderECurrent(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  
  solver->addMagCurrent(
    new CPMLBorderHCurrent(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  solver->addMagCurrent(
    new CPMLBorderHCurrent(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, 1)
  );
  
  initCoefficients(storage);

}

void CPMLBorderOneD::initCoefficients(Storage *storage)
{
  // initialize Kappas here
    
  GridIndex glow  = Globals::instance().gridLow();
  GridIndex ghigh = Globals::instance().gridHigh();

  GridIndex low  = storage->getLow();
  GridIndex high = storage->getHigh();
  
  DataLine *pKappaEdk;
  DataLine *pKappaHdk;

  pKappaEdk = storage->addLine("KappaEdz", 2);
  pKappaHdk = storage->addLine("KappaHdz", 2);
    
  *(storage->addLine("KappaEdy", 1)) = 1; 
  *(storage->addLine("KappaEdx", 0)) = 1;

  *(storage->addLine("KappaHdy", 1)) = 1;
  *(storage->addLine("KappaHdx", 0)) = 1;
    
  std::cerr << "Field Size Edx " << pKappaEdk->getLow()[0] << ", "  
                                 << pKappaEdk->getHigh()[0] <<std::endl;
  
  std::cerr << "Field Size Hdx " << pKappaHdk->getLow()[0] << ", "  
                                 << pKappaHdk->getHigh()[0] <<std::endl;
  
  
    if (low[2]<glow[2]+thickness)
    {
      (*pKappaEdk) = 1.0;
      
      (*pKappaHdk) = 1.0;
      
      for (int i=0; i<thickness; ++i)
      {
        double x  = 1 - double(i)/double(thickness);
        double x3 = x*x*x;
        (*pKappaEdk)(low[2]+i+1) = 1 + (this->kappaMax - 1)*x3;
      }
      
      for (int i=0; i<thickness; ++i)
      {
        double x  = 1 - (double(i)+0.5)/double(thickness);
        double x3 = x*x*x;
        (*pKappaHdk)(low[2]+i) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    if (high[2]>ghigh[2]-thickness)
    {
      for (int i=0; i<thickness; ++i)      
      {
        double x  = 1 - (double(i))/double(thickness);
        double x3 = x*x*x;
        (*pKappaEdk)(high[2]-i-1) = 1 + (this->kappaMax - 1)*x3;
      }
      for (int i=0; i<thickness; ++i)
      {
        double x  = 1 - double(i+0.5)/double(thickness);
        double x3 = x*x*x;
        (*pKappaHdk)(high[2]-i-1) = 1 + (this->kappaMax - 1)*x3;
      }
    }
}

ParameterMap* CPMLBorderOneD::MakeParamMap (ParameterMap* pm)
{
  pm = CurrentFactory::MakeParamMap(pm);
  
  (*pm)["d"] = WParameter(new ParameterValue<int>(&this->thickness,8));
  (*pm)["kappaMax"] = WParameter(new ParameterValue<double>(&this->kappaMax,15));
  (*pm)["aMax"] = WParameter(new ParameterValue<double>(&this->aMax,0.25));
  (*pm)["sigmaMax"] = WParameter(new ParameterValue<double>(&this->sigmaMax,3));
  return pm;
}

//===============================================================
//==========  CPMLBorderCurrent
//===============================================================


CPMLBorderCurrent::CPMLBorderCurrent( int thickness_, Direction dir_, bool isH_,
                                      double kappaMax_, double aMax_, double sigmaMax_, double eps_)
  : thickness(thickness_), dir(dir_), isH(isH_),
    kappaMax(kappaMax_), aMax(aMax_), sigmaMax(sigmaMax_), eps(eps_)
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

void CPMLBorderCurrent::makeCoeff(Storage *storage)
{
  double dt = Globals::instance().dt();
  
  GridIndex low  = pJx->getLow();
  GridIndex high = pJx->getHigh();
  
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
                                        double kappaMax_, double aMax_, double sigmaMax_, double eps_)
  : CPMLBorderCurrent(thickness_,dir_,false,kappaMax_,aMax_,sigmaMax_,eps_)
{}

void CPMLBorderECurrent::initStorage(Storage *storage)
{

  pJx = storage->addBorderLayer("CPMLPsiEx", dir, thickness, 1);
  pJy = storage->addBorderLayer("CPMLPsiEy", dir, thickness, 1);
  pJz = storage->addBorderLayer("CPMLPsiEz", dir, thickness, 1);

  switch (dir)
  {
    case east:  
    case west:
      pPsi[0] = pJy;
      pPsi[1] = pJz;
      pB[0] = &storage->getGrid("By");
      pB[1] = &storage->getGrid("Bz");
      pB[2] = &storage->getGrid("Bx");
      dx = Globals::instance().gridDX();
      break;
    case north: 
    case south:
      pPsi[0] = pJz;
      pPsi[1] = pJx;
      pB[0] = &storage->getGrid("Bz");
      pB[1] = &storage->getGrid("Bx");
      pB[2] = &storage->getGrid("By");
      dx = Globals::instance().gridDY();
      break;
    case up:    
    case down:
      pPsi[0] = pJx;
      pPsi[1] = pJy;
      pB[0] = &storage->getGrid("Bx");
      pB[1] = &storage->getGrid("By");
      pB[2] = &storage->getGrid("Bz");
      dx = Globals::instance().gridDZ();
      break;
  }
  
  if (pJx) makeCoeff(storage);
  else 
  {
    std::cerr << "No E Current here: Direction " << dir << "\n";
  }
}

void CPMLBorderECurrent::stepSchemeInit(double dt)
{}

void CPMLBorderECurrent::stepScheme(double dt)
{
  GridIndex low  = pPsi[0]->getLow();
  GridIndex high = pPsi[0]->getHigh();
    
  DataGrid &Psi0 = *pPsi[0];
  DataGrid &Psi1 = *pPsi[1];
  DataGrid &B0 = *pB[0];
  DataGrid &B1 = *pB[1];
  
  GridIndex ind, indn;
      
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2])
      {
        int j = ind[dim];
        GridIndex indm(ind);
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
                                        double kappaMax_, double aMax_, double sigmaMax_, double eps_)
  : CPMLBorderCurrent(thickness_,dir_,true,kappaMax_,aMax_,sigmaMax_,eps_)
{}

void CPMLBorderHCurrent::initStorage(Storage *storage)
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

  pJx = storage->addBorderLayer("CPMLPsiBx", dir, thickness, distance);
  pJy = storage->addBorderLayer("CPMLPsiBy", dir, thickness, distance);
  pJz = storage->addBorderLayer("CPMLPsiBz", dir, thickness, distance);

  switch (dir)
  {
    case east:  
    case west:
      pPsi[0] = pJy;
      pPsi[1] = pJz;
      pE[0] = &storage->getGrid("Ey");
      pE[1] = &storage->getGrid("Ez");
      pE[2] = &storage->getGrid("Ex");
      dx = Globals::instance().gridDX();
      break;
    case north: 
    case south:
      pPsi[0] = pJz;
      pPsi[1] = pJx;
      pE[0] = &storage->getGrid("Ez");
      pE[1] = &storage->getGrid("Ex");
      pE[2] = &storage->getGrid("Ey");
      dx = Globals::instance().gridDY();
      break;
    case up:    
    case down:
      pPsi[0] = pJx;
      pPsi[1] = pJy;
      pE[0] = &storage->getGrid("Ex");
      pE[1] = &storage->getGrid("Ey");
      pE[2] = &storage->getGrid("Ez");
      dx = Globals::instance().gridDZ();
      break;
  }
  
  if (pJx) makeCoeff(storage);
  else 
  {
    std::cerr << "No H Current here: Direction " << dir << "\n";
  }
}

void CPMLBorderHCurrent::stepSchemeInit(double dt)
{
  stepScheme(0.5*dt);
}

void CPMLBorderHCurrent::stepScheme(double dt)
{
  GridIndex low  = pPsi[0]->getLow();
  GridIndex high = pPsi[0]->getHigh();
    
  DataGrid &Psi0 = *pPsi[0];
  DataGrid &Psi1 = *pPsi[1];
  DataGrid &E0 = *pE[0];
  DataGrid &E1 = *pE[1];
  
  GridIndex ind, indn;
      
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2])
      {
        int j = ind[dim];
        GridIndex indp(ind);
        ++indp[dim];
        
        Psi0(ind[0], ind[1], ind[2]) 
          = bCoeff(j)*Psi0(ind[0], ind[1], ind[2]) 
            + cCoeff(j)*(E1(indp[0], indp[1], indp[2])-E1(ind[0], ind[1], ind[2]))/dx;
            
        Psi1(ind[0], ind[1], ind[2]) 
          = bCoeff(j)*Psi1(ind[0], ind[1], ind[2]) 
            - cCoeff(j)*(E0(indp[0], indp[1], indp[2])-E0(ind[0], ind[1], ind[2]))/dx;
      }
}


