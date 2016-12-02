#include "shortpulseinject.hpp"
#include "specfunc.hpp"

#include <schnek/tools/literature.hpp>

#include <cmath>

//===============================================================
//==========  ShortPulseInject
//===============================================================

void ShortPulseInject::init()
{

  schnek::LiteratureArticle anderBrugge2009("anderBrugge2009", "D. an der Br{\\\"a}ugge and A. Pukhov",
      "Ultrashort focused electromagnetic pulses.",
      "Physical Review E", "2009", "E 79", "016603");

  schnek::LiteratureManager::instance().addReference(
      "Source conditions for ultrashort focused EM pulses",
      anderBrugge2009);
}

pCurrent ShortPulseInject::makeECurrent(int distance_, Direction dir_)
{
  std::cerr << "ShortPulseInjectSourceFunc::makeECurrent\n";
  typedef IncidentSourceECurrent<ShortPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  
  cur->setParam(length, width, om0, TShift, ZShift, Phase, amp, eps, distance_);
  return pCurrent(cur);
}

pCurrent ShortPulseInject::makeHCurrent(int distance_, Direction dir_)
{
  std::cerr << "ShortPulseInjectSourceFunc::makeHCurrent\n";
  typedef IncidentSourceHCurrent<ShortPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(length, width, om0, TShift, ZShift, Phase, amp, eps, distance_);
  return pCurrent(cur);
}

bool ShortPulseInject::needCurrent(Direction dir_)
{
  return (dir_ == down);
}
    
void ShortPulseInject::initParameters(schnek::BlockParameters &blockPars)
{
  IncidentSource::initParameters(blockPars);

  blockPars.addParameter("length", &this->length, 1.0);
  blockPars.addParameter("width", &this->width, 1.0);
  blockPars.addParameter("om0", &this->om0, 2*M_PI);
  blockPars.addParameter("TShift", &this->TShift, 0.0);
  blockPars.addParameter("ZShift", &this->ZShift, 0.0);
  blockPars.addParameter("Phase", &this->Phase, 0.0);

  blockPars.addParameter("amp", &this->amp, 1.0);
  blockPars.addParameter("eps", &this->eps, 1.0);
}


//===============================================================
//==========  ShortPulseInjectSourceFunc
//===============================================================

void ShortPulseInjectSourceFunc::setParam(double length_,
                                          double width_,
                                          double om0_,
                                          double TShift_,
                                          double ZShift_,
                                          double Phase_,
                                          double amp_, 
                                          double eps_,
                                          int distance_)
{
  std::cerr << "ShortPulseInjectSourceFunc::setParam\n";
  
  
  // Grid Spacing and position
  
  DX = MPulse::getDx()[0];
  DY = MPulse::getDx()[1];
  DZ = MPulse::getDx()[2];
  DT = MPulse::getDt();


  Index gridLow  = Index(0);
  Index gridHigh = MPulse::getGlobalMax();
    
  centrex = 0.5*double(gridHigh[0] + gridLow[0]);
  centrey = 0.5*double(gridHigh[1] + gridLow[1]);
  centrez = 0.5*double(gridHigh[2] + gridLow[2]);
  
  // setting most parameters

  length = length_;
  width  = width_;
  om0     = om0_;
  
  eps = eps_;
  dist = distance_;
  
  lightspeed = sqrt(1/eps);
  
  ZRl = 0.5*om0*width*width/lightspeed;
  YComp = Complex(0.0,0.0);

  // Adjusting so that amplitude corresponds to maximum amplitude in the
  // center of the pulse
  
  Phase  = 0.0; //0.5*M_PI;
  TShift  = 0.0; 
  ZShift  = 0.0;

  double Exmax = Efunc(0, 0, 0, 0).real();

  // setting the rest of the parameters
  amp = amp_/Exmax;
  Phase  = 2*M_PI*Phase_;
  
  TShift  = TShift_;
  ZShift  = ZShift_;
}

void ShortPulseInjectSourceFunc
    ::initSourceFunc(pGrid pJx, pGrid pJy, pGrid pJz)
{
}

void ShortPulseInjectSourceFunc::setTime(double time)
{
  // Now calculating on the fly!
  return;
}

Vector ShortPulseInjectSourceFunc::getEField(int i, int j, int k, double time)
{


  double ex=0, ey=0;
  double posxo = (i-centrex)*DX;
  double posxh = (i+0.5-centrex)*DX;
  double posyo = (j-centrey)*DY;
  double posyh = (j+0.5-centrey)*DY;
  double poszo = (k-centrez)*DZ - ZShift;
  double posTime = time - TShift;

  Complex Exc = Efunc(posxh, posyo, poszo, posTime);
  ex = amp*Exc.real();
    
  if (YComp != Complex(0,0))
  {
    Complex Eyc = YComp*Efunc(posxo, posyh, poszo, posTime);
    ey = amp*Eyc.real();
  }
 
  return Vector(ex,ey,0);

}

Vector ShortPulseInjectSourceFunc::getHField(int i, int j, int k, double time)
{

  double bx=0, by=0;
  double posxo = (i-centrex)*DX;
  double posxh = (i+0.5-centrex)*DX;
  double posyo = (j-centrey)*DY;
  double posyh = (j+0.5-centrey)*DY;
  
// THIS IS A HACK
// it should be k+0.5 like in PlaneWaveSource
// I don't know why this is needed here.
  double poszh = (k-0.5-centrez)*DZ - ZShift;
//  double poszh = (k+0.5-centrez)*DZ - ZShift;
  
  double posTime = time-0.5*DT - TShift;

  Complex Bxc = Bfunc(posxo, posyh, poszh, posTime, true);
//  Complex Bxc = Bfunc(posxh, posyo, poszh, posTime, true);
  bx = amp*Bxc.real();
  
  Complex Byc = Bfunc(posxh, posyo, poszh, posTime, false);
//  Complex Byc = Bfunc(posxo, posyh, poszh, posTime, false);
  by = amp*Byc.real();

  return Vector(bx,by,0);

}


