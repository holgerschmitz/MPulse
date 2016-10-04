#include "shortpulseinject.h"
#include <cmath>
#include "specfunc.h"

//===============================================================
//==========  ShortPulseInject
//===============================================================

IncidentSourceCurrent *ShortPulseInject::makeECurrent(int distance_, Direction dir_)
{
  std::cerr << "ShortPulseInjectSourceFunc::makeECurrent\n";
  typedef IncidentSourceECurrent<ShortPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  
  cur->setParam(length, width, om0, TShift, ZShift, Phase, amp, eps, distance_);
  return cur;
}

IncidentSourceCurrent *ShortPulseInject::makeHCurrent(int distance_, Direction dir_)
{
  std::cerr << "ShortPulseInjectSourceFunc::makeHCurrent\n";
  typedef IncidentSourceHCurrent<ShortPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(length, width, om0, TShift, ZShift, Phase, amp, eps, distance_);
  return cur;
}

bool ShortPulseInject::needCurrent(Direction dir_)
{
  return (dir_ == down);
}
    
ParameterMap* ShortPulseInject::MakeParamMap (ParameterMap* pm)
{
  std::cerr << "ShortPulseInjectSourceFunc::MakeParamMap\n";
  pm = IncidentSource::MakeParamMap(pm);
  
  (*pm)["length"] = WParameter(new ParameterValue<double>(&this->length,1.));
  (*pm)["width"] = WParameter(new ParameterValue<double>(&this->width,1.));
  (*pm)["om0"] = WParameter(new ParameterValue<double>(&this->om0,2*M_PI));
  (*pm)["TShift"] = WParameter(new ParameterValue<double>(&this->TShift,0));
  (*pm)["ZShift"] = WParameter(new ParameterValue<double>(&this->ZShift,0));
  (*pm)["Phase"] = WParameter(new ParameterValue<double>(&this->Phase,0));

  (*pm)["amp"] = WParameter(new ParameterValue<double>(&this->amp,1));
  (*pm)["eps"] = WParameter(new ParameterValue<double>(&this->eps,1));
  
  return pm;
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
  
  DX = Globals::instance().gridDX();
  DY = Globals::instance().gridDY();
  DZ = Globals::instance().gridDZ();
  DT = Globals::instance().dt();


  GridIndex gridLow = Globals::instance().gridLow();
  GridIndex gridHigh = Globals::instance().gridHigh();
    
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
    ::initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz)
{
}

void ShortPulseInjectSourceFunc::setTime(int time)
{
  // Now calculating on the fly!
  return;
}

Vector ShortPulseInjectSourceFunc::getEField(int i, int j, int k, int time)
{


  double ex=0, ey=0;
  double posxo = (i-centrex)*DX;
  double posxh = (i+0.5-centrex)*DX;
  double posyo = (j-centrey)*DY;
  double posyh = (j+0.5-centrey)*DY;
  double poszo = (k-centrez)*DZ - ZShift;
  double posTime = time*DT - TShift;

  Complex Exc = Efunc(posxh, posyo, poszo, posTime);
  ex = amp*Exc.real();
    
  if (YComp != Complex(0,0))
  {
    Complex Eyc = YComp*Efunc(posxo, posyh, poszo, posTime);
    ey = amp*Eyc.real();
  }
 
  return Vector(ex,ey,0);

}

Vector ShortPulseInjectSourceFunc::getHField(int i, int j, int k, int time)
{

  double bx=0, by=0;
  double posxo = (i-centrex)*DX;
  double posxh = (i+0.5-centrex)*DX;
  double posyo = (j-centrey)*DY;
  double posyh = (j+0.5-centrey)*DY;
  
//  double poszh = (k+0.5-centrez)*DZ - ZShift;
  double poszh = (k+0.5-centrez)*DZ - ZShift;
//  double posTime = (time+0.5)*DT - TShift;
  double posTime = (time-0.5)*DT - TShift;

  Complex Bxc = Bfunc(posxo, posyh, poszh, posTime, true);
//  Complex Bxc = Bfunc(posxh, posyo, poszh, posTime, true);
  bx = amp*Bxc.real();
  
  Complex Byc = Bfunc(posxh, posyo, poszh, posTime, false);
//  Complex Byc = Bfunc(posxo, posyh, poszh, posTime, false);
  by = amp*Byc.real();

  return Vector(bx,by,0);

}


