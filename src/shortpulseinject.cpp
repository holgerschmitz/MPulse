#include "shortpulseinject.h"
#include <cmath>
#include "specfunc.h"

//===============================================================
//==========  ShortPulseInject
//===============================================================

IncidentSourceCurrent *ShortPulseInject::makeECurrent(int distance_, Direction dir_)
{
  typedef IncidentSourceECurrent<ShortPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  
  cur->setParam(length, width, om0, Shift, Phase, amp, eps, distance_);
  return cur;
}

IncidentSourceCurrent *ShortPulseInject::makeHCurrent(int distance_, Direction dir_)
{
  typedef IncidentSourceHCurrent<ShortPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(length, width, om0, Shift, Phase, amp, eps, distance_);
  return cur;
}

bool ShortPulseInject::needCurrent(Direction dir_)
{
  return (dir_ == down);
}
    
ParameterMap* ShortPulseInject::MakeParamMap (ParameterMap* pm)
{
  pm = IncidentSource::MakeParamMap(pm);
  
  (*pm)["length"] = WParameter(new ParameterValue<double>(&this->length,1.));
  (*pm)["width"] = WParameter(new ParameterValue<double>(&this->width,1.));
  (*pm)["om0"] = WParameter(new ParameterValue<double>(&this->om0,2*M_PI));
  (*pm)["Shift"] = WParameter(new ParameterValue<double>(&this->Shift,0));
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
                                          double Shift_,
                                          double Phase_,
                                          double amp_, 
                                          double eps_,
                                          int distance_)
{
  length = length_;
  width  = width_;
  om0     = om0_;
  Shift  = Shift_;
  Phase  = 2*M_PI*Phase_;
    
  amp = amp_;
  eps = eps_;
  dist = distance_;
  
  lightspeed = sqrt(1/eps);
  ZRl = 0.5*om0*width*width/lightspeed;
  
  
  DX = Globals::instance().gridDX();
  DY = Globals::instance().gridDY();
  DZ = Globals::instance().gridDZ();
  DT = Globals::instance().dt() * lightspeed;

  old_time = -1;
}

void ShortPulseInjectSourceFunc
    ::initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz)
{
  PsiX  = storage->addBorderLayer("IncidentPsiX" , dir, 2, dist, 1);
  PsiY  = storage->addBorderLayer("IncidentPsiY" , dir, 2, dist, 1);
  PsiXp = storage->addBorderLayer("IncidentPsiXp", dir, 2, dist, 1);
  PsiYp = storage->addBorderLayer("IncidentPsiYp", dir, 2, dist, 1);
}

ShortPulseInjectSourceFunc::Complex 
  ShortPulseInjectSourceFunc::calcPsi(double t, double x, double y, double z)
{
  double r2 = x*x+y*y;
  Complex I(0,1);
  Complex q(z,ZRl);
  Complex tp = t - z/lightspeed - r2/(2.*q*lightspeed);
  
  double taubar = 0.5*om0*t;
  double tau = length;
  
  Complex iphi(0,Phase);
    
  Complex wp = Faddeeva_2(tp/tau + I*taubar);
  Complex wn = Faddeeva_2(tp/tau - I*taubar);
  
  Complex psi = 0.5*amp*exp(-taubar*taubar)*(exp(iphi)*wn + exp(-iphi)*wp);
  
  return ZRl*psi/q;
}

void ShortPulseInjectSourceFunc::setTime(int time)
{
  // Only needs to be updated if calculating the electric field
  // This is calculated in the HCurrent
  // Note: All currents share the same Psi fields
  
  if (!isH || (time<=old_time)) return;
  
  old_time = time;
  
  double t = time*DT-Shift/lightspeed;
  // fast copy of the old field
  
  DataGrid &psix  = *PsiX,  &psiy  = *PsiY;
  DataGrid &psixp = *PsiXp, &psiyp = *PsiYp;

  typedef DataGrid::const_storage_iterator CIterator;
  typedef DataGrid::storage_iterator Iterator;
  
  CIterator srcX    = psix.cbegin();
  CIterator srcXEnd = psix.cend();
  Iterator  destX   = psixp.begin();

  CIterator srcY    = psiy.cbegin();
  CIterator srcYEnd = psiy.cend();
  Iterator  destY   = psiyp.begin();
  while (srcX != srcXEnd)
  {
    *destX = *srcX;    ++srcX;    ++destX;
    *destY = *srcY;    ++srcY;    ++destY;
  }
  
  // calculate the new field
  
  GridIndex low  = psix.getLow();
  GridIndex high = psix.getHigh();
  
  int lowx=low[0],   lowy=low[1], lowz=low[2];
  int highx=high[0], highy=high[1], highz=high[2];
  
  GridIndex gridLow = Globals::instance().gridLow();
  GridIndex gridHigh = Globals::instance().gridHigh();
    
  double mx = 0.5*double(gridHigh[0] + gridLow[0]);
  double my = 0.5*double(gridHigh[1] + gridLow[1]);
  double mz = 0.5*double(gridHigh[2] + gridLow[2]);
    
  for (int k=lowz; k<=highz; ++k)
  {
    double poszh = (k+0.5-mz)*DZ;
    
    for (int i=lowx; i<=highx; ++i)
    {
      double posxo = (i-mx)*DX;
      double posxh = (i+0.5-mx)*DX;

      for (int j=lowy; j<=highy; ++j)
      {
          double posyo = (j-my)*DY;        
          double posyh = (j+0.5-my)*DY;
          
          Complex px = calcPsi(t, posxh, posyo, poszh);
          Complex py = YComp*calcPsi(t, posxo, posyh, poszh);
          
          psix(i,j,k) = px.real();
          psiy(i,j,k) = py.real();
      }
    }
  }
  
}

Vector ShortPulseInjectSourceFunc::getEField(int i, int j, int k, int time)
{
  double ex, ey, ez;

  DataGrid &psix  = *PsiX,  &psiy  = *PsiY;
  DataGrid &psixp = *PsiXp, &psiyp = *PsiYp;

  double ax  =   (psix(i,j,k)-psix(i,j,k-1))/DZ;
  double ay  =   (psiy(i,j,k)-psiy(i,j,k-1))/DZ;
  double az  = - (psix(i,j,k)-psix(i-1,j,k))/DX - (psiy(i,j,k)-psiy(i,j-1,k))/DY;
  double axp =   (psixp(i,j,k)-psixp(i,j,k-1))/DZ;
  double ayp =   (psiyp(i,j,k)-psiyp(i,j,k-1))/DZ;
  double azp = - (psixp(i,j,k)-psixp(i-1,j,k))/DX - (psiyp(i,j,k)-psiyp(i,j-1,k))/DY;

  ex = -(ax-axp)/(DT*lightspeed);
  ey = -(ay-ayp)/(DT*lightspeed);
  ez = -(az-azp)/(DT*lightspeed);

  return Vector(ex,ey,ez);
}

Vector ShortPulseInjectSourceFunc::getHField(int i, int j, int k, int time)
{
  double bx, by, bz;

  DataGrid &psix  = *PsiX,  &psiy  = *PsiY;

  double ax  =   (psix(i,j,k)-psix(i,j,k-1))/DZ;
  double ay  =   (psiy(i,j,k)-psiy(i,j,k-1))/DZ;
  double az  = - (psix(i,j,k)-psix(i-1,j,k))/DX - (psiy(i,j,k)-psiy(i,j-1,k))/DY;

  double aypx  =   (psiy(i+1,j,k)-psiy(i+1,j,k-1))/DZ;
  double azpx  = - (psix(i+1,j,k)-psix(i,j,k))/DX - (psiy(i+1,j,k)-psiy(i+1,j-1,k))/DY;
  
  double axpy  =   (psix(i,j+1,k)-psix(i,j+1,k-1))/DZ;
  double azpy  = - (psix(i,j+1,k)-psix(i-1,j+1,k))/DX - (psiy(i,j+1,k)-psiy(i,j,k))/DY;
  
  double axpz  =   (psix(i,j,k+1)-psix(i,j,k))/DZ;
  double aypz  =   (psiy(i,j,k+1)-psiy(i,j,k))/DZ;

  bx = (azpy-az)/DY-(aypz-ay)/DZ;
  by = (axpz-ax)/DZ-(azpx-az)/DX;
  bz = (aypx-ay)/DX-(axpy-ax)/DY;
  
  return Vector(bx, by, bz);
}


