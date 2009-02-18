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
  
  cur->setParam(length, width, om0, Shift, Phase, amp, eps, distance_);
  return cur;
}

IncidentSourceCurrent *ShortPulseInject::makeHCurrent(int distance_, Direction dir_)
{
  std::cerr << "ShortPulseInjectSourceFunc::makeHCurrent\n";
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
  std::cerr << "ShortPulseInjectSourceFunc::MakeParamMap\n";
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
  std::cerr << "ShortPulseInjectSourceFunc::setParam\n";
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
  YComp = Complex(0.0,0.0);

  GridIndex gridLow = Globals::instance().gridLow();
  GridIndex gridHigh = Globals::instance().gridHigh();
    
  centrex = 0.5*double(gridHigh[0] + gridLow[0]);
  centrey = 0.5*double(gridHigh[1] + gridLow[1]);
  centrez = 0.5*double(gridHigh[2] + gridLow[2]);
}

void ShortPulseInjectSourceFunc
    ::initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz)
{
  std::cerr << "ShortPulseInjectSourceFunc::initSourceFunc\n";
  PsiX  = storage->addBorderLayer("IncidentPsiX" , dir, 4, dist-2, 1);
  PsiY  = storage->addBorderLayer("IncidentPsiY" , dir, 4, dist-2, 1);
  PsiXp = storage->addBorderLayer("IncidentPsiXp", dir, 4, dist-2, 1);
  PsiYp = storage->addBorderLayer("IncidentPsiYp", dir, 4, dist-2, 1);
  std::cerr << "ShortPulseInjectSourceFunc::initSourceFunc ... done\n";
}

ShortPulseInjectSourceFunc::Complex 
  ShortPulseInjectSourceFunc::calcPsi(double t, double x, double y, double z)
{
//  return sin(om0*t-z/lightspeed);

  double r2 = x*x+y*y;
  Complex I(0,1);
  Complex q(z,ZRl);
  Complex tp = t - z/lightspeed - r2/(2.*q*lightspeed);
//  Complex tp = t - z/lightspeed;
  
  double tau = length;
  double taubar = 0.5*om0*tau;
  
  return I*(ZRl/q)*exp(I*om0*tp);
  
  Complex iphi(0,Phase);
  
  Complex tptau2 = tp*tp/(tau*tau);
  Complex tautbtp2 = taubar*taubar/tptau2;

// ===============================
// Version 01
  
//  Complex an = tp/tau - I*tau*taubar/tp;
//  Complex ap = tp/tau + I*tau*taubar/tp;
//  
//  Complex e1 = 0.5*exp(-iphi - tptau2 - 2.*I*taubar);
//  Complex e2 = exp(-ap*ap);
//  Complex e3 = exp(2.*iphi + 4.*I*taubar);
//  Complex e4 = exp(-an*an);
  
//  Complex wp = Faddeeva_2(I*tp/tau + tau*taubar/tp);
//  Complex wn = Faddeeva_2(I*tp/tau - tau*taubar/tp);

//  Complex psi = I*(ZRl/q) * e1*(2. - e2*wn + e3*(2.-e4*wp));

// ===============================
// Version 02
  
//  Complex e1 = 0.5*exp(-iphi-2.*tptau2 - tautbtp2);  
//  Complex e2 = exp(2.*iphi + 2.*tptau2);
//  Complex e3 = exp(2.*taubar*(2.*I+taubar/tptau2));

//  Complex wp = Faddeeva_2(I*tp/tau +   tau*taubar/tp);
//  Complex wn = Faddeeva_2(- tp/tau + I*tau*taubar/tp);
 
//  Complex psi = I*(ZRl/q) * e1*(e2*wn - e3*wp);




//  std::cerr << ZRl/q << " " << real(e1) << " " << real(e2) << " " << real(e3) << " "
//    << real(e4) << " " << real(e5) << " " << real(wp) << " "
//    << real(wn)
//    << "\n";
//  return psi;
  
  Complex tpt = tp/tau;
  Complex itb = I*taubar;
  Complex ap = tpt + itb;
  Complex an = tpt + itb;
  
  Complex C = exp(iphi);
  Complex D = exp(-iphi);
  Complex S(0,0);
  
  if (imag(ap)<0)
  {
    S = S + amp*C*exp(-tpt*tpt-2.*itb*tpt);
    C = -C;
    ap = -ap;
  }
  
  if (imag(an)<0)
  {
    S = S + amp*C*exp(-tpt*tpt+2.*itb*tpt);
    D = -D;
    an = -an;
  }
  
  Complex wp = Faddeeva_2(ap);
  Complex wn = Faddeeva_2(an);
  
  Complex psi = S + 0.5*amp*exp(-taubar*taubar)*(C*wn + D*wp);
  
  return ZRl*psi/q;
}

void ShortPulseInjectSourceFunc::setTime(int time)
{
  // Now calculating on the fly!
  return;
  
  // Only needs to be updated if calculating the electric field
  // This is calculated in the HCurrent
  // Note: All currents share the same Psi fields
  
  if (!isH || (time<=old_time)) return;
  
  old_time = time;
  
  double t = time*DT-Shift/lightspeed;
  
  DataGrid &psix  = *PsiX,  &psiy  = *PsiY;
  DataGrid &psixp = *PsiXp, &psiyp = *PsiYp;
   
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
    double poszo = (k-mz)*DZ;
    double poszh = (k+0.5-mz)*DZ;
    
    for (int i=lowx; i<=highx; ++i)
    {
      double posxo = (i-mx)*DX;
      double posxh = (i+0.5-mx)*DX;

      for (int j=lowy; j<=highy; ++j)
      {
          // copy old values first
          
          psixp(i,j,k) = psix(i,j,k);
          psiyp(i,j,k) = psiy(i,j,k);
          
          // calculate new values
          
          double posyo = (j-my)*DY;        
          double posyh = (j+0.5-my)*DY;
          
          // This is the way it should be
          Complex px = calcPsi(t, posxh, posyo, poszh);
          Complex py = YComp*calcPsi(t, posxo, posyh, poszh);
          
          // This is a trial version (on dual grid)
          // Complex px = calcPsi(t, posxo, posyh, poszo);
          // Complex py = YComp*calcPsi(t, posxh, posyo, poszo);
          
          psix(i,j,k) = px.real();
          psiy(i,j,k) = py.real();
      }
    }
  }  
}

Vector ShortPulseInjectSourceFunc::getEField(int i, int j, int k, int time)
{
/*
  double z = k*DZ;
  double t = time*DT;

  return Vector(cos(om0*t-z/lightspeed),0,0);
*/


  double ex=0, ey=0;
  double posxo = (i-centrex)*DX;
  double posxh = (i+0.5-centrex)*DX;
  double posyo = (j-centrey)*DY;
  double posyh = (j+0.5-centrey)*DY;
  double poszo = (k-centrez)*DZ;
  double posTime = time*DT;

  Complex Exc = Efunc(posxh, posyo, poszo, posTime);
  ex = Exc.real()/M_PI;
  
//  std::cerr << i << " "  << j << " "  << k << " " << ex << "\n";
  
  if (YComp != Complex(0,0))
  {
    Complex Eyc = YComp*Efunc(posxo, posyh, poszo, posTime);
    ey = Eyc.real()/M_PI;
  }

//override
/*  double c = cos(0.3), s = sin(0.3);
  double kz = c*om0/lightspeed;
  double kx = s*om0/lightspeed;
  double ky = 0;
  double amp = sin(kx*posxh + ky*posyo + kz*poszo - om0*posTime)/lightspeed;
  
  ex = amp; //*c;
  ey = 0;
*/  
  return Vector(ex,ey,0);

/*
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
*/
/*
  double xo = (k-centrex)*DX;
  double xh = (k+0.5-centrex)*DX;
  double yo = (k-centrey)*DY;
  double yh = (k+0.5-centrey)*DY;
  double zo = (k-centrez)*DZ;
  double t  = time*DT;

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
*/
}

Vector ShortPulseInjectSourceFunc::getHField(int i, int j, int k, int time)
{
/*
  double z = (k-0.5)*DZ;
  double t = time*DT;

  return Vector(0,cos(om0*t-z/lightspeed),0);
*/

  double bx=0, by=0;
  double posxo = (i-centrex)*DX;
  double posxh = (i+0.5-centrex)*DX;
  double posyo = (j-centrey)*DY;
  double posyh = (j+0.5-centrey)*DY;
  
  double poszh = (k+0.5-centrez)*DZ;
  double posTime = (time+0.5)*DT;

  Complex Bxc = Bfunc(posxo, posyh, poszh, posTime, true);
//  Complex Bxc = Bfunc(posxh, posyo, poszh, posTime, true);
  bx = Bxc.real()/M_PI;
  
  Complex Byc = Bfunc(posxh, posyo, poszh, posTime, false);
//  Complex Byc = Bfunc(posxo, posyh, poszh, posTime, false);
  by = Byc.real()/M_PI;

//override
/*  double c = cos(0.3), s = sin(0.3);
  double kz = c*om0/lightspeed;
  double kx = s*om0/lightspeed;
  double ky = 0;
  double amp = sin(kx*posxh + ky*posyo + kz*poszh - om0*posTime);
  bx = 0;
  by = amp;
*/
  return Vector(bx,by,0);


/*
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
*/
}


