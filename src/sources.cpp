#include "sources.h"

bool PlaneWaveSource::needCurrent(Direction dir_)
{
  return (dir_ == down);
}


IncidentSourceCurrent *PlaneWaveSource::makeECurrent(int distance_, Direction dir_)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);
  
  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  
  E[0] *= factor/sqrt(eps);
  E[1] *= factor/sqrt(eps);
  E[2] *= factor/sqrt(eps);
  
  typedef IncidentSourceECurrent<PlaneWaveSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(k, E, H, ramp, eps);
  return cur;
}

IncidentSourceCurrent *PlaneWaveSource::makeHCurrent(int distance_, Direction dir_)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);
  
  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  
  E[0] *= factor/sqrt(eps);
  E[1] *= factor/sqrt(eps);
  E[2] *= factor/sqrt(eps);
    
  typedef IncidentSourceHCurrent<PlaneWaveSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(k, E, H, ramp, eps);
  return cur;
}

ParameterMap* PlaneWaveSource::MakeParamMap (ParameterMap* pm)
{
  pm = IncidentSource::MakeParamMap(pm);
  
  (*pm)["kx"] = WParameter(new ParameterValue<double>(&this->kx,0));
  (*pm)["ky"] = WParameter(new ParameterValue<double>(&this->ky,0));
  (*pm)["kz"] = WParameter(new ParameterValue<double>(&this->kz,1));
  
  (*pm)["Hx"] = WParameter(new ParameterValue<double>(&this->Hx,1));
  (*pm)["Hy"] = WParameter(new ParameterValue<double>(&this->Hy,0));
  (*pm)["Hz"] = WParameter(new ParameterValue<double>(&this->Hz,0));
  
  (*pm)["ramp"] = WParameter(new ParameterValue<double>(&this->ramp,0.5));
  (*pm)["eps"] = WParameter(new ParameterValue<double>(&this->eps,1));
  
  return pm;
}


PlaneWaveSourceEFunc::PlaneWaveSourceEFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void PlaneWaveSourceEFunc::setParam(Vector k_, Vector E_, Vector H_, double ramp_, double eps_)
{
  k = k_;
  E = E_;
  H = H_;
  ramp = ramp_;
  eps = eps_;
  dt = Globals::instance().dt();
  om = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = Globals::instance().gridDX();
  dy = Globals::instance().gridDY();
  dz = Globals::instance().gridDZ();
}


Vector PlaneWaveSourceEFunc
    ::getHField(int i, int j, int l, int time)
{
  double realtime = dt*time;

  double posx = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*l*dz + om*realtime;
  double posy = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;
  double posz = k[0]*i*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
  
  double hx = H[0]*sin(posx);
  double hy = H[1]*sin(posy);
  double hz = H[2]*sin(posz);
  
  if (posx < ramp) hx *= posx/ramp;
  if (posy < ramp) hy *= posy/ramp;
  if (posz < ramp) hz *= posz/ramp;
 
  return Vector(hx, hy, hz);
}


PlaneWaveSourceHFunc::PlaneWaveSourceHFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void PlaneWaveSourceHFunc::setParam(Vector k_, Vector E_, Vector H_, double ramp_, double eps_)
{
  k = k_;
  E = E_;
  H = H_;
  ramp = ramp_;
  eps = eps_;
  dt = Globals::instance().dt();
  om = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = Globals::instance().gridDX();
  dy = Globals::instance().gridDY();
  dz = Globals::instance().gridDZ();
}


Vector PlaneWaveSourceHFunc
    ::getEField(int i, int j, int l, int time)
{
  double realtime = dt*time;

  double posx = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*(l+0.5)*dz + om*realtime;
  double posy = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
  double posz = k[0]*(i+0.5)*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;
  
  double ex = E[0]*sin(posx);
  double ey = E[1]*sin(posy);
  double ez = E[2]*sin(posz);
  
  if (posx < ramp) ex *= posx/ramp;
  if (posy < ramp) ey *= posy/ramp;
  if (posz < ramp) ez *= posz/ramp;
 
  return Vector(ex, ey, ez);
}

