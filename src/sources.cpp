#include "sources.h"
#include "util.h"

//===============================================================
//==========  Plane Wave
//===============================================================

bool PlaneWaveSource::needCurrent(Direction dir_)
{
  return (dir_ == down);
}


IncidentSourceCurrent *PlaneWaveSource::makeECurrent(int distance_, Direction dir_)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);
  Vector H_bg(Hx_bg,Hy_bg,Hz_bg);
  Vector E_bg(Ex_bg,Ey_bg,Ez_bg);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);
  
  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  
  E[0] *= factor/sqrt(eps);
  E[1] *= factor/sqrt(eps);
  E[2] *= factor/sqrt(eps);
  
  typedef IncidentSourceECurrent<PlaneWaveSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(k, E, H, E_bg, H_bg, ramp, eps);
  return cur;
}

IncidentSourceCurrent *PlaneWaveSource::makeHCurrent(int distance_, Direction dir_)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);
  Vector H_bg(Hx_bg,Hy_bg,Hz_bg);
  Vector E_bg(Ex_bg,Ey_bg,Ez_bg);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);
  
  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  
  E[0] *= factor/sqrt(eps);
  E[1] *= factor/sqrt(eps);
  E[2] *= factor/sqrt(eps);
    
  typedef IncidentSourceHCurrent<PlaneWaveSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(k, E, H, E_bg, H_bg, ramp, eps);
  return cur;
}

ParameterMap* PlaneWaveSource::MakeParamMap (ParameterMap* pm)
{
  pm = IncidentSource::MakeParamMap(pm);
  
  (*pm)["kx"] = WParameter(new ParameterValue<double>(&this->kx,0));
  (*pm)["ky"] = WParameter(new ParameterValue<double>(&this->ky,0));
  (*pm)["kz"] = WParameter(new ParameterValue<double>(&this->kz,1));
  
  (*pm)["Hx"] = WParameter(new ParameterValue<double>(&this->Hx,0));
  (*pm)["Hy"] = WParameter(new ParameterValue<double>(&this->Hy,0));
  (*pm)["Hz"] = WParameter(new ParameterValue<double>(&this->Hz,0));
  
  (*pm)["Hx_bg"] = WParameter(new ParameterValue<double>(&this->Hx_bg,0));
  (*pm)["Hy_bg"] = WParameter(new ParameterValue<double>(&this->Hy_bg,0));
  (*pm)["Hz_bg"] = WParameter(new ParameterValue<double>(&this->Hz_bg,0));
  
  (*pm)["Ex_bg"] = WParameter(new ParameterValue<double>(&this->Ex_bg,0));
  (*pm)["Ey_bg"] = WParameter(new ParameterValue<double>(&this->Ey_bg,0));
  (*pm)["Ez_bg"] = WParameter(new ParameterValue<double>(&this->Ez_bg,0));
  
  (*pm)["ramp"] = WParameter(new ParameterValue<double>(&this->ramp,0.5));
  (*pm)["eps"] = WParameter(new ParameterValue<double>(&this->eps,1));
  
  return pm;
}


PlaneWaveSourceEFunc::PlaneWaveSourceEFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void PlaneWaveSourceEFunc::setParam(Vector k_, Vector E_, Vector H_, Vector E_bg_, Vector H_bg_, double ramp_, double eps_)
{
  k = k_;
  E = E_;
  H = H_;
  E_bg = E_bg_;
  H_bg = H_bg_;
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
  double realtime = dt*(time-0.5);

  double posx = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*(l+0.5)*dz + om*realtime;
  double posy = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
  double posz = k[0]*(i+0.5)*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;
  
  double hx = H[0]*sin(posx) + H_bg[0];
  double hy = H[1]*sin(posy) + H_bg[1];
  double hz = H[2]*sin(posz) + H_bg[2];
  
  if (posx < ramp) hx *= posx/ramp;
  if (posy < ramp) hy *= posy/ramp;
  if (posz < ramp) hz *= posz/ramp;
 
  return Vector(hx, hy, hz);
}


PlaneWaveSourceHFunc::PlaneWaveSourceHFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void PlaneWaveSourceHFunc::setParam(Vector k_, Vector E_, Vector H_, Vector E_bg_, Vector H_bg_, double ramp_, double eps_)
{
  k = k_;
  E = E_;
  H = H_;
  E_bg = E_bg_;
  H_bg = H_bg_;
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
//  double realtime = dt*time;
//  double posx = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*(l+0.5)*dz + om*realtime;
//  double posy = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
//  double posz = k[0]*(i+0.5)*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;

  double realtime = dt*time;

  double posx = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*l*dz + om*realtime;
  double posy = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;
  double posz = k[0]*i*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
  
  double ex = E[0]*sin(posx) + E_bg[0];
  double ey = E[1]*sin(posy) + E_bg[1];
  double ez = E[2]*sin(posz) + E_bg[2];
  
  if (posx < ramp) ex *= posx/ramp;
  if (posy < ramp) ey *= posy/ramp;
  if (posz < ramp) ez *= posz/ramp;
 
  return Vector(ex, ey, ez);
}

//===============================================================
//==========  Plane Gaussian Wave Packet
//===============================================================

bool PlaneGaussSource::needCurrent(Direction dir_)
{
  return (dir_ == down);
}


IncidentSourceCurrent *PlaneGaussSource::makeECurrent(int distance_, Direction dir_)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);
  
  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  
  E[0] *= factor/sqrt(eps);
  E[1] *= factor/sqrt(eps);
  E[2] *= factor/sqrt(eps);
  
  typedef IncidentSourceECurrent<PlaneGaussSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(k, E, H, width, offset, eps);
  return cur;
}

IncidentSourceCurrent *PlaneGaussSource::makeHCurrent(int distance_, Direction dir_)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);
  
  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  
  E[0] *= factor/sqrt(eps);
  E[1] *= factor/sqrt(eps);
  E[2] *= factor/sqrt(eps);
    
  typedef IncidentSourceHCurrent<PlaneGaussSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(k, E, H, width, offset, eps);
  return cur;
}

ParameterMap* PlaneGaussSource::MakeParamMap (ParameterMap* pm)
{
  pm = IncidentSource::MakeParamMap(pm);
  
  (*pm)["kx"] = WParameter(new ParameterValue<double>(&this->kx,0));
  (*pm)["ky"] = WParameter(new ParameterValue<double>(&this->ky,0));
  (*pm)["kz"] = WParameter(new ParameterValue<double>(&this->kz,1));
  
  (*pm)["Hx"] = WParameter(new ParameterValue<double>(&this->Hx,0));
  (*pm)["Hy"] = WParameter(new ParameterValue<double>(&this->Hy,0));
  (*pm)["Hz"] = WParameter(new ParameterValue<double>(&this->Hz,0));
  
  (*pm)["width"] = WParameter(new ParameterValue<double>(&this->width,10));
  (*pm)["offset"] = WParameter(new ParameterValue<double>(&this->offset,40));
  (*pm)["eps"] = WParameter(new ParameterValue<double>(&this->eps,1));
  
  return pm;
}


PlaneGaussSourceEFunc::PlaneGaussSourceEFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void PlaneGaussSourceEFunc::setParam(Vector k_, Vector E_, Vector H_, double width_, double offset_, double eps_)
{
  k = k_;
  E = E_;
  H = H_;
  width = width_;
  offset = offset_;
  eps = eps_;
  dt = Globals::instance().dt();
  om = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = Globals::instance().gridDX();
  dy = Globals::instance().gridDY();
  dz = Globals::instance().gridDZ();
}


Vector PlaneGaussSourceEFunc
    ::getHField(int i, int j, int l, int time)
{
  double realtime = dt*(time-0.5);

  double posx = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*(l+0.5)*dz + om*realtime;
  double posy = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
  double posz = k[0]*(i+0.5)*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;
  
  double ampx = H[0]*exp(-sqr( (posx-offset)/width ));
  double ampy = H[1]*exp(-sqr( (posy-offset)/width ));
  double ampz = H[2]*exp(-sqr( (posz-offset)/width ));
  
  double hx = ampx*sin(posx);
  double hy = ampy*sin(posy);
  double hz = ampz*sin(posz);
 
  return Vector(hx, hy, hz);
}


PlaneGaussSourceHFunc::PlaneGaussSourceHFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void PlaneGaussSourceHFunc::setParam(Vector k_, Vector E_, Vector H_, double width_, double offset_, double eps_)
{
  k = k_;
  E = E_;
  H = H_;
  width = width_;
  offset = offset_;
  eps = eps_;
  dt = Globals::instance().dt();
  om = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = Globals::instance().gridDX();
  dy = Globals::instance().gridDY();
  dz = Globals::instance().gridDZ();
}


Vector PlaneGaussSourceHFunc
    ::getEField(int i, int j, int l, int time)
{
//  double realtime = dt*time;
//  double posx = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*(l+0.5)*dz + om*realtime;
//  double posy = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
//  double posz = k[0]*(i+0.5)*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;

  double realtime = dt*time;

  double posx = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*l*dz + om*realtime;
  double posy = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;
  double posz = k[0]*i*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
  
  double ampx = E[0]*exp(-sqr( (posx-offset)/width ));
  double ampy = E[1]*exp(-sqr( (posy-offset)/width ));
  double ampz = E[2]*exp(-sqr( (posz-offset)/width ));

  double ex = ampx*sin(posx);
  double ey = ampy*sin(posy);
  double ez = ampz*sin(posz);
 
  return Vector(ex, ey, ez);
}

