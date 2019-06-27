#include "sources.hpp"

#include <cmath>

//===============================================================
//==========  Plane Wave
//===============================================================

inline double applyPlaneWaveField(double pos, double ramp, double F, double F_bg) {
  if (pos>0) return 0.0;
  double f = F*sin(pos) + F_bg;
  return (pos > -ramp) ? -pos/ramp*f : f;
}

pCurrent PlaneWaveSource::makeECurrent(int distance_, Direction dir_)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);
  Vector H_bg(Hx_bg,Hy_bg,Hz_bg);
  Vector E_bg(Ex_bg,Ey_bg,Ez_bg);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);
  
  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  
  E[0] *= clight*factor/sqrt(eps);
  E[1] *= clight*factor/sqrt(eps);
  E[2] *= clight*factor/sqrt(eps);
  
  typedef IncidentSourceECurrent<PlaneWaveSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_, dir_);
  cur->setParam(k, E, H, E_bg, H_bg, ramp, eps, front);
  return pCurrent(cur);
}

pCurrent PlaneWaveSource::makeHCurrent(int distance_, Direction dir_)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);
  Vector H_bg(Hx_bg,Hy_bg,Hz_bg);
  Vector E_bg(Ex_bg,Ey_bg,Ez_bg);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);
  
  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  
  E[0] *= clight*factor/sqrt(eps);
  E[1] *= clight*factor/sqrt(eps);
  E[2] *= clight*factor/sqrt(eps);
    
  typedef IncidentSourceHCurrent<PlaneWaveSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(k, E, H, E_bg, H_bg, ramp, eps, front);
  return pCurrent(cur);
}

void PlaneWaveSource::initParameters(schnek::BlockParameters &blockPars)
{
  IncidentSource::initParameters(blockPars);

  blockPars.addParameter("kx", &this->kx, 0.0);
  blockPars.addParameter("ky", &this->ky, 0.0);
  blockPars.addParameter("kz", &this->kz, 1.0);

  blockPars.addParameter("Hx", &this->Hx, 0.0);
  blockPars.addParameter("Hy", &this->Hy, 0.0);
  blockPars.addParameter("Hz", &this->Hz, 0.0);

  blockPars.addParameter("Hx_bg", &this->Hx_bg, 0.0);
  blockPars.addParameter("Hy_bg", &this->Hy_bg, 0.0);
  blockPars.addParameter("Hz_bg", &this->Hz_bg, 0.0);

  blockPars.addParameter("Ex_bg", &this->Ex_bg, 0.0);
  blockPars.addParameter("Ey_bg", &this->Ey_bg, 0.0);
  blockPars.addParameter("Ez_bg", &this->Ez_bg, 0.0);

  blockPars.addParameter("ramp", &this->ramp, 0.5);
  blockPars.addParameter("eps", &this->eps, 1.0);
  blockPars.addArrayParameter("front", this->front, Vector(0,0,0));
}


PlaneWaveSourceEFunc::PlaneWaveSourceEFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void PlaneWaveSourceEFunc::setParam(Vector k_, Vector E_, Vector H_, Vector E_bg_, Vector H_bg_, double ramp_, double eps_, const Vector &front_)
{
  k = k_;
  E = E_;
  H = H_;
  E_bg = E_bg_;
  H_bg = H_bg_;
  ramp = ramp_;
  eps = eps_;
  front = front_;
  dt = MPulse::getDt();
  om = clight*sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = MPulse::getDx()[0];
  dy = MPulse::getDx()[1];
  dz = MPulse::getDx()[2];
}


Vector PlaneWaveSourceEFunc::getHField(int i, int j, int l, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx - front[0];
  double y = j*dy - front[1];
  double z = l*dz - front[2];

  double posx = k[0]*x + k[1]*(y+0.5*dy) + k[2]*(z+0.5*dz) - om*realtime;
  double posy = k[0]*(x+0.5*dx) + k[1]*y + k[2]*(z+0.5*dz) - om*realtime;
  double posz = k[0]*(x+0.5*dx) + k[1]*(y+0.5*dy) + k[2]*z - om*realtime;
  
  double hx = applyPlaneWaveField(posx, ramp, H[0], H_bg[0]);
  double hy = applyPlaneWaveField(posy, ramp, H[1], H_bg[1]);
  double hz = applyPlaneWaveField(posz, ramp, H[2], H_bg[2]);
 
//  if ((j==50) && (l==50)) {
//    std::cout << "H " << time << " "  << om << " " << posx << " " << front[0] << " | " << i << " " << hx << " " << hy << " " << hz << std::endl;
//  }

  return Vector(hx, hy, hz);
}


PlaneWaveSourceHFunc::PlaneWaveSourceHFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void PlaneWaveSourceHFunc::setParam(Vector k_, Vector E_, Vector H_, Vector E_bg_, Vector H_bg_, double ramp_, double eps_, const Vector &front_)
{
  k = k_;
  E = E_;
  H = H_;
  E_bg = E_bg_;
  H_bg = H_bg_;
  ramp = ramp_;
  eps = eps_;
  front = front_;
  dt = MPulse::getDt();
  om = clight*sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = MPulse::getDx()[0];
  dy = MPulse::getDx()[1];
  dz = MPulse::getDx()[2];
}


Vector PlaneWaveSourceHFunc::getEField(int i, int j, int l, double time)
{
//  double realtime = dt*time;
//  double posx = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*(l+0.5)*dz + om*realtime;
//  double posy = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
//  double posz = k[0]*(i+0.5)*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;

  double realtime = time;

  double x = i*dx - front[0];
  double y = j*dy - front[1];
  double z = l*dz - front[2];

  double posx = k[0]*(x+0.5*dx) + k[1]*y + k[2]*z - om*realtime;
  double posy = k[0]*x + k[1]*(y+0.5*dy) + k[2]*z - om*realtime;
  double posz = k[0]*x + k[1]*y + k[2]*(z+0.5*dz) - om*realtime;
  
  double ex = applyPlaneWaveField(posx, ramp, E[0], E_bg[0]);
  double ey = applyPlaneWaveField(posy, ramp, E[1], E_bg[1]);
  double ez = applyPlaneWaveField(posz, ramp, E[2], E_bg[2]);
 
//  if ((j==50) && (l==50)) {
//    std::cout << "E " << time << " "  << om << " " << posx << " " << ramp << " | " << i << " " << ex << " " << ey << " " << ez << std::endl;
//  }

  return Vector(ex, ey, ez);
}

//===============================================================
//==========  Plane Gaussian Wave Packet
//===============================================================

pCurrent PlaneGaussSource::makeECurrent(int distance_, Direction dir_)
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
  return pCurrent(cur);
}

pCurrent PlaneGaussSource::makeHCurrent(int distance_, Direction dir_)
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
  return pCurrent(cur);
}

void PlaneGaussSource::initParameters(schnek::BlockParameters &blockPars)
{

  IncidentSource::initParameters(blockPars);

  blockPars.addParameter("kx", &this->kx, 0.0);
  blockPars.addParameter("ky", &this->ky, 0.0);
  blockPars.addParameter("kz", &this->kz, 1.0);

  blockPars.addParameter("Hx", &this->Hx, 0.0);
  blockPars.addParameter("Hy", &this->Hy, 0.0);
  blockPars.addParameter("Hz", &this->Hz, 0.0);

  blockPars.addParameter("width", &this->width, 10.0);
  blockPars.addParameter("offset", &this->offset, 40.0);
  blockPars.addParameter("eps", &this->eps, 1.0);
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
  dt = MPulse::getDt();
  om = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = MPulse::getDx()[0];
  dy = MPulse::getDx()[1];
  dz = MPulse::getDx()[2];
}


Vector PlaneGaussSourceEFunc::getHField(int i, int j, int l, double time)
{
  double realtime = time-0.5*dt;

  double posx = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*(l+0.5)*dz + om*realtime;
  double posy = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
  double posz = k[0]*(i+0.5)*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;
  
  double ampx = H[0]*exp(-std::pow( (posx-offset)/width, 2) );
  double ampy = H[1]*exp(-std::pow( (posy-offset)/width, 2));
  double ampz = H[2]*exp(-std::pow( (posz-offset)/width, 2));
  
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
  dt = MPulse::getDt();
  om = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = MPulse::getDx()[0];
  dy = MPulse::getDx()[1];
  dz = MPulse::getDx()[2];
}


Vector PlaneGaussSourceHFunc::getEField(int i, int j, int l, double time)
{
//  double realtime = dt*time;
//  double posx = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*(l+0.5)*dz + om*realtime;
//  double posy = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
//  double posz = k[0]*(i+0.5)*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;

  double realtime = time;

  double posx = k[0]*(i+0.5)*dx + k[1]*j*dy + k[2]*l*dz + om*realtime;
  double posy = k[0]*i*dx + k[1]*(j+0.5)*dy + k[2]*l*dz + om*realtime;
  double posz = k[0]*i*dx + k[1]*j*dy + k[2]*(l+0.5)*dz + om*realtime;
  
  double ampx = E[0]*exp(-std::pow( (posx-offset)/width, 2));
  double ampy = E[1]*exp(-std::pow( (posy-offset)/width, 2));
  double ampz = E[2]*exp(-std::pow( (posz-offset)/width, 2));

  double ex = ampx*sin(posx);
  double ey = ampy*sin(posy);
  double ez = ampz*sin(posz);
 
  return Vector(ex, ey, ez);
}

