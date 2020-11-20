#include "sources.hpp"

#include "../huerto/constants.hpp"

#include <cmath>

//===============================================================
//==========  Plane Wave
//===============================================================

inline double applyPlaneWaveField(double pos, double ramp, double F, double F_bg) {
  if (pos>0) return 0.0;
  double f = F*sin(pos) + F_bg;
  return (pos > -ramp) ? -pos/ramp*f : f;
}

inline double applyPlaneGaussField(double pos, double width, double F) {
  double r = pos/width;
  return F*exp(-r*r)*sin(pos);
}

pCurrent PlaneWaveSource::makeECurrent(int distance, Direction dir)
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
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setParam(k, E, H, E_bg, H_bg, ramp, eps, front);
  return pCurrent(cur);
}

pCurrent PlaneWaveSource::makeHCurrent(int distance, Direction dir)
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
  CurrentType *cur = new CurrentType(distance, dir, getContext());
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


PlaneWaveSourceEFunc::PlaneWaveSourceEFunc(Direction dir_, bool isH_, SimulationContext &context)
  : dir(dir_), isH(isH_), context(context)
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
  dt = context.getDt();
  om = clight*sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = context.getDx()[0];
  dy = context.getDx()[1];
  dz = context.getDx()[2];
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


PlaneWaveSourceHFunc::PlaneWaveSourceHFunc(Direction dir_, bool isH_, SimulationContext &context)
  : dir(dir_), isH(isH_), context(context)
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
  dt = context.getDt();
  om = clight*sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = context.getDx()[0];
  dy = context.getDx()[1];
  dz = context.getDx()[2];
}


Vector PlaneWaveSourceHFunc::getEField(int i, int j, int l, double time)
{
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

  return Vector(ex, ey, ez);
}

//===============================================================
//==========  Plane Gaussian Wave Packet
//===============================================================

pCurrent PlaneGaussSource::makeECurrent(int distance, Direction dir)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);

  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);

  E[0] *= clight*factor/sqrt(eps);
  E[1] *= clight*factor/sqrt(eps);
  E[2] *= clight*factor/sqrt(eps);

  typedef IncidentSourceECurrent<PlaneGaussSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setParam(k, E, H, width, eps, front);
  return pCurrent(cur);
}

pCurrent PlaneGaussSource::makeHCurrent(int distance, Direction dir)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);

  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);

  E[0] *= clight*factor/sqrt(eps);
  E[1] *= clight*factor/sqrt(eps);
  E[2] *= clight*factor/sqrt(eps);

  typedef IncidentSourceHCurrent<PlaneGaussSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setParam(k, E, H, width, eps, front);
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

  blockPars.addParameter("eps", &this->eps, 1.0);
  blockPars.addArrayParameter("front", this->front, Vector(0,0,0));
}


PlaneGaussSourceEFunc::PlaneGaussSourceEFunc(Direction dir_, bool isH_, SimulationContext &context)
  : dir(dir_), isH(isH_), context(context)
{}

void PlaneGaussSourceEFunc::setParam(Vector k_, Vector E_, Vector H_, double width_, double eps_, const Vector &front_)
{
  k = k_;
  E = E_;
  H = H_;
  double kn = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
  width = width_*kn;
  eps = eps_;
  front = front_;
  dt = context.getDt();
  om = clight*kn/sqrt(eps);

  dx = context.getDx()[0];
  dy = context.getDx()[1];
  dz = context.getDx()[2];
}


Vector PlaneGaussSourceEFunc::getHField(int i, int j, int l, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx - front[0];
  double y = j*dy - front[1];
  double z = l*dz - front[2];

  double posx = k[0]*x + k[1]*(y+0.5*dy) + k[2]*(z+0.5*dz) - om*realtime;
  double posy = k[0]*(x+0.5*dx) + k[1]*y + k[2]*(z+0.5*dz) - om*realtime;
  double posz = k[0]*(x+0.5*dx) + k[1]*(y+0.5*dy) + k[2]*z - om*realtime;

  double hx = applyPlaneGaussField(posx, width, H[0]);
  double hy = applyPlaneGaussField(posy, width, H[1]);
  double hz = applyPlaneGaussField(posz, width, H[2]);

  return Vector(hx, hy, hz);
}


PlaneGaussSourceHFunc::PlaneGaussSourceHFunc(Direction dir_, bool isH_, SimulationContext &context)
  : dir(dir_), isH(isH_), context(context)
{}

void PlaneGaussSourceHFunc::setParam(Vector k_, Vector E_, Vector H_, double width_, double eps_, const Vector &front_)
{
  k = k_;
  E = E_;
  H = H_;
  double kn = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
  width = width_*kn;
  eps = eps_;
  front = front_;
  dt = context.getDt();
  om = clight*kn/sqrt(eps);

  dx = context.getDx()[0];
  dy = context.getDx()[1];
  dz = context.getDx()[2];
}


Vector PlaneGaussSourceHFunc::getEField(int i, int j, int l, double time)
{
  double realtime = time;

  double x = i*dx - front[0];
  double y = j*dy - front[1];
  double z = l*dz - front[2];

  double posx = k[0]*(x+0.5*dx) + k[1]*y + k[2]*z - om*realtime;
  double posy = k[0]*x + k[1]*(y+0.5*dy) + k[2]*z - om*realtime;
  double posz = k[0]*x + k[1]*y + k[2]*(z+0.5*dz) - om*realtime;

  double ex = applyPlaneGaussField(posx, width, E[0]);
  double ey = applyPlaneGaussField(posy, width, E[1]);
  double ez = applyPlaneGaussField(posz, width, E[2]);

  return Vector(ex, ey, ez);
}

