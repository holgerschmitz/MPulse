#ifndef MPULSE_SOURCES_H
#define MPULSE_SOURCES_H

#include "incsource.hpp"

//===============================================================
//==========  Plane Wave
//===============================================================

class PlaneWaveSource : public IncidentSource
{
  public:
    ~PlaneWaveSource() {}
  protected:
    pCurrent makeECurrent(int distance_, Direction dir_);
    pCurrent makeHCurrent(int distance_, Direction dir_);
    
    void initParameters(schnek::BlockParameters &blockPars);
    
    double kx, ky, kz;
    double Hx, Hy, Hz;
    double Hx_bg, Hy_bg, Hz_bg;
    double Ex_bg, Ey_bg, Ez_bg;
    double ramp;
    double eps;
    Vector front;
};

class PlaneWaveSourceEFunc
{
  public:
    PlaneWaveSourceEFunc(Direction dir_, bool isH_);
    void setParam(Vector k_, Vector E_, Vector H_, Vector E_bg_, Vector H_bg_, double ramp_, double eps_, const Vector &front_);

    Vector getHField(int i, int j, int k, double time);

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(int) {}
    
  private:
    Vector k;
    Vector E;
    Vector H;
    Vector E_bg;
    Vector H_bg;
    
    double dt;
    double om;
    double ramp;
    double eps;
    Vector front;
    
    double dx, dy, dz;
    
    Direction dir;
    bool isH;
};

class PlaneWaveSourceHFunc
{
  public:
    PlaneWaveSourceHFunc(Direction dir_, bool isH_);
    void setParam(Vector k_, Vector E_, Vector H_, Vector E_bg_, Vector H_bg_, double ramp_, double eps_, const Vector &front_);

    Vector getEField(int i, int j, int k, double time);

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(int) {}
    
  private:
    Vector k;
    Vector E;
    Vector H;
    Vector E_bg;
    Vector H_bg;
    
    double dt;
    double om;
    double ramp;
    double eps;
    Vector front;
    
    double dx, dy, dz;
    
    Direction dir;
    bool isH;
};


//===============================================================
//==========  Plane Gauss Packet Source
//===============================================================

class PlaneGaussSource : public IncidentSource
{
  public:
    ~PlaneGaussSource() {}
  protected:
    pCurrent makeECurrent(int distance_, Direction dir_);
    pCurrent makeHCurrent(int distance_, Direction dir_);
    
    void initParameters(schnek::BlockParameters &blockPars);
    
    double kx, ky, kz;
    double Hx, Hy, Hz;
    double width;
    double offset;
    double eps;
};

class PlaneGaussSourceEFunc
{
  public:
    PlaneGaussSourceEFunc(Direction dir_, bool isH_);
    void setParam(Vector k_, Vector E_, Vector H_, double width_, double offset_, double eps_);

    Vector getHField(int i, int j, int k, double time);

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(int) {}
    
  private:
    Vector k;
    Vector E;
    Vector H;
    
    double dt;
    double om;
    double width;
    double offset;
    double eps;
    
    double dx, dy, dz;
    
    Direction dir;
    bool isH;
};

class PlaneGaussSourceHFunc
{
  public:
    PlaneGaussSourceHFunc(Direction dir_, bool isH_);
    void setParam(Vector k_, Vector E_, Vector H_, double width_, double offset_, double eps_);

    Vector getEField(int i, int j, int k, double time);

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(int) {}
    
  private:
    Vector k;
    Vector E;
    Vector H;
    
    double dt;
    double om;
    double width;
    double offset;
    double eps;
    
    double dx, dy, dz;
    
    Direction dir;
    bool isH;
};


#endif
