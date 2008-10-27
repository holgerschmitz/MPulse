#ifndef MPULSE_SOURCES_H
#define MPULSE_SOURCES_H

#include "incsource.h"

class PlaneWaveSource : public IncidentSource
{
  public:
    ~PlaneWaveSource() {}
  protected:
    IncidentSourceCurrent *makeECurrent(int distance_, Direction dir_);
    IncidentSourceCurrent *makeHCurrent(int distance_, Direction dir_);
    
    bool needCurrent(Direction dir_);
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
    
    double kx, ky, kz;
    double Hx, Hy, Hz;
    double ramp;
    double eps;
    
};

class PlaneWaveSourceEFunc
{
  public:
    PlaneWaveSourceEFunc(Direction dir_, bool isH_);
    void setParam(Vector k_, Vector E_, Vector H_, double ramp_, double eps_);

    Vector getHField(int i, int j, int k, int time);

    void initSourceFunc(Storage*, DataGrid*, DataGrid*, DataGrid*) {}
    void setTime(int) {}
    
  private:
    Vector k;
    Vector E;
    Vector H;
    
    double dt;
    double om;
    double ramp;
    double eps;
    
    double dx, dy, dz;
    
    Direction dir;
    bool isH;
};

class PlaneWaveSourceHFunc
{
  public:
    PlaneWaveSourceHFunc(Direction dir_, bool isH_);
    void setParam(Vector k_, Vector E_, Vector H_, double ramp_, double eps_);

    Vector getEField(int i, int j, int k, int time);

    void initSourceFunc(Storage*, DataGrid*, DataGrid*, DataGrid*) {}
    void setTime(int) {}
    
  private:
    Vector k;
    Vector E;
    Vector H;
    
    double dt;
    double om;
    double ramp;
    double eps;
    
    double dx, dy, dz;
    
    Direction dir;
    bool isH;
};


#endif
