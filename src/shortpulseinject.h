#ifndef MPULSE_SHORTPULSEINJECT_H
#define MPULSE_SHORTPULSEINJECT_H

#include "incsource.h"

#include <fftw3.h>
#include <complex>
#include <schnek/matrix.h>

/** Calculates a short pulse in the paraxial limit
 */
class ShortPulseInject : public IncidentSource
{
  public:
    virtual ~ShortPulseInject() {}
    
  protected:
    virtual IncidentSourceCurrent *makeECurrent(int distance_, Direction dir_);
    virtual IncidentSourceCurrent *makeHCurrent(int distance_, Direction dir_);
    virtual bool needCurrent(Direction dir_);
    
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
    double lightspeed;

    double length; // corresponds to pulse duration
    double width;
    double om0;
    double TShift;
    double ZShift;
    double Time;
    double Phase;


    double amp;
    double eps;
};

class ShortPulseInjectSourceFunc
{
  public:
    ShortPulseInjectSourceFunc(Direction dir_, bool isH_) : dir(dir_), isH(isH_) {};
    ~ShortPulseInjectSourceFunc() { }
    
    void setParam(double length_, 
                  double width_,
                  double om0_,
                  double TShift_,
                  double ZShift_,
                  double Phase_,
                  double amp_, 
                  double eps_, 
                  int distance_);

    void initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz);

    Vector getEField(int i, int j, int k, int time);
    Vector getHField(int i, int j, int k, int time);

    void setTime(int Time);

    typedef std::complex<double> Complex;
  private:
        
    Complex Efunc(double x, double y, double z, double t);
    Complex Bfunc(double x, double y, double z, double t, bool bx);
           
    double DX;
    double DY;
    double DZ;
    double DT;

    double lightspeed;

    Complex YComp;

    double length; // corresponds to pulse duration
    double width;
    double om0;
    double TShift;
    double ZShift;
    double Phase;

    double ZRl;
    
    int lowx;
    int highx;
    int lowy;
    int highy;
    
    double centrex;
    double centrey;
    double centrez;

    double amp;
    double eps;
    
    int dim;
    int transverse1, transverse2;

    Direction dir;
    bool isH;
    
    int dist;
};

#endif
