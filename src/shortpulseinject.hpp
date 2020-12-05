#ifndef MPULSE_SHORTPULSEINJECT_H
#define MPULSE_SHORTPULSEINJECT_H

#include "../huerto/electromagnetics/source/incsource.hpp"
#include "mpulse.hpp"

#include <fftw3.h>
#include <complex>

/**
 * Calculates a short pulse in the paraxial limit
 */
class ShortPulseInject : public IncidentSource
{
  protected:
    pCurrent makeECurrent(int distance_, Direction dir_);
    pCurrent makeHCurrent(int distance_, Direction dir_);
    bool needCurrent(Direction dir_);

    void initParameters(schnek::BlockParameters &blockPars);
    void init();
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
    ShortPulseInjectSourceFunc(Direction dir_, SimulationContext &context) :
      dir(dir_), context(context) {};

    void setParam(double length_,
                  double width_,
                  double om0_,
                  double TShift_,
                  double ZShift_,
                  double Phase_,
                  double amp_,
                  double eps_,
                  int distance_);

    void initSourceFunc(pGrid pJx, pGrid pJy, pGrid pJz);

    Vector getEField(int i, int j, int k, double time);
    Vector getHField(int i, int j, int k, double time);

    void setTime(double Time);

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

    int dist;
    SimulationContext &context;
};

#endif
