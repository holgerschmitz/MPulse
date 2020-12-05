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
    pCurrent makeECurrent(int distance, Direction dir);
    pCurrent makeHCurrent(int distance, Direction dir);
    bool needCurrent(Direction dir);

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
    ShortPulseInjectSourceFunc(Direction dir, SimulationContext &context) :
      dir(dir), context(context) {};

    void setParam(double length,
                  double width,
                  double om0,
                  double TShift,
                  double ZShift,
                  double Phase,
                  double amp,
                  double eps,
                  int distance);

    void initSourceFunc(pGrid pJx, pGrid pJy, pGrid pJz);

#ifdef HUERTO_ONE_DIM
    Vector3d getEField(int i, double time);
#endif

#ifdef HUERTO_TWO_DIM
    Vector3d getEField(int i, int j, double time);
#endif

#ifdef HUERTO_THREE_DIM
    Vector3d getEField(int i, int j, int k, double time);
#endif

#ifdef HUERTO_ONE_DIM
    Vector3d getHField(int i, double time);
#endif

#ifdef HUERTO_TWO_DIM
    Vector3d getHField(int i, int j, double time);
#endif

#ifdef HUERTO_THREE_DIM
    Vector3d getHField(int i, int j, int k, double time);
#endif

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
