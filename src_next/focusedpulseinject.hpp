#ifndef MPULSE_FOCUSEDPULSEINJECT_H
#define MPULSE_FOCUSEDPULSEINJECT_H

#include "mpulse.hpp"

#include "../huerto/electromagnetics/source/incsource.hpp"

#include <schnek/grid.hpp>

#include <fftw3.h>
#include <complex>
#include <list>

class FocusedPulseDataGenerator
{
  private:
    int currentTime;
    bool initialized;

    SimulationContext &context;

  public:
    FocusedPulseDataGenerator(SimulationContext &context) :
      currentTime(-1),
      initialized(false),
      context(context) {};

    void addEXData(pGrid g);
    void addEYData(pGrid g);
    void addEZData(pGrid g);
    void addBXData(pGrid g);
    void addBYData(pGrid g);
    void addBZData(pGrid g);

    void setTime(int Time);
    void setLow(const Index &low_) { low = low_; }
    void setHigh(const Index &high_) { high = high_; }
    void setSize(int oversample_X_, int oversample_Y_);
    void setShifts(double TShift_, double ZShift_, double Phase_);
    typedef std::complex<double> Complex;
  private:
    typedef schnek::Grid<Complex, 2, HuertoGridChecker> ComplexGrid2d;
    typedef std::list<pGrid> GridList;

    GridList ex_grids;
    GridList ey_grids;
    GridList ez_grids;

    GridList bx_grids;
    GridList by_grids;
    GridList bz_grids;

    double TShift;
    double ZShift;
    double Phase;
    double lightspeed;

    double length; // corresponds to pulse duration
    double width;
    double om0;
    double Time;


    int oversample_X; ///< Oversampling factor
    int oversample_Y; ///< Oversampling factor
    int oversample_Z; ///< Oversampling factor

    /// global size of simulation domain
    int gLowX, gLowY;
    int gHighX, gHighY;

    double centrez;

    double DX;
    double DY;
    double DZ;
    double DT;

    /// size of the grid used in Fourier Transform
    int sizeX;
    int sizeY;

    /// Temporal storage of complex numbers to interface with fftw
    fftw_complex *fft_field_k;
    fftw_complex *fft_field;

    /** fftw plans are set up in the beginning and then used to calculate
     * the Fourier transform
     */
    fftw_plan pfft;

    /** The low index of the fields */
    Index low;
    /** The high index of the fields */
    Index high;

    void init();

    // fill the working grid with data, Fourier transform and copy to the
    // respective grids in the lists
    void calcEx();
    void calcEy();
    void calcEz();
    void calcBx();
    void calcBy();
    void calcBz();

    void calcField(GridList &grids, Complex Shift);

//    Complex ExFunc(double x, double y, double z, double t);
//    Complex Bfunc(double x, double y, double z, double t, bool bx);
    Complex FieldFuncs(double kx, double ky, double z, double t, int fieldid);
};


/**
 * Injects a tightly focused pulse in the long pulse limit
 */
class FocusedPulseInject : public IncidentSource
{
  public:
    /**
     * Overrides #IncidentSource.initCurrents to initialise the data generator
     */
    void initCurrents(CurrentContainer &container);
  protected:

    pCurrent makeECurrent(int distance_, Direction dir_);
    pCurrent makeHCurrent(int distance_, Direction dir_);
    bool needCurrent(Direction dir_);

    void initParameters(schnek::BlockParameters &blockPars);
  private:

    double length; // corresponds to pulse duration
    double width;
    double om0;
    double TShift;
    double ZShift;
    double Time;
    double Phase;

    //Complex YComp;

    int oversample_X;
    int oversample_Y;
    int oversample_Z;

    double amp;
    double eps;
    FocusedPulseDataGenerator generator;
};


class FocusedPulseInjectSourceFunc
{
  public:
    FocusedPulseInjectSourceFunc(Direction dir_, SimulationContext &context)
      : dir(dir_), context(context) {};

    void setParam(double length_,
                  double width_,
                  double om0_,
                  double amp_,
                  double eps_,
                  int distance_,
                  FocusedPulseDataGenerator *generator_);

    void initSourceFunc(pGrid pJx, pGrid pJy, pGrid pJz);

    Vector getEField(int i, int j, int k, int time);
    Vector getHField(int i, int j, int k, int time);

    void setTime(int Time);

  private:


    double ZRl;

    int lowx;
    int highx;
    int lowy;
    int highy;

    double length; // corresponds to pulse duration
    double width;
    double om0;

    double amp;
    double eps;

    int dim;
    int transverse1, transverse2;

    Direction dir;

    // @TODO This needs to be initialised
    bool isH;

    int dist;
    FocusedPulseDataGenerator *generator;

    pGrid x_grid;
    pGrid y_grid;
    pGrid z_grid;

    SimulationContext &context;
};



#endif
