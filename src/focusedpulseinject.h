#ifndef MPULSE_FOCUSEDPULSEINJECT_H
#define MPULSE_FOCUSEDPULSEINJECT_H

#include "incsource.h"

#include <fftw3.h>
#include <complex>
#include <list>
#include <schnek/grid.h>


class FocusedPulseDataGenerator
{
  private:
    int currentTime;
    bool initialized;

  public:
    FocusedPulseDataGenerator() : currentTime(-1), initialized(false) {};

    void addEXData(DataGrid *g);
    void addEYData(DataGrid *g);
    void addEZData(DataGrid *g);
    void addBXData(DataGrid *g);
    void addBYData(DataGrid *g);
    void addBZData(DataGrid *g);

    void setTime(int Time);
    void setLow(const GridIndex &low_) { low = low_; }
    void setHigh(const GridIndex &high_) { high = high_; }
    void setSize(int oversample_X_, int oversample_Y_);
    void setShifts(double TShift_, double ZShift_, double Phase_);
    typedef std::complex<double> Complex;
  private:
    typedef schnek::Grid<Complex, 2, MPulseGridChecker> ComplexGrid2d;
    typedef std::list<DataGrid*> GridList;

    GridList ex_grids;
    GridList ey_grids;
    GridList ez_grids;

    GridList bx_grids;
    GridList by_grids;
    GridList bz_grids;

    double TShift;
    double ZShift;
    double Phase;

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
    GridIndex low;
    /** The high index of the fields */
    GridIndex high;

    void init();

    // fill the working grid with data, fourier transform and copy to the
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


/** Calculates a tightly focused pulse in the long pulse limit
 */
class FocusedPulseInject : public IncidentSource
{
  public:
    virtual ~FocusedPulseInject() {}
    void initCurrents(Storage *storage, FieldSolver *solver);
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

    double lightspeed;

    //Complex YComp;

    double length; // corresponds to pulse duration
    double width;
    double om0;

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
    FocusedPulseInjectSourceFunc(Direction dir_, bool isH_) 
      : dir(dir_), isH(isH_) {};
    ~FocusedPulseInjectSourceFunc() { }
    
    void setParam(double length_, 
                  double width_,
                  double om0_,
                  double amp_, 
                  double eps_, 
                  int distance_,
                  FocusedPulseDataGenerator *generator_);

    void initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz);

    Vector getEField(int i, int j, int k, int time);
    Vector getHField(int i, int j, int k, int time);

    void setTime(int Time);

  private:


    double ZRl;
    
    int lowx;
    int highx;
    int lowy;
    int highy;
    
    double amp;
    double eps;
    
    int dim;
    int transverse1, transverse2;

    Direction dir;
    bool isH;
    
    int dist;
    FocusedPulseDataGenerator *generator;
    
    DataGrid *x_grid;
    DataGrid *y_grid;
    DataGrid *z_grid;
    
};



#endif
