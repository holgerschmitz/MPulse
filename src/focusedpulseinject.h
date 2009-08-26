#ifndef MPULSE_FOCUSEDPULSEINJECT_H
#define MPULSE_FOCUSEDPULSEINJECT_H

#include "incsource.h"

#include <fftw3.h>
#include <complex>
#include <list>
#include <schnek/matrix.h>

/** Calculates a short pulse in the paraxial limit
 */
class FocusedPulseInject : public IncidentSource
{
  public:
    virtual ~FocusedPulseInject() {}
    
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
                  double TShift_,
                  double ZShift_,
                  double Phase_,
                  double amp_, 
                  double eps_, 
                  int distance_,
                  FocusedPulseDataGenerator *generator_);

    void initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz);

    Vector getEField(int i, int j, int k, int time);
    Vector getHField(int i, int j, int k, int time);

    void setTime(int Time);

  private:
        
           
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
    FocusedPulseDataGenerator *generator;
    
    DataGrid x_grid;
    DataGrid y_grid;
    DataGrid z_grid;
    
};


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
    void setLow(GridIndex &low_) { return low = low_; }
    void setHigh(GridIndex &high_) { return high = high_; }

    typedef std::complex<double> Complex;
  private:
    
    std::list<DataGrid*> ex_grids;
    std::list<DataGrid*> ey_grids;
    std::list<DataGrid*> ez_grids;
    
    
    std::list<DataGrid*> bx_grids;
    std::list<DataGrid*> by_grids;
    std::list<DataGrid*> bz_grids;
    
    DataGrid working_grid;
    
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
    
    Complex ExFunc(double x, double y, double z, double t);
    Complex Bfunc(double x, double y, double z, double t, bool bx);
};

#endif
