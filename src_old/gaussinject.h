#ifndef MPULSE_GAUSSINJECT_H
#define MPULSE_GAUSSINJECT_H

#include "incsource.h"

#include <fftw3.h>
#include <complex>
#include <schnek/grid.h>

/** Calculates a gaussian pulse using Fourier Transform
 */
class GaussInject : public IncidentSource
{
  public:
    virtual ~GaussInject() {}
    
  protected:
    virtual IncidentSourceCurrent *makeECurrent(int distance_, Direction dir_);
    virtual IncidentSourceCurrent *makeHCurrent(int distance_, Direction dir_);
    virtual bool needCurrent(Direction dir_);
    
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
    double lightspeed;

    int SizeX;
    int SizeY; 
    int SizeT; 

    int SizeXH;
    int SizeYH;
    int SizeTH;

    double BoxSizeX;
    double BoxSizeY;
    double BoxSizeT;

    double DX;
    double DY;
    double DT;


    double length; // corresponds to pulse duration
    double width;
    double k0;
    double Shift;
    double Time;


    double amp;
    int toff;
    double eps;
};

class GaussInjectSourceFunc
{
  public:
    GaussInjectSourceFunc(Direction dir_, bool isH_) : dir(dir_), isH(isH_) {};
    ~GaussInjectSourceFunc() { clearfftw(); }
    
    void setParam(int    SizeT_, 
                  double length_, 
                  double width_,
                  double k0_,
                  double Shift_, 
                  double Time_,
                  double amp_, 
                  double eps_, 
                  int toff_);

    void initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz);

    Vector getEField(int i, int j, int k, int time);
    Vector getHField(int i, int j, int k, int time);

    void setTime(int Time);
  private:
    
    Vector getField(int i, int j, int k, int time, double factor);
    
    
    DataGrid F1;
    DataGrid F2;
       
    double lightspeed;

    int SizeX;
    int SizeY; 
    int SizeT; 

    int SizeXH;
    int SizeYH;
    int SizeTH;

    double BoxSizeX;
    double BoxSizeY;
    double BoxSizeT;

    double DX;
    double DY;
    double DT;


    double length; // corresponds to pulse duration
    double width;
    double k0;
    double Shift;
    double Time;
    
    int lowx;
    int highx;
    int lowy;
    int highy;

    double amp;
    double eps;
    int toff;
    int Nt;
    int blockCount;
    int blocks;
    int blockToff;
    bool active;
    
    int dim;
    int transverse1, transverse2;

    Direction dir;
    bool isH;
    
  private:
  
    typedef std::complex<double> Complex;
    typedef schnek::Grid<Complex,1> Row;

    Row A;
    DataGrid2d tmpField_TYr, tmpField_TYi;
    DataGrid3d tmpField_TXr, tmpField_TXi;
   
    DataGrid2d debug_om_k, debug_t_k;

    /// Temporal storage of complex numbers to interface with fftw
    fftw_complex *ukx;
    fftw_complex *ux; 
    fftw_complex *uky;
    fftw_complex *uy; 

    fftw_complex *uo;
    fftw_complex *ut; 

    /** fftw plans are set up in the beginning and then used to calculate 
     * the Fourier transform 
     */
    fftw_plan pfftxb, pfftyb, pffttb;

    int sign(double x) { return x>0?1:(x<0?-1:0); }
    
    void makefftw();
    void clearfftw();
    void transformT(int i, int j);
    void transformY(int i);
    void transformX(DataGrid &F);
    void fillEx();
    void fillEy();
    void fillBx();

    struct GaussFunc
    {
      double sx, sy, sz;
      double k0, length, width;
      
      double kp(double x, double y, double z) { return sqrt(x*x + z*z); }
      double kp2(double x, double y, double z) { return x*x + z*z; }
      double kt(double x, double y, double z) { return sqrt(x*x + y*y); }
      double k(double x, double y, double z) { return sqrt(x*x + y*y + z*z); }

      double gauss(double x, double y, double z)
      {
        double zoff = z-k0;
        return exp(-0.25*zoff*zoff*length*length - 0.25*(x*x + y*y)*width*width);
      }
      
      void setShift(double sx_, double sy_, double sz_) { sx=sx_; sy=sy_; sz=sz_; }
      
      double shift(double x, double y, double z) { return -(x*sx + y*sy + z*sz);  }
      
      double excomp(double x, double y, double z) 
      { 
        double D = k(x,y,z)*kp(x,y,z);
        return D==0?0:(-x*y/D); 
      }
      
      double eycomp(double x, double y, double z) 
      { 
        double D = k(x,y,z);
        return D==0?0:(kp(x,y,z)/D); 
      }
      
      double bxcomp(double x, double y, double z) 
      { 
        double D = kp(x,y,z);
        return D==0?0:(-z/D);
      }
    };
    
    GaussFunc gaussfunc;
};

#endif
