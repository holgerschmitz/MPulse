#ifndef MPULSE_INCSOURCE_H
#define MPULSE_INCSOURCE_H

#include "mpulse.h"
#include "rebuild.h"
#include "currentfactory.h"
#include "currents.h"

class Storage;
class IncidentSourceCurrent;

//===============================================================
//==========  Base Classes
//===============================================================

class IncidentSource : public CurrentFactory
{
  public:
    virtual ~IncidentSource() {}
    void initCurrents(Storage *storage, FieldSolver *solver);
    
  protected:
    virtual IncidentSourceCurrent *makeECurrent(int distance_, Direction dir_) = 0;
    virtual IncidentSourceCurrent *makeHCurrent(int distance_, Direction dir_) = 0;
    virtual bool needCurrent(Direction dir_) = 0;
    
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:

    int distance;
};

class IncidentSourceCurrent : public Current
{
  public:
    IncidentSourceCurrent(int distance_, Direction dir_, bool isH_);                       

    virtual ~IncidentSourceCurrent() {}
  protected:
    bool reverse;
    int distance;
    
    int dim;
    int transverse1, transverse2;
    
    Direction dir;
    bool isH;
    int lowOffset;
    int highOffset;

    double dx, dy, dz;
    double dt;
    
    DataGrid *pJ[2];
    double dX[3];
};

template<class SourceFunc>
class IncidentSourceECurrent : public IncidentSourceCurrent, public SourceFunc
{
  public:
    IncidentSourceECurrent(int distance_, Direction dir_);
                       
    virtual ~IncidentSourceECurrent() {}

    void initStorage(Storage *storage);

    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:
    
};

template<class SourceFunc>
class IncidentSourceHCurrent : public IncidentSourceCurrent, public SourceFunc
{
  public:
    IncidentSourceHCurrent(int distance_, Direction dir_);
                       
    virtual ~IncidentSourceHCurrent() {}

    void initStorage(Storage *storage);

    void stepSchemeInit(double dt);
    void stepScheme(double dt);
};

//===============================================================
//==========  Side Inject
//===============================================================

/** Reads the data from a file, maybe scaling the
 *  amplitude by a constant factor
 */
class SideInject : public IncidentSource
{
  public:
    virtual ~SideInject() {}
    
  protected:
    virtual IncidentSourceCurrent *makeECurrent(int distance_, Direction dir_);
    virtual IncidentSourceCurrent *makeHCurrent(int distance_, Direction dir_);
    virtual bool needCurrent(Direction dir_);
    
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
    std::string filename;
    double amp;
    int toff;
    int blocks;
    double eps;
};

class SideInjectSourceFunc
{
  public:
    SideInjectSourceFunc(Direction dir_, bool isH_);
    void setParam(std::string filename_, double amp_, double eps_, int toff_, int blocks_);

    void initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz);

    Vector getEField(int i, int j, int k, int time);
    Vector getHField(int i, int j, int k, int time);

    void setTime(int Time);
  private:
    Vector getField(int i, int j, int k, int time, double factor);
    DataGrid F1;
    DataGrid F2;
    
    std::string filename;
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
};


#include "incsource.t"


#endif
