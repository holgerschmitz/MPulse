#ifndef MPULSE_CPML_BORDER_H
#define MPULSE_CPML_BORDER_H

#include "mpulse.h"
#include "rebuild.h"
#include "currentfactory.h"
#include "currents.h"

class Storage;

class CPMLBorder : public CurrentFactory
{
  public:
    virtual ~CPMLBorder() {}
    void initCurrents(Storage *storage, FieldSolver *solver);
  protected:
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
    
    void initCoefficients(Storage *storage);
    
    int thickness;
    double kappaMax;
    double aMax;
    double sigmaMax;
};

class CPMLBorderOneD : public CurrentFactory
{
  public:
    virtual ~CPMLBorderOneD() {}
    void initCurrents(Storage *storage, FieldSolver *solver);
  protected:
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
    
    void initCoefficients(Storage *storage);
    
    int thickness;
    double kappaMax;
    double aMax;
    double sigmaMax;
};

class CPMLBorderCurrent : public Current
{
  public:
    CPMLBorderCurrent( int thickness_, Direction dir_, bool isH_,
                       double kappaMax_, double aMax_, double sigmaMax_, double eps_);                       

    virtual ~CPMLBorderCurrent() {}
  protected:
    bool reverse;
    int thickness;
    
    int dim;
    int transverse1, transverse2;
    
    Direction dir;
    bool isH;
    int lowOffset;
    int highOffset;
    
    int zerolayer;
    
    double kappaMax;
    double aMax;
    double sigmaMax;
    double eps;
    
    DataLine bCoeff;
    DataLine cCoeff;


    void makeCoeff(Storage *storage);
};

class CPMLBorderECurrent : public CPMLBorderCurrent
{
  public:
    CPMLBorderECurrent( int thickness_, Direction dir_, 
                        double kappaMax_, double aMax_, double sigmaMax_, double eps_);
                       
    virtual ~CPMLBorderECurrent() {}

    void initStorage(Storage *storage);
    
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:
    
    DataGrid *pB[3];
    DataGrid *pPsi[2];
    double dx;
};

class CPMLBorderHCurrent : public CPMLBorderCurrent
{
  public:
    CPMLBorderHCurrent( int thickness_, Direction dir_,
                        double kappaMax_, double aMax_, double sigmaMax_, double eps_);
                       
    virtual ~CPMLBorderHCurrent() {}

    void initStorage(Storage *storage);
    
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:
    DataGrid *pE[3];
    DataGrid *pPsi[2];
    double dx;
};


#endif
