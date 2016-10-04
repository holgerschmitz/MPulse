#ifndef MPULSE_PLASMACURRENT_H
#define MPULSE_PLASMACURRENT_H

#include "rebuild.h"
#include "currents.h"
#include "currentfactory.h"

class Storage;

class PlasmaCurrentFactory : public CurrentFactory
{
  public:    
    virtual ~PlasmaCurrentFactory() {}
    virtual void initCurrents(Storage *storage_, FieldSolver *solver);
  protected:
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
    
    /// e/m for the charge carriers
    double em;
    /// friction coefficient
    double gamma;
    /// constant background density (used mainly for testing)
    double rho_bg;
};

class PlasmaCurrent : public Current, public Rebuildable
{
  protected:
    DataGrid *pEx;
    DataGrid *pEy;
    DataGrid *pEz;

    DataGrid *pRho;
    Storage *storage;
    
    /// e/m for the charge carriers
    double em;
    
    /// friction coefficient
    double gamma;
    
    /// constant background density (used mainly for testing)
    double rho_bg;
  public:
    PlasmaCurrent(double em_, double gamma_, double rho_bg_);
    
    void initStorage(Storage *storage_);
    
    void stepSchemeInit(double dt) {}
    void stepScheme(double dt);
    
};


#endif
