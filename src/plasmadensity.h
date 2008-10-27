#ifndef MPULSE_PLASMADENSITY_H
#define MPULSE_PLASMADENSITY_H

#include "mpulse.h"
#include "rebuild.h"
#include "optfield.h"

class Storage;

class PlasmaDensity : public OptField
{
  protected:
    DataGrid *pEx;
    DataGrid *pEy;
    DataGrid *pEz;
    DataGrid *pRho;
    DataGrid *pSigma;
    Storage *storage;
    
    /// avalanche ionization
    double nu;
    
    /// multi photon ionization
    double mpa;
    
    /// ionization energy
    double Wion;
    
    /// exponent for MPA
    int K;
  public:    
    void initStorage(Storage *storage_);
    
    void stepSchemeInit(double dt) { stepScheme(dt/2); }
    void stepScheme(double dt);
    
    void addToResistivity(DataGrid *pSigma);
    bool addsToResistivity() { return true; }
  protected:
    /// build parametermap
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
};

#endif
