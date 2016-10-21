#ifndef MPULSE_PLASMACURRENT_H
#define MPULSE_PLASMACURRENT_H

#include "current.hpp"


class PlasmaCurrentBlock : public CurrentBlock
{
  public:    
    void initCurrents(CurrentContainer &container);
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    
    /// e/m for the charge carriers
    double em;
    /// ion mass
    double mi;
    /// ion charge number Z
    double Z;
    /// friction coefficient
    double gamma;
};

class PlasmaCurrent : public Current
{
  protected:
    CurrentBlock &plasmaBlock;

    pField pEx;
    pField pEy;
    pField pEz;

    pField pRho;
    
    /// e/m for the charge carriers
    double em;
    /// ion mass
    double mi;
    /// ion charge number Z
    double Z;
    /// friction coefficient
    double gamma;
    
  public:
    PlasmaCurrent(double em_, double mi_, double Z_, double gamma_, CurrentBlock &plasmaBlock_);
    
    void init();
    
    void stepSchemeInit(double dt) {}
    void stepScheme(double dt);
    
};


#endif
