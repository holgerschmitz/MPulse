#ifndef MPULSE_PLASMACURRENT_H
#define MPULSE_PLASMACURRENT_H

#include "../huerto/electromagnetics/current.hpp"

class PlasmaCurrentBlock : public CurrentBlock
{
  public:
    void initCurrents(CurrentContainer &container);
  protected:
    void initParameters(schnek::BlockParameters &blockPars);

    /**
     * charge for the charge carriers
     * 
     * Default: q_e = 1.602176634e-19
     */
    double charge;
    /**
     * mass of the charge carriers
     * 
     * Default: m_e = 9.1093837015e-31
     */
    double mass;
    /// charge number Z
    double Z;
    /// friction coefficient
    double gamma;
};

class PlasmaCurrent : public Current
{
  protected:
    CurrentBlock &plasmaBlock;

    Field Ex;
    Field Ey;
    Field Ez;

    Field Rho;

    /// charge of the charge carriers
    double charge;
    /// mass of the charge carriers
    double mass;
    /// charge number Z
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
