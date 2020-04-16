#ifndef MPULSE_FDTD_PLRC_H
#define MPULSE_FDTD_PLRC_H

#include "fieldsolver.hpp"
#include "current.hpp"
#include "mpulse.hpp"

#include <complex>

class Storage;

/**
 * Core solver for the Finite Difference Time Domain, Piecewise Linear Recursive
 * Convolution (FDTD-PLRC) solver.
 *
 * This core class defines only data members and not any functionality. The PLRC
 * method is assumed to simulate a linear dispersive medium with three
 * Lorentz poles
 *
 * It inherits from #CurrentContainer to hold an arbitrary number of electric
 * currents and "magnetic currents". The latter are current-like fields
 * that enter the Faraday equation much like the electric current enters the
 * Ampere equation. These can be used to provide absorbing boundary layers.
 */
class FDTD_PLRCCore :
        public FieldSolver,
        public CurrentContainer,
        public schnek::BlockContainer<CurrentBlock>
{
  public:
    /**
     * Initialises and registers all internal field grids
     */
    void registerData();

    /**
     * Retrieves the electromagnetic field grids and initialises the child
     * #CurrentBlock instances.
     */
    void init();
  protected:

    /**
     * References to the local grid of the electric field components,
     * \f$\mathbf{E}\f$
     */
    pField pEx, pEy, pEz;

    /**
     * References to the local grid of the magnetic field components,
     * \f$\mathbf{E}\f$
     */
    pField pBx, pBy, pBz;

    /**
     * Electrical conductivity of the medium
     */
    pField pSigma;

    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    pDataLine pKappaEdx, pKappaEdy, pKappaEdz;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    pDataLine pKappaHdx, pKappaHdy, pKappaHdz;

    /**
     * Conductivity factor for the electric fields used by CPML schemes
     */
    pDataLine pCpmlSigmaEx, pCpmlSigmaEy, pCpmlSigmaEz;

    /**
     * Conductivity factor for the magnetic fields used by CPML schemes
     */
    pDataLine pCpmlSigmaHx, pCpmlSigmaHy, pCpmlSigmaHz;

    /**
     * Real part of the accumulator for three poles
     */
    pField pPsiRx[3], pPsiRy[3], pPsiRz[3];

    /**
     * Imaginary part of the accumulator for three poles
     */
    pField pPsiIx[3], pPsiIy[3], pPsiIz[3];

    struct PLRCData
    {
      std::complex<double> dchi0[3];
      std::complex<double> dxi0[3];
      std::complex<double> Crec[3];
      double sumChi0;
      double sumXi0;
    };

    PLRCData plrcData;

    /**
     * value of eps_infty (note that eps_0 = mu_0 = 1)
     */
    double eps;

    /**
     * Value of Delta epsilon_p for the three Lorentz poles
     */
    double LEps[3];

    /**
     * Value of delta_p^2 for the three Lorentz poles
     */
    double LDelta[3];

    /**
     * Value of omega_p^2 for the three Lorentz poles
     */
    double LOm2[3];
};

/**
 * Extension of the core FDTD-PLRC class, #FDTD_PLRCCore, for linear media.
 *
 * It advances the magnetic fields according to the standard FDTD method. The
 * electric fields are integrated according to equations 9.19 and 9.20 in
 * Taflove and Hagness, Computational Electrodynamics 3rd edition, Artech House,
 * Boston, London, 2005
 *
 */
class FDTD_PLRCLinCore : public FDTD_PLRCCore
{
  protected:

    /**
     * Advance the electric fields
     */
    void plrcStepD(double dt,
                   int i, int j, int k,
                   double dx, double dy, double dz,
                   double Jx, double Jy, double Jz);

    /**
     * Advance the magnetic fields
     */
    void plrcStepB(double dt,
                   int i, int j, int k,
                   double dx, double dy, double dz,
                   double Jx, double Jy, double Jz);
};

/**
 * Extension of the core FDTD-PLRC class, #FDTD_PLRCCore, for nonlinear media.
 *
 * It advances the magnetic fields according to the standard FDTD method. The
 * electric fields are integrated according to equations 15 and 20 in
 * JOSA B Vol. 29, Issue 6, pp. 1208-1217 (2012)
 */
class FDTD_PLRCNonlinCore : public FDTD_PLRCCore
{
  protected:

    /**
     * Advance the electric fields
     */
    void plrcStepD(double dt,
                   int i, int j, int k,
                   double dx, double dy, double dz,
                   double Jx, double Jy, double Jz);

    /**
     * Advance the magnetic fields
     */
    void plrcStepB(double dt,
                   int i, int j, int k,
                   double dx, double dy, double dz,
                   double Jx, double Jy, double Jz);

    /**
     * Register the nonlinearity \f$chi\f$
     */
    void initParameters(schnek::BlockParameters &blockPars);

    /**
     * Adds a literature reference
     */
    void init();

    /**
     * Value of the nonlinearity \f$\chi^3\f$
     */
    double chi;

};

/**
 * Generic implementation of the FDTD-PLRC solver
 *
 * This generic class finalises the implementations of #FDTD_PLRCLinCore or
 * #FDTD_PLRCNonlinCore by implementing the #FieldSolver.stepScheme and
 * #FieldSolver.stepSchemeInit methods.
 *
 * The template parameter should be either one of #FDTD_PLRCLinCore or
 * #FDTD_PLRCNonlinCore.
 */
template<class PLRCImplementation>
class FDTD_PLRCSolver : public PLRCImplementation
{
  public:

    /**
     * Initialise the internal fields
     */
    void init();

    /**
     * Perform half a time step on the magnetic fields and call
     * #Current#stepSchemeInit for all electric and magnetic currents.
     */
    void stepSchemeInit(double dt);

    /**
     * Perform a full time step on the electric field and magnetic field and call
     * #Current#stepScheme for all electric and magnetic currents.
     */
    void stepScheme(double dt);
  protected:
    typedef PLRCImplementation Implementation;

    /**
     * Register the parameters to read epsilon, delta and omega for all three
     * Lorentz poles.
     */
    void initParameters(schnek::BlockParameters &blockPars);

  private:

    /**
     * Advance the electric fields
     *
     * Accumulates all magnetic currents and then calls `plrcStepD`
     */
    void stepD(double dt);

    /**
     * Advance the magnetic fields
     *
     * Accumulates all electric currents and then calls `plrcStepB`
     */
    void stepB(double dt);

    /**
     * Initialise the accumulators in #FDTD_PLRCCore at the start of the
     * simulation
     */
    void initAccumulator(double dt);

};

typedef FDTD_PLRCSolver<FDTD_PLRCLinCore> FDTD_PLRCLin;
typedef FDTD_PLRCSolver<FDTD_PLRCNonlinCore> FDTD_PLRCNonlin;

#include "fdtd_plrc.t"

#endif
