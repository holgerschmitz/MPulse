/*
 * fdtd_kerr.hpp
 *
 *  Created on: 21 Dec 2020
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_FDTD_KERR_H
#define MPULSE_FDTD_KERR_H

#include "../huerto/electromagnetics/fieldsolver.hpp"
#include "../huerto/electromagnetics/current.hpp"

#include "../huerto/simulation/simulation_context.hpp"

#include "types.hpp"


/**
 * FDTD field solver with instantaneous Kerr nonlinearity
 *
 * The solver allows any number of electric and "magnetic" currents as well as
 */
class FDTD_Kerr : public FieldSolver,
                  public CurrentContainer,
                  public schnek::BlockContainer<CurrentBlock>
{
  private:

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
     * Value of the nonlinearity \f$\chi^{(3)}\f$ in \f$\mathrm{V}^{-2} \mathrm{m}^2\f$
     */
    double chi;

    /**
     * The relative permittivity \f$\eps_r\f$
     */
    double eps;

#ifdef HUERTO_ONE_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaEdx;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaHdx;
#endif

#ifdef HUERTO_TWO_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaEdx, pKappaEdy;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaHdx, pKappaHdy;
#endif

#ifdef HUERTO_THREE_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaEdx, pKappaEdy, pKappaEdz;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaHdx, pKappaHdy, pKappaHdz;
#endif

  public:

    /**
     * Registers helper grids for sharing
     */
    void registerData();

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);

  protected:

    /**
     * Register the nonlinearity \f$chi\f$
     */
    void initParameters(schnek::BlockParameters &blockPars);

  private:
    void stepD(double dt);
    void stepB(double dt);

};

#endif // MPULSE_FDTD_KERR_H
