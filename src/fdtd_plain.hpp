/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_FDTD_PLAIN_H
#define MPULSE_FDTD_PLAIN_H

#include "fieldsolver.hpp"
#include "current.hpp"
#include "mpulse.hpp"

class Storage;

class FDTD_Plain : public FieldSolver, public CurrentContainer, public schnek::BlockContainer<CurrentBlock>
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
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    pDataLine pKappaEdx, pKappaEdy, pKappaEdz;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    pDataLine pKappaHdx, pKappaHdy, pKappaHdz;
  public:

    /**
     * Registers helper grids for sharing
     */
    void registerData();

    void init();
  
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
    
  private:
    void stepD(double dt);
    void stepB(double dt);
    
};

#endif
