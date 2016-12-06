/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_FDTD_PLAIN_H
#define MPULSE_FDTD_PLAIN_H

#include "fieldsolver.hpp"
#include "mpulse.hpp"

class Storage;

class FDTD_Plain : public FieldSolver
{
  private:
    Field Ex, Ey, Ez;
    Field Bx, By, Bz;
  public:
    void init();
  
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
    
  private:
    void stepD(double dt);
    void stepB(double dt);
    
};

#endif
