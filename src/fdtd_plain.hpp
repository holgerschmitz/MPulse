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
    pField pEx;
    pField pEy;
    pField pEz;
    pField pBx;
    pField pBy;
    pField pBz;
  public:
    void init();
  
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
    
  private:
    void stepD(double dt);
    void stepB(double dt);
    
};

#endif
