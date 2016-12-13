/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_FDTD_H
#define MPULSE_FDTD_H

#include "mpulse.hpp"

class FDTDSolver : public schnek::ChildBlock<FDTDSolver>
{
  private:
    pField pEx, pEy, pEz;
    pField pBx, pBy, pBz;
  protected:
    void registerData();
    void init();

  public:
    void stepScheme(double dt);
    
  private:
    void stepE(double dt);
    void stepB(double dt);
};

typedef boost::shared_ptr<FDTDSolver> pFDTDSolver;

#endif
