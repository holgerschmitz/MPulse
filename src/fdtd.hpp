/*
 * fdtd.hpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef FDTD_HPP_
#define FDTD_HPP_

#include "mpulse.hpp"

class FieldSolver : public schnek::ChildBlock<FieldSolver> {
  private:
    Field Ex, Ey, Ez;
    Field Bx, By, Bz;

    void stepD(double dt);
    void stepB(double dt);
  protected:
    void registerData();
  public:
    void stepSchemeInit(double dt);
    void stepScheme(double dt);

};

typedef boost::shared_ptr<FieldSolver> pFieldSolver;

#endif // FDTD_HPP_
