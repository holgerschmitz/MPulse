/*
 * fdtd.hpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef FDTD_HPP_
#define FDTD_HPP_

#include "mpulse.hpp"

class FieldSolver : public ChildBlock<FieldSolver> {
  private:
    Field<double, 2> Ex, Ey, Ez;
    Field<double, 2> Bx, By, Bz;
    double eps_rel;

    Array<double, 2> x;
    Array<pParameter, 2> x_par;
    Array<pParameter, 3>  E_par, B_par;
    Array<double, 3> initE, initB;

    void stepD(double dt);
    void stepB(double dt);
  protected:
    void initParameters(BlockParameters &parameters);
    void registerData();
    void init();
  public:
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
};

#endif // FDTD_HPP_
