/*
 * fieldsolver.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_FIELDSOLVER_H
#define MPULSE_FIELDSOLVER_H

#include <schnek/variables/blockcontainer.hpp>

class Storage;
class Current;

class FieldSolver : public schnek::ChildBlock<FieldSolver>
{
  public:
    virtual void stepSchemeInit(double dt) = 0;
    virtual void stepScheme(double dt) = 0;    
};

#endif
