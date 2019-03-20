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

/**
 * An interface for field solvers
 *
 * Field solvers form the basis of the EM field simulation. The field solver
 * is a child block of the #MPulse root block.
 *
 * Field solvers should not create their own EM field grids. This is created in
 * the #EMField class. They should instead acquire their copy through the
 * `retrieveData` method.
 */
class FieldSolver : public schnek::ChildBlock<FieldSolver>
{
  public:
    /**
     * Initialise the simulation step.
     *
     * This function will be called once during execution of the simulation.
     * It is not intended to perform any household tasks. Instead is should be used
     * to perform numerical computations that are part of the main scheme but
     * need only be computed before the first time step.
     */
    virtual void stepSchemeInit(double dt) = 0;

    /**
     * Perform a simulation time step
     */
    virtual void stepScheme(double dt) = 0;    
};

typedef boost::shared_ptr<FieldSolver> pFieldSolver;

#endif
