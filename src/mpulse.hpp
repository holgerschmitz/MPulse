/*
 * mpulse.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include "types.hpp"

#include "../huerto/simulation/simulation_context.hpp"
#include "../huerto/simulation/task.hpp"

#include <schnek/grid.hpp>
#include <schnek/variables.hpp>
#include <boost/enable_shared_from_this.hpp>


class FieldSolver;
class EMFields;

/**
 * Base MPulse simulation class
 *
 * This class represents the root block of the simulation
 *
 * MPulse contains two types of child blocks.
 *
 * * #FieldSolver contains the algorithm for solving the EM fields on the grid
 * * #EMFields is a container for the electromagnetic field grids. A child block
 *   of this type not required to be specified through the setup file. If no
 *   #EMField is specified it will be created automatically.
 */
class MPulse : public schnek::Block,
               public schnek::BlockContainer<FieldSolver>,
               public schnek::BlockContainer<EMFields>,
               public boost::enable_shared_from_this<MPulse>,
               public SimulationContext,
               public SimulationTaskRunner
{
  private:
    /**
     * The maximum index of the global grid.
     *
     * This is #gridSize - 1
     */
    Index globalMax;

    /**
     * The CLF factor \f$f_{\mathrm{CLF}}\f$
     *
     * A factor that is multiplied into the timestep #dt of the simulation.
     */
    double cflFactor;

    /**
     * The inner grid range of the local simulation domain
     */
    Range innerRange;

    /**
     * If true, then don't do an inital half time step to bring E and B field in sync.
     */
    bool ignore_initial_time_stagger;
  protected:

    /**
     * Registers the parameters to be read from the setup file
     */
    void initParameters(schnek::BlockParameters &blockPars);
  public:
    /**
     * Create a new MPulse instance
     *
     * The static #instance member is set to `this`
     */
    MPulse();

    /**
     * Initialise some parameters and then calls #initFields and #fillValues to
     * set up the electromagnetic fields
     */
    void init();

    /**
     * Run the simulation.
     *
     * For each field solver child block FieldSolver#stepSchemeInit is called once.
     *
     * The inside the main loop FieldSolver#stepScheme is called repeatedly until
     * the simulation has completed.
     *
     * During execution, the diagnostic manager is executed to produce diagnostic
     * output.
     */
    void execute();
};

#endif // MPULSE_MPULSE_H
