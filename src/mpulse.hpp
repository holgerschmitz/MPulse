/*
 * mpulse.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include "types.hpp"

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
               public boost::enable_shared_from_this<MPulse>
{
  private:
	/**
	 * The reference to the root block instance
	 */
    static MPulse *instance;

    /**
     * The maximum index of the global grid.
     *
     * This is #gridSize - 1
     */
    Index globalMax;

    /**
     * The global grid size, \f$N\f$
     *
     * This quantity is specified through the `N` parameter in the setup file
     */
    Index gridSize;

    /**
     * The physical size of the global grid, \f$L\f$
     *
     * This quantity is specified through the `L` parameter in the setup file
     */
    Vector size;

    /**
     * The grid spacing
     *
     * Calculated as \f$dx_i = L_i / N_i\f$
     */
    Vector dx;

    /**
     * The CLF factor \f$f_{\mathrm{CLF}}\f$
     *
     * A factor that is multiplied into the timestep #dt of the simulation.
     */
    double cflFactor;

    /**
     * The time step, \f$dt\f$, of the simulation
     *
     * The time step is calculated as
     *
     * \f[
     * dt = f_{\mathrm{CLF}} \frac{\min_i(dx_i)}{c}
     * \f]
     */
    double dt;

    /**
     * The total physical simulation time
     *
     * This quantity is specified through the `tMax` parameter in the setup file
     */
    double tMax;

    /**
     * The current physical simulation time
     */
    double time;

    /**
     * The inner grid range of the local simulation domain
     */
    Range innerRange;

    /**
     * The local grid of the electric field components, \f$\mathbf{E}\f$
     *
     * The size of the grid is determined by #innerRange with an additional
     * boundary of two ghost cells.
     *
     * The grids are staggered according to Yee's staggering
     */
    pField pEx, pEy, pEz;

    /**
     * The local grid of the magnetic field components, \f$\mathbf{B}\f$
     *
     * The size of the grid is determined by #innerRange with an additional
     * boundary of two ghost cells.
     *
     * The grids are staggered according to Yee's staggering
     */
    pField pBx, pBy, pBz;

    /**
     * The Cartesian MPI subdivision of the grids containing the electromagnetic field
     */
    schnek::MPICartSubdivision<Field> subdivision;

    /**
     * The free \f$x\f$ coordinate vector
     *
     * This vector is associated with the read-only parameter in the configuration file
     */
    Vector x;

    /**
     * The read-only \f$x\f$ parameter in the configuration file
     */
    schnek::Array<schnek::pParameter, DIMENSION> x_parameters;

    /**
     * The electric field parameter in the configuration file
     *
     * This parameter is used to initialise #pEx, #pEy, and #pEz
     */
    schnek::Array<schnek::pParameter, DIMENSION> E_parameters;

    /**
     * The magnetic field parameter in the configuration file
     *
     * This parameter is used to initialise #pBx, #pBy, and #pBz
     */
    schnek::Array<schnek::pParameter, DIMENSION> B_parameters;

    /**
     * A parameter group containing the \f$x\f$ parameter
     *
     * @todo Check if MPulse.spaceVars is really needed
     */
    schnek::pParametersGroup spaceVars;

    /**
     * The initialisation vector for initialising the electric fields
     */
    Vector initE;

    /**
     * The initialisation vector for initialising the magnetic fields
     */
    Vector initB;
  private:

    /**
     * Initialises the electromagnetic fields.
     *
     * If no #EMField child block was specified, it is created and pre-initialised
     * here.
     *
     * The grids are then obtained from the storage and resized according to the
     * local grid size determined by #subdivision.
     */
    void initFields();

    /**
     * Fills the local grids with values
     *
     * A dependency updater is used to initialise the electromagnetic fields
     * with space dependent values provided through the setup file.
     */
    void fillValues();
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

    /**
     * Get the global maximum grid index, #globalMax
     */
    static Index getGlobalMax() { return instance->globalMax; }

    /**
     * Get the grid spacing #dx
     */
    static Vector getDx() { return instance->dx; }

    /**
     * Get the time step, #dt
     */
    static double getDt() { return instance->dt; }

    /**
     * Get the physical size of the global grid, \f$L\f$
     */
    static Vector getSize() { return instance->size; }

    /**
     * Get the current simulation time
     */
    static double getTime() { return instance->time; }

    /**
     * Get the grid subdivision
     */
    static schnek::DomainSubdivision<Field> &getSubdivision() { return instance->subdivision; };
};

#endif // MPULSE_MPULSE_H
