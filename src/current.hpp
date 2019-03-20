/*
 * current.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: holger
 */


#ifndef MPULSE_CURRENT_H
#define MPULSE_CURRENT_H

#include "mpulse.hpp"

class Current;
typedef boost::shared_ptr<Current> pCurrent;

/**
 * A container for currents that distinguishes between electric and magnetic
 * currents.
 *
 * The magnetic currents are current-like fields that enter the Faraday equation
 * much like the electric current enters the Ampere equation. These can be used
 * to provide absorbing boundary layers. They are not real physical currents but
 * mathematical constructs.
 */
class CurrentContainer
{
  protected:
    /**
     * A list of currents
     */
    typedef std::list<pCurrent> CurrentList;

    /**
     * The list of electric currents
     */
    CurrentList currents;

    /**
     * The list of magnetic currents
     */
    CurrentList magCurrents;
  public:

    /**
     * Add an electric current to the list
     */
    void addCurrent(pCurrent current);

    /**
     * Add a magnetic current to the list
     */
    void addMagCurrent(pCurrent current);
};

/**
 * An abstract block type that calculates currents
 */
class CurrentBlock : public schnek::ChildBlock<CurrentBlock>
{
  public:
    /**
     * Initialise any currents with the current container
     *
     * This usually involves adding the currents to the container.
     *
     * This should be called by the parent block during the `Block.init` phase
     */
    virtual void initCurrents(CurrentContainer &container) = 0;
};

/**
 * A shared pointer to a current block
 */
typedef boost::shared_ptr<CurrentBlock> pCurrentBlock;

/**
 * An abstract base class for currents
 */
class Current
{
  protected:
    /**
     * The grids containing the components of the current
     */
    pGrid pJx, pJy, pJz;
  public:
    /**
     * Virtual destructor
     */
    virtual ~Current() {}
    
    /**
     * Initialise the currents
     *
     * This method will be called during the `Block.init` phase.
     */
    virtual void init() = 0;

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

    /**
     * Check if the current is valid
     *
     * The current is valid if all the components of the current have been
     * initialised
     */
    bool isValid() const { return (pJx) && (pJy) && (pJz); }
    
    /**
     * Get the x-component of the current
     */
    const pGrid getJx() { return pJx; }

    /**
     * Get the y-component of the current
     */
    const pGrid getJy() { return pJy; }

    /**
     * Get the z-component of the current
     */
    const pGrid getJz() { return pJz; }
};


#endif
