/*
 * plasmadensity.hpp
 *
 *  Created on: 22 Nov 2024
 *      Author: Holger Schmitz
 */


#ifndef MPUSE_PLASMADENSITY
#define MPUSE_PLASMADENSITY

#include "types.hpp"
#include "../huerto/simulation/simulation_context.hpp"
#include "../huerto/simulation/initialiser.hpp"

#include <schnek/variables/blockcontainer.hpp>

class MPulse;

/**
 * @brief A container block for plasma density 
 *
 * Multiple sets of densities can be defined
 */
class PlasmaDensity :
        public schnek::ChildBlock<PlasmaDensity>,
        public SimulationEntity
{
  private:
    /// The plasma number density in 1/m^3
    InitialisedField<double> Rho;

    /**
     * Fill the field values from the expressions provided in the setup file
     */
    void fillValues();
  public:

    /**
     * @brief Constructor taking an optional parent block
     *
     * Constructor is compatible with the schnek::Block constructor
     *
     */
    PlasmaDensity(schnek::pBlock parent = schnek::pBlock()) : schnek::ChildBlock<PlasmaDensity>(parent)
    {}

    /**
     * @brief Register the plasma density field
     *
     * The plasma density is created and registered with the Block storage.
     */
    void registerData();

    /**
     * Initialise the parameters available through the setup file
     */
    void initParameters(schnek::BlockParameters &parameters);

    /**
     * Initialise the simulation data
     */
    void init();
};

#endif
