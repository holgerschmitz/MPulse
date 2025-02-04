#ifndef MPULSE_MPI_IONIZATION_HPP
#define MPULSE_MPI_IONIZATION_HPP

#include "types.hpp"

#include "../huerto/simulation/task.hpp"
#include "../huerto/simulation/simulation_context.hpp"

#include <schnek/variables/block.hpp>

class FieldIonization :
        public schnek::Block,
        public SimulationTask,
        public SimulationEntity
{
  private:
    Field Ex;
    Field Ey;
    Field Ez;

    Field Electrons;
    Field Neutrals;

    double Z;
    double Uion;
    double Wion;
    double neff;
    double rat;
  private:
    double computeIonizationRate(double E);
  protected:
    void initParameters(schnek::BlockParameters &blockPars) override;
    void init() override;
  public:
    FieldIonization(schnek::pBlock parent = schnek::pBlock()) : schnek::Block(parent)
    {}
    
    std::string getPhase() override { return "ionization"; }
    void execute() override;
};

#endif