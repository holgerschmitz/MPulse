#ifndef MPULSE_FIELD_IONIZATION_HPP
#define MPULSE_FIELD_IONIZATION_HPP

#include "types.hpp"

#include "../huerto/simulation/task.hpp"
#include "../huerto/simulation/simulation_context.hpp"

#include <schnek/variables/block.hpp>

class MPIIonization :
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

    double rat;
    double mpa;
    double Wion;
    int K;
  protected:
    void initParameters(schnek::BlockParameters &blockPars) override;
    void init() override;
  public:
    MPIIonization(schnek::pBlock parent = schnek::pBlock()) : schnek::Block(parent)
    {}
    
    std::string getPhase() override { return "ionization"; }
    void execute() override;
};

#endif