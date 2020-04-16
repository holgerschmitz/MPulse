/*
 * diagnostic.hpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_DIAGNOSTIC_HPP_
#define MPULSE_DIAGNOSTIC_HPP_

#include "mpulse.hpp"

#include "../huerto/simulation/simulation_context.hpp"

#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

class FieldDiagnostic :
        public schnek::HDFGridDiagnostic<Field, pField, schnek::DeltaTimeDiagnostic>,
        public SimulationEntity
{
  protected:
    typedef HDFGridDiagnostic<Field, pField, schnek::DeltaTimeDiagnostic>::IndexType IndexType;
    IndexType getGlobalMin();
    IndexType getGlobalMax();
    void init();
};

#endif /* SRC_DIAGNOSTIC_HPP_ */
