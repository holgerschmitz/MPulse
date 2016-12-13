/*
 * diagnostic.hpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_DIAGNOSTIC_HPP_
#define MPULSE_DIAGNOSTIC_HPP_

#include "mpulse.hpp"

#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

class FieldDiagnostic : public schnek::HDFGridDiagnostic<Field, pField, schnek::DeltaTimeDiagnostic>
{
  protected:
    typedef HDFGridDiagnostic<Field, pField, schnek::DeltaTimeDiagnostic >::IndexType IndexType;
    IndexType getGlobalMin() { return IndexType(0); }
    IndexType getGlobalMax() { return Simulation::getGlobalMax(); }
};



#endif /* SRC_DIAGNOSTIC_HPP_ */
