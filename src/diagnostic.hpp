/*
 * diagnostic.hpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef DIAGNOSTIC_HPP_
#define DIAGNOSTIC_HPP_

#include "mpulse.hpp"

#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

class FieldDiagnostic : public schnek::HDFGridDiagnostic<Field, Field*, schnek::DeltaTimeDiagnostic> {
  protected:
    typedef schnek::HDFGridDiagnostic<Field, Field*, schnek::DeltaTimeDiagnostic>::IndexType IndexType;
    IndexType getGlobalMin() { return IndexType(0); }
    IndexType getGlobalMax() { return Simulation::getGlobalMax(); }
};

#endif /* DIAGNOSTIC_HPP_ */
