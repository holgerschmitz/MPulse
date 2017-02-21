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

class FieldDiagnostic : public HDFGridDiagnostic<Field<double, 2>, Field<double, 2>*, DeltaTimeDiagnostic> {
  protected:
    IndexType getGlobalMin() { return IndexType(0); }
    IndexType getGlobalMax() { return globalMax; }
};

#endif /* DIAGNOSTIC_HPP_ */
