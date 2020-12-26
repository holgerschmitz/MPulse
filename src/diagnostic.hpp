/*
 * diagnostic.hpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_DIAGNOSTIC_HPP_
#define MPULSE_DIAGNOSTIC_HPP_

#include "mpulse.hpp"

#include "../huerto/diagnostic/field_diagnostic.hpp"
#include "../huerto/diagnostic/slice_diagnostic.hpp"

typedef FieldDiagnostic<Field, pField, schnek::DeltaTimeDiagnostic> MPulseFieldDiagnostic;
typedef GridSliceDiagnostic<Grid, pGrid, schnek::IntervalDiagnostic> SliceDiagnostic;


#endif /* SRC_DIAGNOSTIC_HPP_ */
