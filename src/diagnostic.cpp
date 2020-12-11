/*
 * diagnostic.cpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#include "diagnostic.hpp"
#include "mpulse.hpp"

FieldDiagnostic::IndexType FieldDiagnostic::getGlobalMin()
{
  return getContext().getSubdivision().getGlobalDomain().getLo();
}

FieldDiagnostic::IndexType FieldDiagnostic::getGlobalMax()
{
  return getContext().getSubdivision().getGlobalDomain().getHi();
}

void FieldDiagnostic::init()
{
  SimulationEntity::init(this);
  schnek::HDFGridDiagnostic<Field, pField, schnek::DeltaTimeDiagnostic>::init();
}

