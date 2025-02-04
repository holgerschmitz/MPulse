/*
 * plasmadensity.cpp
 *
 *  Created on: 22 Nov 2024
 *      Author: Holger Schmitz
 */

#include "plasmadensity.hpp"
#include "types.hpp"

#include <schnek/util/logger.hpp>
#include <schnek/grid/domainsubdivision.hpp>
#include <schnek/tools/fieldtools.hpp>

#include <string>
#include <iostream>


void PlasmaDensity::initParameters(schnek::BlockParameters &parameters)
{
  Rho.parameter = parameters.addParameter("Rho", &Rho.value , 0.0);
}

void PlasmaDensity::registerData() {
  addData(fieldName, Rho.field);
}

void PlasmaDensity::fillValues() {
  std::cout << "Filling fields" << std::endl;
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  Vector &x = getContext().getX();
  schnek::Array<schnek::pParameter, DIMENSION> x_parameters = getContext().getXParameter();
  updater.addIndependentArray(x_parameters);

  schnek::fill_field(Rho.field, x, Rho.value, updater, Rho.parameter);
}

void PlasmaDensity::init() {
  schnek::ChildBlock<PlasmaDensity>::init();
  SimulationEntity::init(this);
  const schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();
  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  schnek::Range<double, DIMENSION> domainSize = subdivision.getInnerExtent(getContext().getSize());
  Stagger stagger(false);
  Rho.field.resize(lowIn, highIn, domainSize, stagger, 2);
  fillValues();
}

void NeutralDensity::initParameters(schnek::BlockParameters &parameters)
{
  // Required to call the parent class method to initialize the parameters
  PlasmaDensity::initParameters(parameters);
  parameters.addParameter("A", &A, 1.0);
}

void NeutralDensity::registerData() {
  // Required to call the parent class method to register the data
  PlasmaDensity::registerData();
  addData("A", A);
}