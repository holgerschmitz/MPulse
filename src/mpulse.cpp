/*
 * mpulse.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#include "mpulse.hpp"
#include "diagnostic.hpp"
#include "fdtd.hpp"

#include <schnek/parser.hpp>
#include <schnek/tools/fieldtools.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/util/logger.hpp>

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <mpi.h>

#include <fstream>
#include <string>
#include <unistd.h>

Simulation *Simulation::instance = 0;

Simulation::Simulation() {
  instance = this;
}

void Simulation::initParameters(schnek::BlockParameters &parameters) {
  parameters.addConstant("pi", M_PI);
  parameters.addConstant("clight", clight);
  parameters.addArrayParameter("N", gridSize, 100);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  parameters.addParameter("cflFactor", &cflFactor, 0.99);
  x_parameters = parameters.addArrayParameter("", x, schnek::BlockParameters::readonly);

  E_parameters = parameters.addArrayParameter("E", initE, 0.0);
  B_parameters = parameters.addArrayParameter("B", initB, 0.0);

  spaceVars = schnek::pParametersGroup(new schnek::ParametersGroup());
  spaceVars->addArray(x_parameters);
}


void Simulation::registerData()
{
  dx[0] = size[0] / gridSize[0];
  dx[1] = size[1] / gridSize[1];
  dx[2] = size[2] / gridSize[2];

  subdivision.init(gridSize, 2);

  IndexType low  = subdivision.getLo();
  IndexType high = subdivision.getHi();

  IndexType lowIn  = subdivision.getInnerLo();
  IndexType highIn = subdivision.getInnerHi();

  innerRange = Range(lowIn, highIn);
  schnek::Range<double, 3> domainSize(
      schnek::Array<double, 3>(dx[0]*lowIn[0],dx[1]*lowIn[1],dx[2]*lowIn[2]),
      schnek::Array<double, 3>(dx[0]*highIn[0],dx[1]*highIn[1],dx[2]*highIn[2]));
  schnek::Array<bool, 3> stagger;

  stagger = false;

  Ex = boost::make_shared<Field>(lowIn, highIn, domainSize, exStaggerYee, 2);
  Ey = boost::make_shared<Field>(lowIn, highIn, domainSize, eyStaggerYee, 2);
  Ez = boost::make_shared<Field>(lowIn, highIn, domainSize, ezStaggerYee, 2);

  Bx = boost::make_shared<Field>(lowIn, highIn, domainSize, bxStaggerYee, 2);
  By = boost::make_shared<Field>(lowIn, highIn, domainSize, byStaggerYee, 2);
  Bz = boost::make_shared<Field>(lowIn, highIn, domainSize, bzStaggerYee, 2);

  std::cout << "Register data\n";

  addData("Ex", Ex);
  addData("Ey", Ey);
  addData("Ez", Ez);

  addData("Bx", Bx);
  addData("By", By);
  addData("Bz", Bz);
}

void Simulation::fillValues()
{
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  updater.addIndependentArray(x_parameters);

  schnek::fill_field(*Ex, x, initE[0], updater, E_parameters[0]);
  schnek::fill_field(*Ey, x, initE[1], updater, E_parameters[1]);
  schnek::fill_field(*Ez, x, initE[2], updater, E_parameters[2]);

  schnek::fill_field(*Bx, x, initB[0], updater, B_parameters[0]);
  schnek::fill_field(*By, x, initB[1], updater, B_parameters[1]);
  schnek::fill_field(*Bz, x, initB[2], updater, B_parameters[2]);
}

void Simulation::init() {
//  Field *Ex, *Ey, *Ez;
//  Field *Bx, *By, *Bz;
//
//  retrieveData("Ex", Ex);
//  retrieveData("Ey", Ey);
//  retrieveData("Ez", Ez);
//
//  retrieveData("Bx", Bx);
//  retrieveData("By", By);
//  retrieveData("Bz", Bz);

  globalMax = gridSize - 1;

  subdivision.init(gridSize, 2);

//  for (int i=0; i<3; ++i) dx[i] = size[i] / gridSize[i];
  dt = cflFactor*std::min(dx[0],std::min(dx[1],dx[2]))/clight;

  IndexType lowIn  = subdivision.getInnerLo();
  IndexType highIn = subdivision.getInnerHi();

  fillValues();

//  schnek::Range<double,3> domainSize = subdivision.getInnerExtent(size);
//
//  Ex->resize(lowIn, highIn, domainSize, exStaggerYee, 2);
//  Ey->resize(lowIn, highIn, domainSize, eyStaggerYee, 2);
//  Ez->resize(lowIn, highIn, domainSize, ezStaggerYee, 2);
//
//  Bx->resize(lowIn, highIn, domainSize, bxStaggerYee, 2);
//  By->resize(lowIn, highIn, domainSize, byStaggerYee, 2);
//  Bz->resize(lowIn, highIn, domainSize, bzStaggerYee, 2);
//
//  schnek::pBlockVariables blockVars = getVariables();
//  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));
//
//  schnek::DependencyUpdater updater(depMap);
//
//  updater.addIndependentArray(x_parameters);
//  schnek::fill_field(*Ex, x, initE[0], updater, E_parameters[0]);
//  schnek::fill_field(*Ey, x, initE[1], updater, E_parameters[1]);
//  schnek::fill_field(*Ez, x, initE[2], updater, E_parameters[2]);
//
//  schnek::fill_field(*Bx, x, initB[0], updater, B_parameters[0]);
//  schnek::fill_field(*By, x, initB[1], updater, B_parameters[1]);
//  schnek::fill_field(*Bz, x, initB[2], updater, B_parameters[2]);
}

void Simulation::execute() {
  double time = 0.0;

  schnek::DiagnosticManager::instance().setPhysicalTime(&time);
  schnek::DiagnosticManager::instance().execute();
  while (time<=tMax) {
    double dtAdjust = schnek::DiagnosticManager::instance().adjustDeltaT(dt);

    if (subdivision.master())
      schnek::Logger::instance().out() <<"Time "<< time << std::endl;

    BOOST_FOREACH(pFDTDSolver f, childBlocks()) {
      f->stepScheme(dtAdjust);
    }

    time += dtAdjust;
    schnek::DiagnosticManager::instance().execute();
  }

}

int main (int argc, char** argv) {

  MPI_Init(&argc, &argv);

  try {
    schnek::BlockClasses blocks;

    blocks.registerBlock("simulation").setClass<Simulation>();
    blocks("FieldDiag").setClass<FieldDiagnostic>();
    blocks("FDTDSolver").setClass<FDTDSolver>();

    blocks("simulation").addChildren("FieldDiag")("FDTDSolver");

    std::ifstream in("simulation.setup");
    if (!in) throw std::string("Could not open file 'simulation.setup'");

    schnek::Parser P("simulation", "simulation", blocks);
    registerCMath(P.getFunctionRegistry());

    //P.getFunctionRegistry().registerFunction("random",randomRange);

    schnek::pBlock application = P.parse(in);

    Simulation &simulation = dynamic_cast<Simulation&>(*application);
    simulation.initAll();

    if (simulation.getSubdivision().master()) {
      std::ofstream referencesText("information.tex");
      std::ofstream referencesBib("references.bib");

      schnek::LiteratureManager::instance().writeInformation(referencesText,"references.bib");
      schnek::LiteratureManager::instance().writeBibTex(referencesBib);
      referencesText.close();
      referencesBib.close();
    }

    simulation.execute();
  }
  catch (schnek::ParserError &e) {
    std::cerr << "Parse error in " << e.getFilename() << " at line "
        << e.getLine() << ": " << e.message << std::endl;
    return -1;
  }
  catch (schnek::VariableNotInitialisedException &e) {
    std::cerr << "Variable was not initialised: " << e.getVarName() << std::endl;
    return -1;
  }
  catch (schnek::EvaluationException &e) {
    std::cerr << "Error in evaluation: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (schnek::VariableNotFoundException &e) {
    std::cerr << "Error: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (SchnekException &e) {
    std::cerr << "An error occured" << std::endl;
    return -1;
  }
  catch (std::string &err) {
    std::cerr << "FATAL ERROR: " << err << std::endl;
    return -1;
  }

  MPI_Finalize();

  return 0;
}

// end of main
