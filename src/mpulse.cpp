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
  parameters.addArrayParameter("N", gridSize, 100);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  parameters.addParameter("cflFactor", &cflFactor, 0.5);

  x_par = parameters.addArrayParameter("", x, schnek::BlockParameters::readonly);
  E_par = parameters.addArrayParameter("E", initE, 0.0);
  B_par = parameters.addArrayParameter("B", initB, 0.0);

  spaceVars = schnek::pParametersGroup(new schnek::ParametersGroup());
  spaceVars->addArray(x_par);

  parameters.addConstant("pi", PI);
  parameters.addConstant("clight", clight);
  parameters.addConstant("mu0", mu_0);
  parameters.addConstant("eps0", eps_0);
}

void Simulation::init() {
  globalMax = gridSize-1;
  for (int i=0; i<3; ++i) dx[i] = size[i] / gridSize[i];

  subdivision.init(gridSize, 2);

  IndexType lowIn  = subdivision.getInnerLo();
  IndexType highIn = subdivision.getInnerHi();

  schnek::Range<double, 3> domainSize = subdivision.getInnerExtent(size);

  Field *Ex, *Ey, *Ez;
  Field *Bx, *By, *Bz;

  retrieveData("Ex", Ex);
  retrieveData("Ey", Ey);
  retrieveData("Ez", Ez);

  retrieveData("Bx", Bx);
  retrieveData("By", By);
  retrieveData("Bz", Bz);

  Ex->resize(lowIn, highIn, domainSize, exStaggerYee, 2);
  Ey->resize(lowIn, highIn, domainSize, eyStaggerYee, 2);
  Ez->resize(lowIn, highIn, domainSize, ezStaggerYee, 2);

  Bx->resize(lowIn, highIn, domainSize, bxStaggerYee, 2);
  By->resize(lowIn, highIn, domainSize, byStaggerYee, 2);
  Bz->resize(lowIn, highIn, domainSize, bzStaggerYee, 2);

  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));
  schnek::DependencyUpdater updater(depMap);

  updater.addIndependentArray(x_par);

  schnek::fill_field(*Ex, x, initE[0], updater, E_par[0]);
  schnek::fill_field(*Ey, x, initE[1], updater, E_par[1]);
  schnek::fill_field(*Ez, x, initE[2], updater, E_par[2]);

  schnek::fill_field(*Bx, x, initB[0], updater, B_par[0]);
  schnek::fill_field(*By, x, initB[1], updater, B_par[1]);
  schnek::fill_field(*Bz, x, initB[2], updater, B_par[2]);
}

void Simulation::execute() {
  time = 0.0;

  schnek::DiagnosticManager::instance().setPhysicalTime(&time);
  schnek::DiagnosticManager::instance().execute();

  double minDx = std::min(dx[0], dx[1]);
  double dtBase = cflFactor*minDx/clight;

  dt = schnek::DiagnosticManager::instance().adjustDeltaT(dtBase);
  BOOST_FOREACH(pFieldSolver f, schnek::BlockContainer<FieldSolver>::childBlocks()) {
    f->stepSchemeInit(dt);
  }

  while (time<=tMax) {
    dt = schnek::DiagnosticManager::instance().adjustDeltaT(dtBase);

    if (subdivision.master())
      schnek::Logger::instance().out() <<"Time "<< time << std::endl;

    BOOST_FOREACH(pFieldSolver f, schnek::BlockContainer<FieldSolver>::childBlocks()) {
      f->stepScheme(dt);
    }

    time += dt;
    schnek::DiagnosticManager::instance().execute();
  }

}

int main (int argc, char** argv) {

  MPI_Init(&argc, &argv);

  try {
    schnek::BlockClasses blocks;

    blocks.registerBlock("simulation").setClass<Simulation>();
    blocks("FieldDiag").setClass<FieldDiagnostic>();
    blocks("FDTD").setClass<FieldSolver>();

    blocks("simulation").addChildren("FieldDiag")("FDTD");

    std::ifstream in("simulation.setup");
    if (!in) throw std::string("Could not open file 'simulation.setup'");

    schnek::Parser P("simulation", "simulation", blocks);
    registerCMath(P.getFunctionRegistry());

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
