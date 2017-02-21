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

IndexType globalMax;
MPICartSubdivision<Field<double, 2> > *subdivision;
Array<double, 2> dx;

void SimulationBlock::initParameters(BlockParameters &parameters) {
  parameters.addArrayParameter("N", gridSize, 100);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  parameters.addParameter("cflFactor", &cflFactor, 0.5);

  parameters.addConstant("pi", PI);
  parameters.addConstant("clight", clight);
  parameters.addConstant("mu0", mu_0);
  parameters.addConstant("eps0", eps_0);

}

void SimulationBlock::preInit() {
  globalMax = gridSize-1;
  for (int i=0; i<2; ++i) dx[i] = size[i] / gridSize[i];

  subdivision.init(gridSize, ghostCells);
  ::subdivision = &subdivision;

  Array<int, 2> lo  = subdivision.getInnerLo();
  Array<int, 2> hi = subdivision.getInnerHi();

  Range<double, 2> physRange = subdivision.getInnerExtent(size);

  Field<double, 2> *Ex, *Ey, *Ez;
  Field<double, 2> *Bx, *By, *Bz;

  Array<bool, 2> exStaggerYee(true,  false);
  Array<bool, 2> eyStaggerYee(false, true);
  Array<bool, 2> ezStaggerYee(false, false);

  Array<bool, 2> bxStaggerYee(false, true);
  Array<bool, 2> byStaggerYee(true,  false);
  Array<bool, 2> bzStaggerYee(true,  true);

  retrieveData("Ex", Ex);
  retrieveData("Ey", Ey);
  retrieveData("Ez", Ez);

  retrieveData("Bx", Bx);
  retrieveData("By", By);
  retrieveData("Bz", Bz);

  Ex->resize(lo, hi, physRange, exStaggerYee, ghostCells);
  Ey->resize(lo, hi, physRange, eyStaggerYee, ghostCells);
  Ez->resize(lo, hi, physRange, ezStaggerYee, ghostCells);

  Bx->resize(lo, hi, physRange, bxStaggerYee, ghostCells);
  By->resize(lo, hi, physRange, byStaggerYee, ghostCells);
  Bz->resize(lo, hi, physRange, bzStaggerYee, ghostCells);

}

void SimulationBlock::execute() {
  time = 0.0;

  DiagnosticManager::instance().setPhysicalTime(&time);
  DiagnosticManager::instance().execute();

  double minDx = std::min(dx[0], dx[1]);
  dt = cflFactor*minDx/clight;

  BOOST_FOREACH(boost::shared_ptr<FieldSolver> f, childBlocks()) {
    f->stepSchemeInit(dt);
  }

  while (time<=tMax) {
    if (subdivision.master())
      std::cout <<"Time "<< time << std::endl;

    BOOST_FOREACH(boost::shared_ptr<FieldSolver> f, childBlocks()) {
      f->stepScheme(dt);
    }

    time += dt;
    DiagnosticManager::instance().execute();
  }
}

int main (int argc, char** argv) {

  MPI_Init(&argc, &argv);

  try {
    BlockClasses blocks;

    blocks("simulation").setClass<SimulationBlock>();
    blocks("Diagnostic").setClass<FieldDiagnostic>();
    blocks("FDTD").setClass<FieldSolver>();

    blocks("simulation").addChildren("Diagnostic")("FDTD");

    std::ifstream in("simulation.setup");
    if (!in) throw std::string("Could not open file 'simulation.setup'");

    Parser P("simulation", "simulation", blocks);
    registerCMath(P.getFunctionRegistry());

    pBlock application = P.parse(in);

    SimulationBlock &simulation = dynamic_cast<SimulationBlock&>(*application);
    simulation.initAll();

    if (subdivision->master()) {
      std::ofstream referencesText("information.tex");
      std::ofstream referencesBib("references.bib");

      LiteratureManager::instance().writeInformation(referencesText,"references.bib");
      LiteratureManager::instance().writeBibTex(referencesBib);
      referencesText.close();
      referencesBib.close();
    }

    simulation.execute();
  }
  catch (ParserError &e) {
    std::cerr << "Parse error in " << e.getFilename() << " at line "
        << e.getLine() << ": " << e.message << std::endl;
    return -1;
  }
  catch (VariableNotInitialisedException &e) {
    std::cerr << "Variable was not initialised: " << e.getVarName() << std::endl;
    return -1;
  }
  catch (EvaluationException &e) {
    std::cerr << "Error in evaluation: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (VariableNotFoundException &e) {
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
