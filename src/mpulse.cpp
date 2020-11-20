/*
 * mpulse.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#include "diagnostic.hpp"
#include "fdtd_plrc.hpp"
#include "cpml_border.hpp"
#include "sources.hpp"
#include "shortpulseinject.hpp"
#include "plasmacurrent.hpp"

#include "../huerto/electromagnetics/em_fields.hpp"
#include "../huerto/electromagnetics/fieldsolver.hpp"
#include "../huerto/electromagnetics/fdtd/fdtd_plain.hpp"
#include "../huerto/maths/functions/core.hpp"
#include "../huerto/constants.hpp"

#include <schnek/parser.hpp>
#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/tools/fieldtools.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/util/logger.hpp>

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include <mpi.h>

#include <fstream>
#include <string>
#include <unistd.h>

MPulse::MPulse()
{}

void MPulse::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addArrayParameter("N", gridSize, 100);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  parameters.addParameter("cflFactor", &cflFactor, 0.99);

  registerConstants(parameters);
}

void MPulse::initFields()
{
  if (schnek::BlockContainer<EMFields>::childBlocks().empty())
  {
    boost::shared_ptr<EMFields> fields(new EMFields(shared_from_this()));
    Block::addChild(fields);
    fields->registerData();
    fields->preInit();
  }
}

void MPulse::init()
{
  globalMax = gridSize - 1;
  getSubdivision().init(gridSize, 2);

  for (std::size_t i=0; i<DIMENSION; ++i) dx[i] = size[i] / gridSize[i];
  dt = cflFactor*std::min(dx[0],std::min(dx[1],dx[2]))/clight;

  initFields();
}

void MPulse::execute()
{
  for(pFieldSolver f: schnek::BlockContainer<FieldSolver>::childBlocks())
  {
    f->stepSchemeInit(dt);
  }

  time = 0.0;
  schnek::DiagnosticManager::instance().setPhysicalTime(&time);

  while (time<=tMax)
  {
    schnek::DiagnosticManager::instance().execute();

    if (getSubdivision().master()) {
      schnek::Logger::instance().out() <<"Time "<< time << std::endl;
    }

    for(pFieldSolver f: schnek::BlockContainer<FieldSolver>::childBlocks())
    {
      f->stepScheme(dt);
    }

    time += dt;
  }

  schnek::DiagnosticManager::instance().execute();
}

int main (int argc, char** argv) {

  MPI_Init(&argc, &argv);

  try
  {
    schnek::BlockClasses blocks;

    blocks.registerBlock("mpulse").setClass<MPulse>();
    blocks("EMFields").setClass<EMFields>();
    blocks("FDTD_Plain").setClass<FDTD_Plain>();
    blocks("FDTD_PLRC").setClass<FDTD_PLRCLin>();
    blocks("FDTD_PLRC_Nonlinear").setClass<FDTD_PLRCNonlin>();
    blocks("FieldDiag").setClass<FieldDiagnostic>();
    blocks("CPMLBorder").setClass<CPMLBorder>();
    blocks("ShortPulseInject").setClass<ShortPulseInject>();
    blocks("PlaneWaveSource").setClass<PlaneWaveSource>();
    blocks("PlaneGaussSource").setClass<PlaneGaussSource>();

    blocks("PlasmaCurrent").setClass<PlasmaCurrentBlock>();

    blocks("mpulse").addChildren("EMFields")
        ("FDTD_Plain")("FDTD_PLRC")("FDTD_PLRC_Nonlinear")
        ("FieldDiag");

    blocks("FDTD_Plain").addChildren("CPMLBorder")
        ("PlaneWaveSource")("PlaneGaussSource");
    blocks("FDTD_PLRC").addChildren("CPMLBorder")
        ("ShortPulseInject")("PlaneWaveSource")("PlaneGaussSource")
        ("PlasmaCurrent");
    blocks("FDTD_PLRC_Nonlinear").addChildren("CPMLBorder")
        ("ShortPulseInject")("PlaneWaveSource")("PlaneGaussSource")
        ("PlasmaCurrent");

    std::ifstream in("mpulse.setup");
    if (!in) throw std::string("Could not open file 'mpulse.setup'");

    schnek::Parser P("mpulse", "mpulse", blocks);
    registerCMath(P.getFunctionRegistry());

    //P.getFunctionRegistry().registerFunction("random",randomRange);

    schnek::pBlock application = P.parse(in);

    MPulse &mpulse = dynamic_cast<MPulse&>(*application);
    mpulse.initAll();

    if (mpulse.getSubdivision().master())
    {
      std::ofstream referencesText("information.tex");
      std::ofstream referencesBib("references.bib");

      schnek::LiteratureManager::instance().writeInformation(referencesText,"references.bib");
      schnek::LiteratureManager::instance().writeBibTex(referencesBib);
      referencesText.close();
      referencesBib.close();
    }

    mpulse.execute();
  }
  catch (schnek::ParserError &e)
  {
    std::cerr << "Parse error in " << e.getFilename() << " at line "
        << e.getLine() << ": " << e.message << std::endl;
    return -1;
  }
  catch (schnek::VariableNotInitialisedException &e)
  {
    std::cerr << "Variable was not initialised: " << e.getVarName() << std::endl;
    return -1;
  }
  catch (schnek::EvaluationException &e)
  {
    std::cerr << "Error in evaluation: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (schnek::VariableNotFoundException &e)
  {
    std::cerr << "Error: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (SchnekException &e)
  {
    std::cerr << "An error occured" << std::endl;
    return -1;
  }
  catch (std::string &err)
  {
    std::cerr << "FATAL ERROR: " << err << std::endl;
    return -1;
  }

  MPI_Finalize();

  return 0;
}

// end of main
