/*
 * mpulse.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#include "diagnostic.hpp"
#include "em_fields.hpp"
#include "fieldsolver.hpp"
#include "fdtd_plain.hpp"
#include "fdtd_plrc.hpp"
#include "cpml_border.hpp"
#include "sources.hpp"
#include "shortpulseinject.hpp"
#include "plasmacurrent.hpp"

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


MPulse *MPulse::instance;

MPulse::MPulse()
{
  instance = this;
}

void MPulse::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addArrayParameter("N", gridSize, 100);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  parameters.addParameter("cflFactor", &cflFactor, 0.99);
  x_parameters = parameters.addArrayParameter("", x, schnek::BlockParameters::readonly);

  E_parameters = parameters.addArrayParameter("E", initE, 0.0);
  B_parameters = parameters.addArrayParameter("B", initB, 0.0);

  spaceVars = schnek::pParametersGroup(new schnek::ParametersGroup());
  spaceVars->addArray(x_parameters);

  parameters.addConstant("pi", PI);
  parameters.addConstant("clight", clight);
  parameters.addConstant("me", mass_e);
  parameters.addConstant("mp", mass_p);
  parameters.addConstant("e", unit_charge);
  parameters.addConstant("mu0", mu_0);
  parameters.addConstant("eps0", eps_0);
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

  retrieveData("Ex", pEx);
  retrieveData("Ey", pEy);
  retrieveData("Ez", pEz);

  retrieveData("Bx", pBx);
  retrieveData("By", pBy);
  retrieveData("Bz", pBz);

  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  innerRange = Range(lowIn, highIn);
  schnek::Range<double, DIMENSION> domainSize(
      schnek::Array<double, DIMENSION>(dx[0]*lowIn[0],dx[1]*lowIn[1],dx[2]*lowIn[2]),
      schnek::Array<double, DIMENSION>(dx[0]*(highIn[0]+1),dx[1]*(highIn[1]+1),dx[2]*(highIn[2]+1)));

  pEx->resize(lowIn, highIn, domainSize, exStaggerYee, 2);
  pEy->resize(lowIn, highIn, domainSize, eyStaggerYee, 2);
  pEz->resize(lowIn, highIn, domainSize, ezStaggerYee, 2);

  pBx->resize(lowIn, highIn, domainSize, bxStaggerYee, 2);
  pBy->resize(lowIn, highIn, domainSize, byStaggerYee, 2);
  pBz->resize(lowIn, highIn, domainSize, bzStaggerYee, 2);
}

void MPulse::fillValues()
{
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  updater.addIndependentArray(x_parameters);
  schnek::fill_field(*pEx, x, initE[0], updater, E_parameters[0]);
  schnek::fill_field(*pEy, x, initE[1], updater, E_parameters[1]);
  schnek::fill_field(*pEz, x, initE[2], updater, E_parameters[2]);

  schnek::fill_field(*pBx, x, initB[0], updater, B_parameters[0]);
  schnek::fill_field(*pBy, x, initB[1], updater, B_parameters[1]);
  schnek::fill_field(*pBz, x, initB[2], updater, B_parameters[2]);
}

void MPulse::init()
{
  globalMax = gridSize - 1;
  subdivision.init(gridSize, 2);

  for (std::size_t i=0; i<DIMENSION; ++i) dx[i] = size[i] / gridSize[i];
  dt = cflFactor*std::min(dx[0],std::min(dx[1],dx[2]))/clight;

  initFields();
  fillValues();
}

void MPulse::execute()
{
  BOOST_FOREACH(pFieldSolver f, schnek::BlockContainer<FieldSolver>::childBlocks())
  {
    f->stepSchemeInit(dt);
  }

  time = 0.0;
  schnek::DiagnosticManager::instance().setPhysicalTime(&time);

  while (time<=tMax)
  {
    schnek::DiagnosticManager::instance().execute();

    if (subdivision.master())
      schnek::Logger::instance().out() <<"Time "<< time << std::endl;

      BOOST_FOREACH(pFieldSolver f, schnek::BlockContainer<FieldSolver>::childBlocks())
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
