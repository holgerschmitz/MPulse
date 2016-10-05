/*
 * mpulse.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#include "diagnostic.hpp"

#include <schnek/parser.hpp>
#include <schnek/tools/fieldtools.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/util/logger.hpp>

#include <mpi.h>

#include <fstream>
#include <string>
#include <unistd.h>


void MPulse::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addArrayParameter("N", gridSize);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  x_parameters = parameters.addArrayParameter("", x, schnek::BlockParameters::readonly);

  E_parameters[0] = parameters.addParameter("Ex", &(initE[0]), 0.0);
  E_parameters[1] = parameters.addParameter("Ey", &(initE[1]), 0.0);
  E_parameters[2] = parameters.addParameter("Ez", &(initE[2]), 0.0);

  B_parameters[0] = parameters.addParameter("Bx", &(initB[0]), 0.0);
  B_parameters[1] = parameters.addParameter("By", &(initB[1]), 0.0);
  B_parameters[2] = parameters.addParameter("Bz", &(initB[2]), 0.0);

  spaceVars = schnek::pParametersGroup(new schnek::ParametersGroup());
  spaceVars->addArray(x_parameters);
}

void MPulse::registerData()
{
  addData("Ex", Ex);
  addData("Ey", Ey);
  addData("Ez", Ez);

  addData("Bx", Bx);
  addData("By", By);
  addData("Bz", Bz);
}


void MPulse::fillValues()
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

void MPulse::init()
{
  globalMax = gridSize - 2;

  subdivision.init(gridSize, 2);

  Index low  = subdivision.getLo();
  Index high = subdivision.getHi();

  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  innerRange = Range(lowIn, highIn);
  schnek::Range<double, DIMENSION> domainSize(schnek::Array<double, DIMENSION>(lowIn[0],lowIn[1]), schnek::Array<double, DIMENSION>(highIn[0],highIn[1]));
  schnek::Array<bool, DIMENSION> stagger;

  stagger = false;

  Ex = new Grid(lowIn, highIn, domainSize, exStaggerYee, 2);
  Ey = new Grid(lowIn, highIn, domainSize, eyStaggerYee, 2);
  Ez = new Grid(lowIn, highIn, domainSize, ezStaggerYee, 2);

  Bx = new Grid(lowIn, highIn, domainSize, bxStaggerYee, 2);
  By = new Grid(lowIn, highIn, domainSize, byStaggerYee, 2);
  Bz = new Grid(lowIn, highIn, domainSize, bzStaggerYee, 2);

  fillValues();
}

int main (int argc, char** argv) {

  MPI_Init(&argc, &argv);

  try
  {
    schnek::BlockClasses blocks;

    blocks.registerBlock("mpulse").setClass<MPulse>();
    blocks("FieldDiag").setClass<FieldDiagnostic>();

    blocks("mpulse").addChildren("FieldDiag");

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
