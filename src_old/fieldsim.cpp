#include "fieldsim.h"

#include "fdtd_plain.h"
#include "fdtd_nl.h"
#include "fdtd_disp.h"
#include "fdtd_plrc.h"

#include "globals.h"
#include "process.h"
#include "periodic_bound.h"
#include "mpi_bound.h"

#include "waveinit.h"
#include "gaussinit.h"
#include "pulseinit.h"

#include "freqdiag.h"
#include "followmax.h"



FieldSimulation::FieldSimulation()
{
  
}

FieldSimulation::~FieldSimulation()
{
}

void FieldSimulation::init()
{
  std::cout << "                             *** this->boundary->init() ***" << std::endl;
  this->boundary->init();
      
  dt = Globals::instance().dt();

  std::cout << "                             *** RESIZING ***" << std::endl;
  
  //resize
  resize(this->boundary->RegionLow(),this->boundary->RegionHigh());
  solver->initStorage(this);
    
  //get the global grid spacing
  this->dx = Globals::instance().gridDX();
  this->dy = Globals::instance().gridDY();
  this->dz = Globals::instance().gridDZ();
  
  typedef AllFieldDiag::DiagList::iterator Iter;
  typedef AllFieldDiag::SliceDiagList::iterator SliceIter;
  typedef AllFieldDiag::LineDiagList::iterator LineIter;
  typedef AllFieldDiag::ExtraDiagList::iterator ExtraIter;
  
  for (
    Iter it = fieldDiag.fields.begin();
    it != fieldDiag.fields.end();
    ++it
  )
  {
    (*it)->fetchField(*this);
  }
      
  for (
    SliceIter sit = fieldDiag.slices.begin();
    sit != fieldDiag.slices.end();
    ++sit
  )
  {
    (*sit)->fetchField(*this);
  }
      
  for (
    LineIter lit = fieldDiag.lines.begin();
    lit != fieldDiag.lines.end();
    ++lit
  )
  {
    (*lit)->fetchField(*this);
  }
      
  for (
    ExtraIter xit = fieldDiag.fieldextras.begin();
    xit != fieldDiag.fieldextras.end();
    ++xit
  )
  {
    (*xit)->setStorage(this);
    (*xit)->init();
  }

  if(initializers.empty())
  {
    std::cout << "  NO INITIALIZER SPECIFIED!!!" << std::endl;
    std::cout << "    Starting with zero fields" << std::endl;
  }
  else {
    for (
      FieldSimInitList::iterator iit = initializers.begin();
      iit != initializers.end();
      ++iit
    )
    {
      (*iit)->init(*this);
    }
   
  }
}

void FieldSimulation::execute()
{
  solver->stepScheme(dt);
}

ParameterMap* FieldSimulation::MakeParamMap (ParameterMap* pm) {
  pm = Rebuildable::MakeParamMap(pm);
    
  (*pm)["fdtd-plain"] = WParameter(
      new ParameterRebuild<FDTD_Plain, FieldSolver>(&solver)
  );
  
  (*pm)["fdtd-nl"] = WParameter(
      new ParameterRebuild<FDTD_Nonlinear, FieldSolver>(&solver)
  );
  
  (*pm)["fdtd-disp"] = WParameter(
      new ParameterRebuild<FDTD_Dispersion, FieldSolver>(&solver)
  );
  
  (*pm)["fdtd-plrc-lin"] = WParameter(
      new ParameterRebuild<FDTD_PLRCLin, FieldSolver>(&solver)
  );
  
  (*pm)["fdtd-plrc-nonlin"] = WParameter(
      new ParameterRebuild<FDTD_PLRCNonlin, FieldSolver>(&solver)
  );
  

  (*pm)["single_periodic"] = WParameter(
      new ParameterRebuild<SinglePeriodicBoundary, Boundary>(&boundary)
  );
  
  (*pm)["single_xy_periodic"] = WParameter(
      new ParameterRebuild<SingleXYPeriodicBoundary, Boundary>(&boundary)
  );
  
#ifndef SINGLE_PROCESSOR  
  (*pm)["mpi_slice_periodic"] = WParameter(
      new ParameterRebuild<MPIPeriodicSplitXBoundary, Boundary>(&boundary)
  );
  
  (*pm)["mpi_block_periodic"] = WParameter(
      new ParameterRebuild<MPIPeriodicSplitXYZBoundary, Boundary>(&boundary)
  );
#endif


  (*pm)["wave_init"] = WParameter(
      new ParameterRebuild<PlaneWaveInit, FieldSimInit>(&initializers)
  );
  
  (*pm)["gauss_init"] = WParameter(
      new ParameterRebuild<PlaneGaussPulseInit, FieldSimInit>(&initializers)
  );
  
  (*pm)["pulse_init"] = WParameter(
      new ParameterRebuild<GaussPulseInit, FieldSimInit>(&initializers)
  );
  
  (*pm)["fielddiag"] = WParameter(
      new ParameterRebuild<FieldDiag, FieldDiag>(&fieldDiag.fields)
  );
  
  (*pm)["slicediag"] = WParameter(
      new ParameterRebuild<FieldSliceDiag, FieldSliceDiag>(&fieldDiag.slices)
  );
  
  (*pm)["linediag"] = WParameter(
      new ParameterRebuild<FieldLineDiag, FieldLineDiag>(&fieldDiag.lines)
  );
  
  (*pm)["energy"] = WParameter(
      new ParameterRebuild<FieldEnergyDiag, FieldExtraDiag>(&fieldDiag.fieldextras)
  );
  
  (*pm)["frequency"] = WParameter(
      new ParameterRebuild<FrequencyDiag, FieldExtraDiag>(&fieldDiag.fieldextras)
  );
  
  (*pm)["followmax"] = WParameter(
      new ParameterRebuild<FollowMax, FieldExtraDiag>(&fieldDiag.fieldextras)
  );
  
  return pm;
}
