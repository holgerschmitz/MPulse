#include "process.h"
#include "diagnostic.h"

Process *Process::process;

Process::~Process() 
{
  if (fieldSim) delete fieldSim;
  if (globals) delete globals;
}

void Process::init() 
{
  fieldSim->init();
}

void Process::run() 
{
  int T = Globals::instance().totalTime(); 

  //display time step and call execute for all species +field 
  for (time=0; time<T; ++time) 
  {
    if (Globals::instance().isMaster())
      std::cout << "Cycle " << time << std::endl << std::flush;

    fieldSim->execute();
    
    DiagnosticManager::instance().execute();
  }
}

ParameterMap* Process::MakeParamMap (ParameterMap* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["globals"] 
      = WParameter(new ParameterRebuild<Globals, Globals>(&globals));
      
  (*pm)["fieldsim"]
      = WParameter(new ParameterRebuild<FieldSimulation, FieldSimulation>(&fieldSim));

  return pm;
}


std::string Process::Rebuild(std::istream& in)
{
  std::cerr << "Rebuilding Process ...\n";
  Rebuildable::Rebuild(in);
  
  if (NULL == globals) {
    std::cerr << "No globals specified! Must exit!\n";
    exit(-1);
  }
  if (NULL == fieldSim) {
    std::cerr << "No field simulation specified! Must exit!\n";
    exit(-1);
  }
  return "";
}

