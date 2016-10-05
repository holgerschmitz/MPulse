
#include "diagnostic.h"
#include <sstream>
#include "process.h"
#include "globals.h"
#include "boundary.h"

DiagnosticInterface::DiagnosticInterface() {
  DiagnosticManager::instance().addDiagnostic(this);
  t=0;
}

void DiagnosticInterface::execute() {
  // perform calculations first
  calculate();

  if (singleOut() && !(Globals::instance().isMaster()) ) return;
  
  //if its the first call of execute() and appending is true open the file
  if ((0==t) && appending()) open(parsedFileName());
  
  //if timestep is a multiple of interval perform output
  if ( (t % interval) == 0 )
  { 
    //if not appending open a new file
    if (!appending()) open(parsedFileName());
    write();
    //if not appending close the file
    if (!appending()) close();
  }
  ++t;
}

ParameterMap* DiagnosticInterface::MakeParamMap(ParameterMap* pm)
{
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["file"] = WParameter(
    new ParameterValue<std::string>(&fname, "")
  );
  (*pm)["append"] = WParameter(
    new ParameterValue<std::string>(&append, "n")
  );
  (*pm)["interval"] = WParameter(
    new ParameterValue<int>(&interval, 100)
  );
  return pm;
}

bool DiagnosticInterface::appending() {
  return 'y' == append[0];
}

std::string DiagnosticInterface::parsedFileName() {
  std::string parsed=fname;
  //look up ID of process
  std::ostringstream comrankstr;
  comrankstr << Globals::instance().getUniqueId();
  std::string comrank = comrankstr.str();
  //look up the global time step
  std::ostringstream tstepstr;
  tstepstr << Process::instance().getTime();
  std::string tstep = tstepstr.str();
  //replace placeholders with the appropiate variables
  size_t pos;
  pos = parsed.find("#p");
  if (pos != std::string::npos) parsed.replace(pos,2,comrank);
  pos = parsed.find("#t");
  if (pos != std::string::npos) parsed.replace(pos,2,tstep);
  
  return parsed;
}

DiagnosticManager *DiagnosticManager::theManager = NULL;

DiagnosticManager::DiagnosticManager() {}

void DiagnosticManager::addDiagnostic(DiagnosticInterface *diag)
{
  diags.push_back(pDiagnosticInterface(diag));
}

void DiagnosticManager::execute()
{
  for (
    DiagList::iterator it=diags.begin();
    it != diags.end();
    ++it
  )
  {
    (*it)->execute();
  }
}

DiagnosticManager& DiagnosticManager::instance()
{
  if (theManager==NULL) theManager = new DiagnosticManager;
  return *theManager;
}


/* meyers singelton should work here:
 * DiagnosticManager& DiagnosticManager::instance()
 * {
 *  	static DiagnosticManager theManager ;
 * 	return theManager;
 * }
 */

