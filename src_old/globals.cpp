#include "globals.h"

//-----------------------------------------------------------------------------
Globals *Globals::globals;
int Globals::Argc;
char **Globals::Argv;
//-----------------------------------------------------------------------------
Globals::Globals() {
  globals = this;
  initialized = false;
  IsRestart = false;
  master = true;
  uniqueId = 0;
}
//-----------------------------------------------------------------------------
void Globals::init()
{
  initialized = true;
}
//-----------------------------------------------------------------------------
ParameterMap* Globals::MakeParamMap (ParameterMap* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["grid-x"] = WParameter(new ParameterValue<int>(&GridX,128));
  (*pm)["grid-y"] = WParameter(new ParameterValue<int>(&GridY,128));
  (*pm)["grid-z"] = WParameter(new ParameterValue<int>(&GridZ,128));
  
  (*pm)["grid-dx"] = WParameter(new ParameterValue<double>(&GridDX,0.1));
  (*pm)["grid-dy"] = WParameter(new ParameterValue<double>(&GridDY,0.1));
  (*pm)["grid-dz"] = WParameter(new ParameterValue<double>(&GridDZ,0.1));

  (*pm)["dt"] = WParameter(new ParameterValue<double>(&DT,1e-2));

  (*pm)["T-total"] = WParameter(new ParameterValue<int>(&TotalTime,100000));
  
  return pm;
}   

//-----------------------------------------------------------------------------

std::string Globals::Rebuild(std::istream& in) 
{
  std::string strToken = Rebuildable::Rebuild(in);

  GridLow = GridIndex(0,0,0);
  GridHigh = GridIndex(GridX+1,GridY+1,GridZ+1);

  return strToken;
}
