#include "freqdiag.h"
#include "process.h"
#include "globals.h"
#include "boundary.h"

FrequencyDiag::FrequencyDiag() : storage(0) {}
	  
FrequencyDiag::~FrequencyDiag() {}
	  
void FrequencyDiag::setStorage(Storage *storage_)
{
  storage = storage_;
}

void FrequencyDiag::open(const std::string &fname){
  output.open(fname.c_str());
}

	  
void FrequencyDiag::write()
{
  if (count>lastcount)
  {
    lastcount = count;
    int time = Process::instance().getTime();
    output << time << " " << frequency << std::endl;
//    output << time << " " << lastval << std::endl;
  }
}

void FrequencyDiag::close()
{
  output.close();
}

void FrequencyDiag::init()
{
  std::cerr << "FrequencyDiag::init()\n";
  dt=Globals::instance().dt(); 
  count = -1;
  lastcount = 0;
  lastval = 0;
}

void FrequencyDiag::calculate()
{
  if (!storage->hasGrid(field)) return;
  
  DataGrid &G = storage->getGrid(field);
  double val = 0.15*G(x,y,z); // + 0.85*lastval;
  
  if ((val*lastval<=0) && (fabs(val)*fabs(lastval)>0))
  {
    double d = lastval/(lastval-val);
    double time = (Process::instance().getTime() - 1 + d)*dt;
    if (count<0)
    {
      firstzero = time;
      count = 0;
    }
    else
    {
      ++count;
      frequency = (0.5*count) / (time-firstzero);
    }
    
  }
  lastval = val;

}

ParameterMap* FrequencyDiag::MakeParamMap (ParameterMap* pm) {
  pm = FieldExtraDiag::MakeParamMap(pm);
      
  (*pm)["field"] = WParameter(
    new ParameterValue<std::string>(&field, "Ex")
  );

  (*pm)["x"] = WParameter(
    new ParameterValue<int>(&x, 1)
  );

  (*pm)["y"] = WParameter(
    new ParameterValue<int>(&y, 1)
  );

  (*pm)["z"] = WParameter(
    new ParameterValue<int>(&z, 1)
  );

  return pm;
}
