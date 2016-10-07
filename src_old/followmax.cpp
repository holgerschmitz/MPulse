#include "followmax.h"
#include "process.h"
#include "globals.h"
#include "boundary.h"

FollowMax::FollowMax() : storage(0) {}
	  
FollowMax::~FollowMax() {}
	  
void FollowMax::setStorage(Storage *storage_)
{
  storage = storage_;
}

void FollowMax::open(const std::string &fname){
  output.open(fname.c_str());
  output.precision(10);
}

	  
void FollowMax::write()
{
  int time = Process::instance().getTime();
  if (time<start) return;
  output << time << " " << position << std::endl;
}

void FollowMax::close()
{
  output.close();
}

void FollowMax::init()
{
  std::cerr << "FollowMax::init()\n";

  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();
  GridIndex glow = Globals::instance().gridLow();
  GridIndex ghigh = Globals::instance().gridHigh();
  GridIndex pos(posx, posy, posz);
  
  switch (direction[0]) {
    case 'y':
      dim = 1;
      trans1 = 0;
      trans2 = 2;
      break;
    case 'z':
      dim = 2;
      trans1 = 0;
      trans2 = 1;
      break;
    case 'x':
    default :
      dim = 0;
      trans1 = 1;
      trans2 = 2;
      break;
  }
  
  this->low = low[dim];       
  this->high = high[dim];
  
  iposition = pos[dim];
  position = iposition;
 
  active = (pos[trans1]>low[trans1]) && (pos[trans1]<high[trans1])
    && (pos[trans2]>low[trans2]) && (pos[trans2]<high[trans2]);
  
}

void FollowMax::calculate()
{
  int time = Process::instance().getTime();
  if (time<start) return;
  if (!storage->hasGrid(field)) return;
  
  DataGrid &G = storage->getGrid(field);
  GridIndex pos(posx, posy, posz);
  
  bool maxfound;
  double mid, left, right;
  
  do {
    int pm = iposition - 1;
    int pp = iposition + 1;
    if (pm<=low) pm = high-1;
    if (pp>=high) pp = low+1;
    
    pos[dim] = iposition;
    mid = G(pos[0],pos[1],pos[2]);
    pos[dim] = pm;
    left = G(pos[0],pos[1],pos[2]);
    pos[dim] = pp;
    right = G(pos[0],pos[1],pos[2]);
    
    if (right>mid) {
      iposition = pp;
      maxfound = false;
    } else if (left>mid) {
      iposition = pm;
      maxfound = false;
    } else {
      maxfound = true;
    }
    
  } while (!maxfound);
  
  double offset = (left - right)/(2 * (left - 2 * mid + right)); // maximum of quadradic interpolation
  position = iposition + offset;

}

ParameterMap* FollowMax::MakeParamMap (ParameterMap* pm) {
  pm = FieldExtraDiag::MakeParamMap(pm);
      
  (*pm)["field"] = WParameter(new ParameterValue<std::string>(&field, "Ex"));
  (*pm)["dim"] = WParameter(new ParameterValue<std::string>(&direction,"x"));
  (*pm)["posx"] = WParameter(new ParameterValue<int>(&posx,5));
  (*pm)["posy"] = WParameter(new ParameterValue<int>(&posy,5));
  (*pm)["posz"] = WParameter(new ParameterValue<int>(&posz,5));
  (*pm)["start"] = WParameter(new ParameterValue<int>(&start,0));

  return pm;
}
