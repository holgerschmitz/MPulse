#include "mpulse.h"
#include "incsource.h"
#include "fieldsolver.h"
#include "globals.h"
#include "storage.h"
#include "hdfstream.h"

#include <vector>
#include <sstream>

//===============================================================
//==========  IncidentSource
//===============================================================


void IncidentSource::initCurrents(Storage *storage, FieldSolver *solver)
{ 
  if (needCurrent(north))
  {
    solver->addCurrent(
      makeECurrent(distance, north)
    );
    solver->addMagCurrent(
      makeHCurrent(distance, north)
    );
  }
  
  if (needCurrent(south))
  {
    solver->addCurrent(
      makeECurrent(distance, south)
    );
    solver->addMagCurrent(
      makeHCurrent(distance, south)
    );
  }
  
  if (needCurrent(east))
  {
    solver->addCurrent(
      makeECurrent(distance, east)
    );
    solver->addMagCurrent(
      makeHCurrent(distance, east)
    );
  }
  
  if (needCurrent(west))
  {
    solver->addCurrent(
      makeECurrent(distance, west)
    );
    solver->addMagCurrent(
      makeHCurrent(distance, west)
    );
  }
  
  if (needCurrent(up))
  {
    solver->addCurrent(
      makeECurrent(distance, up)
    );
    solver->addMagCurrent(
      makeHCurrent(distance, up)
    );
  }
  
  if (needCurrent(down))
  {
    solver->addCurrent(
      makeECurrent(distance, down)
    );
    solver->addMagCurrent(
      makeHCurrent(distance, down)
    );
  }
}


ParameterMap* IncidentSource::MakeParamMap (ParameterMap* pm)
{
  pm = CurrentFactory::MakeParamMap(pm);
  
  (*pm)["d"] = WParameter(new ParameterValue<int>(&this->distance,15));
  return pm;
}

//===============================================================
//==========  IncidentSourceCurrent
//===============================================================


IncidentSourceCurrent::IncidentSourceCurrent(int distance_, Direction dir_, bool isH_)
  : distance(distance_), dir(dir_), isH(isH_)
{
  dx = Globals::instance().gridDX();
  dy = Globals::instance().gridDY();
  dz = Globals::instance().gridDZ();
  dt = Globals::instance().dt();

  switch (dir)
  {
    case east:  
    case west:  dim = 0;
                transverse1 = 1;
                transverse2 = 2;
                pJ[0] = pJy;
                pJ[1] = pJz;
                dX[0] = dy;
                dX[1] = dz;
                dX[2] = dx;
                break;
    case north: 
    case south: dim = 1;
                transverse1 = 2;
                transverse2 = 0;
                pJ[0] = pJz;
                pJ[1] = pJx;
                dX[0] = dz;
                dX[2] = dx;
                dX[1] = dy;
                break;
    case up:     
    case down:  dim = 2;
                transverse1 = 0;
                transverse2 = 1;
                pJ[0] = pJx;
                pJ[1] = pJy;
                dX[0] = dx;
                dX[1] = dy;
                dX[2] = dz;
                break;
  }
  
  reverse = ( (dir==east) || (dir==north) || (dir==up) );
}

//===============================================================
//==========  SideInject
//===============================================================


IncidentSourceCurrent *SideInject::makeECurrent(int distance_, Direction dir_)
{
  typedef IncidentSourceECurrent<SideInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(filename, amp, eps, toff, blocks);
  return cur;
}

IncidentSourceCurrent *SideInject::makeHCurrent(int distance_, Direction dir_)
{
  typedef IncidentSourceHCurrent<SideInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(filename, amp, eps, toff, blocks);
  return cur;
}

bool SideInject::needCurrent(Direction dir_)
{
  return (dir_ == down);
}
    
ParameterMap* SideInject::MakeParamMap (ParameterMap* pm)
{
  pm = IncidentSource::MakeParamMap(pm);
  
  (*pm)["file"] = WParameter(new ParameterValue<std::string>(&this->filename,""));
  (*pm)["amp"] = WParameter(new ParameterValue<double>(&this->amp,1));
  (*pm)["eps"] = WParameter(new ParameterValue<double>(&this->eps,1));
  (*pm)["toff"] = WParameter(new ParameterValue<int>(&this->toff,0));
  (*pm)["blocks"] = WParameter(new ParameterValue<int>(&this->blocks,1));
  return pm;
}


SideInjectSourceFunc::SideInjectSourceFunc(Direction dir_, bool isH_)
  : dir(dir_), isH(isH_)
{}

void SideInjectSourceFunc::setParam(std::string filename_, double amp_, double eps_, int toff_, int blocks_)
{
  filename = filename_;
  amp = amp_;
  eps = eps_;
  toff = toff_;
  blocks = blocks_;

  blockCount = 0;
  blockToff = toff;
  Nt = 1;
  active = false;
  std::cerr << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  std::cerr << "AMPLITUDE " << amp << std::endl;
}

void SideInjectSourceFunc
    ::initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz)
{
  switch (dir)
  {
    case east:  
    case west:  dim = 0;
                transverse1 = 1;
                transverse2 = 2;
                break;
    case north: 
    case south: dim = 1;
                transverse1 = 2;
                transverse2 = 0;
                break;
    case up:     
    case down:  dim = 2;
                transverse1 = 0;
                transverse2 = 1;
                break;
  }
}

void SideInjectSourceFunc::setTime(int Time)
{
  if (Time < (toff + blockCount*Nt)) return;
  if (Time >= (toff + blocks*Nt))
  {
    active = false;
    F1.resize(GridIndex(1,1,1));
    F2.resize(GridIndex(1,1,1));
    return;
  }
  
//  int id = Globals::instance().getUniqueId();
  
  std::ostringstream fullfilename;
  fullfilename << filename << blockCount /* << "-" << id */ << ".hdf";
  
  std::ostringstream type;
  type << (isH?"E":"H") ;

  std::string blocks[3];
  blocks[0] = std::string("FieldX") + type.str(); 
  blocks[1] = std::string("FieldY") + type.str();  
  blocks[2] = std::string("FieldZ") + type.str();
  
  std::cerr << "Reading source from file " << filename << "\nBlocks: " 
    << blocks[0] << " " << blocks[1] << " " << blocks[2] << "\n";
  
  HDFistream data(fullfilename.str().c_str());
//  std::cerr << "File open\n";
  if (!data.good())
  {
    std::cerr << "PANIC!!\nCould not open source file " << fullfilename.str() << std::endl;
    exit(-1);
  }
  
//  diagsource << "Block " << (isH?"H":"E") << " " 
//             << blocks[transverse1] << " " << blocks[transverse2] << std::endl;
  
//  std::cerr << "Reading Block " << blocks[transverse1] << std::endl;
  data.setBlockName(blocks[transverse1]);
  
//==========================================================
// The following line is commented because it does not
// agree with the new HDFStream

//  data >> F1;

//  std::cerr << "Reading Block " << blocks[transverse2] << std::endl;
  data.setBlockName(blocks[transverse2]);

//==========================================================
// The following line is commented because it does not
// agree with the new HDFStream

//  data >> F2;
  
//  std::cerr << "Adjusting Parameters\n";
  data.close();
  Nt = F1.getHigh()[2] - F1.getLow()[2] + 1;
  blockToff = toff + blockCount*Nt;
  blockCount++;
  active = true;
//  std::cerr << "DONE Reading source from file " << filename << "\nBlocks: " 
//    << blocks[0] << " " << blocks[1] << " " << blocks[2] << "\n";
}

Vector SideInjectSourceFunc::getField(int i, int j, int k, int time, double factor)
{
  Vector Field(0,0,0);
  if (!active) return Field;
    
  int t = time - blockToff - F1.getLow()[2];
  
  GridIndex index(i,j,k);
    
  Field[transverse1] = factor*amp*F1(index[transverse1], index[transverse2], t);
  Field[transverse2] = factor*amp*F2(index[transverse1], index[transverse2], t);

//  diagsource << (isH?"H":"E") << " " << i << " " << j << " " << k << " "
//             << Field[0] << " " << Field[1] << " " << Field[2] << std::endl;

  return Field;
}

Vector SideInjectSourceFunc::getEField(int i, int j, int k, int time)
{
  return getField(i,j,k,time, 1.0);
}

Vector SideInjectSourceFunc::getHField(int i, int j, int k, int time)
{
  return getField(i,j,k,time, sqrt(eps));
}
