#include "fielddiag.h"
#include "process.h"
#include "globals.h"
#include "boundary.h"
//-----------------------------------------------------------------------------
//---------   FieldDiag
//-----------------------------------------------------------------------------

ParameterMap* FieldDiag::MakeParamMap (ParameterMap* pm)
{
  pm = ParentType::MakeParamMap(pm);

  (*pm)["field"] = WParameter(new ParameterValue<std::string>(&fieldId,""));
  return pm;
}

void FieldDiag::fetchField(Storage &storage)
{
  field.grid = &storage.getGrid(fieldId);
  field.global_min = Globals::instance().gridLow();
  field.global_max = Globals::instance().gridHigh();
  
  this->setField(&field);
}

//-----------------------------------------------------------------------------
//---------   FieldSliceDiag
//-----------------------------------------------------------------------------

void FieldSliceDiag::open(const std::string &fname)
{
  output.open(fname.c_str());
}

FieldSliceDiag::~FieldSliceDiag()
{
  output.close();
}

void FieldSliceDiag::close()
{
  output.close();
}


ParameterMap* FieldSliceDiag::MakeParamMap (ParameterMap* pm)
{
  pm = DiagnosticInterface::MakeParamMap(pm);

  (*pm)["field"] = WParameter(new ParameterValue<std::string>(&fieldId,""));
  (*pm)["plane"] = WParameter(new ParameterValue<std::string>(&plane,"x"));
  (*pm)["pos"] = WParameter(new ParameterValue<int>(&pos,5));
  return pm;
}

void FieldSliceDiag::fetchField(Storage &storage)
{
  this->field = &storage.getGrid(fieldId);
  GridIndex low = storage.getLow();
  GridIndex high = storage.getHigh();
  GridIndex glow = Globals::instance().gridLow();
  GridIndex ghigh = Globals::instance().gridHigh();
  
  switch (plane[0]) {
    case 'y':
      dim1 = 0;
      dim2 = 2;
      normal = 1;
      break;
    case 'z':
      dim1 = 0;
      dim2 = 1;
      normal = 2;
      break;
    case 'x':
    default :
      dim1 = 1;
      dim2 = 2;
      normal = 0;
      break;
  }
  
  low1 = low[dim1];      
  low2 = low[dim2];      
  high1 = high[dim1];    
  high2 = high[dim2];
  
  active = (pos>low[normal]) && (pos<high[normal]);
  output.setActive(active);
    
  slice.resize(GridIndex2d(low1, low2), GridIndex2d(high1, high2));
  sliceContainer.grid = &slice;
  sliceContainer.global_min = GridIndex2d(glow[dim1], glow[dim2]);
  sliceContainer.global_max = GridIndex2d(ghigh[dim1], ghigh[dim2]);
}

void FieldSliceDiag::write()
{
  if (active)
  {
    GridIndex i;
    i[normal] = pos;

    int d1 = dim1;
    int d2 = dim2;

    for (i[d1]=low1; i[d1]<=high1; ++i[d1])
      for (i[d2]=low2; i[d2]<=high2; ++i[d2])
      {
        slice(i[d1],i[d2]) = (*field)(i[0], i[1], i[2]);
      }
  }
    
  output << sliceContainer;
}


//-----------------------------------------------------------------------------
//---------   FieldLineDiag
//-----------------------------------------------------------------------------

void FieldLineDiag::open(const std::string &fname)
{
  if (active) output.open(fname.c_str());
}

FieldLineDiag::~FieldLineDiag()
{
  output.close();
}

void FieldLineDiag::close()
{
  output.close();
}


ParameterMap* FieldLineDiag::MakeParamMap (ParameterMap* pm)
{
  pm = DiagnosticInterface::MakeParamMap(pm);

  (*pm)["field"] = WParameter(new ParameterValue<std::string>(&fieldId,""));
  (*pm)["dim"] = WParameter(new ParameterValue<std::string>(&direction,"x"));
  (*pm)["posx"] = WParameter(new ParameterValue<int>(&posx,5));
  (*pm)["posy"] = WParameter(new ParameterValue<int>(&posy,5));
  (*pm)["posz"] = WParameter(new ParameterValue<int>(&posz,5));
  return pm;
}

void FieldLineDiag::fetchField(Storage &storage)
{
  this->field = &storage.getGrid(fieldId);
  GridIndex low = storage.getLow();
  GridIndex high = storage.getHigh();
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
  
  active = (pos[trans1]>low[trans1]) && (pos[trans1]<high[trans1])
    && (pos[trans2]>low[trans2]) && (pos[trans2]<high[trans2]);
}

void FieldLineDiag::write()
{
  if (active)
  {
    std::string sep("");
    GridIndex i(posx, posy, posz);
    
    for (i[dim]=low; i[dim]<=high; ++i[dim])
    {
      output << sep << (*field)(i[0], i[1], i[2]);
      sep = " ";
    }
    output << std::endl;
  }
}

//-----------------------------------------------------------------------------
//---------   FieldEnergyDiag
//-----------------------------------------------------------------------------

FieldEnergyDiag::FieldEnergyDiag() : storage(0) {}
	  
FieldEnergyDiag::~FieldEnergyDiag() {}
	  
void FieldEnergyDiag::setStorage(Storage *storage_)
{
  storage = storage_;
}

void FieldEnergyDiag::open(const std::string &fname){
  output.open(fname.c_str());
}

	  
void FieldEnergyDiag::write()
{
  int time = Process::instance().getTime();
  
  output << time << " " << energy << std::endl;
}

void FieldEnergyDiag::close()
{
  output.close();
}

void FieldEnergyDiag::calculate()
{
  DataGrid &Ex = storage->getGrid("Ex");
  DataGrid &Ey = storage->getGrid("Ey");
  DataGrid &Ez = storage->getGrid("Ez");
  DataGrid &Bx = storage->getGrid("Bx");
  DataGrid &By = storage->getGrid("By");
  DataGrid &Bz = storage->getGrid("Bz");

  GridIndex low = storage->getLow();
  GridIndex high = storage->getHigh();
  
  energy = 0;
  
  for (int i=low[0]+1; i<high[0]; ++i)
    for (int j=low[1]+1; j<high[1]; ++j)
      for (int k=low[2]+1; k<high[2]; ++k)
      {
        double ex = Ex(i,j,k);
        double ey = Ey(i,j,k);
        double ez = Ez(i,j,k);
        double bx = Bx(i,j,k);
        double by = By(i,j,k);
        double bz = Bz(i,j,k);
        energy += ex*ex + ey*ey + ez*ez + bx*bx + by*by + bz*bz;
      }
  
  const Boundary &bound = storage->getBoundary();
  energy = bound.SumReduce(energy);
  Globals &glob = Globals::instance();
  energy *= glob.gridDX() * glob.gridDY() * glob.gridDZ();
}

