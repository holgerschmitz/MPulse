//===============================================================
//==========  IncidentSourceECurrent
//===============================================================

#include "storage.h"
#include "process.h"

template<class SourceFunc>
IncidentSourceECurrent<SourceFunc>::IncidentSourceECurrent(int distance_, Direction dir_)
  : IncidentSourceCurrent(distance_, dir_, false),
    SourceFunc(dir_, false)
{}

template<class SourceFunc>
void IncidentSourceECurrent<SourceFunc>::initStorage(Storage *storage)
{
  pJx = storage->addBorderLayer("IncidentEJx", IncidentSourceCurrent::dir, 1, distance);
  pJy = storage->addBorderLayer("IncidentEJy", IncidentSourceCurrent::dir, 1, distance);
  pJz = storage->addBorderLayer("IncidentEJz", IncidentSourceCurrent::dir, 1, distance);
  
  DataGrid *allJ[3];
  allJ[0] = pJx;
  allJ[1] = pJy;
  allJ[2] = pJz;
  
  pJ[0] = allJ[IncidentSourceCurrent::transverse1];
  pJ[1] = allJ[IncidentSourceCurrent::transverse2];
  
  if ((pJx!=0) && (pJy!=0) && (pJz!=0))
     this->initSourceFunc(storage, pJx, pJy, pJz);
}

template<class SourceFunc>
void IncidentSourceECurrent<SourceFunc>::stepSchemeInit(double dt)
{
  stepScheme(0.5*dt);
}

template<class SourceFunc>
void IncidentSourceECurrent<SourceFunc>::stepScheme(double dt)
{
  DataGrid &J0 = *pJ[0];
  DataGrid &J1 = *pJ[1];

  GridIndex low  = J0.getLow();
  GridIndex high = J0.getHigh();
  
  GridIndex ind, indn;
  
  int Time = Process::instance().getTime();
  this->setTime(Time);
  
  int off[3] = {0,0,0};
  off[IncidentSourceCurrent::dim] = reverse?-1:-1;
  
  double DX = dX[2];
  
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
  {
    int x = ind[0]+off[0];
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
    {
      int y = ind[1]+off[1];
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2])
      {
        int z = ind[2]+off[2];
        Vector H = this->getHField(x,y,z,Time);
        J0(ind[0], ind[1], ind[2]) = H[IncidentSourceCurrent::transverse2]/DX;
        J1(ind[0], ind[1], ind[2]) = H[IncidentSourceCurrent::transverse1]/DX;
      }
    }
  }
}


//===============================================================
//==========  IncidentSourceHCurrent
//===============================================================

template<class SourceFunc>
IncidentSourceHCurrent<SourceFunc>::IncidentSourceHCurrent(int distance_, Direction dir_)
  : IncidentSourceCurrent(distance_, dir_, true),
    SourceFunc(dir_, true)
{}

template<class SourceFunc>
void IncidentSourceHCurrent<SourceFunc>::initStorage(Storage *storage)
{
  if (!reverse) distance = distance-1;
  pJx = storage->addBorderLayer("IncidentHJx", IncidentSourceCurrent::dir, 1, distance);
  pJy = storage->addBorderLayer("IncidentHJy", IncidentSourceCurrent::dir, 1, distance);
  pJz = storage->addBorderLayer("IncidentHJz", IncidentSourceCurrent::dir, 1, distance);

  DataGrid *allJ[3];
  allJ[0] = pJx;
  allJ[1] = pJy;
  allJ[2] = pJz;
  
  pJ[0] = allJ[IncidentSourceCurrent::transverse1];
  pJ[1] = allJ[IncidentSourceCurrent::transverse2];

  this->initSourceFunc(storage, pJx, pJy, pJz);
}

template<class SourceFunc>
void IncidentSourceHCurrent<SourceFunc>::stepSchemeInit(double dt)
{
  stepScheme(0.5*dt);
}

template<class SourceFunc>
void IncidentSourceHCurrent<SourceFunc>::stepScheme(double dt)
{
  DataGrid &J0 = *pJ[0];
  DataGrid &J1 = *pJ[1];
  
  GridIndex low  = J0.getLow();
  GridIndex high = J0.getHigh();
  
  GridIndex ind, indn;
  
  int Time = Process::instance().getTime();
  this->setTime(Time);
  
  int off[3] = {0,0,0};
  off[IncidentSourceCurrent::dim] = reverse?+1:0;
  
  double DX = dX[2];
  
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
  {
    int x = ind[0]+off[0];
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
    {
      int y = ind[1]+off[1];
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2])
      {
        int z = ind[2]+off[2];
        Vector E = this->getEField(x,y,z,Time);
        J0(ind[0], ind[1], ind[2]) = E[IncidentSourceCurrent::transverse2]/DX;
        J1(ind[0], ind[1], ind[2]) = E[IncidentSourceCurrent::transverse1]/DX;
      }
    }
  }
}


