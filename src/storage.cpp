#include "storage.h"
#include "globals.h"
#include "boundary.h"



Storage::Storage()
  : low(0,0,0),
    high(1,1,1),
    boundary(0),
    grids(),
    gridsN(),
    gridsS(),
    gridsE(),
    gridsW(),
    groups()
{}

Storage::~Storage()
{
  delete boundary;
  GridDeleter deleteGrid;
  forAllGrids(grids, deleteGrid);
  forAllGrids(gridsN, deleteGrid);
  forAllGrids(gridsS, deleteGrid);
  forAllGrids(gridsE, deleteGrid);
  forAllGrids(gridsW, deleteGrid);
}

void Storage::resize(GridIndex low_, GridIndex high_)
{
  low = low_;
  high = high_;
  GridResizer resizer(low,high);
  forAllGrids(grids, resizer);
  if (
     (gridsN.size() + gridsS.size() 
    + gridsE.size() + gridsW.size() 
    + gridsU.size() + gridsD.size() 
    + lines.size() ) > 0
  )
  {
    std::cerr << "Resizing not implemented\n";
    exit(-1);
  }
}
    
DataGrid &Storage::getGrid(const std::string &gridid)
{
  return *(addGrid(gridid));
}

bool Storage::hasGrid(const std::string &gridid)
{
  return (grids.count(gridid) > 0);
}

DataGrid *Storage::addGrid(const std::string &gridid)
{
  return addGrid(gridid, low, high);
}

DataGrid *Storage::addGrid(const std::string &gridid, GridIndex lowg, GridIndex highg)
{
  DataGrid *g = 0;
  if (grids.count(gridid) == 0)
  {
    g = new DataGrid(lowg, highg);
    (*g) = 0;
    grids[gridid] = g;
  } else
  {
    g = grids[gridid];
  }
  return g;
}

DataGrid &Storage::getBorderLayer(const std::string &gridid, Direction dir)
{
  switch (dir)
  {
    case north: return *(gridsN[gridid]); break;
    case south: return *(gridsS[gridid]); break;
    case east:  return *(gridsE[gridid]); break;
    case west:  return *(gridsW[gridid]); break;
    case up:    return *(gridsU[gridid]); break;
    case down:  return *(gridsD[gridid]); break;
  }
  return *(gridsN[gridid]);
}

bool Storage::hasBorderLayer(const std::string &gridid, Direction dir)
{
  switch (dir)
  {
    case north: return (gridsN.count(gridid) > 0); break;
    case south: return (gridsS.count(gridid) > 0); break;
    case east:  return (gridsE.count(gridid) > 0); break;
    case west:  return (gridsW.count(gridid) > 0); break;
    case up:    return (gridsU.count(gridid) > 0); break;
    case down:  return (gridsD.count(gridid) > 0); break;
  }
  return false;
}

DataGrid *Storage::addBorderLayer(const std::string &gridid, 
                                  Direction dir, 
                                  int thickness, 
                                  int distance, 
                                  int ghostcells)
{
  DataGrid *g = 0;
  GridMap *gm;
  
//  std::cerr << "Adding border layer " << gridid << " " << dir << std::endl;
  switch (dir)
  {
    case north: gm = &gridsN; break;
    case south: gm = &gridsS; break;
    case east:  gm = &gridsE; break;
    case west:  gm = &gridsW; break;
    case up:    gm = &gridsU; break;
    case down:  
    default:    gm = &gridsD; break;
  }
  
  if (gm->count(gridid) == 0)
  {
    GridIndex b_low, b_high;
    if (getBorderExtent(dir, thickness, distance, ghostcells, b_low, b_high))
    {
      g = new DataGrid(b_low, b_high);
      (*g) = 0;
      gm->operator[](gridid) = g;
    }
  } else
  {
    g = grids[gridid];
  }
  return g;
}

DataLine &Storage::getLine(const std::string &lineid)
{
  return *(lines[lineid]);
}

bool Storage::hasLine(const std::string &lineid)
{
  return (lines.count(lineid) > 0);
}


DataLine *Storage::addLine(const std::string &lineid, int orientation)
{
  DataLine *ln = 0;
  if (lines.count(lineid) == 0)
  {
    ln = new DataLine(low[orientation], high[orientation]);
    (*ln) = 0;
    lines[lineid] = ln;
  } else
  {
    ln = lines[lineid];
  }
  return ln;
}

template<class Func>
void Storage::forAllGrids(const Func &func)
{
  forAllGrids(grids, func);
}

template<class Func>
void Storage::forAllGrids(GridMap &gm, Func &func)
{
  for (GridMap::iterator it=gm.begin(); it!=gm.end(); ++it)
  {
    func(it->first, it->second);
  }
}


void Storage::addToGroup(
    const std::string &groupid,
    const std::string &gridid
  )
{
  groups[groupid].push_back(gridid);
}
        
void Storage::applyBoundary(const std::string &groupid)
{
  IdList &gr = groups[groupid];
  for (IdList::iterator it=gr.begin(); it!=gr.end(); ++it)
  {
    DataGrid &g = *(grids[*it]);
    boundary->exchangeX(g);
    boundary->exchangeY(g);
    boundary->exchangeZ(g);
  }
}



bool Storage::getBorderExtent
  (
    Direction dir,
    int thickness,
    int distance,
    int ghostcells,
    GridIndex &blow,
    GridIndex &bhigh
  )
{
  bool haveBorder = false;
  GridIndex glow  = Globals::instance().gridLow();
  GridIndex ghigh = Globals::instance().gridHigh();

  blow[0] = low[0]-ghostcells;
  blow[1] = low[1]-ghostcells;
  blow[2] = low[2]-ghostcells;
  
  bhigh[0] = high[0]+ghostcells;
  bhigh[1] = high[1]+ghostcells;
  bhigh[2] = high[2]+ghostcells;
  
  switch (dir)
  {
    case west:
      if (low[0]<glow[0]+thickness+distance)
      {
        bhigh[0] = glow[0]+thickness-1+distance;
        blow[0] = glow[0]+distance;
        haveBorder = true;
      }
      break;
    case south:
      if (low[1]<glow[1]+thickness+distance)
      {
        bhigh[1] = glow[1]+thickness-1+distance;
        blow[1] = glow[1]+distance;
        haveBorder = true;
      }
      break;
    case down:
      if (low[2]<glow[2]+thickness+distance)
      {
        bhigh[2] = glow[2]+thickness-1+distance;
        blow[2] = glow[2]+distance;
        haveBorder = true;
      }
      break;
    case east:
      if (high[0]>ghigh[0]-thickness-distance)
      {
        blow[0] = ghigh[0]-thickness+1-distance;
        bhigh[0] = ghigh[0]-distance;
        haveBorder = true;
      }
      break;
    case north:
      if (high[1]>ghigh[1]-thickness-distance)
      {
        blow[1] = ghigh[1]-thickness+1-distance;
        bhigh[1] = ghigh[1]-distance;
        haveBorder = true;
      }
      break;
    case up:
      if (high[2]>ghigh[2]-thickness-distance)
      {
        blow[2] = ghigh[2]-thickness+1-distance;
        bhigh[2] = ghigh[2]-distance;
        haveBorder = true;
      }
      break;
  }
  /*
  std::cerr << "Function getBorderExtent: " << dir << " " << thickness << std::endl;
  std::cerr << "Grid Low " << glow[0] << " " << glow[1] << " " << glow[2] << std::endl;
  std::cerr << "Grid High " << ghigh[0] << " " << ghigh[1] << " " << ghigh[2] << std::endl;
  std::cerr << "Border Low " << blow[0] << " " << blow[1] << " " << blow[2] << std::endl;
  std::cerr << "Border High " << bhigh[0] << " " << bhigh[1] << " " << bhigh[2] << std::endl;
  */
  return haveBorder;
}

void Storage::GridDeleter::operator()(std::string, DataGrid* grid)
{
  delete grid;
}

Storage::GridResizer::GridResizer(const GridIndex &low_, const GridIndex &high_)
  : low(low_), high(high_) {}
  
void Storage::GridResizer::operator()(std::string, DataGrid* grid)
{
  grid->resize(low,high);
}

