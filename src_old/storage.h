#ifndef MPULSE_STORAGE_H
#define MPULSE_STORAGE_H

#include "mpulse.h"
#include "boundary.h"
#include <map>
#include <list>
#include <string>

class Boundary;


/** Storage Base class for storing the fields
 */
class Storage
{
  public:
    /** Construct the fields with zero size*/
    Storage();
    
    /** Destructor (empty) */
    virtual ~Storage();

    /** Resizes the fields to span from low to high */
    void resize(GridIndex low_, GridIndex high_);

  protected:
    /** The low index of the fields */
    GridIndex low;
    /** The high index of the fields */
    GridIndex high;
    
    /** The grid delta x */
    double dx;
    /** The grid delta y */
    double dy;
    /** The grid delta z */
    double dz;
    
  public:
    
    DataGrid &getGrid(const std::string &gridid);
    bool hasGrid(const std::string &gridid);
    
    template<class Func>
    void forAllGrids(const Func &func);
    
    DataGrid *addGrid(const std::string &gridid);
    DataGrid *addGrid(const std::string &gridid, GridIndex lowg, GridIndex highg);
    void addToGroup(const std::string &groupid, const std::string &gridid);

    DataGrid &getBorderLayer(const std::string &gridid, Direction dir);
    bool hasBorderLayer(const std::string &gridid, Direction dir);
    DataGrid *addBorderLayer(const std::string &gridid, 
                             Direction dir, 
                             int thickness, 
                             int distance=0, 
                             int bordercells=0);
    
    DataLine &getLine(const std::string &lineid);
    bool hasLine(const std::string &lineid);
    DataLine *addLine(const std::string &lineid, int orientation);

    /** Apply the boundary condition to the electric field*/
    void applyBoundary(const std::string &groupid);
    
    const GridIndex &getLow() { return low; }
    const GridIndex &getHigh() { return high; }
    
    double getDx() { return dx; }
    double getDy() { return dy; }
    double getDz() { return dz; }

  protected:
    Boundary *boundary;
    
  private:
    typedef std::map<std::string, DataGrid*> GridMap;
    typedef std::map<std::string, DataLine*> LineMap;
    typedef std::list<std::string> IdList;
    GridMap grids;
    GridMap gridsN;
    GridMap gridsS;
    GridMap gridsE;
    GridMap gridsW;
    GridMap gridsU;
    GridMap gridsD;
    
    LineMap lines;
    
    std::map<std::string,IdList> groups;

    template<class Func>
    void forAllGrids(GridMap &gm, Func &func);
    
    bool getBorderExtent(Direction dir, 
                         int thickness, 
                         int distance, 
                         int bordercells, 
                         GridIndex &blow, 
                         GridIndex &bhigh);
    
    struct GridDeleter 
    {
      void operator()(std::string, DataGrid* grid); 
    };
    
    struct GridResizer 
    {
      const GridIndex &low, &high;
      GridResizer(const GridIndex &low_, const GridIndex &high_);
      void operator()(std::string, DataGrid* grid); 
    };
  public: 
    const Boundary& getBoundary() const { return *boundary; }
};




#endif // MPULSE_STORAGE_H
