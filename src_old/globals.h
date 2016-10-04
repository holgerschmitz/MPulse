#ifndef MPULSE_GLOBALS_H
#define MPULSE_GLOBALS_H

#include "mpulse.h"
#include "rebuild.h"

#include <cmath>

/** @file globals.h
 * @brief global parameter class
 *
 * Defines a class for storing global Parameters
 */

/** A class for storing parameters, inherits Rebuildable.
 * Not all parameters might be needed for all types of simulation.
 */
class Globals : public Rebuildable {
  private:
    bool master;
    int uniqueId;
    /// global cell count x
    int GridX;
    /// global cell count y
    int GridY;
    /// global cell count z
    int GridZ;
    /// Size of a single grid cell in x-direction
    double GridDX;
    /// Size of a single grid cell in y-direction
    double GridDY;
    /// Size of a single grid cell in z-direction
    double GridDZ;
    /// time step
    double DT;
    /// total number of time steps
    int TotalTime;

    /// used for initilaizing the Fields
    GridIndex GridLow;
    /// used for initilaizing the Fields
    GridIndex GridHigh;

    ///will be set "true" if its a restart
    bool IsRestart;

    ///pointer to global parameters, will be set to "this" in constructor
    static Globals *globals;

    /// set to "true" if parameters class is initialized
    bool initialized;

    ///static member, number of the parameters main was called with
    static int Argc;
    /// static member, first argument main was called with
    static char **Argv;
  protected:
      ///Create the parameter map
    virtual ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  public:  
    
    /// constructor, sets up a not initialized parameters object, restart is "false", sets globals to "this"
    Globals();
    /// initialize and set initialized to "true"
    void init();
    
    /// accessor method, returns GridX
    int gridX() { return GridX; }
    /// accessor method, returns GridY
    int gridY() { return GridY; }
    /// accessor method, returns GridZ
    int gridZ() { return GridZ; }
    
    /// accessor method, returns grid spacing in x-direction
    double gridDX() { return GridDX; }
    /// accessor method, returns grid spacing in y-direction
    double gridDY() { return GridDY; }
    /// accessor method, returns grid spacing in z-direction
    double gridDZ() { return GridDZ; }

    /// accessor method, returns cell volume
    double volumeQuant() { return fabs(GridDX*GridDY*GridDZ); }
    
    /// accessor method, returns time step
    double dt() { return DT; }
    /// accessor method, returns total number of time steps
    int totalTime() { return TotalTime; }
    
    ///accessor method, returns isRestart
    bool isRestart() { return IsRestart; }
    ///sets isRestart to its argument
    void setRestart(bool IsRestart_) { IsRestart = IsRestart_; }

    ///accessor method, returns master
    bool isMaster() { return master; }
    ///sets master to its argument
    void setMaster(bool master_) { master = master_; }
    
    ///accessor method, returns master
    int getUniqueId() { return uniqueId; }
    ///sets master to its argument
    void setUniqueId(int uniqueId_) { uniqueId = uniqueId_; }
    
    ///accessor method, returns GridLow
    const GridIndex& gridLow() { return GridLow; }
    ///accessor method, returns GridHigh
    const GridIndex& gridHigh() { return GridHigh; }
    
    ///accessor method, eturns argc
    static int getArgc() { return Argc; }
    ///sets argc to its argument
    static void setArgc(int Argc_) { Argc = Argc_; }
    ///accessor method, returns argv
    static char **getArgv() { return Argv; }
    ///stes argv to its argument
    static void setArgv(char** Argv_) { Argv = Argv_; }
    
    ///if parameters is not initialized, initialize. in any case: return parameters
    static Globals &instance() { 
      if ( !(globals->initialized) ) globals->init();
      return *globals; 
    }
      
    ///rebuild
    std::string Rebuild(std::istream& in);
};


#endif // MPULSE_GLOBALS_H
