#ifndef MPULSE_PROCESS_H
#define MPULSE_PROCESS_H
//-----------------------------------------------------------------------------
#include "rebuild.h"
#include "parameter.h"
#include "globals.h"
#include "fieldsim.h"

//-----------------------------------------------------------------------------
/** @file process.h
  * @brief contains the process classes
  *
  * Here the process classes are defined, which are the global objects of this code.
  */

/** @brief The process class
  *
  * This class has only one instance and controls setup, rebuilding, initializing and
  * executing of all objects of the programm. First it calls all the Rebuild methods, 
  * builds the correct objects and then distributes the execution to the execxte() members
  * of the appropiate classes. It also prints global (error) message on the screen.
  */
class Process : public Rebuildable {
  private:
    ///parameters
    Globals *globals;
    
    /// pointer to process
    static Process *process;
    
    /// pointer to the field simulation class
    FieldSimulation *fieldSim;
    
    /// time step
    int time;
  public:
	//-----------------------------------------------------------------------
	//constructor & destructor
    /** @brief constructor, 
      *
      * sets static member process to "this" and time step to 0. Since process is
      * static, the lifetime of the dynamical allocated Process object is prolonged 
      * indefinetly. Which means, it has to be deleted in finalize()
      */
    Process()
    {
      process = this;
      time = 0;
    }
    
    ///destructor, deletes boundary
    ~Process();
	
    /// returns time
    int getTime() { return time; }
    
    /// returns reference to current instance of the process class
    static Process& instance() { return *process; }
    
    /** @brief Initialize the fields and species
     *
     * Call the initialize members of each items in the SpeciesList species and 
     * the field from parameters.
     */
    void init();
    
    /** @brief Run the process.
     *
     * Which means: For the desired number of timesteps call the execute() members of
     * the objects stored in the speciesList
     */
    void run();
        
// ------------- Rebuildable -------------

  protected:
    ///build parametermap, containing boundary, species and general parameters
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  public:
    ///rebuild from file, call Rebuildable::Rebuild and makeProcess
    std::string Rebuild (std::istream& in);
};

#endif //MPULSE_PROCESS_H
