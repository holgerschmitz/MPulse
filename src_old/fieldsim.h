#ifndef MPULSE_FIELDSIM_H
#define MPULSE_FIELDSIM_H

#include "rebuild.h"
#include "init.h"
#include "storage.h"
#include "fielddiag.h"
#include "optfield.h"

class FieldSolver;

class FieldSimulation : public Storage, public Rebuildable
{
  private:
    double dt;
    
  public:
    /** @brief Default constructor */
    FieldSimulation ();

    /** @brief Destructor, deletes boundary */
    virtual ~FieldSimulation ();

    /** @brief Initializes the FieldSimulation
     */
    void init ();

    /** @brief Perform one timestep.
     */
    void execute();

// ------------- Rebuildable -------------

  protected:
    ///build parametermap, containing boundary, species and general parameters
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);

    AllFieldDiag fieldDiag;
    
    FieldSolver *solver;
    
    FieldSimInit *initializer;
    
};

#endif
