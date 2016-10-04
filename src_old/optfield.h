#ifndef MPULSE_OPTFIELDS_H
#define MPULSE_OPTFIELDS_H

#include "rebuild.h"

class Storage;

class OptField : public Rebuildable
{
  public:    
    virtual void initStorage(Storage *storage_) = 0;
    
    virtual void stepSchemeInit(double dt) = 0;
    virtual void stepScheme(double dt) = 0;
    
    /** should the field be excecuted before the magnetic field
     * time step (true) or just before the electric field timestep (false)
     * Note that the corresponding currents are executed AFTER the optional fields.
     */
    virtual bool isH() { return false; }
};

typedef std::list<OptField*> OptFieldList;

#endif
