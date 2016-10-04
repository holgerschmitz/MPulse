#ifndef MPULSE_CURRENTFACTORY_H
#define MPULSE_CURRENTFACTORY_H

#include "mpulse.h"
#include "rebuild.h"

class Storage;
class FieldSolver;

class CurrentFactory : public Rebuildable
{
  public:    
    virtual ~CurrentFactory() {}
    virtual void initCurrents(Storage *storage_, FieldSolver *solver) = 0;
};

#endif
