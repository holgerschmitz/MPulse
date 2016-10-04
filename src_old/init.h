#ifndef MPULSE_INIT_H
#define MPULSE_INIT_H

#include "rebuild.h"

class Storage;

class FieldSimInit : public Rebuildable
{
  public:
    virtual void init(Storage &fields) = 0;
};

#endif
