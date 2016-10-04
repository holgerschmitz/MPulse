#ifndef MPULSE_CURRENTS_H
#define MPULSE_CURRENTS_H

#include "mpulse.h"
#include "rebuild.h"

class Storage;

class Current
{
  protected:
    DataGrid *pJx;
    DataGrid *pJy;
    DataGrid *pJz;
  public:
    virtual ~Current() {}
    virtual void initStorage(Storage *storage_) = 0;
    
    virtual void stepSchemeInit(double dt) = 0;
    virtual void stepScheme(double dt) = 0;

    bool isValid() const { return (pJx!=0) && (pJy!=0) && (pJz!=0); }
    
    const DataGrid *getJx() { return pJx; }
    const DataGrid *getJy() { return pJy; }
    const DataGrid *getJz() { return pJz; }
    
    struct CurrentInvalidPredicate
    {
      bool operator()(const Current *bc) { return !(bc->isValid()); }
    };
};

#endif
