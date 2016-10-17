#ifndef MPULSE_CURRENTS_H
#define MPULSE_CURRENTS_H

#include "mpulse.hpp"

class Storage;

class Current : public schnek::ChildBlock<Current>
{
  protected:
    pField pJx;
    pField pJy;
    pField pJz;
  public:
    virtual ~Current() {}
    
    virtual void stepSchemeInit(double dt) = 0;
    virtual void stepScheme(double dt) = 0;

    bool isValid() const { return (pJx) && (pJy) && (pJz); }
    
    bool isMagneticCurrent() { return false; }

    const pField getJx() { return pJx; }
    const pField getJy() { return pJy; }
    const pField getJz() { return pJz; }
};

typedef boost::shared_ptr<Current> pCurrent;

#endif
