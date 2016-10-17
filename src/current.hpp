#ifndef MPULSE_CURRENT_H
#define MPULSE_CURRENT_H

#include "mpulse.hpp"

class CurrentBlock : public schnek::ChildBlock<CurrentBlock>
{
  public:
    virtual ~CurrentBlock() {}
    virtual void initCurrents(FieldSolver *solver) = 0;
};

class Current
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
    
    const pField getJx() { return pJx; }
    const pField getJy() { return pJy; }
    const pField getJz() { return pJz; }
};

typedef boost::shared_ptr<Current> pCurrent;

#endif
