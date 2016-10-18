/*
 * current.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: holger
 */


#ifndef MPULSE_CURRENT_H
#define MPULSE_CURRENT_H

#include "mpulse.hpp"

class Current;
typedef boost::shared_ptr<Current> pCurrent;

class CurrentContainer
{
  protected:
    typedef std::list<pCurrent> CurrentList;

    CurrentList currents;
    CurrentList magCurrents;
  public:
    void addCurrent(pCurrent current);
    void addMagCurrent(pCurrent current);
};

class CurrentBlock : public schnek::ChildBlock<CurrentBlock>
{
  public:
    virtual ~CurrentBlock() {}
    virtual void initCurrents(CurrentContainer &container) = 0;
};

typedef boost::shared_ptr<CurrentBlock> pCurrentBlock;

class Current
{
  protected:
    pGrid pJx;
    pGrid pJy;
    pGrid pJz;
  public:
    virtual ~Current() {}
    
    virtual void init() = 0;
    virtual void stepSchemeInit(double dt) = 0;
    virtual void stepScheme(double dt) = 0;

    bool isValid() const { return (pJx) && (pJy) && (pJz); }
    
    const pGrid getJx() { return pJx; }
    const pGrid getJy() { return pJy; }
    const pGrid getJz() { return pJz; }
};


#endif
