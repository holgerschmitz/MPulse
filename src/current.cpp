/*
 * current.cpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "current.hpp"

#include <boost/make_shared.hpp>

void CurrentContainer::addCurrent(pCurrent current)
{
  current->init();
  if (current->isValid())
    this->currents.push_back(current);
}

void CurrentContainer::addMagCurrent(pCurrent current)
{
  current->init();
  if (current->isValid())
    this->magCurrents.push_back(current);
}

void CurrentContainer::sumCurrents()
{
  Field &jxT = *this->pJx;
  Field &jyT = *this->pJy;
  Field &jzT = *this->pJz;

  jxT = 0;
  jyT = 0;
  jzT = 0;

  BOOST_FOREACH(pCurrent current, this->currents)
  {
    Grid &jx = *current->getJx();
    Grid &jy = *current->getJy();
    Grid &jz = *current->getJz();

    Index low = jx.getLo();
    Index high = jx.getHi();
    for (int i=low[0]; i<=high[0]; ++i)
      for (int j=low[1]; j<=high[1]; ++j)
        for (int k=low[2]; k<=high[2]; ++k)
        {
          jxT(i,j,k) += jx(i,j,k);
          jyT(i,j,k) += jy(i,j,k);
          jzT(i,j,k) += jz(i,j,k);
        }
  }
}

void CurrentContainer::sumMagCurrents()
{
  Field &jxT = *this->pMx;
  Field &jyT = *this->pMy;
  Field &jzT = *this->pMz;

  jxT = 0;
  jyT = 0;
  jzT = 0;

  BOOST_FOREACH(pCurrent current, this->magCurrents)
  {
    Grid &jx = *current->getJx();
    Grid &jy = *current->getJy();
    Grid &jz = *current->getJz();

    Index low = jx.getLo();
    Index high = jx.getHi();
    for (int i=low[0]; i<=high[0]; ++i)
      for (int j=low[1]; j<=high[1]; ++j)
        for (int k=low[2]; k<=high[2]; ++k)
        {
          jxT(i,j,k) += jx(i,j,k);
          jyT(i,j,k) += jy(i,j,k);
          jzT(i,j,k) += jz(i,j,k);
        }
  }
}

void CurrentContainer::init()
{

  schnek::DomainSubdivision<Field> &subdivision = MPulse::getSubdivision();
  Index lowIn = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  schnek::Range<double, DIMENSION> domainSize(schnek::Array<double, DIMENSION>(0,0,0), MPulse::getSize());

  pJx = boost::make_shared<Field>(lowIn, highIn, domainSize, exStaggerYee, 2);
  pJy = boost::make_shared<Field>(lowIn, highIn, domainSize, eyStaggerYee, 2);
  pJz = boost::make_shared<Field>(lowIn, highIn, domainSize, ezStaggerYee, 2);

  pMx = boost::make_shared<Field>(lowIn, highIn, domainSize, bxStaggerYee, 2);
  pMy = boost::make_shared<Field>(lowIn, highIn, domainSize, byStaggerYee, 2);
  pMz = boost::make_shared<Field>(lowIn, highIn, domainSize, bzStaggerYee, 2);
}


