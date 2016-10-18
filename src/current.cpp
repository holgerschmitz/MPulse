/*
 * current.cpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "current.hpp"

void CurrentContainer::addCurrent(pCurrent current)
{
  this->currents.push_back(current);
}

void CurrentContainer::addMagCurrent(pCurrent current)
{
  this->magCurrents.push_back(current);
}

