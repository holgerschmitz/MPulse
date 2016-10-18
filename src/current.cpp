/*
 * current.cpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "current.hpp"

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

