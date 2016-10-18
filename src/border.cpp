/*
 * border.cpp
 *
 *  Created on: 18 Oct 2016
 *      Author: Holger Schmitz
 */

#include "border.hpp"

bool getBorderExtent(Direction dir,
                     int thickness,
                     int distance,
                     Index &blow,
                     Index &bhigh)
{
  bool haveBorder = false;

  schnek::DomainSubdivision<Field> &subdivision = MPulse::getSubdivision();
  Index glow  = Index(0);
  Index ghigh = MPulse::getGlobalMax();

  Index low  = subdivision.getInnerLo();
  Index high = subdivision.getInnerHi();

  blow[0] = low[0];
  blow[1] = low[1];
  blow[2] = low[2];

  bhigh[0] = high[0];
  bhigh[1] = high[1];
  bhigh[2] = high[2];

//  int borderfit = 0;

  switch (dir)
  {
    case west:
      if (low[0]<glow[0]+thickness+distance)
      {
        bhigh[0] = glow[0]+thickness-1+distance;
        blow[0] = glow[0]+distance;

        haveBorder = true;
      }
      break;
    case south:
      if (low[1]<glow[1]+thickness+distance)
      {
        bhigh[1] = glow[1]+thickness-1+distance;
        blow[1] = glow[1]+distance;

        haveBorder = true;
      }
      break;
    case down:
      if ((low[2]<glow[2]+thickness+distance) && (high[2]>=glow[2]+thickness))
      {
        bhigh[2] = glow[2]+thickness-1+distance;
        blow[2] = glow[2]+distance;

        haveBorder = true;
      }
      break;
    case east:
      if (high[0]>ghigh[0]-thickness-distance)
      {
        blow[0] = ghigh[0]-thickness+1-distance;
        bhigh[0] = ghigh[0]-distance;
        haveBorder = true;
      }
      break;
    case north:
      if (high[1]>ghigh[1]-thickness-distance)
      {
        blow[1] = ghigh[1]-thickness+1-distance;
        bhigh[1] = ghigh[1]-distance;

        haveBorder = true;
      }
      break;
    case up:
      if (high[2]>ghigh[2]-thickness-distance)
      {
        blow[2] = ghigh[2]-thickness+1-distance;
        bhigh[2] = ghigh[2]-distance;

        haveBorder = true;
      }
      break;
  }
  return haveBorder;
}

