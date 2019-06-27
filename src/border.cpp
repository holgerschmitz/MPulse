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
                     Index &bhigh,
                     bool isH)
{
  bool haveBorder = false;
  int distanceLow = distance - (isH?1:0);
  int distanceHigh = distance - 0; //(isH?0:1);

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

  Index coords(0,1,2);

//  int borderfit = 0;

  switch (dir)
  {
    case west:
    case east:
      coords = Index(0,1,2);
      break;
    case south:
    case north:
      coords = Index(1,2,0);
      break;
    case down:
    case up:
      coords = Index(2,0,1);
      break;
  }

  int normal = coords[0];
  int t1 = coords[1];
  int t2 = coords[2];

  switch (dir)
  {
    case west:
    case south:
    case down:
      bhigh[normal] = std::min(glow[normal]+thickness-1+distanceLow, high[normal]);
      blow[normal] = std::max(glow[normal]+distanceLow, low[normal]);

      break;
    case east:
    case north:
    case up:
      blow[normal] = std::max(ghigh[normal]-thickness+2-distance, low[normal]);
      bhigh[normal] = std::min(ghigh[normal]-distance+1, high[normal]);

      break;
  }

  haveBorder = (blow[normal] <= high[normal]) && (bhigh[normal] >= low[normal]);

  if (haveBorder) {
    bhigh[t1] = std::min(ghigh[t1]-distanceHigh, high[t1]);
    blow[t1] = std::max(glow[t1]+distance, low[t1]);

    bhigh[t2] = std::min(ghigh[t2]-distanceHigh, high[t2]);
    blow[t2] = std::max(glow[t2]+distance, low[t2]);

    haveBorder = haveBorder
        && (blow[t1] <= high[t1]) && (bhigh[t1] >= low[t1])
        && (blow[t2] <= high[t2]) && (bhigh[t2] >= low[t2]);
  }

  std::cout << "Border " << haveBorder << " " << blow << " " << bhigh << std::endl;

  return haveBorder;
}

