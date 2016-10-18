/*
 * border.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef BORDER_HPP_
#define BORDER_HPP_

#include "mpulse.hpp"

bool getBorderExtent(Direction dir,
                     int thickness,
                     int distance,
                     Index &blow,
                     Index &bhigh);


#endif /* BORDER_HPP_ */
