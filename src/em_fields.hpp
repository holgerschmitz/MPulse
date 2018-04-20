/*
 * em_fields.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_EM_FIELDS_H
#define MPULSE_EM_FIELDS_H

#include "types.hpp"

#include <schnek/variables/blockcontainer.hpp>

class MPulse;

/** @brief A container block for electromagnetic fields
 *
 */
class EMFields : public schnek::ChildBlock<EMFields>
{
  private:
    friend class MPulse;
    pField pEx;
    pField pEy;
    pField pEz;
    pField pBx;
    pField pBy;
    pField pBz;
  public:
    EMFields(schnek::pBlock parent = schnek::pBlock()) : schnek::ChildBlock<EMFields>(parent)
    {}
    void registerData();
};

#endif
