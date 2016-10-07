/*
 * diagnostic.hpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_DIAGNOSTIC_HPP_
#define MPULSE_DIAGNOSTIC_HPP_

#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

class FieldDiagnostic : public schnek::HDFGridDiagnostic<Field, Field* >
{
  protected:
    typedef HDFGridDiagnostic<Field, Field* >::IndexType IndexType;
    IndexType getGlobalMin();
    IndexType getGlobalMax();
};



#endif /* SRC_DIAGNOSTIC_HPP_ */
