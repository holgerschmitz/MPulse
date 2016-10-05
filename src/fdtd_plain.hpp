/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */


#ifndef MPULSE_FDTD_PLAIN_H
#define MPULSE_FDTD_PLAIN_H

#include "fieldsolver.hpp"
#include "mpulse.hpp"

class Storage;

class FDTD_Plain : public FieldSolver
{
  private:
    DataGrid *pEx;
    DataGrid *pEy;
    DataGrid *pEz;
    DataGrid *pBx;
    DataGrid *pBy;
    DataGrid *pBz;

    Storage *storage;
  public:
    void init(Storage *storage_);
  
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
    
  private:
    void stepD(double dt);
    void stepB(double dt);
    
};

#endif
