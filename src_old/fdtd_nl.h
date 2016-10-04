#ifndef MPULSE_FDTD_NL_H
#define MPULSE_FDTD_NL_H

#include "fieldsolver.h"
#include "fielddiag.h"
#include "mpulse.h"

class Storage;

class FDTD_Nonlinear : public FieldSolver
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
    void initStorage(Storage *storage_);
  
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:
    /// build parametermap
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);

  private:
    void stepD(double dt);
    void stepB(double dt);
    
    /// value of eps_0*eps_infty (note that mu_0 = 0)
    double eps;
    
    /// value of eps_0*chi^(3)
    double chi;
    
};

#endif
