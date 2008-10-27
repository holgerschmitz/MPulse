#ifndef MPULSE_FDTD_DISP_H
#define MPULSE_FDTD_DISP_H

#include "fieldsolver.h"
#include "fielddiag.h"
#include "mpulse.h"
#include "currents.h"

class Storage;

class FDTD_Dispersion : public FieldSolver
{
  private:
    DataGrid *pEx;
    DataGrid *pEy;
    DataGrid *pEz;
    DataGrid *pBx;
    DataGrid *pBy;
    DataGrid *pBz;

    // Lorentz polarization for three poles
    DataGrid *pPx[3];
    DataGrid *pPy[3];
    DataGrid *pPz[3];
    
    // Lorentz polarization from the previous timestep
    DataGrid *pPxp[3];
    DataGrid *pPyp[3];
    DataGrid *pPzp[3];
    
    Storage *storage;
    
    typedef std::list<Current*> CurrentList;
    CurrentList currents;
    
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
    
    /// value of Delta epsilon_p^2 for the three Lorentz poles
    double LEps2[3];
    /// value of omega_p^2 for the three Lorentz poles
    double LOm2[3];
};

#endif
