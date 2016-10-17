#ifndef MPULSE_FDTD_PLRC_H
#define MPULSE_FDTD_PLRC_H

#include "fieldsolver.hpp"
#include "mpulse.hpp"

#include <complex>

class Storage;

class FDTD_PLRCCore : public FieldSolver, public schnek::BlockContainer<Current>
{
  public:
    void registerData();
    void init();
  protected:
    pField pEx;
    pField pEy;
    pField pEz;
    pField pBx;
    pField pBy;
    pField pBz;
    
    // electrical conductivity of the medium
    pField pSigma;

    pDataLine pKappaEdx;
    pDataLine pKappaEdy;
    pDataLine pKappaEdz;

    pDataLine pKappaHdx;
    pDataLine pKappaHdy;
    pDataLine pKappaHdz;

    pDataLine pCpmlSigmaEx;
    pDataLine pCpmlSigmaEy;
    pDataLine pCpmlSigmaEz;

    pDataLine pCpmlSigmaHx;
    pDataLine pCpmlSigmaHy;
    pDataLine pCpmlSigmaHz;
    
    // Real part of the accumulator for three poles
    pField pPsiRx[3];
    pField pPsiRy[3];
    pField pPsiRz[3];
    
    // Imaginary part of the accumulator for three poles
    pField pPsiIx[3];
    pField pPsiIy[3];
    pField pPsiIz[3];
    
    typedef std::list<Current*> CurrentList;
    CurrentList currents;
    CurrentList magCurrents;

    OptFieldList optfields;
    OptFieldList optfieldsE;
    OptFieldList optfieldsH;
    
    OptFieldList optfieldsRes;
    
    typedef std::list<CurrentFactory*> CurrentFactoryList;
    CurrentFactoryList currentFactories;
    
    struct PLRCData
    {
      std::complex<double> dchi0[3];
      std::complex<double> dxi0[3];
      std::complex<double> Crec[3];
      double sumChi0;
      double sumXi0;
    };

    PLRCData plrcData;
    
    /// value of eps_infty (note that eps_0 = mu_0 = 1)
    double eps;
    
    /// value of Delta epsilon_p for the three Lorentz poles
    double LEps[3];
    /// value of delta_p^2 for the three Lorentz poles
    double LDelta[3];
    /// value of omega_p^2 for the three Lorentz poles
    double LOm2[3];
  private:
  
    template<class StorageHolder>
    struct InitStorageFunctor
    {
      Storage *storage;
      InitStorageFunctor(Storage *storage_) : storage(storage_) {}
      void operator()(StorageHolder *holder) { holder->initStorage(storage); }
    };

    
    
};

class FDTD_PLRCLinCore : public FDTD_PLRCCore
{
  protected: 
    void plrcStepD(double dt, 
                   int i, int j, int k, 
                   double dx, double dy, double dz,
                   double Jx, double Jy, double Jz);
                   
    void plrcStepB(double dt, 
                   int i, int j, int k, 
                   double dx, double dy, double dz,
                   double Jx, double Jy, double Jz);
                   
    ParameterMap* CustomParamMap (ParameterMap* pm = NULL);
};

class FDTD_PLRCNonlinCore : public FDTD_PLRCCore
{
  protected: 
    void plrcStepD(double dt, 
                   int i, int j, int k, 
                   double dx, double dy, double dz,
                   double Jx, double Jy, double Jz);
                   
    void plrcStepB(double dt, 
                   int i, int j, int k, 
                   double dx, double dy, double dz,
                   double Jx, double Jy, double Jz);
                   
    ParameterMap* CustomParamMap (ParameterMap* pm = NULL);
      
    /// value of chi^(3)
    double chi;
    
};

template<class PLRCImplementation>
class FDTD_PLRCSolver : public PLRCImplementation
{
  public:
    void initStorage(Storage *storage_);
    void addSigma(Current *current);
    void addCurrent(Current *current);
    void addMagCurrent(Current *current);
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:
    typedef PLRCImplementation Implementation;
    /// build parametermap
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);

  private:
    pField pJx;
    pField pJy;
    pField pJz;
    
    pField pMx;
    pField pMy;
    pField pMz;

    void stepD(double dt);
    void stepB(double dt);
    
    void initAccumulator(double dt);
    void sumCurrents();
    void sumMagCurrents();
    
    typedef typename Implementation::CurrentList CurrentList;
    typedef typename Implementation::CurrentFactoryList CurrentFactoryList;
    
};

typedef FDTD_PLRCSolver<FDTD_PLRCLinCore> FDTD_PLRCLin;
typedef FDTD_PLRCSolver<FDTD_PLRCNonlinCore> FDTD_PLRCNonlin;

#include "fdtd_plrc.t"

#endif
