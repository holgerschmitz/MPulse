#include "plasmadensity.h"
#include "plasmacurrent.h"
#include "cpml_border.h"
#include "incsource.h"
#include "sources.h"
#include "gaussinject.h"
#include "shortpulseinject.h"
#include <boost/function.hpp>
#include <boost/bind.hpp>

//===============================================================
//==========  FDTD_PLRCSolver
//===============================================================


template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::initStorage(Storage *storage_)
{
  if ((this->LOm2[0] == 0) || (this->LOm2[1] == 0) || (this->LOm2[2] == 0))
  {
    std::cerr << "FDTD_PLRCSolver: Omega should not be zero!\n";
    exit(-1);
  }
  
  for (typename OptFieldList::iterator it=this->optfields.begin();
        it != this->optfields.end();
        ++it
      )
  {
    OptField *field = (*it);
    field->initStorage(storage_);
    if (field->isH())
      this->optfieldsH.push_back(field);
    else
      this->optfieldsE.push_back(field);
  }
  
  for (typename CurrentFactoryList::iterator it=this->currentFactories.begin();
        it != this->currentFactories.end();
        ++it
      )
  {
    (*it)->initCurrents(storage_, this);
  }
  
  this->coreInitStorage(storage_); 
  
  if (!(this->currents.empty()))
  {
    pJx = this->storage->addGrid("JxTotal");
    pJy = this->storage->addGrid("JyTotal");
    pJz = this->storage->addGrid("JzTotal");
    std::cerr << "Added Currents\n";
  }
  else
  {
    pJx = 0;
    pJy = 0;
    pJz = 0;
  }
  
  if (!(this->magCurrents.empty()))
  {
    pMx = this->storage->addGrid("MxTotal");
    pMy = this->storage->addGrid("MyTotal");
    pMz = this->storage->addGrid("MzTotal");
    std::cerr << "Added Magnetic Currents\n";
  }
  else
  {
    pMx = 0;
    pMy = 0;
    pMz = 0;
  }
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::addCurrent(Current *current)
{
  this->currents.push_back(current);
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::addMagCurrent(Current *current)
{
  this->magCurrents.push_back(current);
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepSchemeInit(double dt)
{
  initAccumulator(dt);

  stepB(0.5*dt);
  
  for ( typename OptFieldList::iterator it=this->optfieldsE.begin();
        it != this->optfieldsE.end();
        ++it )
  {
    (*it)->stepSchemeInit(dt);
  }

  for ( typename CurrentList::iterator it = this->currents.begin(); 
        it != this->currents.end(); 
        ++it )
  {
    (*it)->stepSchemeInit(dt);
  }

  for ( typename OptFieldList::iterator it=this->optfieldsH.begin();
        it != this->optfieldsH.end();
        ++it )
  {
    (*it)->stepSchemeInit(dt);
  }

  for ( typename CurrentList::iterator it = this->magCurrents.begin(); 
        it != this->magCurrents.end(); 
        ++it )
  {
    (*it)->stepSchemeInit(dt);
  }
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepScheme(double dt)
{
  (*this->pSigma) = 0;
  for ( typename OptFieldList::iterator it=this->optfieldsE.begin();
        it != this->optfieldsE.end();
        ++it )
  {
    (*it)->stepScheme(dt);
  }

  for 
  (
    typename CurrentList::iterator it = this->currents.begin(); 
    it != this->currents.end(); 
    ++it
  )
  {
    (*it)->stepScheme(dt);
  }
  stepD(dt);
  
  for 
  ( 
    typename OptFieldList::iterator it=this->optfieldsH.begin();
    it != this->optfieldsH.end();
    ++it
  )
  {
    (*it)->stepScheme(dt);
  }

  for 
  (
    typename CurrentList::iterator it = this->magCurrents.begin(); 
    it != this->magCurrents.end(); 
    ++it
  )
  {
    (*it)->stepScheme(dt);
  }
  stepB(dt);
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepD(double dt)
{
  GridIndex low = this->storage->getLow();
  GridIndex high = this->storage->getHigh();

  double dx = this->storage->getDx();
  double dy = this->storage->getDy();
  double dz = this->storage->getDz();

  
  /// value of beta_p for the three Lorentz poles
  double beta[3];
  /// value of gamma_p for the three Lorentz poles
  double gamma[3];
  
  std::complex<double> I(0,1);
  std::complex<double> phasor[3];
  
  std::complex<double> chi0[3];
  std::complex<double> xi0[3];
  
  std::complex<double> sumChi = 0.0;
  std::complex<double> sumXi = 0.0;

  this->plrcData.sumChi0 = 0;
  this->plrcData.sumXi0  = 0;

  for (int n=0;n<3;++n)
  {
    beta[n] = sqrt(this->LOm2[n] - this->LDelta[n]*this->LDelta[n]);
    gamma[n] = this->LEps[n]*this->LOm2[n]/beta[n];
    phasor[n] = this->LDelta[n] - I*beta[n];
    this->plrcData.Crec[n] = exp(-phasor[n]*dt);
    chi0[n] = (I*gamma[n]*(this->plrcData.Crec[n]-1.0)/phasor[n] );
    xi0[n] = I*gamma[n] *(this->plrcData.Crec[n]*(phasor[n]*dt + 1.0) - 1.0)
              / (phasor[n]*phasor[n]*dt);
    
    this->plrcData.dchi0[n] = chi0[n]*(1.0-this->plrcData.Crec[n]);
    this->plrcData.dxi0[n]  = xi0[n] *(1.0-this->plrcData.Crec[n]);
    
    sumChi += chi0[n];
    sumXi  += xi0[n];
  }

  this->plrcData.sumChi0 = std::real(sumChi);
  this->plrcData.sumXi0  = std::real(sumXi);

  double jx(0), jy(0), jz(0);
  if (this->pJx != 0) sumCurrents();

  for (int i=low[0]+1; i<high[0]; ++i)
    for (int j=low[1]+1; j<high[1]; ++j)
      for (int k=low[2]+1; k<high[2]; ++k)
      {
        if (this->pJx != 0)
        {
          jx = (*this->pJx)(i,j,k);
          jy = (*this->pJy)(i,j,k);
          jz = (*this->pJz)(i,j,k);
        }
        
        this->plrcStepD(dt, i, j, k, dx, dy, dz, jx, jy, jz);
      }
        
        
  this->storage->applyBoundary("E");
}


template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::initAccumulator(double dt)
{
  DataGrid &Ex = *this->pEx;
  DataGrid &Ey = *this->pEy;
  DataGrid &Ez = *this->pEz;

  GridIndex low = this->storage->getLow();
  GridIndex high = this->storage->getHigh();
  
  /// value of beta_p for the three Lorentz poles
  double beta[3];
  /// value of gamma_p for the three Lorentz poles
  double gamma[3];
  
  std::complex<double> I(0,1);
  std::complex<double> phasor[3];
  std::complex<double> xi0[3];

  for (int n=0;n<3;++n)
  {
    beta[n] = sqrt(this->LOm2[n] - this->LDelta[n]*this->LDelta[n]);
    gamma[n] = this->LEps[n]*this->LOm2[n]/beta[n];
    phasor[n] = this->LDelta[n] - I*beta[n];
    this->plrcData.Crec[n] = exp(-phasor[n]*dt);

    xi0[n] = I*gamma[n] *(this->plrcData.Crec[n]*(phasor[n]*dt + 1.0) - 1.0)
              / (phasor[n]*phasor[n]*dt);
    
    this->plrcData.dxi0[n]  = xi0[n] *(1.0-this->plrcData.Crec[n]);
    
  }

  for (int i=low[0]+1; i<high[0]; ++i)
    for (int j=low[1]+1; j<high[1]; ++j)
      for (int k=low[2]+1; k<high[2]; ++k)
  {
    double ex = Ex(i,j,k);
    double ey = Ey(i,j,k);
    double ez = Ez(i,j,k);
            
    for (int n=0;n<3;++n)
    {
      REAL &pxr = this->pPsiRx[n]->operator()(i,j,k);
      REAL &pyr = this->pPsiRy[n]->operator()(i,j,k);
      REAL &pzr = this->pPsiRz[n]->operator()(i,j,k);
      
      REAL &pxi = this->pPsiIx[n]->operator()(i,j,k);
      REAL &pyi = this->pPsiIy[n]->operator()(i,j,k);
      REAL &pzi = this->pPsiIz[n]->operator()(i,j,k);
      
      std::complex<double> px = std::complex<double>(pxr,pxi);
      std::complex<double> py = std::complex<double>(pyr,pyi);
      std::complex<double> pz = std::complex<double>(pzr,pzi);

      px = this->plrcData.dxi0[n]*ex + this->plrcData.Crec[n]*px;
      py = this->plrcData.dxi0[n]*ey + this->plrcData.Crec[n]*py;
      pz = this->plrcData.dxi0[n]*ez + this->plrcData.Crec[n]*pz;
      
      pxr = std::real(px);
      pyr = std::real(py);
      pzr = std::real(pz);

      pxi = std::imag(px);
      pyi = std::imag(py);
      pzi = std::imag(pz);

    }
  }
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepB(double dt)
{
  GridIndex low = this->storage->getLow();
  GridIndex high = this->storage->getHigh();

  double dx = this->storage->getDx();
  double dy = this->storage->getDy();
  double dz = this->storage->getDz();

  double jx(0), jy(0), jz(0);
  if (this->pMx != 0) sumMagCurrents();

  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      for (int k=low[2]; k<high[2]; ++k)
      {
        if (this->pMx != 0)
        {
          jx = (*this->pMx)(i,j,k);
          jy = (*this->pMy)(i,j,k);
          jz = (*this->pMz)(i,j,k);
        }
        
        this->plrcStepB(dt, i, j, k, dx, dy, dz, jx, jy, jz);
      }
        
  this->storage->applyBoundary("B");
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::sumCurrents()
{
  DataGrid &jxT = *this->pJx;
  DataGrid &jyT = *this->pJy;
  DataGrid &jzT = *this->pJz;
  
  jxT = 0;
  jyT = 0;
  jzT = 0;
  
  for 
  (
    typename CurrentList::iterator it = this->currents.begin(); 
    it != this->currents.end(); 
    ++it
  )
  {
    const DataGrid &jx = *(*it)->getJx();
    const DataGrid &jy = *(*it)->getJy();
    const DataGrid &jz = *(*it)->getJz();
    
    GridIndex low = jx.getLow();
    GridIndex high = jx.getHigh();
    for (int i=low[0]; i<=high[0]; ++i)
      for (int j=low[1]; j<=high[1]; ++j)
        for (int k=low[2]; k<=high[2]; ++k)
        {
          jxT(i,j,k) += jx(i,j,k);
          jyT(i,j,k) += jy(i,j,k);
          jzT(i,j,k) += jz(i,j,k);
        }
  } 
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::sumMagCurrents()
{
  DataGrid &jxT = *this->pMx;
  DataGrid &jyT = *this->pMy;
  DataGrid &jzT = *this->pMz;
  
  jxT = 0;
  jyT = 0;
  jzT = 0;
  
  for 
  (
    typename CurrentList::iterator it = this->magCurrents.begin(); 
    it != this->magCurrents.end(); 
    ++it
  )
  {
    const DataGrid &jx = *(*it)->getJx();
    const DataGrid &jy = *(*it)->getJy();
    const DataGrid &jz = *(*it)->getJz();
    
    GridIndex low = jx.getLow();
    GridIndex high = jx.getHigh();
    for (int i=low[0]; i<=high[0]; ++i)
      for (int j=low[1]; j<=high[1]; ++j)
        for (int k=low[2]; k<=high[2]; ++k)
        {
          jxT(i,j,k) += jx(i,j,k);
          jyT(i,j,k) += jy(i,j,k);
          jzT(i,j,k) += jz(i,j,k);
        }
  } 
}

template<class PLRCImplementation>
ParameterMap* FDTD_PLRCSolver<PLRCImplementation>::MakeParamMap (ParameterMap* pm)
{
  pm = FieldSolver::MakeParamMap(pm);
  pm = Implementation::CustomParamMap(pm);

  (*pm)["eps"] = WParameter(new ParameterValue<double>(&this->eps,1.0));

  (*pm)["L1eps"] = WParameter(new ParameterValue<double>(&this->LEps[0],0.0));
  (*pm)["L2eps"] = WParameter(new ParameterValue<double>(&this->LEps[1],0.0));
  (*pm)["L3eps"] = WParameter(new ParameterValue<double>(&this->LEps[2],0.0));
  
  (*pm)["L1delta"] = WParameter(new ParameterValue<double>(&this->LDelta[0],0.0));
  (*pm)["L2delta"] = WParameter(new ParameterValue<double>(&this->LDelta[1],0.0));
  (*pm)["L3delta"] = WParameter(new ParameterValue<double>(&this->LDelta[2],0.0));
  
  (*pm)["L1Om2"] = WParameter(new ParameterValue<double>(&this->LOm2[0],1.0));
  (*pm)["L2Om2"] = WParameter(new ParameterValue<double>(&this->LOm2[1],1.0));
  (*pm)["L3Om2"] = WParameter(new ParameterValue<double>(&this->LOm2[2],1.0));
  
  (*pm)["plasma"] = WParameter(
      new ParameterRebuild<PlasmaDensity, OptField>(&this->optfields)
  );
  
  (*pm)["plasma_current"] = WParameter(
      new ParameterRebuild<PlasmaCurrentFactory, CurrentFactory>(&this->currentFactories)
  );
  
  (*pm)["cpml_border"] = WParameter(
      new ParameterRebuild<CPMLBorder, CurrentFactory>(&this->currentFactories)
  );

  (*pm)["cpml_border_1d"] = WParameter(
      new ParameterRebuild<CPMLBorderOneD, CurrentFactory>(&this->currentFactories)
  );
  
  (*pm)["side_inject"] = WParameter(
      new ParameterRebuild<SideInject, CurrentFactory>(&this->currentFactories)
  );

  (*pm)["gauss_inject"] = WParameter(
      new ParameterRebuild<GaussInject, CurrentFactory>(&this->currentFactories)
  );

  (*pm)["short_pulse_inject"] = WParameter(
      new ParameterRebuild<ShortPulseInject, CurrentFactory>(&this->currentFactories)
  );

  (*pm)["wave_inject"] = WParameter(
      new ParameterRebuild<PlaneWaveSource, CurrentFactory>(&this->currentFactories)
  );

  return pm;
}
