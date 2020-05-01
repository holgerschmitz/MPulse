

//===============================================================
//==========  FDTD_PLRCSolver
//===============================================================


template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::init()
{
  PLRCImplementation::init();

  if ((this->LOm2[0] == 0) || (this->LOm2[1] == 0) || (this->LOm2[2] == 0))
  {
    std::cerr << "FDTD_PLRCSolver: Omega should not be zero!\n";
    exit(-1);
  }
}


template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepSchemeInit(double dt)
{
  initAccumulator(dt);

  stepB(0.5*dt);

  BOOST_FOREACH(pCurrent current, this->currents)
  {
    current->stepSchemeInit(dt);
  }

  BOOST_FOREACH(pCurrent current, this->magCurrents)
  {
    current->stepSchemeInit(dt);
  }
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepScheme(double dt)
{
  (*this->pSigma) = 0;


  BOOST_FOREACH(pCurrent current, this->currents)
  {
    current->stepScheme(dt);
  }

  stepD(dt);


  BOOST_FOREACH(pCurrent current, this->magCurrents)
  {
    current->stepScheme(dt);
  }

  stepB(dt);
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepD(double dt)
{
  Index low = this->pEx->getInnerLo();
  Index high = this->pEx->getInnerHi();

  double dx = getContext().getDx()[0];
  double dy = getContext().getDx()[1];
  double dz = getContext().getDx()[2];

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

  double jx, jy, jz;
  this->sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
      {
        jx = (*this->pJx)(i,j,k);
        jy = (*this->pJy)(i,j,k);
        jz = (*this->pJz)(i,j,k);

        this->plrcStepD(dt, i, j, k, dx, dy, dz, jx, jy, jz);
      }

  getContext().getSubdivision().exchange(*this->pEx);
  getContext().getSubdivision().exchange(*this->pEy);
  getContext().getSubdivision().exchange(*this->pEz);
}


template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::initAccumulator(double dt)
{
  Field &Ex = *this->pEx;
  Field &Ey = *this->pEy;
  Field &Ez = *this->pEz;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  // value of beta_p for the three Lorentz poles
  double beta[3];
  // value of gamma_p for the three Lorentz poles
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

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double ex = Ex(i,j,k);
    double ey = Ey(i,j,k);
    double ez = Ez(i,j,k);

    for (int n=0;n<3;++n)
    {
      double &pxr = this->pPsiRx[n]->operator()(i,j,k);
      double &pyr = this->pPsiRy[n]->operator()(i,j,k);
      double &pzr = this->pPsiRz[n]->operator()(i,j,k);

      double &pxi = this->pPsiIx[n]->operator()(i,j,k);
      double &pyi = this->pPsiIy[n]->operator()(i,j,k);
      double &pzi = this->pPsiIz[n]->operator()(i,j,k);

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
  Index low = this->pBx->getInnerLo();
  Index high = this->pBx->getInnerHi();

  double dx = getContext().getDx()[0];
  double dy = getContext().getDx()[1];
  double dz = getContext().getDx()[2];

  double jx, jy, jz;
  this->sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
      {
        jx = (*this->pMx)(i,j,k);
        jy = (*this->pMy)(i,j,k);
        jz = (*this->pMz)(i,j,k);

        this->plrcStepB(dt, i, j, k, dx, dy, dz, jx, jy, jz);
      }

  getContext().getSubdivision().exchange(*this->pBx);
  getContext().getSubdivision().exchange(*this->pBy);
  getContext().getSubdivision().exchange(*this->pBz);
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::initParameters(schnek::BlockParameters &blockPars)
{
  PLRCImplementation::initParameters(blockPars);

  blockPars.addParameter("eps", &this->eps, 1.0);

  blockPars.addParameter("L1eps", &this->LEps[0], 0.0);
  blockPars.addParameter("L2eps", &this->LEps[1], 0.0);
  blockPars.addParameter("L3eps", &this->LEps[2], 0.0);

  blockPars.addParameter("L1delta", &this->LDelta[0], 0.0);
  blockPars.addParameter("L2delta", &this->LDelta[1], 0.0);
  blockPars.addParameter("L3delta", &this->LDelta[2], 0.0);

  blockPars.addParameter("L1Om2", &this->LOm2[0], 1.0);
  blockPars.addParameter("L2Om2", &this->LOm2[1], 1.0);
  blockPars.addParameter("L3Om2", &this->LOm2[2], 1.0);

//  (*pm)["plasma"] = WParameter(
//      new ParameterRebuild<PlasmaDensity, OptField>(&this->optfields)
//  );
//
//  (*pm)["plasma_current"] = WParameter(
//      new ParameterRebuild<PlasmaCurrentFactory, CurrentFactory>(&this->currentFactories)
//  );
//
//  (*pm)["cpml_border"] = WParameter(
//      new ParameterRebuild<CPMLBorder, CurrentFactory>(&this->currentFactories)
//  );
//
//  (*pm)["cpml_border_1d"] = WParameter(
//      new ParameterRebuild<CPMLBorderOneD, CurrentFactory>(&this->currentFactories)
//  );
//
//  (*pm)["side_inject"] = WParameter(
//      new ParameterRebuild<SideInject, CurrentFactory>(&this->currentFactories)
//  );
//
//  (*pm)["gauss_inject"] = WParameter(
//      new ParameterRebuild<GaussInject, CurrentFactory>(&this->currentFactories)
//  );
//
//  (*pm)["short_pulse_inject"] = WParameter(
//      new ParameterRebuild<ShortPulseInject, CurrentFactory>(&this->currentFactories)
//  );
//
//  (*pm)["focused_pulse_inject"] = WParameter(
//      new ParameterRebuild<FocusedPulseInject, CurrentFactory>(&this->currentFactories)
//  );
//
//  (*pm)["wave_inject"] = WParameter(
//      new ParameterRebuild<PlaneWaveSource, CurrentFactory>(&this->currentFactories)
//  );
//
//  (*pm)["plane_gauss_inject"] = WParameter(
//      new ParameterRebuild<PlaneGaussSource, CurrentFactory>(&this->currentFactories)
//  );

}
