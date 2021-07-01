

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


  BOOST_FOREACH(pCurrent current, this->currents) {
    current->stepScheme(dt);
  }

  stepD(dt);


  BOOST_FOREACH(pCurrent current, this->magCurrents) {
    current->stepScheme(dt);
  }

  stepB(dt);
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepD(double dt)
{
  Index low = this->pEx->getInnerLo();
  Index high = this->pEx->getInnerHi();

  Vector dx = this->getContext().getDx();

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

  Vector j;
  this->sumCurrents();

  Index pos;

  for (pos[0]=low[0]; pos[0]<=high[0]; ++pos[0]) {
#ifndef HUERTO_ONE_DIM
    for (pos[1]=low[1]; pos[1]<=high[1]; ++pos[1]) {
#endif
#ifdef HUERTO_THREE_DIM
      for (pos[2]=low[2]; pos[2]<=high[2]; ++pos[2]) {
#endif
        double jx = (*this->pJx)[pos];
        double jy = (*this->pJy)[pos];
        double jz = (*this->pJz)[pos];

        this->plrcStepD(dt, pos, dx, jx, jy, jz);
#ifdef HUERTO_THREE_DIM
      }
#endif
#ifndef HUERTO_ONE_DIM
    }
#endif
  }

  this->getContext().getSubdivision().exchange(*this->pEx);
  this->getContext().getSubdivision().exchange(*this->pEy);
  this->getContext().getSubdivision().exchange(*this->pEz);
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
    // Kelley & Luebbers Eq (21)
    beta[n] = sqrt(this->LOm2[n] - this->LDelta[n]*this->LDelta[n]);
    gamma[n] = this->LEps[n]*this->LOm2[n]/beta[n];
    phasor[n] = this->LDelta[n] - I*beta[n];
    this->plrcData.Crec[n] = exp(-phasor[n]*dt);

    // Kelley & Luebbers Eq (23)
    xi0[n] = I*gamma[n] *(this->plrcData.Crec[n]*(phasor[n]*dt + 1.0) - 1.0)
              / (phasor[n]*phasor[n]*dt);

    this->plrcData.dxi0[n]  = xi0[n] *(1.0-this->plrcData.Crec[n]);

  }

  Index pos;
  for (pos[0]=low[0]; pos[0]<=high[0]; ++pos[0]) {
#ifndef HUERTO_ONE_DIM
    for (pos[1]=low[1]; pos[1]<=high[1]; ++pos[1]) {
#endif
#ifdef HUERTO_THREE_DIM
      for (pos[2]=low[2]; pos[2]<=high[2]; ++pos[2]) {
#endif
        double ex = Ex[pos];
        double ey = Ey[pos];
        double ez = Ez[pos];

        for (int n=0;n<3;++n)
        {
          double &pxr = (*this->pPsiRx[n])[pos];
          double &pyr = (*this->pPsiRy[n])[pos];
          double &pzr = (*this->pPsiRz[n])[pos];

          double &pxi = (*this->pPsiIx[n])[pos];
          double &pyi = (*this->pPsiIy[n])[pos];
          double &pzi = (*this->pPsiIz[n])[pos];

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
#ifdef HUERTO_THREE_DIM
      }
#endif
#ifndef HUERTO_ONE_DIM
    }
#endif
  }
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::stepB(double dt)
{
  Index low = this->pBx->getInnerLo();
  Index high = this->pBx->getInnerHi();

  Vector dx = this->getContext().getDx();

  double jx, jy, jz;
  this->sumMagCurrents();

  Index pos;
  for (pos[0]=low[0]; pos[0]<=high[0]; ++pos[0]) {
#ifndef HUERTO_ONE_DIM
    for (pos[1]=low[1]; pos[1]<=high[1]; ++pos[1]) {
#endif
#ifdef HUERTO_THREE_DIM
      for (pos[2]=low[2]; pos[2]<=high[2]; ++pos[2]) {
#endif
        jx = (*this->pMx)[pos];
        jy = (*this->pMy)[pos];
        jz = (*this->pMz)[pos];

        this->plrcStepB(dt, pos, dx, jx, jy, jz);
#ifdef HUERTO_THREE_DIM
      }
#endif
#ifndef HUERTO_ONE_DIM
    }
#endif
  }

  this->getContext().getSubdivision().exchange(*this->pBx);
  this->getContext().getSubdivision().exchange(*this->pBy);
  this->getContext().getSubdivision().exchange(*this->pBz);
}

template<class PLRCImplementation>
void FDTD_PLRCSolver<PLRCImplementation>::initParameters(schnek::BlockParameters &blockPars) {
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
}
