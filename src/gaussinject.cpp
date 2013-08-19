#include "gaussinject.h"

bool testNaN(double x)
{
  return !((x<1) || (x>0));
}

//===============================================================
//==========  GaussInject
//===============================================================

IncidentSourceCurrent *GaussInject::makeECurrent(int distance_, Direction dir_)
{
  typedef IncidentSourceECurrent<GaussInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(SizeT, length, width, k0, Shift, Time, amp, eps, toff);
  return cur;
}

IncidentSourceCurrent *GaussInject::makeHCurrent(int distance_, Direction dir_)
{
  typedef IncidentSourceHCurrent<GaussInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(SizeT, length, width, k0, Shift, Time, amp, eps, toff);
  return cur;
}

bool GaussInject::needCurrent(Direction dir_)
{
  return (dir_ == down);
}
    
ParameterMap* GaussInject::MakeParamMap (ParameterMap* pm)
{
  pm = IncidentSource::MakeParamMap(pm);
  
  (*pm)["SizeT"] = WParameter(new ParameterValue<int>(&this->SizeT,4096));
  (*pm)["length"] = WParameter(new ParameterValue<double>(&this->length,1.));
  (*pm)["width"] = WParameter(new ParameterValue<double>(&this->width,1.));
  (*pm)["k0"] = WParameter(new ParameterValue<double>(&this->k0,2*M_PI));
  (*pm)["Shift"] = WParameter(new ParameterValue<double>(&this->Shift,0));
  (*pm)["Time"] = WParameter(new ParameterValue<double>(&this->Time,0));

  (*pm)["amp"] = WParameter(new ParameterValue<double>(&this->amp,1));
  (*pm)["eps"] = WParameter(new ParameterValue<double>(&this->eps,1));
  
  (*pm)["toff"] = WParameter(new ParameterValue<int>(&this->toff,0));
  return pm;
}


//===============================================================
//==========  GaussInjectSourceFunc
//===============================================================

void GaussInjectSourceFunc::setParam(int    SizeT_, 
                                    double length_, 
                                    double width_,
                                    double k0_,
                                    double Shift_, 
                                    double Time_,
                                    double amp_, 
                                    double eps_, 
                                    int toff_)
{
  SizeT  = SizeT_;
  length = length_;
  width  = width_;
  k0     = k0_;
  Shift  = Shift_;
  Time   = Time_;
  
  gaussfunc.k0 = k0;
  gaussfunc.width = width;
  gaussfunc.length = length;
  
  amp = amp_;
  eps = eps_;
  toff = toff_;
  
  lightspeed = sqrt(1/eps);
  
  SizeX = Globals::instance().gridX()+2;
  SizeY = Globals::instance().gridY()+2;
  
  SizeTH = SizeT/2;
  SizeXH = SizeX/2;
  SizeYH = SizeY/2;

  DX = Globals::instance().gridDX();
  DY = Globals::instance().gridDY();
  DT = Globals::instance().dt() * lightspeed;

  BoxSizeX = SizeX*DX;
  BoxSizeY = SizeY*DY;
  BoxSizeT = SizeT*DT;

  A.resize(Row::IndexType(SizeT));
  
  makefftw();

  blockCount = 0;
  blockToff = toff;
  Nt = 1;
  active = false;
}

void GaussInjectSourceFunc
    ::initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz)
{
  switch (dir)
  {
    case east:  
    case west:  
    case north: 
    case south: std::cerr << "INTERNAL ERROR!\n";
                exit(-1);
                break;
    case up:     
    case down:  dim = 2;
                transverse1 = 0;
                transverse2 = 1;
                break;
  }
  
  GridIndex low  = storage->getLow();
  GridIndex high = storage->getHigh();
  
  lowx = low[0];
  lowy = low[1];
  highx = high[0];
  highy = high[1];
  
  Nt = 0x1 << int(round(log(high[2]-low[2])/log(2) + 1e-6));
//  Nt = SizeT;
  
  std::cerr << "Block Size is " << Nt << "\n";
    
  tmpField_TYr.resize(GridIndex2d(SizeY, Nt));
  tmpField_TYi.resize(GridIndex2d(SizeY, Nt));

  tmpField_TXr.resize(GridIndex3d(0, lowy, 0), GridIndex3d(SizeX-1, highy, Nt-1));
  tmpField_TXi.resize(GridIndex3d(0, lowy, 0), GridIndex3d(SizeX-1, highy, Nt-1));
  
  F1.resize(GridIndex3d(lowx, lowy, 0), GridIndex3d(highx, highy, Nt-1));
  F2.resize(GridIndex3d(lowx, lowy, 0), GridIndex3d(highx, highy, Nt-1));
  
//  debug_om_k.resize(GridIndex2d(SizeY, SizeT));
//  debug_t_k.resize(GridIndex2d(SizeY, SizeT));
}

void GaussInjectSourceFunc::setTime(int time_)
{
  if (time_ < (toff + blockCount*Nt)) return;
  if (time_ >= (toff + SizeT))
  {
    active = false;
    F1.resize(GridIndex(1,1,1));
    F2.resize(GridIndex(1,1,1));
    return;
  }
  
// here we have to generate the fields
  std::cerr << "Calculating Source ...";
  if (!isH)
  {
    gaussfunc.setShift(0, 0.5*DX, Shift+0.5*DX);
    fillBx();
  }
  else
  {
    Time -= 0.5*DT;
    gaussfunc.setShift(0.5*DX, 0, Shift);
    fillEx();
    gaussfunc.setShift(0, 0.5*DX, Shift);
    fillEy();
    Time += 0.5*DT;
  }
  
  blockToff = toff + blockCount*Nt;
  blockCount++;
  active = true;
  std::cerr << " done\n";
}

Vector GaussInjectSourceFunc::getField(int i, int j, int k, int time, double factor)
{
  Vector Field(0,0,0);
  if (!active) return Field;
    
  int t = time - blockToff - F1.getLo()[2];
  
  GridIndex index(i,j,k);
    
  Field[transverse1] = factor*amp*F1(index[transverse1], index[transverse2], t);
  Field[transverse2] = factor*amp*F2(index[transverse1], index[transverse2], t);

//  diagsource << (isH?"H":"E") << " " << i << " " << j << " " << k << " "
//             << Field[0] << " " << Field[1] << " " << Field[2] << std::endl;

  return Field;
}

Vector GaussInjectSourceFunc::getEField(int i, int j, int k, int time)
{
  if (!active) return Vector(0,0,0);

  int t = time - blockToff - F1.getLo()[2];
    
  return Vector(amp*F1(i, j, t), amp*F2(i, j, t), 0);
}

Vector GaussInjectSourceFunc::getHField(int i, int j, int k, int time)
{
  
  if (!active) return Vector(0,0,0);

  int t = time - blockToff - F1.getLo()[2];
  
  return Vector(sqrt(eps)*amp*F1(i, j, t), 0, 0);
}

//===============================================================
//==========  GaussInjectSourceFunc: computation


void GaussInjectSourceFunc::makefftw()
{
  ukx = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SizeX);
  ux  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SizeX);

  uky = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SizeY);
  uy  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SizeY);
  
  uo = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SizeT);
  ut = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SizeT);
  
  pfftxb = fftw_plan_dft_1d(SizeX, ukx, ux, FFTW_FORWARD, FFTW_ESTIMATE);
  pfftyb = fftw_plan_dft_1d(SizeY, uky, uy, FFTW_FORWARD, FFTW_ESTIMATE);
  pffttb = fftw_plan_dft_1d(SizeT, uo, ut, FFTW_BACKWARD, FFTW_ESTIMATE);
}

void GaussInjectSourceFunc::clearfftw()
{
  std::cerr << "========DESTROYING PLANS!!!!========\n";
  fftw_destroy_plan(pfftxb);
  fftw_destroy_plan(pfftyb);
  fftw_destroy_plan(pffttb);
  fftw_free(ukx);
  fftw_free(ux);
  fftw_free(uky);
  fftw_free(uy);
  fftw_free(uo);
  fftw_free(ut);
}

void GaussInjectSourceFunc::transformT(int i, int j)
{
  double f =(2*M_PI*length/BoxSizeT);
  for (int k=0; k<SizeT; ++k)
  {
    int kh = (k+SizeTH)%SizeT;
    Complex c = A[kh];
    uo[k][0] = std::real(c);
    uo[k][1] = std::imag(c);
    
  }
  
  fftw_execute(pffttb);
  
  int blockmin = blockCount*Nt;
  for (int kh=0; kh<Nt; ++kh)
  {
    int k = (kh+blockmin+SizeTH)%SizeT;
    tmpField_TYr(j,kh) = f*ut[k][0];
    tmpField_TYi(j,kh) = f*ut[k][1];
  }
}

void GaussInjectSourceFunc::transformY(int i)
{
  double f = (width/BoxSizeY);

  for (int k=0; k<Nt; ++k)
  {
    for (int j=0; j<SizeY; ++j)
    {
      int jh = (j+SizeYH)%SizeY;
      uky[jh][0] = tmpField_TYr(j,k);
      uky[jh][1] = tmpField_TYi(j,k);
      
    }
    fftw_execute(pfftyb);
    for (int j=lowy; j<=highy; ++j)
    {
      int jh = (j+SizeYH)%SizeY;
      tmpField_TXr(i, j, k) = f*uy[jh][0];
      tmpField_TXi(i, j, k) = f*uy[jh][1];
    }
  }
}

void GaussInjectSourceFunc::transformX(DataGrid &F)
{
  double f = (width/BoxSizeX);
  
  for (int j=lowy; j<=highy; ++j)
  {
    for (int k=0; k<Nt; ++k)
    {
      for (int i=0; i<SizeX; ++i)
      {
        int ih = (i+SizeXH)%SizeX;
        ukx[ih][0] = tmpField_TXr(i,j,k);
        ukx[ih][1] = tmpField_TXi(i,j,k);
      }
      
      fftw_execute(pfftxb);
            
      for (int i=lowx; i<=highx; ++i)
      {
        int ih = (i+SizeXH)%SizeX;
        F(i, j, k) = f*ux[ih][0];
      }
    }
  }
}

void GaussInjectSourceFunc::fillEx()
{
  Complex I(0,1);
  for (int i=0; i<SizeX; ++i)
  {  
    for (int j=0; j<SizeY; ++j)
    {  
      for (int k=0; k<SizeT; ++k)
      {  
        double kx = 2*M_PI*(i-SizeXH)/(BoxSizeX);
        double ky = 2*M_PI*(j-SizeYH)/(BoxSizeY);
        double om = 2*M_PI*(k-SizeTH)/(BoxSizeT);

//        if ((om*om <= (ky*ky + kx*kx)) || (om==0)) A[k] = 0;
//        else
//        {
//          double kz = sqrt(om*om - ky*ky - kx*kx)*sign(om);
//          if (testNaN(kz)) kz=0;

        double K = om*om - ky*ky - kx*kx;
        if (fabs(K)<1e-9) K=0;
        if (K<0) A[k] = Complex(0.0,0.0);
        else
        {
          double kz = sqrt(K)*sign(om);
          
          double a = gaussfunc.excomp(kx,ky,kz)*gaussfunc.gauss(kx,ky,kz);
          double ph = gaussfunc.shift(kx,ky,kz) + 0.5*M_PI;

          A[k] = Complex(a,0)*exp(-I*(ph+Time*om));
          
//          if (i==SizeXH) debug_om_k(j,k) = std::real(A[k]);
        }
        
      }

      transformT(i,j);
//      if (i==SizeXH) {
//        for (int k=0; k<SizeT; ++k) debug_t_k(j,k) = tmpField_TYr(j,k);
//      }
    }
    transformY(i);
  }
  transformX(F1);
  
//  std::ofstream fdebug_om_k("debug_om_k.out");
//  std::ofstream fdebug_t_k ("debug_t_k.out");
//  std::ofstream fdebug_t_y ("debug_t_y.out");
//  std::ofstream fdebug_t_x ("debug_t_x.out");
//  
//  for (int j=0; j<SizeY; ++j)
//  {  
//    for (int k=0; k<SizeT; ++k)
//    {
//      fdebug_om_k << j << " " << k << " " << debug_om_k(j,k) << "\n";
//      fdebug_t_k  << j << " " << k << " " << debug_t_k (j,k) << "\n";
//      fdebug_t_y  << j << " " << k << " " << tmpField_TXr(SizeXH,j,k) << "\n";
//      fdebug_t_x  << j << " " << k << " " << F1(SizeXH,j,k) << "\n";
//    }
//    fdebug_om_k << "\n";
//    fdebug_t_k  << "\n";
//    fdebug_t_y  << "\n";
//    fdebug_t_x  << "\n";
//  }
//  fdebug_om_k.close();
//  fdebug_t_k.close();
//  fdebug_t_y.close();
//  fdebug_t_x.close();
//  exit(0);

}

void GaussInjectSourceFunc::fillEy()
{  
  Complex I(0,1);
  for (int i=0; i<SizeX; ++i)
  {  
    for (int j=0; j<SizeY; ++j)
    {  
      for (int k=0; k<SizeT; ++k)
      {  
        double kx = 2*M_PI*(i-SizeXH)/(BoxSizeX);
        double ky = 2*M_PI*(j-SizeYH)/(BoxSizeY);
        double om = 2*M_PI*(k-SizeTH)/(BoxSizeT);

//        if ((om*om <= (ky*ky + kx*kx)) || (om==0)) A[k] = 0;
//        else
//        {
//          double kz = sqrt(om*om - ky*ky - kx*kx)*sign(om);
//          if (testNaN(kz)) kz=0;

        double K = om*om - ky*ky - kx*kx;
        if (fabs(K)<1e-9) K=0;
        if (K<0) A[k] = Complex(0.0,0.0);
        else
        {
          double kz = sqrt(K)*sign(om);

          double a = gaussfunc.eycomp(kx,ky,kz)*gaussfunc.gauss(kx,ky,kz);
          double ph = gaussfunc.shift(kx,ky,kz) + 0.5*M_PI;

          A[k] = Complex(a,0)*exp(-I*(ph+Time*om));
        }
      }

      transformT(i,j);
    }
    transformY(i);
  }
  transformX(F2);
}

void GaussInjectSourceFunc::fillBx()
{
  Complex I(0,1);
  for (int i=0; i<SizeX; ++i)
  {  
    for (int j=0; j<SizeY; ++j)
    {  
      for (int k=0; k<SizeT; ++k)
      {  
        double kx = 2*M_PI*(i-SizeXH)/(BoxSizeX);
        double ky = 2*M_PI*(j-SizeYH)/(BoxSizeY);
        double om = 2*M_PI*(k-SizeTH)/(BoxSizeT);

//        if ((om*om <= (ky*ky + kx*kx)) || (om==0)) A[k] = 0;
//        else
//        {
//          double kz = sqrt(om*om - ky*ky - kx*kx)*sign(om);
//          if (testNaN(kz)) kz=0;

        double K = om*om - ky*ky - kx*kx;
        if (fabs(K)<1e-9) K=0;
        if (K<0) A[k] = Complex(0.0,0.0);
        else
        {
          double kz = sqrt(K)*sign(om);

          double a = gaussfunc.bxcomp(kx,ky,kz)*gaussfunc.gauss(kx,ky,kz);
          double ph = gaussfunc.shift(kx,ky,kz) + 0.5*M_PI;

          A[k] = Complex(a,0)*exp(-I*(ph+Time*om));        
        }
      }

      transformT(i,j);
    }
    transformY(i);
  }
  transformX(F1);
  

  
}
