#include "focusedpulseinject.hpp"

#include "specfunc.hpp"

#include <cmath>

//===============================================================
//==========  FocusedPulseInject
//===============================================================

void FocusedPulseInject::initCurrents(CurrentContainer &container)
{
  generator.setSize(oversample_X, oversample_Y);
  generator.setShifts(TShift, ZShift, 2*M_PI*Phase);
  IncidentSource::initCurrents(container);
}

pCurrent FocusedPulseInject::makeECurrent(int distance_, Direction dir_)
{
  std::cerr << "FocusedPulseInject::makeECurrent\n";
  typedef IncidentSourceECurrent<FocusedPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  
  cur->setParam(length, width, om0, amp, eps, distance_, &generator);
  return pCurrent(cur);
}

pCurrent FocusedPulseInject::makeHCurrent(int distance_, Direction dir_)
{
  std::cerr << "FocusedPulseInject::makeHCurrent\n";
  typedef IncidentSourceHCurrent<FocusedPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(length, width, om0, amp, eps, distance_, &generator);
  return pCurrent(cur);
}

bool FocusedPulseInject::needCurrent(Direction dir_)
{
  return true;
}

void FocusedPulseInject::initParameters(schnek::BlockParameters &blockPars)
{
  IncidentSource::initParameters(blockPars);
  
  blockPars.addParameter("length", &this->length,1.);
  blockPars.addParameter("width", &this->width,1.);
  blockPars.addParameter("om0", &this->om0,2*M_PI);
  blockPars.addParameter("TShift", &this->TShift,0.);
  blockPars.addParameter("ZShift", &this->ZShift,0.);
  blockPars.addParameter("Phase", &this->Phase,0.);
  blockPars.addParameter("OversampleX", &this->oversample_X,1);
  blockPars.addParameter("OversampleY", &this->oversample_Y,1);

  blockPars.addParameter("amp", &this->amp,1.);
  blockPars.addParameter("eps", &this->eps,1.);
}


//===============================================================
//==========  FocusedPulseInjectSourceFunc
//===============================================================

void FocusedPulseInjectSourceFunc::setParam(double length_,
                                          double width_,
                                          double om0_,
                                          double amp_, 
                                          double eps_,
                                          int distance_, 
                                          FocusedPulseDataGenerator *generator_)
{
  std::cerr << "FocusedPulseInjectSourceFunc::setParam\n";
  
  generator = generator_;
  
  // Grid Spacing and position
  Index gridLow = MPulse::getSubdivision().getLo();
  Index gridHigh = MPulse::getSubdivision().getHi();
  
  // setting most parameters

  length = length_;
  width  = width_;
  om0     = om0_;
  
  eps = eps_;
  dist = distance_;
  
  const int lightspeed = sqrt(1/eps);
  
  ZRl = 0.5*om0*width*width/lightspeed;
//  YComp = Complex(0.0,0.0);

//  double Exmax = Efunc(0, 0, 0, 0).real();

  // setting the rest of the parameters
//  amp = amp_/Exmax;
}

void FocusedPulseInjectSourceFunc::initSourceFunc(pGrid pJx, pGrid pJy, pGrid pJz)
{
  std::cerr << "FocusedPulseInjectSourceFunc::initSourceFunc\n";

  const schnek::DomainSubdivision<Field> &subdivision = MPulse::getSubdivision();
  Index lo = subdivision.getInnerLo();
  Index hi = subdivision.getInnerHi();
  generator->setLow(lo);
  generator->setHigh(hi);
  
  Index blow, bhigh;
  if (!getBorderExtent(dir, 1, dist, blow, bhigh, isH)) return;
  x_grid = boost::make_shared<Grid>(blow, bhigh);
  y_grid = boost::make_shared<Grid>(blow, bhigh);
  z_grid = boost::make_shared<Grid>(blow, bhigh);

  if (isH)
  {
    generator->addBXData(x_grid);
    generator->addBYData(y_grid);
    generator->addBZData(z_grid);
  }
  else
  {
    generator->addEXData(x_grid);
    generator->addEYData(y_grid);
    generator->addEZData(z_grid);
  }
}

void FocusedPulseInjectSourceFunc::setTime(int time)
{
  generator->setTime(time);
}

Vector FocusedPulseInjectSourceFunc::getEField(int i, int j, int k, int time)
{
  return Vector((*x_grid)(i,j,k), (*y_grid)(i,j,k), (*z_grid)(i,j,k));
}

Vector FocusedPulseInjectSourceFunc::getHField(int i, int j, int k, int time)
{
  return Vector((*x_grid)(i,j,k), (*y_grid)(i,j,k), (*z_grid)(i,j,k));
}

//===============================================================
//==========  FocusedPulseDataGenerator
//===============================================================

void FocusedPulseDataGenerator::addEXData(pGrid g)
{
  ex_grids.push_back(g);
}

void FocusedPulseDataGenerator::addEYData(pGrid g)
{
  ey_grids.push_back(g);
}

void FocusedPulseDataGenerator::addEZData(pGrid g)
{
  ez_grids.push_back(g);
}

void FocusedPulseDataGenerator::addBXData(pGrid g)
{
  bx_grids.push_back(g);
}

void FocusedPulseDataGenerator::addBYData(pGrid g)
{
  by_grids.push_back(g);
}

void FocusedPulseDataGenerator::addBZData(pGrid g)
{
  bz_grids.push_back(g);
}

void FocusedPulseDataGenerator::setSize(int oversample_X_, int oversample_Y_)
{
  oversample_X = oversample_X_;
  oversample_Y = oversample_Y_;
  Index low = MPulse::getSubdivision().getLo();
  Index high = MPulse::getSubdivision().getHi();

  DX = MPulse::getDx()[0];
  DY = MPulse::getDx()[1];
  DZ = MPulse::getDx()[2];
  DT = MPulse::getDt();
  
  sizeX = oversample_X*(high[0]-low[0]+1);
  sizeY = oversample_Y*(high[1]-low[1]+1);
  
  gLowX = low[0];
  gLowY = low[1];
  gHighX = high[0];
  gHighY = high[1];

  centrez = 0.5*double(high[2] + low[2]);

//  working_grid.resize(GridIndex2d(sizeX, sizeY));
//  working_grid_tmp.resize(GridIndex2d(sizeX, sizeY));

  fft_field_k = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY);
  fft_field   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY);
  
  pfft = fftw_plan_dft_2d(sizeX, sizeY, fft_field_k, fft_field, FFTW_FORWARD, FFTW_ESTIMATE);
}

void FocusedPulseDataGenerator::setShifts(double TShift_, double ZShift_, double Phase_)
{
  TShift = TShift_;
  ZShift = ZShift_;
  Phase = Phase_;
}

void FocusedPulseDataGenerator::setTime(int Time)
{
  if (!initialized) this->init();
  if (Time==currentTime) return;
  
  currentTime = Time;
  
  calcEx();
  calcEy();
  calcEz();
  calcBx();
  calcBy();
  calcBz();

}

void FocusedPulseDataGenerator::init()
{
  throw std::string("Something needs to happen here");
}

void FocusedPulseDataGenerator::calcEx()
{

  double DX = MPulse::getDx()[0];
  double DY = MPulse::getDx()[1];

  double dKx = 1./(sizeX*DX);
  double dKy = 1./(sizeY*DY);

  int sizeXh = sizeX/2;
  int sizeYh = sizeY/2;

  Complex I(0.0, 1.0);

  for (int z=low[2]; z<=high[2]; ++z)
  {
    // 1. Fill working_grid with analytic values

    //start at -1 because we use prefix increment
    int fftPos = -1;

    for (int x=0; x<sizeX; ++x)
    {
      double kx = ((x+sizeXh)%sizeX)*dKx;
      Complex Shift = Phase*exp(-I*0.5*DX*kx);

      // coordinates are shifted 0.5DX to the right,
      // so in local coord we have to shift function 0.5DX to the left

      for (int y=0; y<sizeY; ++y)
      {
        double ky = ((y+sizeYh)%sizeY)*dKy;

        double poszo = (z-centrez)*DZ - ZShift;
        double posTime = currentTime*DT - TShift;

        Complex Ex_k = Shift*FieldFuncs(kx, ky, poszo, posTime, 0); // 0 for Ex component
        fft_field_k[++fftPos][0] = Ex_k.real();
        fft_field_k[  fftPos][1] = Ex_k.imag();
      }
    }

    // 2. perform 2d FFT in x and y
    fftw_execute(pfft);

    // 3. copy values from working_grid into respective boundary grids

    // some values needed for the index transform
    int GhX = (gHighX - gLowX + 1)/2;
    int GhY = (gHighY - gLowY + 1)/2;
    int offsetX = sizeX-GhX-gLowX;
    int offsetY = sizeY-GhY-gLowY;


    for (GridList::iterator git=ex_grids.begin(); git!=ex_grids.end(); ++git)
    {
      Grid &grid = *(*git);

      Index low = grid.getLo();
      Index high = grid.getHi();

      if ((z<low[2]) || (z>high[2])) continue;

      // 3.1 What is the corresponding domain in configuration space. Remember,
      // oversampling will produce a domain that is larger than the simulation domain.

      // iterate over the boundary grid's domain
      for (int x=low[0]; x<=high[0]; ++x)
      {
        // coordinate in the fft field
        int xf = (offsetX+x) % sizeX;
        for (int y=low[1]; y<=high[1]; ++y)
        {
          // coordinate in the fft field
          int yf = (offsetY+y) % sizeY;
          // getting the real part
          grid(x,y,z) = fft_field[sizeY*xf+yf][0];
        }
      }

    }
  }
}

void FocusedPulseDataGenerator::calcEy()
{

  double DX = MPulse::getDx()[0];
  double DY = MPulse::getDx()[1];

  double dKx = 1./(sizeX*DX);
  double dKy = 1./(sizeY*DY);

  int sizeXh = sizeX/2;
  int sizeYh = sizeY/2;

  Complex I(0.0, 1.0);

  for (int z=low[2]; z<=high[2]; ++z)
  {
    // 1. Fill working_grid with analytic values

    //start at -1 because we use prefix increment
    int fftPos = -1;

    for (int x=0; x<sizeX; ++x)
    {
      double kx = ((x+sizeXh)%sizeX)*dKx;


      for (int y=0; y<sizeY; ++y)
      {
        double ky = ((y+sizeYh)%sizeY)*dKy;

        // coordinates are shifted 0.5DY to the right,
        // so in local coord we have to shift function 0.5DY to the left
        Complex Shift = Phase*exp(-I*0.5*DY*ky);

        double poszo = (z-centrez)*DZ - ZShift;
        double posTime = currentTime*DT - TShift;

        Complex Ey_k = Shift*FieldFuncs(kx, ky, poszo, posTime, 1); // 1 for Ey component
        fft_field_k[++fftPos][0] = Ey_k.real();
        fft_field_k[  fftPos][1] = Ey_k.imag();
      }
    }

    // 2. perform 2d FFT in x and y
    fftw_execute(pfft);

    // 3. copy values from working_grid into respective boundary grids

    // some values needed for the index transform
    int GhX = (gHighX - gLowX + 1)/2;
    int GhY = (gHighY - gLowY + 1)/2;
    int offsetX = sizeX-GhX-gLowX;
    int offsetY = sizeY-GhY-gLowY;


    for (GridList::iterator git=ey_grids.begin(); git!=ey_grids.end(); ++git)
    {
      Grid &grid = *(*git);

      Index low = grid.getLo();
      Index high = grid.getHi();

      if ((z<low[2]) || (z>high[2])) continue;

      // 3.1 What is the corresponding domain in configuration space. Remember,
      // oversampling will produce a domain that is larger than the simulation domain.

      // iterate over the boundary grid's domain
      for (int x=low[0]; x<=high[0]; ++x)
      {
        // coordinate in the fft field
        int xf = (offsetX+x) % sizeX;
        for (int y=low[1]; y<=high[1]; ++y)
        {
          // coordinate in the fft field
          int yf = (offsetY+y) % sizeY;
          // getting the real part
          grid(x,y,z) = fft_field[sizeY*xf+yf][0];
        }
      }

    }
  }
}

void FocusedPulseDataGenerator::calcEz()
{

  double DX = MPulse::getDx()[0];
  double DY = MPulse::getDx()[1];

  double dKx = 1./(sizeX*DX);
  double dKy = 1./(sizeY*DY);

  int sizeXh = sizeX/2;
  int sizeYh = sizeY/2;

  Complex I(0.0, 1.0);
  Complex Shift = Phase;

  for (int z=low[2]; z<=high[2]; ++z)
  {
    // 1. Fill working_grid with analytic values

    //start at -1 because we use prefix increment
    int fftPos = -1;

    for (int x=0; x<sizeX; ++x)
    {
      double kx = ((x+sizeXh)%sizeX)*dKx;

      for (int y=0; y<sizeY; ++y)
      {
        double ky = ((y+sizeYh)%sizeY)*dKy;

        double poszo = (z+0.5-centrez)*DZ - ZShift;
        double posTime = currentTime*DT - TShift;

        Complex Ez_k = Shift*FieldFuncs(kx, ky, poszo, posTime, 2); // 2 for Ez component
        fft_field_k[++fftPos][0] = Ez_k.real();
        fft_field_k[  fftPos][1] = Ez_k.imag();
      }
    }

    // 2. perform 2d FFT in x and y
    fftw_execute(pfft);

    // 3. copy values from working_grid into respective boundary grids

    // some values needed for the index transform
    int GhX = (gHighX - gLowX + 1)/2;
    int GhY = (gHighY - gLowY + 1)/2;
    int offsetX = sizeX-GhX-gLowX;
    int offsetY = sizeY-GhY-gLowY;


    for (GridList::iterator git=ez_grids.begin(); git!=ez_grids.end(); ++git)
    {
      Grid &grid = *(*git);

      Index low = grid.getLo();
      Index high = grid.getHi();

      if ((z<low[2]) || (z>high[2])) continue;

      // 3.1 What is the corresponding domain in configuration space. Remember,
      // oversampling will produce a domain that is larger than the simulation domain.

      // iterate over the boundary grid's domain
      for (int x=low[0]; x<=high[0]; ++x)
      {
        // coordinate in the fft field
        int xf = (offsetX+x) % sizeX;
        for (int y=low[1]; y<=high[1]; ++y)
        {
          // coordinate in the fft field
          int yf = (offsetY+y) % sizeY;
          // getting the real part
          grid(x,y,z) = fft_field[sizeY*xf+yf][0];
        }
      }

    }
  }
}

void FocusedPulseDataGenerator::calcBx()
{

  double DX = MPulse::getDx()[0];
  double DY = MPulse::getDx()[1];

  double dKx = 1./(sizeX*DX);
  double dKy = 1./(sizeY*DY);

  int sizeXh = sizeX/2;
  int sizeYh = sizeY/2;

  Complex I(0.0, 1.0);

  for (int z=low[2]; z<=high[2]; ++z)
  {
    // 1. Fill working_grid with analytic values

    //start at -1 because we use prefix increment
    int fftPos = -1;

    for (int x=0; x<sizeX; ++x)
    {
      double kx = ((x+sizeXh)%sizeX)*dKx;

      for (int y=0; y<sizeY; ++y)
      {
        double ky = ((y+sizeYh)%sizeY)*dKy;

        // coordinates are shifted 0.5DX to the right,
        // so in local coord we have to shift function 0.5DX to the left
        Complex Shift = Phase*exp(-I*0.5*DY*ky);

        double poszo = (z+0.5-centrez)*DZ - ZShift;
        double posTime = currentTime*DT - TShift;

        Complex Bx_k = Shift*FieldFuncs(kx, ky, poszo, posTime, 3); // 3 for Bx component
        fft_field_k[++fftPos][0] = Bx_k.real();
        fft_field_k[  fftPos][1] = Bx_k.imag();
      }
    }

    // 2. perform 2d FFT in x and y
    fftw_execute(pfft);

    // 3. copy values from working_grid into respective boundary grids

    // some values needed for the index transform
    int GhX = (gHighX - gLowX + 1)/2;
    int GhY = (gHighY - gLowY + 1)/2;
    int offsetX = sizeX-GhX-gLowX;
    int offsetY = sizeY-GhY-gLowY;


    for (GridList::iterator git=bx_grids.begin(); git!=bx_grids.end(); ++git)
    {
      Grid &grid = *(*git);

      Index low = grid.getLo();
      Index high = grid.getHi();

      if ((z<low[2]) || (z>high[2])) continue;

      // 3.1 What is the corresponding domain in configuration space. Remember,
      // oversampling will produce a domain that is larger than the simulation domain.

      // iterate over the boundary grid's domain
      for (int x=low[0]; x<=high[0]; ++x)
      {
        // coordinate in the fft field
        int xf = (offsetX+x) % sizeX;
        for (int y=low[1]; y<=high[1]; ++y)
        {
          // coordinate in the fft field
          int yf = (offsetY+y) % sizeY;
          // getting the real part
          grid(x,y,z) = fft_field[sizeY*xf+yf][0];
        }
      }

    }
  }
}

void FocusedPulseDataGenerator::calcBy()
{

  double DX = MPulse::getDx()[0];
  double DY = MPulse::getDx()[1];

  double dKx = 1./(sizeX*DX);
  double dKy = 1./(sizeY*DY);

  int sizeXh = sizeX/2;
  int sizeYh = sizeY/2;

  Complex I(0.0, 1.0);

  for (int z=low[2]; z<=high[2]; ++z)
  {
    // 1. Fill working_grid with analytic values

    //start at -1 because we use prefix increment
    int fftPos = -1;

    for (int x=0; x<sizeX; ++x)
    {
      double kx = ((x+sizeXh)%sizeX)*dKx;

      // coordinates are shifted 0.5DX to the right,
      // so in local coord we have to shift function 0.5DX to the left
      Complex Shift = Phase*exp(-I*0.5*DX*kx);

      for (int y=0; y<sizeY; ++y)
      {
        double ky = ((y+sizeYh)%sizeY)*dKy;

        double poszo = (z+0.5-centrez)*DZ - ZShift;
        double posTime = currentTime*DT - TShift;

        Complex By_k = Shift*FieldFuncs(kx, ky, poszo, posTime, 4); // 4 for By component
        fft_field_k[++fftPos][0] = By_k.real();
        fft_field_k[  fftPos][1] = By_k.imag();
      }
    }

    // 2. perform 2d FFT in x and y
    fftw_execute(pfft);

    // 3. copy values from working_grid into respective boundary grids

    // some values needed for the index transform
    int GhX = (gHighX - gLowX + 1)/2;
    int GhY = (gHighY - gLowY + 1)/2;
    int offsetX = sizeX-GhX-gLowX;
    int offsetY = sizeY-GhY-gLowY;


    for (GridList::iterator git=by_grids.begin(); git!=by_grids.end(); ++git)
    {
      Grid &grid = *(*git);

      Index low = grid.getLo();
      Index high = grid.getHi();

      if ((z<low[2]) || (z>high[2])) continue;

      // 3.1 What is the corresponding domain in configuration space. Remember,
      // oversampling will produce a domain that is larger than the simulation domain.

      // iterate over the boundary grid's domain
      for (int x=low[0]; x<=high[0]; ++x)
      {
        // coordinate in the fft field
        int xf = (offsetX+x) % sizeX;
        for (int y=low[1]; y<=high[1]; ++y)
        {
          // coordinate in the fft field
          int yf = (offsetY+y) % sizeY;
          // getting the real part
          grid(x,y,z) = fft_field[sizeY*xf+yf][0];
        }
      }

    }
  }
}

void FocusedPulseDataGenerator::calcBz()
{

  double DX = MPulse::getDx()[0];
  double DY = MPulse::getDx()[1];

  double dKx = 1./(sizeX*DX);
  double dKy = 1./(sizeY*DY);

  int sizeXh = sizeX/2;
  int sizeYh = sizeY/2;

  Complex I(0.0, 1.0);

  for (int z=low[2]; z<=high[2]; ++z)
  {
    // 1. Fill working_grid with analytic values

    //start at -1 because we use prefix increment
    int fftPos = -1;

    for (int x=0; x<sizeX; ++x)
    {
      double kx = ((x+sizeXh)%sizeX)*dKx;

      for (int y=0; y<sizeY; ++y)
      {
        double ky = ((y+sizeYh)%sizeY)*dKy;

        // coordinates are shifted 0.5DX to the right,
        // so in local coord we have to shift function 0.5DX to the left
        Complex Shift = Phase*exp(-I*0.5*(DX*kx + DY*ky));

        double poszo = (z-centrez)*DZ - ZShift;
        double posTime = currentTime*DT - TShift;

        Complex Bz_k = Shift*FieldFuncs(kx, ky, poszo, posTime, 5); // 5 for Bz component
        fft_field_k[++fftPos][0] = Bz_k.real();
        fft_field_k[  fftPos][1] = Bz_k.imag();
      }
    }

    // 2. perform 2d FFT in x and y
    fftw_execute(pfft);

    // 3. copy values from working_grid into respective boundary grids

    // some values needed for the index transform
    int GhX = (gHighX - gLowX + 1)/2;
    int GhY = (gHighY - gLowY + 1)/2;
    int offsetX = sizeX-GhX-gLowX;
    int offsetY = sizeY-GhY-gLowY;


    for (GridList::iterator git=bz_grids.begin(); git!=bz_grids.end(); ++git)
    {
      Grid &grid = *(*git);

      Index low = grid.getLo();
      Index high = grid.getHi();

      if ((z<low[2]) || (z>high[2])) continue;

      // 3.1 What is the corresponding domain in configuration space. Remember,
      // oversampling will produce a domain that is larger than the simulation domain.

      // iterate over the boundary grid's domain
      for (int x=low[0]; x<=high[0]; ++x)
      {
        // coordinate in the fft field
        int xf = (offsetX+x) % sizeX;
        for (int y=low[1]; y<=high[1]; ++y)
        {
          // coordinate in the fft field
          int yf = (offsetY+y) % sizeY;
          // getting the real part
          grid(x,y,z) = fft_field[sizeY*xf+yf][0];
        }
      }

    }
  }
}
