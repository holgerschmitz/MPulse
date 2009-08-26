#include "shortpulseinject.h"
#include <cmath>
#include "specfunc.h"

//===============================================================
//==========  FocusedPulseInject
//===============================================================

IncidentSourceCurrent *FocusedPulseInject::makeECurrent(int distance_, Direction dir_)
{
  std::cerr << "FocusedPulseInject::makeECurrent\n";
  typedef IncidentSourceECurrent<FocusedPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  
  cur->setParam(length, width, om0, TShift, ZShift, Phase, amp, eps, distance_, &generator);
  return cur;
}

IncidentSourceCurrent *FocusedPulseInject::makeHCurrent(int distance_, Direction dir_)
{
  std::cerr << "FocusedPulseInject::makeHCurrent\n";
  typedef IncidentSourceHCurrent<FocusedPulseInjectSourceFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance_,dir_);
  cur->setParam(length, width, om0, TShift, ZShift, Phase, amp, eps, distance_, &generator);
  return cur;
}

bool FocusedPulseInject::needCurrent(Direction dir_)
{
  return true;
}
    
ParameterMap* FocusedPulseInject::MakeParamMap (ParameterMap* pm)
{
  std::cerr << "FocusedPulseInject::MakeParamMap\n";
  pm = IncidentSource::MakeParamMap(pm);
  
  (*pm)["length"] = WParameter(new ParameterValue<double>(&this->length,1.));
  (*pm)["width"] = WParameter(new ParameterValue<double>(&this->width,1.));
  (*pm)["om0"] = WParameter(new ParameterValue<double>(&this->om0,2*M_PI));
  (*pm)["TShift"] = WParameter(new ParameterValue<double>(&this->TShift,0));
  (*pm)["ZShift"] = WParameter(new ParameterValue<double>(&this->ZShift,0));
  (*pm)["Phase"] = WParameter(new ParameterValue<double>(&this->Phase,0));

  (*pm)["amp"] = WParameter(new ParameterValue<double>(&this->amp,1));
  (*pm)["eps"] = WParameter(new ParameterValue<double>(&this->eps,1));
  
  return pm;
}


//===============================================================
//==========  FocusedPulseInjectSourceFunc
//===============================================================

void FocusedPulseInjectSourceFunc::setParam(double length_,
                                          double width_,
                                          double om0_,
                                          double TShift_,
                                          double ZShift_,
                                          double Phase_,
                                          double amp_, 
                                          double eps_,
                                          int distance_, 
                                          FocusedPulseDataGenerator *generator_)
{
  std::cerr << "FocusedPulseInjectSourceFunc::setParam\n";
  
  generator = generator_;
  
  // Grid Spacing and position
  
  DX = Globals::instance().gridDX();
  DY = Globals::instance().gridDY();
  DZ = Globals::instance().gridDZ();
  DT = Globals::instance().dt();


  GridIndex gridLow = Globals::instance().gridLow();
  GridIndex gridHigh = Globals::instance().gridHigh();
    
  centrex = 0.5*double(gridHigh[0] + gridLow[0]);
  centrey = 0.5*double(gridHigh[1] + gridLow[1]);
  centrez = 0.5*double(gridHigh[2] + gridLow[2]);
  
  // setting most parameters

  length = length_;
  width  = width_;
  om0     = om0_;
  
  eps = eps_;
  dist = distance_;
  
  lightspeed = sqrt(1/eps);
  
  ZRl = 0.5*om0*width*width/lightspeed;
  YComp = Complex(0.0,0.0);

  // Adjusting so that amplitude corresponds to maximum amplitude in the
  // center of the pulse
  
  Phase  = 0.0; //0.5*M_PI;
  TShift  = 0.0; 
  ZShift  = 0.0;

  double Exmax = Efunc(0, 0, 0, 0).real();

  // setting the rest of the parameters
  amp = amp_/Exmax;
  Phase  = 2*M_PI*Phase_;
  
  TShift  = TShift_;
  ZShift  = ZShift_;
}

void FocusedPulseInjectSourceFunc
    ::initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz)
{
  generator->setLow(storage->getLow());
  generator->setHigh(storage->getHigh());
  
  std::cerr << "FocusedPulseInjectSourceFunc::initSourceFunc\n";
  if (isH)
  {
    x_grid  = storage->addBorderLayer("IncidentBX" , dir, 1, dist, dist);
    y_grid  = storage->addBorderLayer("IncidentBY" , dir, 1, dist, dist);
    z_grid  = storage->addBorderLayer("IncidentBZ" , dir, 1, dist, dist);
    
    generator->addBXData(x_grid);
    generator->addBYData(y_grid);
    generator->addBZData(z_grid);
  }
  else
  {
    x_grid  = storage->addBorderLayer("IncidentEX" , dir, 1, dist, dist);
    y_grid  = storage->addBorderLayer("IncidentEY" , dir, 1, dist, dist);
    z_grid  = storage->addBorderLayer("IncidentEZ" , dir, 1, dist, dist);
    
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
  return Vector(x_grid(i,j,k), y_grid(i,j,k), z_grid(i,j,k));
}

Vector FocusedPulseInjectSourceFunc::getHField(int i, int j, int k, int time)
{
  return Vector(x_grid(i,j,k), y_grid(i,j,k), z_grid(i,j,k));
}

//===============================================================
//==========  FocusedPulseDataGenerator
//===============================================================

void FocusedPulseDataGenerator::addEXData(DataGrid *g)
{
  ex_grids.push_back(g);
}

void FocusedPulseDataGenerator::addEYData(DataGrid *g)
{
  ey_grids.push_back(g);
}

void FocusedPulseDataGenerator::addEZData(DataGrid *g)
{
  ez_grids.push_back(g);
}

void FocusedPulseDataGenerator::addBXData(DataGrid *g)
{
  bx_grids.push_back(g);
}

void FocusedPulseDataGenerator::addBYData(DataGrid *g)
{
  by_grids.push_back(g);
}

void FocusedPulseDataGenerator::addBZData(DataGrid *g)
{
  bz_grids.push_back(g);
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

void FocusedPulseDataGenerator::calcEx()
{
  for (int z=low[2]; z<=high[2]; ++z)
  {
    for (int x=low[0]; x<=high[0]; ++x)
    {
      for (int y=low[1]; y<=high[1]; ++y)
      {
        double posxh = (i+0.5-centrex)*DX;
        double posyo = (j-centrey)*DY;
        double poszo = (k-centrez)*DZ - ZShift;
        double posTime = time*DT - TShift;

        working_grid(x,y) = ExFunc(posxh, posyo, poszo, posTime);
        // 1. Fill working_grid with analytic values
        // 2. perform 2d FFT in x and y
        // 3. copy values from working_grid into respective boundary grids
      }
    }
  }
}

void FocusedPulseDataGenerator::calcEy();
void FocusedPulseDataGenerator::calcEz();
void FocusedPulseDataGenerator::calcBx();
void FocusedPulseDataGenerator::calcBy();
void FocusedPulseDataGenerator::calcBz();
