#include "plasmacurrent.hpp"

#include "../huerto/electromagnetics/fieldsolver.hpp"

#include <memory>
#include <chrono>

std::chrono::duration<double> total_time1(0);
std::chrono::duration<double> total_time2(0);
int iteration_count = 0;

void PlasmaCurrentBlock::initParameters(schnek::BlockParameters &blockPars)
{
  CurrentBlock::initParameters(blockPars);

  blockPars.addParameter("charge", &charge, 1.602176634e-19);
  blockPars.addParameter("mass", &mass, 9.1093837015e-31);
  blockPars.addParameter("gamma", &gamma, 0.01);
  blockPars.addParameter("Z", &Z, 1.0);
}

void PlasmaCurrentBlock::registerData() {
  plasmaCurrent = std::make_shared<PlasmaCurrent>(charge, mass, Z, gamma, boost::ref(*this));
  plasmaCurrent->registerData();
}

void PlasmaCurrentBlock::initCurrents(CurrentContainer &container)
{
  container.addCurrent(plasmaCurrent);
}

PlasmaCurrent::PlasmaCurrent(double charge_, double mass_, double Z_, double gamma_, CurrentBlock &plasmaBlock_)
  : plasmaBlock(plasmaBlock_), charge(charge_), mass(mass_), Z(Z_), gamma(gamma_)
{}

void PlasmaCurrent::registerData() {
  plasmaBlock.addData("PlasmaJx", Jx);
  plasmaBlock.addData("PlasmaJy", Jy);
  plasmaBlock.addData("PlasmaJz", Jz);
}

void PlasmaCurrent::init()
{
  schnek::DomainSubdivision<Field> &subdivision = plasmaBlock.getContext().getSubdivision();

  std::cout << "PlasmaCurrent::init called" << std::endl;

  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  plasmaBlock.retrieveData("Ex", Ex);
  plasmaBlock.retrieveData("Ey", Ey);
  plasmaBlock.retrieveData("Ez", Ez);

  plasmaBlock.retrieveData("Rho", Rho);


  Jx.resize(lowIn, highIn);
  Jy.resize(lowIn, highIn);
  Jz.resize(lowIn, highIn);
}

void PlasmaCurrent::stepScheme(double dt)
{
  Index low = Jx.getLo();
  Index high = Jx.getHi();

  const double gdtn = 1-0.5*gamma*dt;
  const double gdtd = 1+0.5*gamma*dt;
  const double emdt = dt*Z*Z*charge*charge/mass;

#ifdef HUERTO_ONE_DIM
  Range1d::LimitType low1D(low[0]), high1D(high[0]);
  Range1d range(low1D, high1D);

  // std::cout<<"PlasmaCurrent:stepScheme 1D Range"<<std::endl;

  for (auto &pos : range) {
    double &jx = Jx[pos];
    double &jy = Jy[pos];
    double &jz = Jz[pos];
    double rho = Rho[pos];

    jx = (jx*gdtn + emdt*Ex[pos]*rho)/gdtd;
    jy = (jy*gdtn + emdt*Ey[pos]*rho)/gdtd;
    jz = (jz*gdtn + emdt*Ez[pos]*rho)/gdtd;
  }
#endif

#ifdef HUERTO_TWO_DIM
  Range2d::LimitType low2D(low[0], low[1]), high2D(high[0], high[1]);
  Range2d range(low2D, high2D);

  iteration_count++;

  // Case 1
  auto start = std::chrono::high_resolution_clock::now();

  for (int i=low[0]; i<high[0]; ++i) {
    for (int j=low[1]; j<high[1]; ++j) {
      double &jx = Jx(i,j);
      double &jy = Jy(i,j);
      double &jz = Jz(i,j);
      double rho = Rho(i,j);

      jx = (jx*gdtn + emdt*Ex(i,j)*rho)/gdtd;
      jy = (jy*gdtn + emdt*Ey(i,j)*rho)/gdtd;
      jz = (jz*gdtn + emdt*Ez(i,j)*rho)/gdtd;
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  
  std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time1 += duration;

  // Case 2
  start = std::chrono::high_resolution_clock::now();

  for (auto &pos : range) {
      double &jx = Jx[pos];
      double &jy = Jy[pos];
      double &jz = Jz[pos];
      double rho = Rho[pos];

      jx = (jx*gdtn + emdt*Ex[pos]*rho)/gdtd;
      jy = (jy*gdtn + emdt*Ey[pos]*rho)/gdtd;
      jz = (jz*gdtn + emdt*Ez[pos]*rho)/gdtd;
  }

  end = std::chrono::high_resolution_clock::now();
  
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time2 += duration;

  if (iteration_count % 100 == 0){
    std::cout << "Total execution time (1): " << total_time1.count() << " ms" << std::endl;
    std::cout << "Total execution time (2): " << total_time2.count() << " ms" << std::endl;
    iteration_count = 0;
  }

#endif

#ifdef HUERTO_THREE_DIM
  std::cout<<"PlasmaCurrent:stepScheme 3D Range"<<std::endl;

  for (int i=low[0]; i<high[0]; ++i) {
    for (int j=low[1]; j<high[1]; ++j) {
      for (int k=low[2]; k<high[2]; ++k) {
        double &jx = Jx(i,j,k);
        double &jy = Jy(i,j,k);
        double &jz = Jz(i,j,k);
        double rho = Rho(i,j,k);

        jx = (jx*gdtn + emdt*Ex(i,j,k)*rho)/gdtd;
        jy = (jy*gdtn + emdt*Ey(i,j,k)*rho)/gdtd;
        jz = (jz*gdtn + emdt*Ez(i,j,k)*rho)/gdtd;
      }
    }
  }
#endif

}

