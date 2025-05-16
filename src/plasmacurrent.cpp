#include "plasmacurrent.hpp"

#include "../huerto/electromagnetics/fieldsolver.hpp"

#include <schnek/grid/gridstorage/kokkos-storage.hpp>
#include <schnek/grid/array.hpp>
#include <schnek/grid/grid.hpp>
#include <schnek/grid/field.hpp>


#include <memory>
#include <chrono>
#include <cstdlib>

std::chrono::duration<double> total_time1(0);
std::chrono::duration<double> total_time2(0);
std::chrono::duration<double> total_time3(0);
std::chrono::duration<double> total_time4(0);
double total_time41(0);
std::chrono::duration<double> total_time5(0);
double total_time51(0);
std::chrono::duration<double> total_time6(0);
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
  std::cout << "PlasmaCurrent::registerData called" << std::endl;

  plasmaBlock.addData("PlasmaJx", Jx);
  plasmaBlock.addData("PlasmaJy", Jy);
  plasmaBlock.addData("PlasmaJz", Jz);

  // plasmaBlock.addData("Ex", Ex);
  // plasmaBlock.addData("Ey", Ey);
  // plasmaBlock.addData("Ez", Ez);

  // plasmaBlock.addData("Rho", Rho);
}

void PlasmaCurrent::init()
{
  schnek::DomainSubdivision<Field> &subdivision = plasmaBlock.getContext().getSubdivision();

  // std::cout << "PlasmaCurrent::init called" << std::endl;

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
  
  // std::cout<<"PlasmaCurrent:stepScheme"<<std::endl;

#ifdef HUERTO_ONE_DIM
  Range1d::LimitType low1D(low[0]), high1D(high[0]);
  Range1d range(low1D, high1D);


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

  // std::cout << std::endl;
  auto dims = Jx.getDims();
  // std::cout << "Jx Dims : " << dims[0] << ", " << dims[1] << std::endl;
  dims = Jy.getDims();
  // std::cout << "Jy Dims : " << dims[0] << ", " << dims[1] << std::endl;
  dims = Jz.getDims();
  // std::cout << "Jz Dims : " << dims[0] << ", " << dims[1] << std::endl;
  dims = Rho.getDims();
  // std::cout << "Rho Dims : " << dims[0] << ", " << dims[1] << std::endl;
  
  srand(99);
  int rn = low[0] + (rand() % (high[0] - low[0] + 1));
  // std::cout << "(" << low[0] << ", " << high[0] << ") Random : " << rn << std::endl;
  rn = low[0];

  Grid Jx1, Jy1, Jz1;
  Field Ex1, Ey1, Ez1, Rho1;

  // Creating temporary vars
  Jx1.resize(Jx.getLo(), Jx.getHi());
  Jy1.resize(Jy.getLo(), Jy.getHi());
  Jz1.resize(Jz.getLo(), Jz.getHi());
  Ex1.resize(Ex.getLo(), Ex.getHi(), Ex.getDomain(), Ex.getStagger(), Ex.getghostCells());
  Ey1.resize(Ey.getLo(), Ey.getHi(), Ey.getDomain(), Ey.getStagger(), Ey.getghostCells());
  Ez1.resize(Ez.getLo(), Ez.getHi(), Ez.getDomain(), Ez.getStagger(), Ez.getghostCells());
  Rho1.resize(Rho.getLo(), Rho.getHi(), Rho.getDomain(), Rho.getStagger(), Rho.getghostCells());

  Grid Jx2, Jy2, Jz2;
  Field Ex2, Ey2, Ez2, Rho2;

  // Creating temporary vars
  Jx2.resize(Jx.getLo(), Jx.getHi());
  Jy2.resize(Jy.getLo(), Jy.getHi());
  Jz2.resize(Jz.getLo(), Jz.getHi());
  Ex2.resize(Ex.getLo(), Ex.getHi(), Ex.getDomain(), Ex.getStagger(), Ex.getghostCells());
  Ey2.resize(Ey.getLo(), Ey.getHi(), Ey.getDomain(), Ey.getStagger(), Ey.getghostCells());
  Ez2.resize(Ez.getLo(), Ez.getHi(), Ez.getDomain(), Ez.getStagger(), Ez.getghostCells());
  Rho2.resize(Rho.getLo(), Rho.getHi(), Rho.getDomain(), Rho.getStagger(), Rho.getghostCells());
  
  for (auto &pos : range) {
    Jx1[pos] = Jx[pos]; Jy1[pos] = Jy[pos]; Jz1[pos] = Jz[pos];
    Ex1[pos] = Ex[pos]; Ey1[pos] = Ey[pos]; Ez1[pos] = Ez[pos];
    Rho1[pos] = Rho[pos];
  }
  // Creating temporary vars

  if (Jx(rn, rn) != 0) std::cout << "(Before 1) Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(Before 1) Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(Before 1) Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;

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

  for (auto &pos : range) {
    Jx1[pos] = Jx[pos]; Jy1[pos] = Jy[pos]; Jz1[pos] = Jz[pos];
    Ex1[pos] = Ex[pos]; Ey1[pos] = Ey[pos]; Ez1[pos] = Ez[pos];
    Rho1[pos] = Rho[pos];
  }

  for (auto &pos : range) {
    Jx2[pos] = Jx[pos]; Jy2[pos] = Jy[pos]; Jz2[pos] = Jz[pos];
    Ex2[pos] = Ex[pos]; Ey2[pos] = Ey[pos]; Ez2[pos] = Ez[pos];
    Rho2[pos] = Rho[pos];
  }

  if (Jx(rn, rn) != 0) std::cout << "(After 1)  Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(After 1)  Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(After 1)  Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;

  // Re-assignment
  for (auto &pos : range) {
    Jx[pos] = Jx1[pos];
    Jy[pos] = Jy1[pos];
    Jz[pos] = Jz1[pos];
  }
  // Re-assignment
  
  std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time1 += duration;
  
  if (Jx(rn, rn) != 0) std::cout << "(Before 2) Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(Before 2) Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(Before 2) Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;
  
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
  
  if (Jx(rn, rn) != 0) std::cout << "(After 2)  Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(After 2)  Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(After 2)  Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;
  
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time2 += duration;
  
  // perform check
  bool equalFlag = true;

  for (auto &pos : range) {
    if ((Jx[pos] != Jx2[pos] && Jx[pos] != 0) || (Jy[pos] != Jy2[pos] && Jy[pos] != 0) || (Jy[pos] != Jy2[pos] && Jy[pos] != 0)){
      equalFlag = false;
      break;
    }
  }
  
  if (!equalFlag) std::cout << "For Case 2, the data equal check is : " << equalFlag << std::endl;
  
  // Re-assignment
  for (auto &pos : range) {
    Jx[pos] = Jx1[pos];
    Jy[pos] = Jy1[pos];
    Jz[pos] = Jz1[pos];
  }
  // Re-assignment
  
  if (Jx(rn, rn) != 0) std::cout << "(Before 3) Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(Before 3) Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(Before 3) Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;
  
  // Case 3
  start = std::chrono::high_resolution_clock::now();
  
  Jx.parallel_func(low, high, SCHNEK_DEVICE_LAMBDA(const auto& pos) {
    double jx = Jx[pos];
    double jy = Jy[pos];
    double jz = Jz[pos];
    double rho = Rho[pos];
    
    Jx.set(pos, (jx*gdtn + emdt*Ex[pos]*rho)/gdtd);
    Jy.set(pos, (jy*gdtn + emdt*Ey[pos]*rho)/gdtd);
    Jz.set(pos, (jz*gdtn + emdt*Ez[pos]*rho)/gdtd);
  });
  
  end = std::chrono::high_resolution_clock::now();
  
  if (Jx(rn, rn) != 0) std::cout << "(After 3)  Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(After 3)  Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(After 3)  Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;
  
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time3 += duration;
  
  // Re-assignment
  for (auto &pos : range) {
    Jx[pos] = Jx1[pos];
    Jy[pos] = Jy1[pos];
    Jz[pos] = Jz1[pos];
  }
  // Re-assignment

  if (Jx(rn, rn) != 0) std::cout << "(Before 4) Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(Before 4) Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(Before 4) Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;

  // Case 4
  start = std::chrono::high_resolution_clock::now();

  schnek::kokkos_utils::parallel_kokkos_parallel_for<Field>(low, high, SCHNEK_DEVICE_LAMBDA(const auto& pos) {
    double jx = Jx[pos];
    double jy = Jy[pos];
    double jz = Jz[pos];  
    double rho = Rho[pos];
    
    // Kokkos::Timer timer;
    
    Jx.set(pos, (jx*gdtn + emdt*Ex[pos]*rho)/gdtd);
    Jy.set(pos, (jy*gdtn + emdt*Ey[pos]*rho)/gdtd);
    Jz.set(pos, (jz*gdtn + emdt*Ez[pos]*rho)/gdtd);

    // double elapsed = timer.seconds(); // Convert to ms

    // Kokkos::atomic_add(&total_time41, elapsed);
  });

  end = std::chrono::high_resolution_clock::now();

  if (Jx(rn, rn) != 0) std::cout << "(After 4)  Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(After 4)  Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(After 4)  Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;
  
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time4 += duration;

  // Re-assignment
  for (auto &pos : range) {
    Jx[pos] = Jx1[pos];
    Jy[pos] = Jy1[pos];
    Jz[pos] = Jz1[pos];
  }
  // Re-assignment

  if (Jx(rn, rn) != 0) std::cout << "(Before 5) Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(Before 5) Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(Before 5) Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;

  // Case 5
  start = std::chrono::high_resolution_clock::now();

  schnek::kokkos_utils::parallel_kokkos_iteration<Field>(low, high, SCHNEK_DEVICE_LAMBDA(const auto& pos) {
    double jx = Jx[pos];
    double jy = Jy[pos];
    double jz = Jz[pos];  
    double rho = Rho[pos];

    // Kokkos::Timer timer;
    
    Jx.set(pos, (jx*gdtn + emdt*Ex[pos]*rho)/gdtd);
    Jy.set(pos, (jy*gdtn + emdt*Ey[pos]*rho)/gdtd);
    Jz.set(pos, (jz*gdtn + emdt*Ez[pos]*rho)/gdtd);
    
    // double elapsed = timer.seconds(); // Convert to ms

    // Kokkos::atomic_add(&total_time51, elapsed);
  });

  end = std::chrono::high_resolution_clock::now();

  if (Jx(rn, rn) != 0) std::cout << "(After 5)  Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(After 5)  Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(After 5)  Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;

  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time5 += duration;

  // Re-assignment
  for (auto &pos : range) {
    Jx[pos] = Jx1[pos];
    Jy[pos] = Jy1[pos];
    Jz[pos] = Jz1[pos];
  }
  // Re-assignment

  if (Jx(rn, rn) != 0) std::cout << "(Before 6) Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(Before 6) Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(Before 6) Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;

  // Case 6
  start = std::chrono::high_resolution_clock::now();

  schnek::kokkos_utils::parallel_kokkos_parallel_for_v1<Field>(low, high, SCHNEK_DEVICE_LAMBDA(const auto& pos) {
    double jx = Jx[pos];
    double jy = Jy[pos];
    double jz = Jz[pos];  
    double rho = Rho[pos];
    
    Jx.set(pos, (jx*gdtn + emdt*Ex[pos]*rho)/gdtd);
    Jy.set(pos, (jy*gdtn + emdt*Ey[pos]*rho)/gdtd);
    Jz.set(pos, (jz*gdtn + emdt*Ez[pos]*rho)/gdtd);
  });

  end = std::chrono::high_resolution_clock::now();

  if (Jx(rn, rn) != 0) std::cout << "(After 6)  Jx("<<rn<<", "<<rn<<") : "<<Jx(rn, rn)<<" - "<<Jx1(rn, rn)<<std::endl;
  if (Jy(rn, rn) != 0) std::cout << "(After 6)  Jy("<<rn<<", "<<rn<<") : "<<Jy(rn, rn)<<" - "<<Jy1(rn, rn)<<std::endl;
  if (Jz(rn, rn) != 0) std::cout << "(After 6)  Jz("<<rn<<", "<<rn<<") : "<<Jz(rn, rn)<<" - "<<Jz1(rn, rn)<<std::endl;
  
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time6 += duration;
  
  if (iteration_count % 10 == 0){
    std::cout << "Total execution time (1): " << total_time1.count() << " ms" << std::endl;
    std::cout << "Total execution time (2): " << total_time2.count() << " ms" << std::endl;
    std::cout << "Total execution time (3): " << total_time3.count() << " ms" << std::endl;
    std::cout << "Total execution time (4): " << total_time4.count() << " ms " << total_time41 << std::endl;
    std::cout << "Total execution time (5): " << total_time5.count() << " ms " << total_time51 << std::endl;
    std::cout << "Total execution time (6): " << total_time6.count() << " ms" << std::endl;
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