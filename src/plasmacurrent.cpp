#include "plasmacurrent.hpp"

#include "../huerto/electromagnetics/fieldsolver.hpp"

#include <schnek/grid/gridstorage/kokkos-storage.hpp>
#include <schnek/grid/array.hpp>
#include <schnek/grid/grid.hpp>
#include <schnek/grid/field.hpp>
// #include <schnek/grid/grid_utils.hpp>

// #include <Kokkos_Core.hpp>

#include <memory>
#include <chrono>

std::chrono::duration<double> total_time1(0);
std::chrono::duration<double> total_time2(0);
std::chrono::duration<double> total_time3(0);
std::chrono::duration<double> total_time4(0);
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

  // Testing
  Grid2d grid1(low2D, high2D);
  Grid2d grid2(low2D, high2D);

  // g(low2D[0], high2D[1]) = 1.0;
  fill_kokkos_grid(grid1, 0.5);

  grid2 = grid1;

  // auto dims = Jx.getDims();
  dims = Jx.getDims();

  // std::cout << "Jx Dims : " << dims[0] << ", " << dims[1] << std::endl;

  double grid1Sum = grid1.reduce(std::plus<double>(), 0.0);
  double grid2Sum = grid2.reduce(std::plus<double>(), 0.0);
  
  // std::cout << "Sum of grid1 elements: " << grid1Sum << std::endl;
  // std::cout << "Sum of grid2 elements: " << grid2Sum << std::endl;
  
  auto grid1Size = grid1.getSize();
  // std::cout << "Size of grid1 : " << grid1Size << std::endl;
  
  typedef schnek::Field<double, 2, HuertoGridChecker, schnek::KokkosDefaultGridStorage> Field2N;
  typedef schnek::Grid<double, 2, HuertoGridChecker, schnek::KokkosDefaultGridStorage> Grid2N;
  // int n = 100;
  
  Grid2N Jx_1, Jy_1, Jz_1;
  Field2N Ex_1, Ey_1, Ez_1;
  
  Jx_1.resize(Jx.getLo(), Jx.getHi());
  Jy_1.resize(Jy.getLo(), Jy.getHi());
  Jz_1.resize(Jz.getLo(), Jz.getHi());
  
  for (auto &pos : range) {
    Jx_1[pos] = Jx[pos];
    Jy_1[pos] = Jy[pos];
    Jz_1[pos] = Jz[pos];

    // if (Jx[pos] != 0) std::cout << "Jx("<<pos[0]<<" ,"<<pos[1]<<") : "<<Jx_1[pos]<<" - "<<Jx[pos]<<std::endl;
    // if (Jy[pos] != 0) std::cout << "Jy("<<pos[0]<<" ,"<<pos[1]<<") : "<<Jy_1[pos]<<" - "<<Jy[pos]<<std::endl;
    // if (Jz[pos] != 0) std::cout << "Jz("<<pos[0]<<" ,"<<pos[1]<<") : "<<Jz_1[pos]<<" - "<<Jz[pos]<<std::endl;
  }

  Ex_1.resize(Ex.getLo(), Ex.getHi(), Ex.getDomain(), Ex.getStagger(), Ex.getghostCells());
  Ey_1.resize(Ey.getLo(), Ey.getHi(), Ey.getDomain(), Ey.getStagger(), Ey.getghostCells());
  Ez_1.resize(Ez.getLo(), Ez.getHi(), Ez.getDomain(), Ez.getStagger(), Ez.getghostCells());
  
  for (auto &pos : range) {
    Ex_1[pos] = Ex[pos];
    Ey_1[pos] = Ey[pos];
    Ez_1[pos] = Ez[pos];

    // if (Ex[pos] != 0) std::cout << "Ex("<<pos[0]<<" ,"<<pos[1]<<") : "<<Ex_1[pos]<<" - "<<Ex[pos]<<std::endl;
    // if (Ey[pos] != 0) std::cout << "Ey("<<pos[0]<<" ,"<<pos[1]<<") : "<<Ey_1[pos]<<" - "<<Ey[pos]<<std::endl;
    // if (Ez[pos] != 0) std::cout << "Ez("<<pos[0]<<" ,"<<pos[1]<<") : "<<Ez_1[pos]<<" - "<<Ez[pos]<<std::endl;
  }

  Field2N f1;

  auto f1dims = f1.getDims();
  auto f1Size = f1.getSize();
  double f1Sum = f1.reduce(std::plus<double>(), 0.0);

  // std::cout << "f1 Dims : " << f1dims[0] << ", " << f1dims[1] << std::endl;
  // std::cout << "f1 Size : " << f1Size << std::endl;
  // std::cout << "f1 Sum  : " << f1Sum << std::endl;
  
  auto JxSize = Jx.getSize();
  // std::cout << "Jx Size : " << JxSize << std::endl;

  Field2N Rho_1;

  Rho_1.resize(
    Rho.getLo(), 
    Rho.getHi(),     
    Rho.getDomain(), 
    Rho.getStagger(),
    Rho.getghostCells()
    // 0
  );

  // std::cout << "Rho ghostCells : "<<Rho.getghostCells()<<std::endl;

  for (auto &pos : range) {
    Rho_1[pos] = Rho[pos];
    // std::cout << "Rho("<<pos[0]<<" ,"<<pos[1]<<") : "<<Rho_1[pos]<<" - "<<Rho[pos]<<std::endl;
  }

  // throws error
  // double JxSum = Jx.reduce(std::plus<double>(), 0.0);
  // std::cout << "Sum of Jx elements: " << JxSum << std::endl;

  double Jz1Sum = Jz_1.reduce(std::plus<double>(), 0.0);
  std::cout << "Sum of Jz_1 elements: " << Jz1Sum << std::endl;

  // Case 3
  start = std::chrono::high_resolution_clock::now();

  // WORKS using redefined variables
  Kokkos::parallel_for("plasmacurrent_2D",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
      {low[0], low[1]}, 
      {high[0], high[1]}
    ),
    KOKKOS_LAMBDA (const int i, const int j) {
      // Create index for the current position
      typename Field2N::IndexType pos;
      pos[0] = i;
      pos[1] = j;
      
      // Access field values
      double jx = Jx_1[pos];
      double jy = Jy_1[pos];
      double jz = Jz_1[pos];
      double rho = Rho_1[pos];
      
      // THESE ASSIGNMENT OPERATIONS DO NOT WORK - error
      // error occurs because when using [] inside a KOKKOS_LAMBDA, you don't get a reference that can be modified
      // it is designed to execute on devices (not the host) - so references to the host are not allowed
      // Jx_1[pos] = (jx*gdtn + emdt*Ex_1[pos]*rho)/gdtd;
      // Jy_1[pos] = (jy*gdtn + emdt*Ey_1[pos]*rho)/gdtd;
      // Jz_1[pos] = (jz*gdtn + emdt*Ez_1[pos]*rho)/gdtd;

      Jx_1.set(pos, (jx*gdtn + emdt*Ex_1[pos]*rho)/gdtd);
      Jy_1.set(pos, (jy*gdtn + emdt*Ey_1[pos]*rho)/gdtd);
      Jz_1.set(pos, (jz*gdtn + emdt*Ez_1[pos]*rho)/gdtd);
    }
  );

  Kokkos::fence();

  end = std::chrono::high_resolution_clock::now();
  
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time3 += duration;

  // Case 4
  start = std::chrono::high_resolution_clock::now();

  Jx_1.parallel_func(low, high, [&](const auto& pos) {
    // Access field values
    double jx = Jx_1[pos];
    double jy = Jy_1[pos];
    double jz = Jz_1[pos];
    double rho = Rho_1[pos];
    
    // Calculate new values using set method
    Jx_1.set(pos, (jx*gdtn + emdt*Ex_1[pos]*rho)/gdtd);
    Jy_1.set(pos, (jy*gdtn + emdt*Ey_1[pos]*rho)/gdtd);
    Jz_1.set(pos, (jz*gdtn + emdt*Ez_1[pos]*rho)/gdtd);
  });

  // Jx.parallel_func(low, high, [&](const auto& pos) {
  //   // Access field values
  //   double jx = Jx[pos];
  //   double jy = Jy[pos];
  //   double jz = Jz[pos];
  //   double rho = Rho[pos];
    
  //   // Calculate new values using set method
  //   Jx.set(pos, (jx*gdtn + emdt*Ex[pos]*rho)/gdtd);
  //   Jy.set(pos, (jy*gdtn + emdt*Ey[pos]*rho)/gdtd);
  //   Jz.set(pos, (jz*gdtn + emdt*Ez[pos]*rho)/gdtd);
  // });

  Kokkos::fence();

  end = std::chrono::high_resolution_clock::now();
  
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  total_time4 += duration;

  if (iteration_count % 10 == 0){
    std::cout << "Total execution time (1): " << total_time1.count() << " ms" << std::endl;
    std::cout << "Total execution time (2): " << total_time2.count() << " ms" << std::endl;
    std::cout << "Total execution time (3): " << total_time3.count() << " ms" << std::endl;
    std::cout << "Total execution time (4): " << total_time4.count() << " ms" << std::endl;
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

