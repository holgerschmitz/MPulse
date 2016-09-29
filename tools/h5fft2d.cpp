#include <schnek/grid.h>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include "../src/hdfstream.h"
#include "../src/mpulse.h"

#include <sstream>
#include <fstream>
#include <complex>
#include <fftw3.h>

namespace po = boost::program_options;

typedef std::complex<double> Complex;
typedef schnek::Grid<Complex, 2, MPulseGridChecker> CDataGrid2d;


void performFFT(const DataGrid2d &inGrid, CDataGrid2d &outGrid)
{
  fftw_complex *in, *out;
  fftw_plan p;

  int dimx = inGrid.getDims(0);
  int dimy = inGrid.getDims(1);
  int lowx = inGrid.getLow(0);
  int lowy = inGrid.getLow(1);
  int highx = inGrid.getHigh(0);
  int highy = inGrid.getHigh(1);

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimx*dimy);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimx*dimy);

  for (int x=lowx, nx=0; x<=highx; ++x, ++nx)
  {
    for (int y=lowy, ny=0; y<=highy; ++y, ++ny)
    {
      in[nx*dimy+ny][0] = inGrid(x,y);
      in[nx*dimy+ny][1] = 0;
//      Complex z(inGrid(x,y),0.0);
//      in[nx*dimy+ny] = *(reinterpret_cast<fftw_complex*>(&z));
    }
  }

  p = fftw_plan_dft_2d(dimx, dimy, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


  fftw_execute(p); 
  
  for (int x=lowx, nx=0; x<=highx; ++x, ++nx)
  {
    for (int y=lowy, ny=0; y<=highy; ++y, ++ny)
    {
      int px = (dimx/2 + nx) % dimx;
      int py = (dimy/2 + ny) % dimy;
      outGrid(px,py) = Complex(out[nx*dimy+ny][0], out[nx*dimy+ny][1]);
    }
  }
  
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);

}

void filterCos2(DataGrid2d &grid)
{
  std::cerr << "Cos^2 filter\n";
  int lowx = grid.getLow(0);
  int lowy = grid.getLow(1);
  int highx = grid.getHigh(0);
  int highy = grid.getHigh(1);
  
  double mx = 0.5*(lowx+highx);
  double my = 0.5*(lowy+highy);
  
  double wx = M_PI/(highx-lowx+1);
  double wy = M_PI/(highy-lowy+1);
  
  for (int nx=lowx; nx<=highx; ++nx)
  {
    for (int ny=lowy; ny<=highy; ++ny)
    {
      double cx = cos(wx*(nx-mx));
      double cy = cos(wy*(ny-my));
      grid(nx,ny) = grid(nx,ny)*cx*cx*cy*cy;
    }
  }
}

void calcAbs(const CDataGrid2d &inGrid, DataGrid2d &outGrid)
{
  int lowx = inGrid.getLow(0);
  int lowy = inGrid.getLow(1);
  int highx = inGrid.getHigh(0);
  int highy = inGrid.getHigh(1);
  
  for (int nx=lowx; nx<=highx; ++nx)
  {
    for (int ny=lowy; ny<=highy; ++ny)
    {
      outGrid(nx,ny) = abs(inGrid(nx,ny));
    }
  }
}

void calcNorm(const CDataGrid2d &inGrid, DataGrid2d &outGrid)
{
  int lowx = inGrid.getLow(0);
  int lowy = inGrid.getLow(1);
  int highx = inGrid.getHigh(0);
  int highy = inGrid.getHigh(1);
  
  for (int nx=lowx; nx<=highx; ++nx)
  {
    for (int ny=lowy; ny<=highy; ++ny)
    {
      outGrid(nx,ny) = norm(inGrid(nx,ny));
    }
  }
}

void calcPhase(const CDataGrid2d &inGrid, DataGrid2d &outGrid)
{
  int lowx = inGrid.getLow(0);
  int lowy = inGrid.getLow(1);
  int highx = inGrid.getHigh(0);
  int highy = inGrid.getHigh(1);
  
  for (int nx=lowx; nx<=highx; ++nx)
  {
    for (int ny=lowy; ny<=highy; ++ny)
    {
      outGrid(nx,ny) = arg(inGrid(nx,ny));
    }
  }
}

void calcImag(const CDataGrid2d &inGrid, DataGrid2d &outGrid)
{
  int lowx = inGrid.getLow(0);
  int lowy = inGrid.getLow(1);
  int highx = inGrid.getHigh(0);
  int highy = inGrid.getHigh(1);
  
  for (int nx=lowx; nx<=highx; ++nx)
  {
    for (int ny=lowy; ny<=highy; ++ny)
    {
      outGrid(nx,ny) = imag(inGrid(nx,ny));
    }
  }
}

void calcReal(const CDataGrid2d &inGrid, DataGrid2d &outGrid)
{
  int lowx = inGrid.getLow(0);
  int lowy = inGrid.getLow(1);
  int highx = inGrid.getHigh(0);
  int highy = inGrid.getHigh(1);
  
  for (int nx=lowx; nx<=highx; ++nx)
  {
    for (int ny=lowy; ny<=highy; ++ny)
    {
      outGrid(nx,ny) = real(inGrid(nx,ny));
    }
  }
}

void calcLog(const CDataGrid2d &inGrid, DataGrid2d &outGrid)
{
  int lowx = inGrid.getLow(0);
  int lowy = inGrid.getLow(1);
  int highx = inGrid.getHigh(0);
  int highy = inGrid.getHigh(1);
  
  for (int nx=lowx; nx<=highx; ++nx)
  {
    for (int ny=lowy; ny<=highy; ++ny)
    {
      
      outGrid(nx,ny) = log(abs(inGrid(nx,ny))+1e-20);
    }
  }
}


int main(int argc, char**argv) {
  std::string infile;
  std::string outfile;
  std::string operation;
  std::string filter;
  std::string blockname;
  int dimx;
  int dimy;
  
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("input,i", po::value<std::string>(&infile),
          "set input file name base (use %p for the processor number)")
      ("output,o", po::value<std::string>(&outfile),
          "set output file name")
      ("xdim,x", po::value<int>(&dimx),
          "x dimension of the field")
      ("ydim,y", po::value<int>(&dimy),
          "y dimension of the field")
      ("block,b", po::value<std::string>(&blockname),
          "name of the data block in the hdf file")
      ("operation,p", po::value<std::string>(&operation)->default_value("abs"),
          "perform operation on complex fourier field (abs, phase, norm, real, imag, log)")
      ("filter,f", po::value<std::string>(&filter)->default_value("none"),
          "filter input (none, cos2)");


  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if ((vm.count("help")>0) ||
      (vm.count("input")==0) ||
      (vm.count("output")==0)) {
      std::cerr << desc << "\n";
      return 1;
  }
  
  DataGrid2dContainer inContainer;
  DataGrid2d inGrid(GridIndex2d(dimx,dimy));
  inContainer.grid = &inGrid;
  inContainer.active = true;

  std::cerr << "opening " << infile.c_str() << "\n";
  HDFistream input(infile.c_str());
  if (vm.count("block")>0)
  {
    input.setBlockName(blockname);
  }
  std::cerr << "reading " << infile.c_str() << "\n";
  input >> inContainer; 
  input.close();

  if (filter=="cos2")
    filterCos2(inGrid);


  CDataGrid2d fftGrid(inGrid.getDims());
  
  performFFT(inGrid, fftGrid);
  
  DataGrid2dContainer outContainer;
  DataGrid2d outGrid(inGrid.getDims());
  outContainer.grid = &outGrid;
  outContainer.active = true;
  outContainer.global_min = GridIndex2d(0,0);
  outContainer.global_max = GridIndex2d(dimx,dimy);
  
  if (operation=="abs")
    calcAbs(fftGrid, outGrid);
  else if (operation=="phase")
    calcPhase(fftGrid, outGrid);
  else if (operation=="norm")
    calcNorm(fftGrid, outGrid);
  else if (operation=="real")
    calcReal(fftGrid, outGrid);
  else if (operation=="imag")
    calcImag(fftGrid, outGrid); 
  else if (operation=="log")
    calcLog(fftGrid, outGrid);
  else
    calcAbs(fftGrid, outGrid);
  

  HDFostream output(outfile.c_str());
  output << outContainer;
  output.close();


}
