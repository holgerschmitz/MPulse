#include "mpi_bound.h"

#include "mpulse.h"
#include "globals.h"
#include "factor.h"

#include <sstream>
#include <fstream>

#ifndef SINGLE_PROCESSOR

/* **************************************************************
 *                 MPIPeriodicSplitXBoundary                    *
 ****************************************************************/

MPIPeriodicSplitXBoundary::MPIPeriodicSplitXBoundary()
{} 

void MPIPeriodicSplitXBoundary::init()
{
  Low = Globals::instance().gridLow();
  High = Globals::instance().gridHigh();
  
  MPI_Comm_size(MPI_COMM_WORLD,&ComSize);

  int periodic = true;

  MPI_Cart_create(MPI_COMM_WORLD,1,&ComSize,&periodic,true,&comm); 
  MPI_Comm_rank(comm,&ComRank);

  MPI_Cart_coords(comm,ComRank,1,&mycoord); 

  MPI_Cart_shift(comm,0,1,&leftcoord,&rightcoord); 

  double width = (High[0]-2.)/double(ComSize);

  if (ComRank>0) 
      Low[0] = int(width*mycoord);

  if (ComRank<(ComSize-1))
      High[0] = int(width*(mycoord+1))+1;

  exchSize = (High[1]-Low[1]+1)*(High[2]-Low[2]+1);

  sendarr = new double[exchSize];
  recvarr = new double[exchSize];
  
  Globals::instance().setMaster(this->master());
  Globals::instance().setUniqueId(this->getUniqueId());

  std::ostringstream S;
  S << "boundary"<<ComRank<<".dat"<<char(0);
  std::ofstream O(S.str().c_str());
  std::cout << "Coord: " << mycoord << "\n";
  std::cout  << "Low: " << Low[0] << " " << Low[1] << "\n";
  std::cout  << "High: " << High[0] << " " << High[1] << "\n";
  O.close();
}

MPIPeriodicSplitXBoundary::~MPIPeriodicSplitXBoundary() {
    MPI_Finalize();
    delete[] sendarr;
    delete[] recvarr;
}

void MPIPeriodicSplitXBoundary::exchangeX(DataGrid3d &field)
{
  int xi;
  int yi;
  int zi;
  MPI_Status stat; 
    
  int arr_ind = 0;
  xi = Low[0] + 1;
  
  for (yi = Low[1]; yi <= High[1]; ++yi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      sendarr[arr_ind++] = field(xi, yi, zi);

  MPI_Sendrecv(sendarr, exchSize, MPI_DOUBLE, leftcoord, 0, 
               recvarr, exchSize, MPI_DOUBLE, rightcoord, 0, 
               comm, &stat); 

  arr_ind = 0;
  xi = High[0];
  
  for (yi = Low[1]; yi <= High[1]; ++yi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      field(xi,yi,zi) = recvarr[arr_ind++]; 

  arr_ind = 0;
  xi = High[0] - 1;
  
  for (yi = Low[1]; yi <= High[1]; ++yi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      sendarr[arr_ind++] = field(xi, yi, zi);

  MPI_Sendrecv(sendarr, exchSize, MPI_DOUBLE, rightcoord, 0, 
               recvarr, exchSize, MPI_DOUBLE, leftcoord, 0, 
               comm, &stat); 

  arr_ind = 0;
  xi = Low[0];
  
  for (yi = Low[1]; yi <= High[1]; ++yi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      field(xi,yi,zi) = recvarr[arr_ind++]; 

}

void MPIPeriodicSplitXBoundary::exchangeY(DataGrid3d &field) 
{
  const GridIndex3d &UBound = field.getHi();
  const GridIndex3d &LBound = field.getLo();
  
  int xi;
  int zi;

  int my0=UBound[1], my1=my0-1;
  int ly0=LBound[1], ly1=ly0+1;
  
  for (xi = LBound[0]; xi <= UBound[0]; ++xi)
  {
    for (zi = LBound[2]; zi <= UBound[2]; ++zi)
    {
      field(xi, ly0, zi) = field(xi, my1, zi);
      field(xi, my0, zi) = field(xi, ly1, zi);
    }
  }
}

void MPIPeriodicSplitXBoundary::exchangeZ(DataGrid3d &field) 
{
  const GridIndex3d &UBound = field.getHi();
  const GridIndex3d &LBound = field.getLo();
  
  int xi;
  int yi;

  int mz0=UBound[2], mz1=mz0-1;
  int lz0=LBound[2], lz1=lz0+1;
  
  for (xi = LBound[0]; xi <= UBound[0]; ++xi)
  {
    for (yi = LBound[1]; yi <= UBound[1]; ++yi)
    {
      field(xi, yi, lz0) = field(xi, yi, mz1);
      field(xi, yi, mz0) = field(xi, yi, lz1);
    }
  }
}


double MPIPeriodicSplitXBoundary::AvgReduce(double val) const {
  double result;
    //this collects results from all nodes, MPI_SUM returns the sum of all results
  MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
  return result/double(ComSize);
}

double MPIPeriodicSplitXBoundary::MaxReduce(double val) const {
  double result;
    //this collects results from all nodes, MPI_MAX returns the maximum of all results
  MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_MAX, comm);
  return result;
}

double MPIPeriodicSplitXBoundary::SumReduce(double val) const {
  double result;
    //this collects results from all nodes, MPI_MAX returns the maximum of all results
  MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
  return result;
}


const GridIndex &MPIPeriodicSplitXBoundary::RegionLow() const {
    return Low;
}

const GridIndex &MPIPeriodicSplitXBoundary::RegionHigh() const {
    return High;
}

/* **************************************************************
 *                 MPIPeriodicSplitXYBoundary                    *
 ****************************************************************/

MPIPeriodicSplitXYZBoundary::MPIPeriodicSplitXYZBoundary()
{} 


void MPIPeriodicSplitXYZBoundary::init()
{ 
  Low = Globals::instance().gridLow();
  High = Globals::instance().gridHigh();
  
  MPI_Comm_size(MPI_COMM_WORLD,&ComSize);

  int periodic[2] = {true, true};

  std::vector<int> box(3);
  for (int i=0; i<3; ++i)
  {
    box[i] = High[i]-Low[i]-1;
  }

  std::vector<int> eqDims;

  equalFactors(ComSize, 3, eqDims, box);

  std::copy(eqDims.begin(), eqDims.end(), dims);

  MPI_Cart_create(MPI_COMM_WORLD,3,dims,periodic,true,&comm); 
  MPI_Comm_rank(comm,&ComRank);

  MPI_Cart_coords(comm,ComRank,3,mycoord); 

  MPI_Cart_shift(comm,0,1,&xprevcoord,&xnextcoord); 
  MPI_Cart_shift(comm,1,1,&yprevcoord,&ynextcoord); 
  MPI_Cart_shift(comm,2,1,&zprevcoord,&znextcoord);

  double width[3];
  width[0] = (High[0]-1.)/double(dims[0]);
  width[1] = (High[1]-1.)/double(dims[1]);
  width[2] = (High[2]-1.)/double(dims[2]);

  for (int i=0; i<3; ++i)
  {
    if (mycoord[i]>0) 
      Low[i] = int(width[i]*mycoord[i]);

    if (mycoord[i]<(dims[i]-1))
      High[i] = int(width[i]*(mycoord[i]+1))+1;

  }


  exchSize[0] = (High[1]-Low[1]+1)*(High[2]-Low[2]+1);
  exchSize[1] = (High[0]-Low[0]+1)*(High[2]-Low[2]+1);
  exchSize[2] = (High[0]-Low[0]+1)*(High[1]-Low[1]+1);

  sendarrx = new double[exchSize[0]];
  recvarrx = new double[exchSize[0]];
  sendarry = new double[exchSize[1]];
  recvarry = new double[exchSize[1]];
  sendarrz = new double[exchSize[2]];
  recvarrz = new double[exchSize[2]];
 
  Globals::instance().setMaster(this->master());
  Globals::instance().setUniqueId(this->getUniqueId());

  std::ostringstream S;
  S << "boundary"<<ComRank<<".dat"<<char(0);
  std::ofstream O(S.str().c_str());
  O << "Coord: " << mycoord[0] << " " << mycoord[1] << " " << mycoord[2] << "\n";
  O << "Dims: " << dims[0] << " " << dims[1] << " " << dims[2] << "\n";
  O << "Low: " << Low[0] << " " << Low[1] << " " << Low[2] << "\n";
  O << "High: " << High[0] << " " << High[1] << " " << High[2] << "\n";
  O.close();
}

MPIPeriodicSplitXYZBoundary::~MPIPeriodicSplitXYZBoundary() 
{
    MPI_Finalize();
    delete[] sendarrx;
    delete[] recvarrx;
    delete[] sendarry;
    delete[] recvarry;
    delete[] sendarrz;
    delete[] recvarrz;
}

void MPIPeriodicSplitXYZBoundary::exchangeX(DataGrid3d &field)
{
  int xi;
  int yi;
  int zi;
  MPI_Status stat; 
    
  int arr_ind = 0;
  xi = Low[0] + 1;
  
  for (yi = Low[1]; yi <= High[1]; ++yi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      sendarrx[arr_ind++] = field(xi, yi, zi);

  if (arr_ind!=exchSize[0]) std::cerr << "Error X 1\n";

  MPI_Sendrecv(sendarrx, exchSize[0], MPI_DOUBLE, xprevcoord, 0, 
               recvarrx, exchSize[0], MPI_DOUBLE, xnextcoord, 0, 
               comm, &stat); 
  arr_ind = 0;
  xi = High[0];
  
  for (yi = Low[1]; yi <= High[1]; ++yi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      field(xi,yi,zi) = recvarrx[arr_ind++]; 

  if (arr_ind!=exchSize[0]) std::cerr << "Error X 2\n";

  arr_ind = 0;
  xi = High[0] - 1;
  
  for (yi = Low[1]; yi <= High[1]; ++yi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      sendarrx[arr_ind++] = field(xi, yi, zi);

  if (arr_ind!=exchSize[0]) std::cerr << "Error X 3\n";

  MPI_Sendrecv(sendarrx, exchSize[0], MPI_DOUBLE, xnextcoord, 0, 
               recvarrx, exchSize[0], MPI_DOUBLE, xprevcoord, 0, 
               comm, &stat); 

  arr_ind = 0;
  xi = Low[0];
  
  for (yi = Low[1]; yi <= High[1]; ++yi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      field(xi,yi,zi) = recvarrx[arr_ind++]; 

  if (arr_ind!=exchSize[0]) std::cerr << "Error X 4\n";

}

void MPIPeriodicSplitXYZBoundary::exchangeY(DataGrid3d &field)
{
  int xi;
  int yi;
  int zi;
  MPI_Status stat; 
    
  int arr_ind = 0;
  yi = Low[1] + 1;
  
  for (xi = Low[0]; xi <= High[0]; ++xi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      sendarry[arr_ind++] = field(xi, yi, zi);

  if (arr_ind!=exchSize[1]) std::cerr << "Error Y 1\n";

  MPI_Sendrecv(sendarry, exchSize[1], MPI_DOUBLE, yprevcoord, 0, 
               recvarry, exchSize[1], MPI_DOUBLE, ynextcoord, 0, 
               comm, &stat); 

  arr_ind = 0;
  yi = High[1];
  
  for (xi = Low[0]; xi <= High[0]; ++xi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      field(xi,yi,zi) = recvarry[arr_ind++]; 

  if (arr_ind!=exchSize[1]) std::cerr << "Error Y 2\n";

  arr_ind = 0;
  yi = High[1] - 1;
  
  for (xi = Low[0]; xi <= High[0]; ++xi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      sendarry[arr_ind++] = field(xi, yi, zi);

  if (arr_ind!=exchSize[1]) std::cerr << "Error Y 3\n";

  MPI_Sendrecv(sendarry, exchSize[1], MPI_DOUBLE, ynextcoord, 0, 
               recvarry, exchSize[1], MPI_DOUBLE, yprevcoord, 0, 
               comm, &stat); 

  arr_ind = 0;
  yi = Low[1];
  
  for (xi = Low[0]; xi <= High[0]; ++xi)
    for (zi = Low[2]; zi <= High[2]; ++zi)
      field(xi,yi,zi) = recvarry[arr_ind++]; 

  if (arr_ind!=exchSize[1]) std::cerr << "Error Y 4\n";

}

void MPIPeriodicSplitXYZBoundary::exchangeZ(DataGrid3d &field)
{
  int xi;
  int yi;
  int zi;
  MPI_Status stat; 
    
  int arr_ind = 0;
  zi = Low[2] + 1;
  
  for (xi = Low[0]; xi <= High[0]; ++xi)
    for (yi = Low[1]; yi <= High[1]; ++yi)
      sendarrz[arr_ind++] = field(xi, yi, zi);

  if (arr_ind!=exchSize[2]) std::cerr << "Error Z 1\n";

  MPI_Sendrecv(sendarrz, exchSize[2], MPI_DOUBLE, zprevcoord, 0, 
               recvarrz, exchSize[2], MPI_DOUBLE, znextcoord, 0, 
               comm, &stat); 

  arr_ind = 0;
  zi = High[2];
  
  for (xi = Low[0]; xi <= High[0]; ++xi)
    for (yi = Low[1]; yi <= High[1]; ++yi)
      field(xi,yi,zi) = recvarrz[arr_ind++]; 

  if (arr_ind!=exchSize[2]) std::cerr << "Error Z 2\n";

  arr_ind = 0;
  zi = High[2] - 1;
  
  for (xi = Low[0]; xi <= High[0]; ++xi)
    for (yi = Low[1]; yi <= High[1]; ++yi)
      sendarrz[arr_ind++] = field(xi, yi, zi);

  if (arr_ind!=exchSize[2]) std::cerr << "Error Z 3\n";

  MPI_Sendrecv(sendarrz, exchSize[2], MPI_DOUBLE, znextcoord, 0, 
               recvarrz, exchSize[2], MPI_DOUBLE, zprevcoord, 0, 
               comm, &stat); 

  arr_ind = 0;
  zi = Low[2];
  
  for (xi = Low[0]; xi <= High[0]; ++xi)
    for (yi = Low[1]; yi <= High[1]; ++yi)
      field(xi,yi,zi) = recvarrz[arr_ind++]; 

  if (arr_ind!=exchSize[2]) std::cerr << "Error Z 4\n";
}


double MPIPeriodicSplitXYZBoundary::AvgReduce(double val) const {
  double result;
  MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
  return result/double(ComSize);
}

double MPIPeriodicSplitXYZBoundary::MaxReduce(double val) const {
  double result;
  MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_MAX, comm);
  return result;
}

double MPIPeriodicSplitXYZBoundary::SumReduce(double val) const {
  double result;
  MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
  return result;
}

const GridIndex &MPIPeriodicSplitXYZBoundary::RegionLow() const {
  return Low;
}

const GridIndex &MPIPeriodicSplitXYZBoundary::RegionHigh() const {
  return High;
}


#endif
