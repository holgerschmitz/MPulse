#include "hdfstream.h"

#ifdef USE_HDF_PARALLEL
#include <mpi.h>
#endif

HDFstream::HDFstream()
  : file_id(-1),
    status(0),
    blockname("data"),
    sets_count(0),
    active(true),
    activeModified(false)
#ifdef USE_HDF_PARALLEL
    , commSet(false)
#endif
{}

HDFstream::HDFstream(const HDFstream& hdf)
  : file_id(hdf.file_id),
    status(hdf.status),
    blockname(hdf.blockname),
    sets_count(hdf.sets_count),
    active(true),
    activeModified(false)
#ifdef USE_HDF_PARALLEL
    , commSet(false)
#endif
{}

HDFstream &HDFstream::operator=(const HDFstream& hdf)
{
  file_id = hdf.file_id;
  status = hdf.status;
  sets_count = hdf.sets_count;
  blockname = hdf.blockname;
  active = hdf.active;
  activeModified = hdf.activeModified;
#ifdef USE_HDF_PARALLEL
  mpiComm = hdf.mpiComm;
  commSet = hdf.commSet;
#endif
  return *this;
}

HDFstream::~HDFstream()
{
  close();
}

void HDFstream::close()
{
  if (file_id >= 0) {
    H5Fclose (file_id);
  }
  file_id = -1;
}

bool HDFstream::good() const
{
  return (file_id>=0);
}

void HDFstream::setBlockName(std::string blockname_)
{
  blockname = blockname_;
  sets_count = -1;
}

std::string HDFstream::getNextBlockName()
{
  std::ostringstream bname;
  bname << blockname;

  if (sets_count<0) sets_count = 1;
  else
  {
    bname << sets_count++;
  }
  return bname.str();
}

#ifdef USE_HDF_PARALLEL
void HDFstream::makeMPIGroup()
{
  if (!activeModified) {
    if (!commSet)
    {
      mpiComm = MPI_COMM_WORLD;
      commSet = true;
    }
    return;
  }
  
  int rank, size;
  MPI_Group worldGroup, activeGroup;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
  
  int *inputArr = new int[size];
  int *activeArr = new int[size];
  
  for (int i=0; i<size; ++i) inputArr[i] = 0;
  if (active) inputArr[rank] = 1;

  MPI_Allreduce(inputArr, activeArr, size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int count = 0;  
  for (int i=0; i<size; ++i)
  {
    if (activeArr[i]>0) {
      inputArr[count] = i;
      ++count;
    }
  }
  
  MPI_Group_incl(worldGroup, count, inputArr, &activeGroup);
  MPI_Comm_create(MPI_COMM_WORLD, activeGroup, &mpiComm);

  delete[] activeArr;
  delete[] inputArr;
  
  activeModified = false;
}
#endif


// ----------------------------------------------------------------------

HDFistream::HDFistream() 
  : HDFstream() {}

HDFistream::HDFistream(const char* fname)
  : HDFstream()
{
  open(fname);
}

HDFistream::HDFistream(const HDFistream& hdf)
  : HDFstream(hdf)
{}

int HDFistream::open(const char* fname)
{
  close();

#ifdef USE_HDF_PARALLEL
  makeMPIGroup();
  if (active)
  {
    /* setup file access template */
    hid_t plist_id = H5Pcreate (H5P_FILE_ACCESS);

    /* set Parallel access with communicator */
//    H5Pset_fapl_mpio(plist_id, mpiComm, MPI_INFO_NULL);   
    H5Pset_fapl_mpiposix(plist_id, mpiComm, 0);   
    /* open the file collectively */
    file_id = H5Fopen (fname, H5F_ACC_RDONLY, plist_id);
    /* Release file-access template */
    H5Pclose(plist_id);
  }
#else
  if (active)
    file_id = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif  
  sets_count = 0;
  return 1;
}


// ----------------------------------------------------------------------

HDFostream::HDFostream()
   : HDFstream()
{}

HDFostream::HDFostream(const HDFostream& hdf)
  : HDFstream(hdf)
{}

HDFostream::HDFostream(const char* fname)
   : HDFstream()
{
  open(fname);
}

int HDFostream::open(const char* fname)
{
  sets_count = 0;
  
#ifdef USE_HDF_PARALLEL
  makeMPIGroup();
  if (active)
  {
    /* setup file access template */
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    /* set Parallel access with communicator */
//    H5Pset_fapl_mpio(plist_id, mpiComm, MPI_INFO_NULL);   
    H5Pset_fapl_mpiposix(plist_id, mpiComm, 0);   
    file_id = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    /* Release file-access template */
    H5Pclose(plist_id);
  }
#else
  if (active)
    file_id = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif  

  return file_id;
}

// ----------------------------------------------------------------------

template<>
const hid_t H5DataType<int>::type = H5T_NATIVE_INT;

template<>
const hid_t H5DataType<float>::type = H5T_NATIVE_FLOAT;

template<>
const hid_t H5DataType<double>::type = H5T_NATIVE_DOUBLE;

// ----------------------------------------------------------------------

std::ostream &operator<<(std::ostream& out, const DataGridContainer &data)
{
  DataGrid &grid = *(data.grid);
  GridIndex low = grid.getLow();
  GridIndex high = grid.getHigh();
  
  for (int i=low[0]; i<=high[0]; ++i)
  {
    for (int j=low[1]; j<=high[1]; ++j)
    {
      for (int k=low[2]; k<=high[2]; ++k)
      {
        out << i << " " << j << " " << k << " " << grid(i,j,k) << "\n";
      }
      out << "\n";
    }
  }
  return out;
}

