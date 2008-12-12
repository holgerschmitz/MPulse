#include <sstream>

#ifndef SINGLE_PROCESSOR
#include <mpi.h>
#endif

template<typename TYPE, int RANK, template<int> class Checking>
HDFistream& HDFistream::operator>>(MatrixContainer<TYPE, RANK, Checking>& m)
{
#ifndef SINGLE_PROCESSOR
  if (!active) return *this;
#endif

  std::string dset_name = getNextBlockName();

  typedef typename schnek::Matrix<TYPE, RANK, Checking>::IndexType IndexType;
  
  IndexType mdims = m.grid->getDims();
  IndexType mlow = m.grid->getLow();
  IndexType mhigh = m.grid->getHigh();
  
  hsize_t dims[RANK];
  hsize_t locdims[RANK];
  hsize_t start[RANK];
    
  for (int i=0; i<RANK; ++i) 
  {
    int gmin = m.global_min[i];
    dims[i] = 1 + m.global_max[i] - gmin;
    locdims[i] = mdims[i];
    start[i] = mlow[i] - gmin;
  }
  
  TYPE *data = &(*m.grid->begin());

  /* open the dataset collectively */
  hid_t dataset = H5Dopen(file_id, dset_name.c_str());
  assert(dataset != -1);
  
#ifndef SINGLE_PROCESSOR
  /* create a file dataspace independently */
  hid_t file_dataspace = H5Dget_space(dataset);
  assert(file_dataspace != -1);

  hid_t ret=H5Sselect_hyperslab(file_dataspace, 
                                H5S_SELECT_SET, 
                                start, 
                                NULL,
	                              locdims, 
                                NULL);
  assert(ret != -1);

  /* create a memory dataspace independently */
  hid_t mem_dataspace = H5Screate_simple(RANK, locdims, NULL);
  assert (mem_dataspace != -1);


  /* read the data independently */
  ret = H5Dread(dataset,
                H5DataType<TYPE>::type, 
                mem_dataspace, 
                file_dataspace,
	              H5P_DEFAULT,
                data);
  assert(ret != -1);

  /* release all IDs created */
  H5Sclose(mem_dataspace);
  H5Sclose(file_dataspace);
#else
  /* read the data on single processor */
  hid_t ret = H5Dread(dataset,
                      H5DataType<TYPE>::type, 
                      H5S_ALL, 
                      H5S_ALL,
	                    H5P_DEFAULT,
                      data);
  assert(ret != -1);
#endif

  /* close dataset collectively */
  ret=H5Dclose(dataset);
  assert(ret != -1);

  
  return *this;
}

template<typename TYPE, int RANK, template<int> class Checking>
HDFostream& HDFostream::operator<< (const MatrixContainer<TYPE, RANK, Checking>& m)
{
#ifndef SINGLE_PROCESSOR
  if (!active) return *this;
#endif

  std::string dset_name = getNextBlockName();
  
  typedef typename schnek::Matrix<TYPE, RANK, Checking>::IndexType IndexType;
  
  IndexType mdims = m.grid->getDims();
  IndexType mlow = m.grid->getLow();
  IndexType mhigh = m.grid->getHigh();
  
  hsize_t dims[RANK];
  hsize_t locdims[RANK];
  hsize_t start[RANK];
  
  bool empty = false;
  
  for (int i=0; i<RANK; ++i) 
  {
    int gmin = m.global_min[i];
    dims[i] = 1 + m.global_max[i] - gmin;
    locdims[i] = mdims[i];
    start[i] = mlow[i] - gmin;
    
    if (locdims[i]<=0) empty = true;
    
    if (dims[i]<(start[i]+locdims[i]))
    {
      std::cerr << "FATAL ERROR!\n"
        << "Dimension " << i << ":\n  global size: " << dims[i]
        << "\n  local start: " << start[i]
        << "\n  local size: " << locdims[i]
        << "\n  global min: " << gmin << "\n";
      exit(-1);
    }
  }
  
  const TYPE *data = &(*m.grid->cbegin());
  hid_t ret;
  
  /* setup dimensionality object */
  hid_t sid = H5Screate_simple (RANK, dims, NULL);  
  
  /* create a dataset collectively */
  hid_t dataset = H5Dcreate(file_id, 
                            dset_name.c_str(), 
                            H5DataType<TYPE>::type, 
                            sid, 
                            H5P_DEFAULT);

#ifndef SINGLE_PROCESSOR
  /* create a file dataspace independently */
  hid_t file_dataspace = H5Dget_space(dataset);

  ret = H5Sselect_hyperslab(file_dataspace,  H5S_SELECT_SET, 
                              start, NULL, locdims, NULL);
  assert(ret != -1);

  /* create a memory dataspace independently */
  hid_t mem_dataspace = H5Screate_simple (RANK, locdims, NULL);
//  hid_t mem_dataspace = H5Dget_space(dataset);

  /* write data independently */
  ret = H5Dwrite(dataset, 
                 H5DataType<TYPE>::type, 
                 mem_dataspace, 
                 file_dataspace,	    
	               H5P_DEFAULT, 
                 data);

  assert(ret != -1);

  /* release dataspace ID */
  H5Sclose(mem_dataspace);
  H5Sclose(file_dataspace);
#else
  /* write data on single processor */
  hid_t ret = H5Dwrite(dataset, 
                       H5DataType<TYPE>::type, 
                       H5S_ALL, 
                       H5S_ALL,	    
	                     H5P_DEFAULT, 
                       data);
  assert(ret != -1);
#endif			    

  /* close dataset collectively */					    
  ret=H5Dclose(dataset);
  assert(ret != -1);
    
  /* release all IDs created */
  H5Sclose(sid);
    
  return *this;
}

