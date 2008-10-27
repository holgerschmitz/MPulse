#include <sstream>



template<typename TYPE, int RANK, template<int> class Checking>
HDFistream& HDFistream::operator>>(schnek::Matrix<TYPE, RANK, Checking>& m)
{
  std::string dset_name = getNextBlockName();

  long extent[2*RANK];
  
  typename schnek::Matrix<TYPE, RANK, Checking>::IndexType mlow; 
  typename schnek::Matrix<TYPE, RANK, Checking>::IndexType mhigh; 
  
  H5LTget_attribute_long(file_id,dset_name.c_str(),"extent",extent);
  
  for (int i=0; i<RANK; ++i) 
  {
    mlow[i] = extent[2*i];
    mhigh[i] = extent[2*i+1];
//    std::cerr << "Resizing " << i << ": " << mlow[i] << " " << mhigh[i] << "\n";	
  }
  
  m.resize(mlow,mhigh);
  
  TYPE *data = &(*m.begin());

  H5LTread_dataset (file_id, dset_name.c_str(), H5DataType<TYPE>::type, data);
  
  return *this;
}

template<typename TYPE, int RANK, template<int> class Checking>
HDFostream& HDFostream::operator<< (const schnek::Matrix<TYPE, RANK, Checking>& m)
{
  std::string dset_name = getNextBlockName();
  
  typename schnek::Matrix<TYPE, RANK, Checking>::IndexType mdims = m.getDims();
  typename schnek::Matrix<TYPE, RANK, Checking>::IndexType mlow = m.getLow();
  typename schnek::Matrix<TYPE, RANK, Checking>::IndexType mhigh = m.getHigh();
  
  hsize_t dims[RANK];
  long extent[2*RANK];
  
  for (int i=0; i<RANK; ++i) 
  {
    dims[i] = mdims[i];
    extent[2*i] = mlow[i];
    extent[2*i + 1] = mhigh[i];
  }
  
  const TYPE *data = &(*m.cbegin());
    
  H5LTmake_dataset(file_id, dset_name.c_str(), RANK, dims, H5DataType<TYPE>::type, data);
  H5LTset_attribute_long (file_id, dset_name.c_str(), "extent", extent, 2*RANK);
  
  return *this;
}

