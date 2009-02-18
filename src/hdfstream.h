#ifndef HDFSTREAM_H
#define HDFSTREAM_H
//-----------------------------------------------------------------------------

#include <hdf5.h>
#include <iostream>

#include <schnek/matrix.h>
#include "mpulse.h"


//-----------------------------------------------------------------------------
  /** @file hdfstream.h
    * @brief IO class for HDF
    *
    * IO classes for handling HDF files and streams
    */
//-----------------------------------------------------------------------------
//HDFstream

template<typename TYPE, int RANK, template<int> class Checking>
struct MatrixContainer
{
  bool active;
  schnek::Matrix<TYPE, RANK, Checking> *grid;
  typename schnek::Matrix<TYPE, RANK, Checking>::IndexType global_min;
  typename schnek::Matrix<TYPE, RANK, Checking>::IndexType global_max;
};

typedef MatrixContainer<double, 3, MPulseGridChecker> DataGridContainer;
typedef MatrixContainer<double, 2, MPulseGridChecker> DataGrid2dContainer;

std::ostream &operator<<(std::ostream& out, const DataGridContainer &);


/** @brief IO class for handling HDF files
  *  
  * This is the abstract base class for HDF-IO- classes.
  * Implements the basic operations on HDFstreams as virtual methods.
  */
class HDFstream {
  protected:
    /// HDF5 File id
    hid_t       file_id; 
    
    /// HDF5 Error status
    herr_t      status;
    
    /// name of the datablock to be read or written
    std::string blockname;
    /// counter for the sets with a given blockname read from or written to the file
    int sets_count;  
    
    /// Specifies if the stream is active in this process (in case of parallel execution)
    bool active;
    bool activeModified;
    
  public:
    /// constructor 
    HDFstream();
    ///copy constructur 
    HDFstream(const HDFstream&);
    /// destructor
    virtual ~HDFstream();

    /// open  file 
    virtual int open(const char*)=0;
    
    /// close file 
    virtual void close();
    
    /// return true=1 if data are still available 
    virtual bool good() const;
    
    void setBlockName(std::string blockname_);
    /// assign 
    HDFstream& operator = (const HDFstream&);
    
    void setActive(bool active_) { active = active_; activeModified = true; }
    
  protected:
    std::string getNextBlockName();
    
#ifndef SINGLE_PROCESSOR
    void makeMPIGroup();
  
    MPI_Comm mpiComm;
    bool commSet;
#endif

};

//HDFstream
//---------------------------------------------------------------
//HDFistream
/** @brief Input stream for HDF files */
class HDFistream : public HDFstream {
  public:
    /// constructor 
    HDFistream();
    
    /// copy constructor */
    HDFistream(const HDFistream&);
    
    /// constructor, opens HDF file "fname", selects first dataset 
    HDFistream(const char* fname);

    /// opens HDF file "fname", selects first dataset 
    virtual int open(const char*);
       
    /// stream input operator for a schnek::Matrix 
    template<typename TYPE, int RANK, template<int> class Checking>
    HDFistream& operator>>(MatrixContainer<TYPE, RANK, Checking>& m);
};
//HDFistream
//-----------------------------------------------------------------------------
//HDFostream
/** @brief output stream for HDF files */
class HDFostream : public HDFstream {
  public:
    /// constructor 
    HDFostream();
    
    /// copy constructor 
    HDFostream(const HDFostream&);
    
    /// constructor, opens HDF file "fname" 
    HDFostream(const char* fname);

    /// open file 
    virtual int open(const char*);
    
    /// stream output operator for a matrix
    template<typename TYPE, int RANK, template<int> class Checking>
    HDFostream& operator<< (const MatrixContainer<TYPE, RANK, Checking>&);    
};

template<typename TYPE>
struct H5DataType{
  static const hid_t type;
};



//HDFostream
//-----------------------------------------------------------------------------

#include "hdfstream.t"

//-----------------------------------------------------------------------------
#endif // HDFSTREAM_H
