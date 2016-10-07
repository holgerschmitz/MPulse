#ifndef MPULSE_FOLLOWMAX_H
#define MPULSE_FOLLOWMAX_H

#include "mpulse.h"
#include "storage.h"
#include "fieldsim.h"

#include <fstream>
#include <string>

class FollowMax : public FieldExtraDiag {
  private:
	  /// Type of field
    Storage *storage;
    
    double position;
    int iposition;
    
    int start;
    
    int posx, posy, posz;
    std::string direction;
    
    int dim;
    int trans1, trans2;
    int low, high;
    
    bool active;
	      
    std::string field;
    
    /// Stream to write to
    std::ofstream output;
  public:
	  ///default constructor
    FollowMax();
	  
    ///destrcutor
    ~FollowMax();
	  
    void init();
    
    ///set the pointer to the field
    void setStorage(Storage *storage);
  protected:
	  ///open a file stream, 
    void open(const std::string &);
	  
    ///write diagnostics
    void write();
	  
    ///close the file 
    void close();
	  
    ///returns the single_out member
    bool singleOut() { return true; }
    
    /// calculate the field energy
    void calculate();
    
    ParameterMap* MakeParamMap (ParameterMap* pm);
};

#endif
