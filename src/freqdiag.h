#ifndef MPULSE_FREQDIAG_H
#define MPULSE_FREQDIAG_H

#include "mpulse.h"
#include "storage.h"
#include "fieldsim.h"

#include <fstream>
#include <string>

class FrequencyDiag : public FieldExtraDiag {
  private:
	  /// Type of field
    Storage *storage;
    
    double frequency;
    int count;
    int lastcount;
        
    int x, y, z;
    double dt;
    
    double lastval;
    double firstzero;
	      
    std::string field;
    
    /// Stream to write to
    std::ofstream output;
	  
  public:
	  ///default constructor
    FrequencyDiag();
	  
    ///destrcutor
    ~FrequencyDiag();
	  
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
