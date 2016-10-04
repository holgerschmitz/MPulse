#ifndef MPULSE_FIELDDIAG_H
#define MPULSE_FIELDDIAG_H

#include "mpulse.h"
#include "hdfstream.h"
#include "diagnostic.h"
#include "storage.h"

#include <fstream>

class FieldDiag : public SimpleDiagnostic<DataGridContainer,HDFostream>
//class FieldDiag : public SimpleDiagnostic<DataGridContainer,std::ofstream>
{
  public:
    void fetchField(Storage &storage);
  private:
    std::string fieldId;
    DataGridContainer field;
  protected:
    typedef SimpleDiagnostic<DataGridContainer,HDFostream> ParentType;
//    typedef SimpleDiagnostic<DataGridContainer,std::ofstream> ParentType;
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);  
};

class FieldSliceDiag : public DiagnosticInterface
// : public SimpleDiagnostic<DataGrid2dContainer,HDFostream>
{
// methods and properties formerly implemented in SimpleDiagnostic

  private:	  
    /// Stream to write to
    HDFostream output;
	  
    ///"true" if diagnostic is performed globaly
    bool single_out;
  public:
	  ///default constructor
    FieldSliceDiag() { single_out=false; }
	  
    ///destrcutor
    ~FieldSliceDiag();
    
  protected:
	  ///open a file stream, 
    void open(const std::string &);
	  
    ///write diagnostics
    void write();
	  
    ///close the file 
    void close();
	  
    ///returns the single_out member
    bool singleOut() { return single_out; }
  public:
	  ///set the single_out member
    void setSingleOut(bool single_out_) { single_out = single_out_; }

// methods and properties inherently part of FieldSliceDiag

  public:
    void fetchField(Storage &storage);
  private:
    std::string fieldId;
    DataGrid2d slice;
    DataGrid2dContainer sliceContainer;
    DataGrid *field;
    
    int pos;
    std::string plane;
    
    int normal;
    int dim1, dim2;
    int low1, low2;
    int high1, high2;
    
    bool active;
  protected:
    typedef SimpleDiagnostic<DataGrid2dContainer,HDFostream> ParentType;
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);


};


class FieldExtraDiag : public DiagnosticInterface {
  public:
    virtual void setStorage(Storage *storage) = 0;
    virtual void init() {}
};

class FieldEnergyDiag : public FieldExtraDiag {
  private:
	  /// Type of field
    Storage *storage;
    
    double energy;
	  
    /// Stream to write to
    std::ofstream output;
	  
  public:
	  ///default constructor
    FieldEnergyDiag();
	  
    ///destrcutor
    ~FieldEnergyDiag();
	  
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
};


struct AllFieldDiag
{
  AllFieldDiag() {}
  typedef std::list<FieldDiag*> DiagList;
  typedef std::list<FieldSliceDiag*> SliceDiagList;
  typedef std::list<FieldExtraDiag*> ExtraDiagList;
  DiagList fields;
  SliceDiagList slices;
  ExtraDiagList fieldextras;
};

#endif
