#ifndef MPULSE_DIAGNOSTIC_H
#define MPULSE_DIAGNOSTIC_H

#include <iostream>
#include <list>
#include <string>
#include "stlpwrapper.h"
#include "rebuild.h"
//-----------------------------------------------------------------------------
/** @file diagnostic.h 
 * @brief Interface for diagnostic tasks.
 * 
 * This interface can be used to implement different types of diagnostics.
 * The simplest form is implemented as SimpleDiagnostic, which writes the fields
 * to an output stream of arbitrary type.
 * The DiagnosticInterface is closely related to the DiagnosticManager. When an
 * instance of the interface is created it will register itself with the 
 * DiagnosticManager. This then takes the responsibility of calling the execute
 * method of the DiagnosticInterface.
 */
 
//-----------------------------------------------------------------------------
//---------   DiagnosticInterface
//-----------------------------------------------------------------------------

/** @brief Interface for diagnostic tasks.
 * 
 * This interface can be used to implement different types of diagnostics.
 * The DiagnosticInterface is closely related to the DiagnosticManager. When an
 * instance of the interface is created it will register itself with the 
 * DiagnosticManager. This then takes the responsibility of calling the execute
 * method of the DiagnosticInterface.
 * This the abstract base class for all Diagnostic Interface Types.
 */
class DiagnosticInterface : public Rebuildable {
  private:
    /// The file name into which to write
    std::string fname;
      
    /** A string specifying whether to append to the file or whether to
     * write into a new file at each turn.
     * Only the first letter is checked for equality to 'y'
     */
    std::string append;
	  
    /// each interval timesteps the actual output is performed
    int interval;
	  
    ///timestep
    int t;
  public:
	  ///constructor
    DiagnosticInterface();
	  
    ///execute, calls the open, write and close method to perform the output
    void execute();
  protected:
	  ///prototype of the open method
    virtual void open(const std::string &)=0;
	  
    ///prototype of the write method
    virtual void write()=0;
	  
    ///prototype of the close method
    virtual void close()=0;
	  
    ///returns "false" if not overwritten by the derived class
    virtual bool singleOut() { return false; }
    
    /** can be overwritten do perform some calculation on all
     *  processes, independent on singleOut
     */
    virtual void calculate() {}
	  
    ///create the parameter map
    virtual ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
	  ///if "true" append the output to the stream, if "false" each  write to a new file
    bool appending();
    
	  ///returns the name of the output file in a string
    std::string parsedFileName();
};
//DiagnosticINterface
//-----------------------------------------------------------------------------
///wrapped pointer to DignosticInterfaces
typedef PtrWrapper<DiagnosticInterface> pDiagnosticInterface;
///list of wrapped pointers to DiagnosticInterfaces which are registered in the manager
typedef std::list<pDiagnosticInterface> DiagList;
//-----------------------------------------------------------------------------
//DiagnosticManager
///@todo meyers singleton?! see the .cpp for this idea
/** @brief The diagnostic manager class
  *
  * Every instance of DiagnosticInterface registers here. 
  * DiagnosticManager then takes care of calling the execute members of each interface.
  */
class DiagnosticManager {
  private:
	//should be deleted if meyers singleton is put to work
	///the one and only instance of diagnosticManager 
      static DiagnosticManager *theManager;
	///list of the registered interfaces (pointers to these)
      DiagList diags;
  public:
	///returns a reference to this instance
      static DiagnosticManager& instance();
	///add a interface to the list
      void addDiagnostic(DiagnosticInterface*);
	///call execute for each item of the list
      void execute();
  private:
	///default constructor, is private to enforce singleton behavior
      DiagnosticManager();
	///copy constructor, is private to enforce singleton behavior
      DiagnosticManager(DiagnosticManager&);
};
//DiagnosticManager
//-----------------------------------------------------------------------------
//SimpleDiagnsotic for writing to streams

/** @brief a simple diagnostic interface, derived from DiagnosticInterface
  *
  * This can be used for writing results to files or streams for further use.
  */
template<class Type, class StreamType>
class SimpleDiagnostic : public DiagnosticInterface {
  private:
	  /// Type of field
    Type *field;
	  
    /// Stream to write to
    StreamType output;
	  
    ///"true" if diagnostic is performed globaly
    bool single_out;
  public:
	  ///default constructor
    SimpleDiagnostic() { single_out=false; }
	  
    ///destrcutor
    ~SimpleDiagnostic();
	  
    ///set the pointer to the field
    void setField(Type*);
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
};
//SimpleDiagnostic
//-----------------------------------------------------------------------------
// including the template members
#include "diagnostic.t"
//-----------------------------------------------------------------------------
#endif //DIAGNOSTIC_H
