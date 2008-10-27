#ifndef MPULSE_PARAMETER_H
#define MPULSE_PARAMETER_H

/** @file parameter.h
 * @brief Classes reading input files, storing parameters
 *
 *  The classes in this file are used for reading parameters from
 *  an input file. The parameters are stored in a parameter map
 */
//----------------------------------------------------------------------------
#include <map>
#include <string>
#include <iostream>
#include "stlpwrapper.h"
//----------------------------------------------------------------------------
//forward declarations
class Rebuildable;
class Parameter;
//----------------------------------------------------------------------------
///wrapped pointer to parameter
typedef PtrWrapper<Parameter> WParameter;
/// map strings against wrapped pointers to parameter ojects (used by all derivations of of Rebuildable)
typedef std::map<std::string,WParameter> ParameterMap;
//----------------------------------------------------------------------------
//Parameter
/** @brief An abstract base class for parameters.
 *
 *  There is only one abstract method: Rebuild. This method has to rebuild 
 *  the data from the input file and return the next token of the
 *  in a string.
 *  @todo Is the < operator really needed? 
 */
class Parameter {
  public:
    /// Default constructor
    Parameter () {};
    /// Virtual destructor
    virtual ~Parameter () {};
   
    /** @brief This method has to rebuild the data
     *
     *  from the input file and return the next token of th input file
     *  in a string.
     *
     * This method needs to be overwritten by derived classes.
     */
    virtual std::string Rebuild (std::istream& in) = 0; 

    /// Always returns true since sorting is irrelevant
    friend bool operator< (const Parameter& left, const Parameter& right) {
      return true; 
    };
};

//-----------------------------------------------------------------------------
//ParameterValue
/** @brief Single value parameter of some type.
 *
 *  This reads in a single value of some type TYPE from the input file
 *  and stores it under a given address. If the parameter is not encountered
 *  in the input, a default value is set.
 */
template<class TYPE>
class ParameterValue : public Parameter {
  private:
    TYPE* pValue;   ///< Pointer to the variable to set
    TYPE Default;   ///< A default value
  public:
    /// Constructor, takes a pointer to the variable and a default value
    ParameterValue (TYPE* _pValue, TYPE _Default) : pValue(_pValue) { 
      SetDefault(_Default);
    };
    /// Virtual destructor
    virtual ~ParameterValue () {};
    
    /// Sets the default value for TYPE
    void SetDefault (TYPE _Default) { 
      Default = _Default;
      *pValue = Default; 
    };
    
    /// Reads the value for TYPE from the stream and returns the next token
    virtual std::string Rebuild(std::istream& in);
}; 

//----------------------------------------------------------------------------
//ParameterRebuild
/** @brief Base class for reading a subclass of the Rebuildable class from the input. 
 *
 *  A Rebuildable class is read from the input file by invoking the creating the
 *  object and calling its Rebuildable::Rebuild method. The parameters of the newly 
 *  built Rebuildable have to be enclosed in curly brackets.
 */
template<class Type, class BaseType>
class ParameterRebuild : public Parameter {
  public:
  protected:
    BaseType **value; ///< Pointer to the pointer of the object
    typedef std::list<BaseType*> BaseList;
    BaseList *values;
  public:
    /// Constructor takes a pointer to the parent Rebuildable
    ParameterRebuild (BaseType **value_) 
      : value(value_),
        values(NULL)
    { 
      *value = NULL; 
    }

    /// Constructor takes a pointer to the parent Rebuildable
    ParameterRebuild (BaseList *values_) 
      : value(NULL),
        values(values_) {}
    
    /// Destructor
    virtual ~ParameterRebuild () {}

    /** @brief Creates a new Rebuildable, calls its 
     * Rebuildable::Rebuild method 
     */
    virtual std::string Rebuild (std::istream& in);
    /// Abstract method that should return a pointer to a new task
    virtual BaseType* NewInstance ()  { return new Type; }
}; 

//-----------------------------------------------------------------------------
//including the implementations for template members
#include "parameter.t"
//-----------------------------------------------------------------------------
#endif // PARAMETER_H
