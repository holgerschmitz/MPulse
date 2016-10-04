#ifndef MPULSE_REBUILD_H
#define MPULSE_REBUILD_H

/** @file task.h
 *  @brief the fundametal class Rebuildable
 *
 *  Contains the declarations of the fundamental Rebuildable class 
 */
//-----------------------------------------------------------------------------
#include <list>
#include <vector>
#include <string>
#include "parameter.h"
//-----------------------------------------------------------------------------

/** @brief fundamental base class for all the rebuildable objects
 *
 * All rebuildable objects inherit from this base class. 
 * 
 */
class Rebuildable {
  public:
    ///default constructor
    Rebuildable () {}
    /// virtual destructor
    virtual ~Rebuildable () {} 

  protected:
    /** @brief Register the parameters needed for the Task
     *
     *  A derived class should overwrite this whenever it needs additional 
     *  parameters from the setup file. It should then ALWAYS call the 
     *  MakeParamMap of its superclass. 
     *  If none present allocate new parameter map and return pointer to it.
     *  If a existing map is passed, do nothing.
     */
    virtual ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  public:
    /** @brief Rebuilds the task from the setup. This normally does not need 
     *  to be overwritten.
     *
     *  Rebuild always returns the next token that does not belong to the
     *  object setup. Calls the rebuild methods of the objects stored in the parameter map.
     *  Will also remove all comments from the instream.
     */
    virtual std::string Rebuild (std::istream& in);
    ///virtual method,  does nothing unless overwritten by derived classes
    virtual void finalize() {};
}; 

#endif // MPULSE_REBUILD_H

