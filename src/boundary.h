/** @file boundary.h
 * @brief Interface for wrapping and exchanging boundaries
 *
 *  This interface is used to exchange the boundaries of distribution
 *  functions and scalar fields. It can be implemented to define
 *  periodic boundaries or exchange data with other processes.
 *
 *  See the page @ref indices for a discussion on the numerical
 *  ranges of the fields
 */
//-----------------------------------------------------------------------------
#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "mpulse.h"
#include "globals.h"
#include "rebuild.h"

//-----------------------------------------------------------------------------
//Boundary

/** @brief Interface for wrapping and exchanging boundaries .
 *
 *  This interface is used to exchange the boundaries of distribution
 *  functions and scalar fields. It can be implemented to define
 *  periodic boundaries or exchange data with other processes.
 *  This is the (abstract) base class of all boundary classes.
 *
 *  See the page @ref indices for a discussion on the numerical
 *  ranges of the fields
 */
class Boundary : public Rebuildable {
  public:
    /// Default constructor, 
    Boundary() {}
    /** @brief Virtual destructor
     *  
     *  We need a virtual destructor because the class has 
     *  virtual methods
     */
    virtual ~Boundary() {}
    
    /** Initialize the boundary
     */
    virtual void init() = 0;

    /** @brief Exchange the boundaries of a field function
     *  in the x-direction
     */
    virtual void exchangeX(DataGrid3d &field) = 0;

    /** @brief Exchange the boundaries of a field function
     *  in the y-direction
     */
    virtual void exchangeY(DataGrid3d &field) = 0;

    /** @brief Exchange the boundaries of a field function
     *  in the z-direction
     */
    virtual void exchangeZ(DataGrid3d &field) = 0;

    /// Return the average of a single value over all the processes
    virtual double AvgReduce(double) const = 0;

    /// Return the average of a single value over all the processes
    virtual double SumReduce(double) const = 0;

    /// Return the maximum of a single value over all the processes
    virtual double MaxReduce(double) const = 0;

    /// Return the lower bound of the field
    virtual const GridIndex &RegionLow() const = 0;

    /// Return the upper bound of the field
    virtual const GridIndex &RegionHigh() const = 0;

    /// Return true if this is the master process and false otherwise
    virtual bool master() const = 0;

    /// Return the process number
    virtual int procnum() const = 0;

    ///get a unique Id
    virtual int getUniqueId() const = 0;
};

//Boundary



#endif //BOUNDARY_H
