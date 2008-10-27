#ifndef PERIODIC_BOUND_H
#define PERIODIC_BOUND_H

#include "boundary.h"

//-----------------------------------------------------------------------------
//SinglePeriodicBoundary

/** @brief Implements Boudary to supply a periodic system running on
 *  a single processor
 */
class SinglePeriodicBoundary : public Boundary {
  public:   
    /// No initialization needs to be done
    void init() {};
    /// Wraps the boundaries in x-direction
    void exchangeX(DataGrid3d &field);
    /// Wraps the boundaries in y-direction
    void exchangeY(DataGrid3d &field);
    /// Wraps the boundaries in z-direction
    void exchangeZ(DataGrid3d &field);
        
    /// There is no average to be calculated
    double AvgReduce(double val) const { return val; }
    /// There is no sum to be calculated
    double SumReduce(double val) const { return val; }
    /// There is no maximum to be calculated
    double MaxReduce(double val) const { return val; }
    /// Returns the global lower bound of the distribution function
    const GridIndex &RegionLow() const;
    /// Returns the global upper bound of the distribution function
    const GridIndex &RegionHigh() const;
    /// There is only one process, so master always returns true
    bool master() const { return true; }
    /// The process number is always zero
    int procnum() const { return 0; }
    /// The unique id number is always zero
    int getUniqueId() const { return 0; }
};
//SinglePeriodicBoundary

/** @brief Implements Boudary to supply a periodic system running on
 *  a single processor
 */
class SingleXYPeriodicBoundary : public Boundary {
  public:   
    /// No initialization needs to be done
    void init() {};
    /// Wraps the boundaries in x-direction
    void exchangeX(DataGrid3d &field);
    /// Wraps the boundaries in y-direction
    void exchangeY(DataGrid3d &field);
    /// Wraps the boundaries in z-direction
    void exchangeZ(DataGrid3d &field);
        
    /// There is no average to be calculated
    double AvgReduce(double val) const { return val; }
    /// There is no sum to be calculated
    double SumReduce(double val) const { return val; }
    /// There is no maximum to be calculated
    double MaxReduce(double val) const { return val; }
    /// Returns the global lower bound of the distribution function
    const GridIndex &RegionLow() const;
    /// Returns the global upper bound of the distribution function
    const GridIndex &RegionHigh() const;
    /// There is only one process, so master always returns true
    bool master() const { return true; }
    /// The process number is always zero
    int procnum() const { return 0; }
    /// The unique id number is always zero
    int getUniqueId() const { return 0; }
};
//SinglePeriodicBoundary

#endif
