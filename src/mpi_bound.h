#ifndef MPI_BOUND_H
#define MPI_BOUND_H

#include "boundary.h"

#ifndef SINGLE_PROCESSOR
#define MPICH_SKIP_MPICXX
// extern "C" {
#include <mpi.h>
// }
//-----------------------------------------------------------------------------
//MPIPeriodicSplitBoundary

/** @brief a boundary class for mutliple processor runs
  *
  * Is design to be exchanged via the MPI protocol
  * This implementation splits only the x-axis into rectangles.
  */
class MPIPeriodicSplitXBoundary : public Boundary {
  protected:
    /// The number of processes
    int ComSize;

    /// The rank of the current process
    int ComRank;

    /// The Comm object referring to the cartesian process grid
    MPI_Comm comm;

    /// The coordinates of this process
    int mycoord;

    int leftcoord;  ///< The rank of the left neighbour process
    int rightcoord; ///< The rank of the right  neighbour process

    /** @brief The size of the array that needs to be exchanged, 
     *  when the exchangeX method is called
     */
    int exchSize;

    double *sendarr; ///< buffer holding the data to be send (size: exchSize)
    double *recvarr; ///< buffer holding the received data (size: exchSize)

    /// The size of the scalar fields when reducing
    int scalarSize;
    
    ///The position of the lower corner of the local piece of the grid
    GridIndex Low;
    ///<The position of the lower corner of the local piece of the dgrid
    GridIndex High;
  public:
    ///default constructor
    MPIPeriodicSplitXBoundary();
      
    /// Virtual destructor deleting all the allocated arrays
    ~MPIPeriodicSplitXBoundary();

    ///initialize
    void init();

    /** @brief Exchanges the boundaries in x-direction.
     *  The two outmost simulated cells are sent and the surrounding 
     *  two ghost cells are filled with values
     */
    void exchangeX(DataGrid3d &field);

    /** @brief Exchanges the boundaries in y-direction.
     *  The two outmost simulated cells are sent and the surrounding 
     *  two ghost cells are filled with values
     */
    void exchangeY(DataGrid3d &field);

    /** @brief Exchanges the boundaries in z-direction.
     *  The two outmost simulated cells are sent and the surrounding 
     *  two ghost cells are filled with values
     */
    void exchangeZ(DataGrid3d &field);

    /** @brief Use MPIALLReduce to calculate the sum and then divide
     *  by the number of processes.
     */
    double AvgReduce(double val) const;

    /** @brief Use MPIALLReduce to calculate the maximum
     */
    double MaxReduce(double val) const;
      
    /// Use MPIALLReduce to calculate the maximum
    double SumReduce(double val) const;

    /// Returns the global lower bound of the distribution function
    const GridIndex &RegionLow()const;

    /// Returns the global upper bound of the distribution function
    const GridIndex &RegionHigh()const;

    /// The process with the rank zero is designated master process
    bool master() const { return ComRank==0; }

    /// Returns the comm rank as given by mpi
    int procnum() const { return ComRank; }

    ///returns an ID, which is identical with the coordinates
    int getUniqueId() const { return mycoord; }
};
//MPIPeriodicSplitXBoundary
//-----------------------------------------------------------------------------
//MPIPeriodicSplitXYBoundary

/** @brief a boundary class for mutliple processor runs
 *
 * Is designed to be exchanged via the MPI protocol.
 * Here splitting is performed in both spatial directions.
 */
class MPIPeriodicSplitXYZBoundary : public Boundary {
  protected:
    /// The number of processes
    int ComSize;
  
    /// The rank of the current process
    int ComRank;
    
    /// The Comm object referring to the cartesian process grid
    MPI_Comm comm;

      
    int xprevcoord; ///< The rank of the left neighbour process in x
    int xnextcoord; ///< The rank of the right neighbour process in x
      
    int yprevcoord; ///< The rank of the left neighbour process in y
    int ynextcoord; ///< The rank of the right neighbour process in y

    int zprevcoord; ///< The rank of the left neighbour process in z
    int znextcoord; ///< The rank of the right neighbour process in z
    
    ///dimensions
    int dims[3];
    /// The cartesian coordinates of this process
    int mycoord[3];

    /** @brief The size of the array that needs to be exchanged, 
     *  when the exchangeX or the exchangeY method is called
     */
    int exchSize[3];
      
    double *sendarrx; ///< send Buffer for exchanging data in x-direction (size: exchSize[0])
    double *recvarrx; ///< receive Buffer for exchanging data in x-direction (size: exchSize[0])
      
    double *sendarry; ///< send Buffer for exchanging data in y-direction (size: exchSize[1])
    double *recvarry; ///< receive Buffer for exchanging data in y-direction (size: exchSize[1])
      
    double *sendarrz; ///< send Buffer for exchanging data in z-direction (size: exchSize[2])
    double *recvarrz; ///< receive Buffer for exchanging data in z-direction (size: exchSize[2])

    /// The size of the scalar fields when reducing
    int scalarSize;
      
    /// The positions of the lower corner of the local piece of the grid
    GridIndex Low;
    ///<The positions of the upper corner of the local piece of the grid
    GridIndex High;
    
  public:
    ///default constructor
    MPIPeriodicSplitXYZBoundary();
      
    /// Virtual destructor deleting all the allocated arrays
    ~MPIPeriodicSplitXYZBoundary();
      
    ///initialize
    void init();

    /** @brief Exchanges the boundaries in x-direction.
     *
     *  The two outmost simulated cells are sent and the surrounding 
     *  two ghost cells are filled with values
     */
    void exchangeX(DataGrid3d &field);
      
    /** @brief Exchanges the boundaries in y-direction.
     *
     *  The two outmost simulated cells are sent and the surrounding 
     *  two ghost cells are filled with values
     */
    void exchangeY(DataGrid3d &field);
      
    /** @brief Exchanges the boundaries in z-direction.
     *
     *  The two outmost simulated cells are sent and the surrounding 
     *  two ghost cells are filled with values
     */
    void exchangeZ(DataGrid3d &field);

    /// Use MPIALLReduce to calculate the sum and then divide by the number of processes. 
    double AvgReduce(double val) const;
      
    /// Use MPIALLReduce to calculate the maximum
    double MaxReduce(double val) const;
      
    /// Use MPIALLReduce to calculate the maximum
    double SumReduce(double val) const;
      
    /// Returns the global lower bound of the distribution function
    const GridIndex &RegionLow()const;
      
    /// Returns the global upper bound of the distribution function
    const GridIndex &RegionHigh()const;
      
    /// The process with the rank zero is designated master process
    bool master() const { return ComRank==0; }
      
    /// Returns the comm rank as given by mpi
    int procnum() const { return ComRank; }
    ///returns an ID, which consists of the Dimensions and coordinates
    int getUniqueId() const { 
      return dims[2]*(dims[1]*mycoord[0] + mycoord[1]) + mycoord[2]; 
    }
    /// returns "true" since this is periodic in the x-direction
    virtual bool periodicX() { return true; } 
    /// returns "true" since this is periodic in the y-direction
    virtual bool periodicY() { return true; }
};
//MPIPeriodicSplitXYBoundary
//-----------------------------------------------------------------------------
#endif // multiple processor

#endif
