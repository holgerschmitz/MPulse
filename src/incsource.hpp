#ifndef MPULSE_INCSOURCE_H
#define MPULSE_INCSOURCE_H

#include "mpulse.hpp"
#include "current.hpp"

class IncidentSourceCurrent;


//===============================================================
//==========  Base Classes
//===============================================================

/**
 * An abstract base class for sources on the simulation boundary
 *
 * Sources for time-dependent electromagnetic fields are represented by
 * transverse currents on the boundary.
 */
class IncidentSource : public CurrentBlock
{
  public:
    /**
     * Initialise any currents
     *
     * This will check all 6 faces of the boundary box. If #needCurrent indicates
     * that currents are needed on that boundary, #makeECurrent and #makeHCurrent
     * are used to create the transverse currents.
     *
     * The current created in this way will usually extend #IncidentSourceCurrent
     */
    void initCurrents(CurrentContainer &container);

  protected:

    /**
     * Create an #IncidentSourceECurrent for a single boundary at a given distance
     */
    virtual pCurrent makeECurrent(int distance_, Direction dir_) = 0;

    /**
     * Create an #IncidentSourceBCurrent for a single boundary at a given distance
     */
    virtual pCurrent makeHCurrent(int distance_, Direction dir_) = 0;

    /**
     * Initialise the setup parameters
     */
    void initParameters(schnek::BlockParameters &blockPars);
  private:

    /**
     * The distance of the currents from the outer boundary of the simulation box
     */
    int distance;
};

/**
 * An incident source current
 *
 * Incident source currents are intended to be used with #IncidentSource current
 * blocks. They define a transverse current on a boundary plane.
 */
class IncidentSourceCurrent : public Current
{
  public:

    /**
     * Initialise the incident source current
     *
     * @param distance_  the distance of the injection plane from the simulation
     *                   border
     * @param dir_       the direction of the border from which the wave should
     *                   be injected
     * @param isH_       a flag indicating whether the current is a magnetic
     *                   current.
     */
    IncidentSourceCurrent(int distance_, Direction dir_, bool isH_, SimulationContext &context);

  protected:
    /**
     * A flag indicating whether the wave should propagate in the negative
     * coordinate direction
     */
    bool reverse;

    /**
     * The distance of the injection plane from the simulation
     * border
     */
    int distance;

    /**
     * The propagation axis
     */
    int dim;

    /**
     * The transverse axes
     */
    int transverse1, transverse2;

    /**
     * The direction of the border from which the wave should be injected
     */
    Direction dir;

    /**
     * A flag indicating whether the current is a magnetic current.
     */
    bool isH;
    int lowOffset;
    int highOffset;

    /**
     * The grid spacing
     */
    double dx, dy, dz;

    /**
     * The time step
     */
    double dt;

    /**
     * References to the transverse components of the current from #Current.pJx,
     * #Current.pJy, and #Current.pJz
     */
    pGrid pJ[2];

    /**
     * A re-ordered vector of grid spacings in the local coordinate system in
     * which the propagation direction is along the positive z-axis
     */
    double dX[3];

    /**
     * The simulation context
     */
    SimulationContext &context;
};

/**
 * Specialisation of #IncidentSourceCurrent for electric currents
 *
 * This function is templated with a source function that should provide the
 * value of the **magnetic field**. `SourceFunc` must expose a method
 * `getHField` to obtain this field and a method `initSourceFunc` to initialise
 * the source function.
 */
template<class SourceFunc>
class IncidentSourceECurrent : public IncidentSourceCurrent, public SourceFunc
{
  public:

    /**
     * Create a new IncidentSourceECurrent instance
     *
     * @param distance_  the distance of the injection plane from the simulation
     *                   border
     * @param dir_       the direction of the border from which the wave should
     *                   be injected
     */
    IncidentSourceECurrent(int distance_, Direction dir_, SimulationContext &context);

    void init();
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:

};

/**
 * Specialisation of #IncidentSourceCurrent for magnetic currents
 *
 * This function is templated with a source function that should provide the
 * value of the **electric field**. `SourceFunc` must expose a method
 * `getEField` to obtain this field and a method `initSourceFunc` to initialise
 * the source function.
 */
template<class SourceFunc>
class IncidentSourceHCurrent : public IncidentSourceCurrent, public SourceFunc
{
  public:

    /**
     * Create a new IncidentSourceBCurrent instance
     *
     * @param distance_  the distance of the injection plane from the simulation
     *                   border
     * @param dir_       the direction of the border from which the wave should
     *                   be injected
     */
    IncidentSourceHCurrent(int distance_, Direction dir_, SimulationContext &context);

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);
};

//===============================================================
//==========  Side Inject
//===============================================================

///** Reads the data from a file, maybe scaling the
// *  amplitude by a constant factor
// */
//class SideInject : public IncidentSource
//{
//  public:
//    virtual ~SideInject() {}
//
//  protected:
//    virtual pCurrent makeECurrent(int distance_, Direction dir_);
//    virtual pCurrent makeHCurrent(int distance_, Direction dir_);
//    virtual bool needCurrent(Direction dir_);
//
//    void initParameters(schnek::BlockParameters &blockPars);
//  private:
//    std::string filename;
//    double amp;
//    int toff;
//    int blocks;
//    double eps;
//};
//
//class SideInjectSourceFunc
//{
//  public:
//    SideInjectSourceFunc(Direction dir_, bool isH_);
//    void setParam(std::string filename_, double amp_, double eps_, int toff_, int blocks_);
//
//    void initSourceFunc(Storage *storage, DataGrid *pJx, DataGrid *pJy, DataGrid *pJz);
//
//    Vector getEField(int i, int j, int k, int time);
//    Vector getHField(int i, int j, int k, int time);
//
//    void setTime(int Time);
//  private:
//    Vector getField(int i, int j, int k, int time, double factor);
//    DataGrid F1;
//    DataGrid F2;
//
//    std::string filename;
//    double amp;
//    double eps;
//    int toff;
//    int Nt;
//    int blockCount;
//    int blocks;
//    int blockToff;
//    bool active;
//
//    int dim;
//    int transverse1, transverse2;
//
//    Direction dir;
//    bool isH;
//};


#include "incsource.t"


#endif
