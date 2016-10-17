/*
 * mpulse.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include <schnek/grid.hpp>

#ifdef NDEBUG
#define MPulseGridChecker schnek::GridNoArgCheck
#else
#define MPulseGridChecker schnek::GridAssertCheck
#endif

static const size_t DIMENSION = 3;

typedef schnek::Array<int, DIMENSION> Index;
typedef schnek::Array<double, DIMENSION> Vector;

typedef schnek::Field<double, DIMENSION> Field;
typedef boost::shared_ptr<Field> pField;

typedef schnek::Field<double, 1> DataLine;
typedef boost::shared_ptr<DataLine> pDataLine;

typedef schnek::Range<int, DIMENSION> Range;
typedef schnek::Array<bool, DIMENSION> Stagger;

static const double clight = 299792458;
static const double clight2 = clight*clight;

static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true );

static const Stagger bxStaggerYee(false, true,  true );
static const Stagger byStaggerYee(true,  false, true );
static const Stagger bzStaggerYee(true,  true,  false);

class MPulse : public schnek::Block, public schnek::BlockContainer<FieldSolver>
{
  private:
    static MPulse *instance;
    Index globalMax;
    Index gridSize;
    Vector size;
    Vector dx;

    double cflFactor;
    double dt;


    double tMax;
    Range innerRange;
    pField Ex, Ey, Ez;
    pField Bx, By, Bz;

    schnek::MPICartSubdivision<Field> subdivision;

    Vector x;
    schnek::Array<schnek::pParameter, DIMENSION> x_parameters;
    schnek::Array<schnek::pParameter, DIMENSION> E_parameters;
    schnek::Array<schnek::pParameter, DIMENSION> B_parameters;
    schnek::pParametersGroup spaceVars;

    Vector initE;
    Vector initB;
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void registerData();
    void fillValues();
  public:
    MPulse();
    void init();
    void execute();

    static Index getGlobalMax() { return instance->globalMax; }
    static Vector getDx() { return instance->dx; }
    static Vector getSize() { return instance->size; }
    static schnek::DomainSubdivision<Field> &getSubdivision() { return instance->subdivision; };
};

#endif // MPULSE_MPULSE_H
